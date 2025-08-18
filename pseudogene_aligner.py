#!/usr/bin/env python3
import argparse, sys, re
from time import sleep
from Bio import Entrez, SeqIO, pairwise2
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ------------- GRCh38 chromosome -> RefSeq accession -------------
CHR_TO_REFSEQ_GRCh38 = {
    "1":"NC_000001.11","2":"NC_000002.12","3":"NC_000003.12","4":"NC_000004.12",
    "5":"NC_000005.10","6":"NC_000006.12","7":"NC_000007.14","8":"NC_000008.11",
    "9":"NC_000009.12","10":"NC_000010.11","11":"NC_000011.10","12":"NC_000012.12",
    "13":"NC_000013.11","14":"NC_000014.9","15":"NC_000015.10","16":"NC_000016.10",
    "17":"NC_000017.11","18":"NC_000018.10","19":"NC_000019.10","20":"NC_000020.11",
    "21":"NC_000021.9","22":"NC_000022.11","X":"NC_000023.11","Y":"NC_000024.10","MT":"NC_012920.1"
}

# ---------------- helpers ----------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Align parent gene mRNA to a pseudogene (FASTA or genomic coords), report per-exon homology, and plot gapped tracks."
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--accession", help="Parent RefSeq mRNA accession (e.g., NM_020436.5)")
    src.add_argument("--gene", help="Parent gene symbol (resolves to a human RefSeq mRNA)")

    pg = p.add_mutually_exclusive_group(required=True)
    pg.add_argument("--pseudo-fasta", help="Path to pseudogene FASTA (e.g., pseudogene.fasta)")
    pg.add_argument("--pseudo-coords", help="GRCh38 genomic coords, e.g., chr5:37209364-37212378[:strand]")

    p.add_argument("--build", default="GRCh38", help="Genome build for --pseudo-coords (default: GRCh38)")
    p.add_argument("--prefix", required=True, help="Output prefix (filenames will start with this)")
    p.add_argument("--email", required=True, help="Your email for NCBI Entrez")
    p.add_argument("--ssl-ignore", action="store_true", help="TEMP: ignore SSL cert errors if needed")
    return p.parse_args()

def resolve_accession_from_gene(gene_symbol, email):
    Entrez.email = email
    term = f'{gene_symbol}[Gene] AND Homo sapiens[Organism] AND refseq[filter] AND mRNA[Filter]'
    h = Entrez.esearch(db="nucleotide", term=term, retmax=50)
    ids = Entrez.read(h)["IdList"]
    h.close()
    if not ids:
        sys.exit(f"[ERROR] No RefSeq mRNA found for gene '{gene_symbol}'. Try --accession.")
    h = Entrez.efetch(db="nucleotide", id=",".join(ids), rettype="acc", retmode="text")
    accs = [line.strip() for line in h if line.strip()]
    h.close()
    nm = [a for a in accs if a.startswith("NM_")]
    return nm[0] if nm else accs[0]

def fetch_transcript_and_exons(accession, email):
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    try:
        record = SeqIO.read(handle, "gb")
    except Exception as e:
        sys.exit(f"[ERROR] Failed to fetch/parse GenBank for {accession}: {e}")
    finally:
        handle.close()
    seq = str(record.seq)
    exons = []
    exon_count = 1
    for feature in record.features:
        if feature.type == "exon":
            start = int(feature.location.start)
            end = int(feature.location.end)
            exon_number = feature.qualifiers.get("number", [str(exon_count)])[0]
            exons.append({"exon": exon_number, "start": start, "end": end})
            exon_count += 1
    if not exons:
        sys.exit(f"[ERROR] No exon features found in {accession}. Try a different transcript.")
    return seq, exons

def parse_coords(s):
    s = s.strip().replace(",", "")
    m = re.match(r"^chr?([0-9XYMT]+):(\d+)-(\d+)(?::([+\-12]))?$", s, re.IGNORECASE)
    if not m:
        sys.exit(f"[ERROR] Could not parse --pseudo-coords: '{s}'. Expected 'chr2:123-456[:strand]'.")
    chrom = m.group(1).upper()
    start = int(m.group(2)); end = int(m.group(3))
    if start > end: start, end = end, start
    strand_token = (m.group(4) or '+')
    strand = 1 if strand_token in ['+', '1', None] else 2  # NCBI: 1=forward, 2=reverse
    return chrom, start, end, strand

def fetch_pseudo_from_coords(coords, build, email):
    # normalize build name: accept GRCh38 / grch38 / hg38
    build_norm = (build or "").replace(".", "").upper()
    if build_norm not in {"GRCH38", "HG38"}:
        sys.exit("[ERROR] Only GRCh38/hg38 supported for --pseudo-coords right now.")

    chrom, start, end, strand = parse_coords(coords)  # handles chr17:.. or 17:.. and optional strand
    acc = CHR_TO_REFSEQ_GRCh38.get(chrom)
    if not acc:
        sys.exit(f"[ERROR] No RefSeq accession for chromosome '{chrom}' in GRCh38 map.")

    Entrez.email = email
    handle = Entrez.efetch(
        db="nucleotide",
        id=acc,
        rettype="fasta",
        retmode="text",
        seq_start=start,
        seq_stop=end,
        strand=strand  # 1=forward, 2=reverse (NCBI auto reverse-complements)
    )
    try:
        rec = SeqIO.read(handle, "fasta")
    except Exception as e:
        sys.exit(f"[ERROR] Failed to fetch pseudogene sequence for {coords}: {e}")
    finally:
        handle.close()

    strand_char = "-" if strand == 2 else "+"
    return str(rec.seq), f"{acc}:{start}-{end}:{strand_char}"

def dedup_legend(ax, loc='upper right'):
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if by_label:
        ax.legend(by_label.values(), by_label.keys(), loc=loc, fontsize=9)

# ---------------- main ----------------
def main():
    args = parse_args()

    # optional SSL workaround
    if args.ssl_ignore:
        import ssl
        ssl._create_default_https_context = ssl._create_unverified_context

    # resolve transcript
    transcript_id = args.accession or resolve_accession_from_gene(args.gene, args.email)

    # Step 1
    print(f"Step 1: Fetching {transcript_id} record and exon annotations...")
    parent_seq, exons = fetch_transcript_and_exons(transcript_id, args.email)
    print(f"Fetched mRNA length: {len(parent_seq)}")
    print(f"Number of exons found: {len(exons)}")
    sleep(0.2)

    # Step 2: pseudogene input
    if args.pseudo_fasta:
        print("\nStep 2: Loading pseudogene sequence from FASTA...")
        try:
            pg_rec = SeqIO.read(args.pseudo_fasta, "fasta")
        except Exception as e:
            sys.exit(f"[ERROR] Could not read FASTA '{args.pseudo_fasta}': {e}")
        pseudo_seq = str(pg_rec.seq)
        pseudo_label = args.pseudo_fasta
    else:
        print("\nStep 2: Fetching pseudogene sequence from genomic coordinates...")
        pseudo_seq, pseudo_label = fetch_pseudo_from_coords(args.pseudo_coords, args.build, args.email)

    print(f"Pseudogene sequence length: {len(pseudo_seq)} ({pseudo_label})")
    if not parent_seq or not pseudo_seq:
        sys.exit("[ERROR] One of the sequences is empty.")

    # Step 3: alignment
    print("\nStep 3: Performing local pairwise alignment...")
    aligns = pairwise2.align.localms(parent_seq, pseudo_seq, 2, -1, -0.5, -0.1)
    if not aligns:
        sys.exit("[ERROR] No alignment found.")
    aln_parent, aln_pseudo, score, begin, end = aligns[0]
    print(f"Alignment score: {score}")
    print(f"Alignment region in parent cDNA: {begin} to {end} (length {end - begin})")

    # Step 4: overlapping exons
    print("\nStep 4: Identifying overlapping exons...")
    overlapping = [ex for ex in exons if not (ex["end"] <= begin or ex["start"] >= end)]
    print(f"Number of overlapping exons: {len(overlapping)}")
    for ex in overlapping:
        print(f" - Exon {ex['exon']}: {ex['start']}–{ex['end']}")

    # Step 5: overall identity (gapless)
    print("\nStep 5: Calculating overall percent identity...")
    matches = valid = 0
    for a, b in zip(aln_parent, aln_pseudo):
        if a != '-' and b != '-':
            valid += 1
            if a.upper() == b.upper():
                matches += 1
    overall_identity = 100 * matches / valid if valid else 0.0
    print(f"Overall percent identity (gapless positions only): {overall_identity:.2f}%")

    print("\nCalculating detailed per-exon identity and coverage...")
    per_exon_stats = []
    for ex in overlapping:
        ovl_start = max(begin, ex["start"])
        ovl_end   = min(end,   ex["end"])
        if ovl_start >= ovl_end:
            continue
        pa = aln_parent[ovl_start - begin : ovl_end - begin]
        pb = aln_pseudo[ovl_start - begin : ovl_end - begin]

        m = gp = gq = v = 0
        for c1, c2 in zip(pa, pb):
            if c1 == '-': gp += 1
            if c2 == '-': gq += 1
            if c1 != '-' and c2 != '-':
                v += 1
                if c1.upper() == c2.upper():
                    m += 1

        exon_len = ex["end"] - ex["start"]
        aligned_len = len(pa)  # includes gaps
        coverage = 100.0 * aligned_len / exon_len if exon_len else 0.0
        identity_per_exon = 100.0 * m / exon_len if exon_len else 0.0

        per_exon_stats.append({
            "exon": ex["exon"], "start": ovl_start, "end": ovl_end,
            "matches": m, "gaps_parent": gp, "gaps_pseudo": gq,
            "valid_positions": v, "exon_length": exon_len, "aligned_length": aligned_len,
            "coverage_percent": coverage, "identity_percent": identity_per_exon
        })

        print(f"Exon {ex['exon']}: {ovl_start}–{ovl_end} | len {exon_len} | "
              f"aligned {aligned_len} | matches {m} | gaps(parent) {gp} gaps(pseudo) {gq} | "
              f"valid {v} | cov {coverage:.2f}% | ID {identity_per_exon:.2f}%")

    # Step 6: global plot (unchanged)
    print("\nStep 6: Generating plot...")
    plt.figure(figsize=(14, 3))
    plt.hlines(1, 0, len(parent_seq), colors='lightgray', linewidth=6, label="mRNA")
    cmap = cm.get_cmap('tab20', len(exons))
    for i, ex in enumerate(exons):
        plt.hlines(1, ex["start"], ex["end"], colors=cmap(i), linewidth=10, label=f"Exon {ex['exon']}")
    plt.hlines(1.1, begin, end, colors='red', linewidth=10, label="Pseudogene Alignment")
    plt.yticks([])
    plt.xlabel("mRNA cDNA position (nt)")
    plt.title("Exon Structure with Pseudogene Homology Region")
    dedup_legend(plt.gca(), loc='upper right')
    plt.tight_layout()
    global_png = f"{args.prefix}_pseudogene_alignment.png"
    plt.savefig(global_png); plt.close()
    print(f"Plot saved to {global_png}")

    # Step 6B: gapped exon track (unchanged, with labels)
    print("\nStep 6B: Generating exon-level gapped alignment track...")
    fig, ax = plt.subplots(figsize=(14, 2 + 0.3 * max(1, len(overlapping))))
    y_pos = 1.0
    for ex in overlapping:
        st = next((e for e in per_exon_stats if e['exon'] == ex['exon']), None)
        if not st: continue

        pa = aln_parent[st['start'] - begin : st['end'] - begin]
        pb = aln_pseudo[st['start'] - begin : st['end'] - begin]

        ax.hlines(y_pos, st['start'], st['end'], colors='lightgray', linewidth=6,
                  label='Exon Region' if y_pos == 1.0 else "")

        had_tick = False
        for i, (b1, b2) in enumerate(zip(pa, pb)):
            if b1 != '-' and b2 != '-':
                pos = st['start'] + i
                ax.hlines(y_pos, pos, pos+1, colors='blue', linewidth=6,
                          label='Aligned Bases' if y_pos == 1.0 and not had_tick else "")
                had_tick = True

        pid = st.get('identity_percent', 0.0)
        cov = st.get('coverage_percent', 0.0)
        gap_only = (st.get('valid_positions', 0) == 0)
        label = f"Exon {ex['exon']} ({pid:.1f}% ID, {cov:.1f}% cov" + (", gap-only" if gap_only else "") + ")"
        ax.text(st['start'] - 50, y_pos, label, va='center', fontsize=9)

        y_pos -= 0.4

    ax.set_xlim(0, len(parent_seq))
    ax.set_yticks([])
    ax.set_xlabel("mRNA cDNA position (nt)")
    ax.set_title("Gapped Alignment per Exon between mRNA and Pseudogene")
    dedup_legend(ax, loc='upper right')
    gapped_png = f"{args.prefix}_pseudogene_gapped_exons.png"
    plt.tight_layout(); plt.savefig(gapped_png); plt.close()
    print(f"Gapped exon alignment plot saved to {gapped_png}")

    # Step 7: stats file (unchanged)
    print("\nStep 7: Writing stats to output file...")
    stats_txt = f"{args.prefix}_pseudogene_alignment_stats.txt"
    with open(stats_txt, "w") as f:
        f.write(f"Alignment Score: {score}\n")
        f.write(f"Alignment Region (cDNA coords): {begin}–{end}\n")
        f.write(f"Alignment Length: {end - begin}\n")
        f.write(f"Percent Identity (gapless positions only): {overall_identity:.2f}%\n\n")
        f.write("Overlapping Exons:\n")
        for ex in overlapping:
            f.write(f" - Exon {ex['exon']}: {ex['start']}–{ex['end']} (cDNA coords)\n")
        f.write("\nPer-Exon Identity and Coverage:\n")
        f.write(f"{'Exon':<6}{'Start':>8}{'End':>8}{'Matches':>10}{'Gaps_parent':>12}{'Gaps_Pseudo':>12}"
                f"{'Valid_Pos':>12}{'Exon_Len':>10}{'Aligned_Len':>12}{'Coverage(%)':>14}{'Identity(%)':>14}\n")
        for s in per_exon_stats:
            f.write(f"{s['exon']:<6}{s['start']:>8}{s['end']:>8}"
                    f"{s['matches']:>10}{s['gaps_parent']:>12}{s['gaps_pseudo']:>12}"
                    f"{s['valid_positions']:>12}{s['exon_length']:>10}{s['aligned_length']:>12}"
                    f"{s['coverage_percent']:>14.2f}{s['identity_percent']:>14.2f}\n")
    print(f"Alignment stats saved to {stats_txt}")
    print("\nDone!")

if __name__ == "__main__":
    main()
