# Pseudogene Alignment & Visualization

This script aligns a **parent gene transcript (RefSeq accession)** against a **pseudogene sequence** and generates:
- Exon-level alignment statistics (coverage, identity, gaps).
- Visualization of aligned regions (gapped alignment tracks).
- Tabular outputs for downstream inspection.

---

## Requirements

- Python ≥ 3.9  
- Biopython ≥ 1.81  
- Matplotlib ≥ 3.7  
- NCBI Entrez access (requires an email address)  

Install requirements:

```bash
pip install biopython matplotlib
```

---

## Usage

Run the script directly from the command line.

### **Option 1: Provide a pseudogene FASTA file**
```bash
python pseudogene_aligner.py \
  --accession NM_020436.5 \
  --pseudo-fasta pseudogene.fasta \
  --prefix SALL4 \
  --email your_email@domain.com
```

- `--accession` : RefSeq mRNA accession of the parent gene (e.g., NM_020436.5).  
- `--pseudo-fasta` : FASTA file containing the pseudogene sequence.  
- `--prefix` : Prefix for output file names (plots, text files).  
- `--email` : Your email address (required by NCBI Entrez).  

---

### **Option 2: Fetch pseudogene sequence by genomic coordinates**
```bash
python pseudogene_aligner.py \
  --accession NM_000285.4 \
  --pseudo-coords 17:55560906-55561577:+ \
  --prefix PEPD \
  --email your_email@domain.com
```

- `--pseudo-coords` : Genomic coordinates of pseudogene (format: `chrom:start-end:strand`).  
  - **Example:** `17:55560906-55561577:+`  
  - Strand can be `+` or `-`.  
- Currently only **GRCh38** is supported for coordinate fetching.  

---

## Outputs

1. **Alignment summary (text file)**  
   - Overall alignment score, region, identity.  
   - Per-exon coverage and identity table.  

2. **Visualizations (PNG files)**  
   - Exon-level gapped alignment tracks.  
   - Color-coded exon alignment maps.  

3. **Console log**  
   - Step-by-step progress with exon statistics.  

---

## Example

```bash
python pseudogene_aligner.py \
  --accession NM_001159.5 \
  --pseudo-fasta MY_PSEUDO.fasta \
  --prefix MYGENE \
  --email me@example.com
```

Output:
- `MYGENE_alignment.txt` – alignment summary  
- `MYGENE_exon_alignment.png` – visual alignment track  
- `MYGENE_gapped_exons.png` – gapped exon visualization  
