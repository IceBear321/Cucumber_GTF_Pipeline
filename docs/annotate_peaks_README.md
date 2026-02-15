# Peak Annotation Script

## Usage

```bash
python3 /data2/czh/TEL/cucumber/MeRIP_Seq_1/annotate_peaks_cucumber.py
```

## Input Files

The script automatically uses the following files:
- Peak file (forward): `exomePeak2_fwd/exomePeak2_fwd_whole_genome/peaks.csv`
- Peak file (reverse): `exomePeak2_rev/exomePeak2_rev_whole_genome/peaks.csv`
- GTF file (forward strand): `/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.fwd.gtf`
- GTF file (reverse strand): `/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.rev.gtf`

## Output File

- `exomePeak2_annotated_peaks_cucumber_strand_corrected.tsv`

## Output Format

| Column | Description |
|--------|-------------|
| chr | Chromosome |
| peak_start | Peak start position |
| peak_end | Peak end position |
| strand | Strand direction |
| geneid | Annotated gene ID |
| feature | Feature type |
| log2FC | log2 fold change |
| pvalue | p-value |
| fdr | FDR-corrected p-value |

## Feature Priority

From highest to lowest priority:
1. three_prime_utr
2. stop_codon (CDS +/-10bp)
3. exon
4. start_codon (CDS +/-10bp)
5. five_prime_utr
6. intron
7. intergenic

## Customization

To modify input files, edit the variables in the script:
- `PEAK_FWD`: Forward peak CSV file path
- `PEAK_REV`: Reverse peak CSV file path
- `GTF_FWD`: Forward strand GTF file path
- `GTF_REV`: Reverse strand GTF file path
- `OUTPUT`: Output file path
