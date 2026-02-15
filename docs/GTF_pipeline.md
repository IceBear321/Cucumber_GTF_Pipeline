# ChineseLong v3 GTF Complete Processing Pipeline

## Background

The original GTF file (`ChineseLong_v3.gtf`) has mRNA 5'UTR and 3'UTR regions that are much shorter than those detected by RNA-Seq data, leading to inaccurate transcript information.

This pipeline uses MeRIP-Seq Input sequencing data to reassemble the genome annotation to obtain more accurate UTR regions, with quality control (filtering low-quality transcripts, correcting strand orientation).

## Input Data

| File | Description |
|------|-------------|
| `ChineseLong_v3.gtf` | Reference genome annotation (only provides CDS positions) |
| `ChineseLong_genome_v3.fa` | Cucumber genome sequence |
| `cs_9300_Input_fwd.bam` | Forward strand RNA-Seq data |
| `cs_9300_Input_rev.bam` | Reverse strand RNA-Seq data |
| `ChineseLong_CDS_v3.fa.gz` | Reference CDS sequences (for validation) |

## Pipeline Steps

### Step 1: Create CDS-only Reference

```bash
python3 create_cds_ref.py
```

- Extract gene and CDS from `ChineseLong_v3.gtf`, remove mRNA/exon information
- Prevent StringTie from over-relying on old annotation
- Output: `ChineseLong_v3_CDS_only.gtf`

### Step 2: StringTie Assembly (Forward and Reverse separately)

```bash
# Forward strand assembly
stringtie -G ChineseLong_v3_CDS_only.gtf \
    -o stringtie_fwd.gtf -l fwd -p 8 -f 0.01 -g 50 \
    cs_9300_Input_fwd.bam

# Reverse strand assembly
stringtie -G ChineseLong_v3_CDS_only.gtf \
    -o stringtie_rev.gtf -l rev -p 8 -f 0.01 -g 50 \
    cs_9300_Input_rev.bam
```

Parameter description:
- `-f 0.01`: Minimum FPKM threshold, allows discovery of low-expression transcripts
- `-g 50`: Maximum gap length allowed, for discovering new exons

### Step 3: Merge Assemblies

```bash
echo -e "stringtie_fwd.gtf\nstringtie_rev.gtf" > merge_list.txt
stringtie --merge -G ChineseLong_v3_CDS_only.gtf \
    -o stringtie_merged.gtf merge_list.txt
```

### Step 4: Compare with Reference

```bash
gffcompare -r ChineseLong -o gff_v3.gtfcompare_out \
    stringtie_merged.gtf
```

### Step 5: Generate Final Annotation (with complete UTR)

```bash
python3 create_final_v4.py
```

- Output: `ChineseLong_v3.final.gtf`
- Contains complete mRNA/exon/CDS/5'UTR/3'UTR information

### Step 6: Filter and Correct Strand

```bash
python3 filter_and_correct_strand.py
```

This step automatically:
1. Remove transcripts starting with MSTRG (newly discovered unknown genes)
2. Remove transcripts <500bp (possibly incomplete transcripts)
3. Analyze strand direction using `featureCounts -s 0`
4. Flip transcripts with incorrect strand orientation:
   - + strand genes: if rev/fwd > 10 → flip to - strand
   - - strand genes: if fwd/rev > 10 → flip to + strand

- Output: `ChineseLong_v3.final.strand_corrected.gtf`

## Results Summary

### GTF Feature Statistics (Final Output)

| Feature | Count |
|---------|-------|
| mRNA | ~21,257 |
| exon | ~280,899 |
| CDS | ~121,072 |
| five_prime_utr | ~14,758 |
| three_prime_utr | ~14,957 |

### Filtering Statistics

| Filter | Count |
|--------|-------|
| Original transcripts | 25,787 |
| Removed MSTRG | 1,470 |
| Removed <500bp | 3,358 |
| Total removed | 4,530 |
| Remaining transcripts | 21,257 |

### Strand Correction Statistics

| Type | Total | Pass | Need Flip |
|------|-------|------|-----------|
| + strand genes | 10,590 | 10,117 | 473 |
| - strand genes | 10,667 | 10,222 | 445 |
| Total | 21,257 | 20,339 (95.7%) | 918 (4.3%) |

### Transcript Naming Statistics

| Type | Count | Description |
|------|-------|-------------|
| CsaV3* | ~24,317 | Mapped to original gene ID |
| MSTRG* | 1,470 | Newly discovered genes / no match |

## Core Scripts

| File | Function |
|------|----------|
| `create_cds_ref.py` | Create CDS-only reference GTF |
| `create_final_v4.py` | Generate final annotation (with complete mRNA/exon/CDS/UTR) |
| `filter_and_correct_strand.py` | Filter low-quality transcripts + correct strand |
| `GTF_pipeline.md` | Complete pipeline documentation |

## Output Files

| File | Description |
|------|-------------|
| `ChineseLong_v3_CDS_only.gtf` | CDS-only reference (Step 1 output) |
| `ChineseLong_v3.final.gtf` | Reassembled GTF (Step 5 output) |
| `ChineseLong_v3.final.strand_corrected.gtf` | Final corrected GTF (Step 6 output) |

## Notes

1. **Reverse strand data verification**: Reverse strand data (`cs_9300_Input_rev.bam`) needs strand-specific alignment confirmation first
2. **class_code explanation in assembly results**:
   - `=`: Exact match to reference transcript
   - `s`: Strand matches but structure differs
   - `u`: Unknown gene (newly discovered)
   - `n`: Partial match
3. **Strand direction logic**: Use `featureCounts -s 0` to count fwd and rev BAM separately, determine if strand direction is correct based on strand-specific library construction mode

## UTR Length Comparison

| Metric | Original Annotation | New Annotation |
|--------|-------------------|----------------|
| 5'UTR average length | 502.4 bp | ~950 bp |
| 3'UTR average length | 680.1 bp | ~940 bp |
| 5'UTR count | 15,483 | 14,758 |
| 3'UTR count | 15,449 | 14,957 |
