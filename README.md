# Cucumber (ChineseLong v3) GTF Processing Pipeline

A complete pipeline for processing cucumber genome annotation using RNA-Seq data, including transcript assembly, quality filtering, strand direction correction, and peak annotation.

## Overview

This pipeline improves genome annotation by:
1. Creating CDS-only reference from existing annotation
2. Assembling transcripts using StringTie with RNA-Seq data
3. Filtering out low-quality transcripts (MSTRG, short <500bp)
4. Correcting strand orientation based on strand-specific sequencing
5. Annotating peaks to genomic features

For peak annotation, see the separate `annotate_peaks_cucumber.sh` script in the parent directory.

## Pipeline Steps

### Step 1: Create CDS-only Reference
```bash
python3 scripts/create_cds_ref.py
```

### Step 2: Transcript Assembly (requires StringTie)
```bash
# Forward strand
stringtie -G data/ChineseLong_v3_CDS_only.gtf \
    -o data/stringtie_fwd.gtf -l fwd -p 8 -f 0.01 -g 50 \
    /path/to/test_fwd.bam

# Reverse strand
stringtie -G data/ChineseLong_v3_CDS_only.gtf \
    -o data/stringtie_rev.gtf -l rev -p 8 -f 0.01 -g 50 \
    /path/to/test_rev.bam
```

### Step 3: Merge Assemblies
```bash
echo -e "data/stringtie_fwd.gtf\ndata/stringtie_rev.gtf" > data/merge_list.txt
stringtie --merge -G data/ChineseLong_v3_CDS_only.gtf \
    -o data/stringtie_merged.gtf data/merge_list.txt
```

### Step 4: Compare with Reference
```bash
gffcompare -r data/ChineseLong_v3.gtf -o data/gffcompare_out \
    data/stringtie_merged.gtf
```

### Step 5: Generate Final GTF
```bash
python3 scripts/create_final_v4.py
```

### Step 6: Filter and Correct Strand
```bash
python3 scripts/filter_and_correct_strand.py
```

## Required Data Files

Download from: http://cucurbitgenomics.org/v2/ftp/genome/cucumber/Chinese_long/v3/

| File | Description |
|------|-------------|
| `ChineseLong_v3.gtf` | Reference genome annotation |
| `ChineseLong_CDS_v3.fa.gz` | Reference CDS sequences |
| `test_fwd.bam` | Forward strand RNA-Seq |
| `test_rev.bam` | Reverse strand RNA-Seq |

## Output Files

| File | Description |
|------|-------------|
| `ChineseLong_v3.final.strand_corrected.gtf` | Final processed GTF |
| `ChineseLong_v3.final.gtf` | Intermediate GTF (before strand correction) |

## Statistics

### Filtering Results
- Original transcripts: 25,787
- Removed MSTRG: 1,470
- Removed <500bp: 3,358
- Remaining: 21,257

### Strand Correction Results
- + strand OK: 10,117 (95.5%)
- - strand OK: 10,222 (95.9%)
- Total flipped: 918 (4.3%)

## Requirements

- Python 3.8+
- StringTie
- GFFCompare
- featureCounts (subread)
- pandas
- bedtools (for peak annotation)

## Installation

```bash
# Install required tools
conda install -c bioconda stringtie gffcompare subread bedtools
```

## Citation

If you use this pipeline, please cite:
- ChineseLong v3 genome annotation
- StringTie paper
- RNA-Seq methodology

## License

MIT License

## Author

Zhihe Cai

## Contact

zh.cai@pku.edu.cn
