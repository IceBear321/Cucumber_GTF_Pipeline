# Peak Annotation Script

## 使用方法

```bash
python3 /data2/czh/TEL/cucumber/MeRIP_Seq_1/annotate_peaks_cucumber.py
```

## 输入文件

脚本会自动使用以下文件：
- Peak文件 (正向): `exomePeak2_fwd/exomePeak2_fwd_whole_genome/peaks.csv`
- Peak文件 (反向): `exomePeak2_rev/exomePeak2_rev_whole_genome/peaks.csv`
- GTF文件 (正链): `/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.fwd.gtf`
- GTF文件 (负链): `/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.rev.gtf`

## 输出文件

- `exomePeak2_annotated_peaks_cucumber_strand_corrected.tsv`

## 输出格式

| 列名 | 说明 |
|------|------|
| chr | 染色体 |
| peak_start | peak起始位置 |
| peak_end | peak结束位置 |
| strand | 链方向 |
| geneid | 注释到的基因ID |
| feature | 区域类型 |
| log2FC | log2 fold change |
| pvalue | p值 |
| fdr | FDR校正后的p值 |

## Feature优先级

按优先级从高到低：
1. three_prime_utr
2. stop_codon (CDS ±10bp)
3. exon
4. start_codon (CDS ±10bp)
5. five_prime_utr
6. intron
7. intergenic

## 自定义修改

如需修改输入文件，编辑脚本中的变量：
- `PEAK_FWD`: 正向peak CSV文件路径
- `PEAK_REV`: 反向peak CSV文件路径
- `GTF_FWD`: 正链GTF文件路径
- `GTF_REV`: 负链GTF文件路径
- `OUTPUT`: 输出文件路径
