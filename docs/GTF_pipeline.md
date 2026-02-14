# 黄瓜基因组注释文件 (ChineseLong_v3) 完整处理流程

## 背景

原始GTF文件 (`ChineseLong_v3.gtf`) 中的mRNA的5'UTR和3'UTR区域比RNA-Seq测序数据检测到的短很多，导致转录本信息不准确。

本流程使用MeRIP-Seq的Input测序数据对基因组注释进行重新组装，以获得更准确的UTR区域，并进行质量控制（剔除低质量转录本、校正链方向）。

## 使用的数据

| 文件 | 描述 |
|------|------|
| `ChineseLong_v3.gtf` | 参考基因组注释 (仅用于提供CDS位置) |
| `ChineseLong_genome_v3.fa` | 黄瓜基因组序列 |
| `cs_9300_Input_fwd.bam` | 正链RNA-Seq测序数据 |
| `cs_9300_Input_rev.bam` | 负链RNA-Seq测序数据 |
| `ChineseLong_CDS_v3.fa.gz` | 参考CDS序列 (用于验证) |

## 流程步骤

### Step 1: 创建CDS-only参考

```bash
python3 create_cds_ref.py
```

- 从`ChineseLong_v3.gtf`提取gene和CDS，去除mRNA/exon信息
- 避免StringTie过度依赖旧注释
- 输出: `ChineseLong_v3_CDS_only.gtf`

### Step 2: StringTie组装 (正链和负链分开)

```bash
# 正链组装
stringtie -G ChineseLong_v3_CDS_only.gtf \
    -o stringtie_fwd.gtf -l fwd -p 8 -f 0.01 -g 50 \
    cs_9300_Input_fwd.bam

# 负链组装
stringtie -G ChineseLong_v3_CDS_only.gtf \
    -o stringtie_rev.gtf -l rev -p 8 -f 0.01 -g 50 \
    cs_9300_Input_rev.bam
```

参数说明:
- `-f 0.01`: 最低FPKM阈值，允许发现低表达转录本
- `-g 50`: 允许的最短gap，用于发现新外显子

### Step 3: 合并组装结果

```bash
echo -e "stringtie_fwd.gtf\nstringtie_rev.gtf" > merge_list.txt
stringtie --merge -G ChineseLong_v3_CDS_only.gtf \
    -o stringtie_merged.gtf merge_list.txt
```

### Step 4: 与参考比较

```bash
gffcompare -r ChineseLong_v3.gtf -o gffcompare_out \
    stringtie_merged.gtf
```

### Step 5: 生成最终注释 (含完整UTR)

```bash
python3 create_final_v4.py
```

- 输出: `ChineseLong_v3.final.gtf`
- 包含完整的mRNA/exon/CDS/5'UTR/3'UTR信息

### Step 6: 过滤和链方向校正

```bash
python3 filter_and_correct_strand.py
```

此步骤自动完成:
1. 剔除MSTRG开头的转录本 (新发现的未知基因)
2. 剔除长度<500bp的转录本 (可能是不完整的转录本)
3. 使用featureCounts -s 0分析链方向
4. 翻转错误链方向的转录本:
   - +链基因：如果rev/fwd > 10 → 翻转成-链
   - -链基因：如果fwd/rev > 10 → 翻转成+链

- 输出: `ChineseLong_v3.final.strand_corrected.gtf`

## 结果统计

### GTF特征统计 (最终输出)

| 特征 | 数量 |
|------|------|
| mRNA | ~21,257 |
| exon | ~280,899 |
| CDS | ~121,072 |
| five_prime_utr | ~14,758 |
| three_prime_utr | ~14,957 |

### 过滤统计

| 过滤项 | 数量 |
|--------|------|
| 原始转录本 | 25,787 |
| 删除MSTRG | 1,470 |
| 删除<500bp | 3,358 |
| 总删除 | 4,530 |
| 剩余转录本 | 21,257 |

### 链方向校正统计

| 类型 | 总数 | 满足要求 | 需要翻转 |
|------|------|----------|----------|
| +链基因 | 10,590 | 10,117 | 473 |
| -链基因 | 10,667 | 10,222 | 445 |
| 总计 | 21,257 | 20,339 (95.7%) | 918 (4.3%) |

### 转录本命名统计

| 类型 | 数量 | 说明 |
|------|------|------|
| CsaV3* | ~24,317 | 映射到原始基因ID |
| MSTRG* | 1,470 | 新发现的基因/无匹配 |

## 核心脚本

| 文件 | 功能 |
|------|------|
| `create_cds_ref.py` | 创建CDS-only参考GTF |
| `create_final_v4.py` | 生成最终注释文件(含完整mRNA/exon/CDS/UTR) |
| `filter_and_correct_strand.py` | 过滤低质量转录本 + 校正链方向 |
| `GTF_pipeline.md` | 完整流程文档 |

## 输出文件

| 文件 | 描述 |
|------|------|
| `ChineseLong_v3_CDS_only.gtf` | CDS-only参考 (Step1输出) |
| `ChineseLong_v3.final.gtf` | 重组装后的GTF (Step5输出) |
| `ChineseLong_v3.final.strand_corrected.gtf` | 最终校正后的GTF (Step6输出) |

## 注意事项

1. **负链数据验证**: 负链数据 (`cs_9300_Input_rev.bam`) 需要先进行strand-specific的比对确认
2. **组装结果中的class_code说明**:
   - `=` : 完全匹配参考转录本
   - `s` : 链匹配但结构不同
   - `u` : 未知基因(新发现)
   - `n` : 部分匹配
3. **链方向判断逻辑**: 使用featureCounts -s 0对fwd和rev BAM分别计数，根据链特异性建库模式判断链方向是否正确

## UTR长度对比

| 指标 | 原始注释 | 新注释 |
|------|----------|--------|
| 5'UTR平均长度 | 502.4 bp | ~950 bp |
| 3'UTR平均长度 | 680.1 bp | ~940 bp |
| 5'UTR数量 | 15,483 | 14,758 |
| 3'UTR数量 | 15,449 | 14,957 |
