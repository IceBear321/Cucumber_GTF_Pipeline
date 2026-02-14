#!/usr/bin/env python3
"""
GTF过滤和链方向校正脚本
功能：
1. 剔除MSTRG开头的转录本
2. 剔除长度<500bp的转录本
3. 使用featureCounts分析链方向
4. 翻转错误链方向的转录本
"""
import os
import re
import subprocess
import pandas as pd

class GTFFilterAndCorrector:
    def __init__(self, input_gtf, bam_fwd, bam_rev, output_gtf, ratio_threshold=10):
        self.input_gtf = input_gtf
        self.bam_fwd = bam_fwd
        self.bam_rev = bam_rev
        self.output_gtf = output_gtf
        self.ratio_threshold = ratio_threshold
        self.featurecounts = '/home/czh/miniconda3/bin/featureCounts'
        
    def get_transcript_info(self, gtf_file):
        """从GTF提取转录本信息"""
        transcripts = {}
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                if fields[2] != 'mRNA':
                    continue
                
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                tid = None
                gid = None
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        tid = attr.split('"')[1]
                    elif attr.startswith('gene_id'):
                        gid = attr.split('"')[1]
                
                if tid:
                    length = end - start + 1
                    transcripts[tid] = {
                        'gene_id': gid,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'length': length
                    }
        return transcripts
    
    def filter_mstrg_and_short(self, transcripts):
        """剔除MSTRG开头和长度<500bp的转录本"""
        keep_transcripts = {}
        remove_mstrg = 0
        remove_short = 0
        
        for tid, info in transcripts.items():
            if tid.startswith('MSTRG'):
                remove_mstrg += 1
                continue
            if info['length'] < 500:
                remove_short += 1
                continue
            keep_transcripts[tid] = info
        
        print(f"MSTRG转录本: 删除 {remove_mstrg}")
        print(f"短转录本(<500bp): 删除 {remove_short}")
        print(f"保留转录本: {len(keep_transcripts)}")
        
        return set(keep_transcripts.keys())
    
    def run_featurecounts(self, keep_transcripts):
        """运行featureCounts"""
        work_dir = os.path.dirname(self.output_gtf)
        prefix = os.path.splitext(os.path.basename(self.output_gtf))[0]
        
        fwd_count = os.path.join(work_dir, f'{prefix}_fwd.txt')
        rev_count = os.path.join(work_dir, f'{prefix}_rev.txt')
        
        cmd_fwd = [self.featurecounts, '-s', '0', '-p', '-a', self.input_gtf,
                   '-o', fwd_count, self.bam_fwd, '-T', '8']
        cmd_rev = [self.featurecounts, '-s', '0', '-p', '-a', self.input_gtf,
                   '-o', rev_count, self.bam_rev, '-T', '8']
        
        print("运行featureCounts (fwd)...")
        subprocess.run(cmd_fwd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        print("运行featureCounts (rev)...")
        subprocess.run(cmd_rev, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        return fwd_count, rev_count
    
    def analyze_strand(self, fwd_count, rev_count):
        """分析链方向，确定需要翻转的转录本"""
        fwd = pd.read_csv(fwd_count, sep='\t', comment='#')
        rev = pd.read_csv(rev_count, sep='\t', comment='#')
        
        fwd_col = fwd.columns[-1]
        rev_col = rev.columns[-1]
        
        def get_main_strand(s):
            s = str(s)
            if '+' in s:
                return '+'
            elif '-' in s:
                return '-'
            return '.'
        
        merged = pd.DataFrame({
            'Geneid': fwd['Geneid'],
            'Strand': fwd['Strand'].apply(get_main_strand),
            'fwd': pd.to_numeric(fwd[fwd_col]),
            'rev': pd.to_numeric(rev[rev_col])
        })
        
        merged['fwd_rev_ratio'] = (merged['fwd'] + 1) / (merged['rev'] + 1)
        merged['rev_fwd_ratio'] = (merged['rev'] + 1) / (merged['fwd'] + 1)
        
        plus = merged[merged['Strand'] == '+']
        minus = merged[merged['Strand'] == '-']
        
        # +链需要翻转: rev/fwd > threshold
        # -链需要翻转: fwd/rev > threshold
        plus_need_flip = set(plus[plus['rev_fwd_ratio'] > self.ratio_threshold]['Geneid'])
        minus_need_flip = set(minus[minus['fwd_rev_ratio'] > self.ratio_threshold]['Geneid'])
        
        flip_dict = {}
        for tid in plus_need_flip:
            flip_dict[tid + '.1'] = '-'
        for tid in minus_need_flip:
            flip_dict[tid + '.1'] = '+'
        
        print(f"+链需要翻转: {len(plus_need_flip)}")
        print(f"-链需要翻转: {len(minus_need_flip)}")
        print(f"总共需要翻转: {len(flip_dict)}")
        
        return flip_dict
    
    def apply_filter_and_correction(self, keep_transcripts, flip_dict):
        """应用过滤和链方向校正"""
        count_mrna = 0
        count_exon = 0
        count_cds = 0
        count_utr5 = 0
        count_utr3 = 0
        count_other = 0
        flip_count = 0
        
        with open(self.input_gtf, 'r') as f, open(self.output_gtf, 'w') as out:
            for line in f:
                if line.startswith('#'):
                    out.write(line)
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    out.write(line)
                    continue
                
                feature = fields[2]
                strand = fields[6]
                attributes = fields[8]
                
                tid = None
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        tid = attr.split('"')[1]
                        break
                
                # 检查是否在保留列表中
                tid_base = tid.replace('.1', '') if tid else None
                if tid_base not in keep_transcripts:
                    continue
                
                # 检查是否需要翻转链方向
                if tid in flip_dict:
                    fields[6] = flip_dict[tid]
                    flip_count += 1
                
                out.write('\t'.join(fields) + '\n')
                
                # 统计
                if feature == 'mRNA':
                    count_mrna += 1
                elif feature == 'exon':
                    count_exon += 1
                elif feature == 'CDS':
                    count_cds += 1
                elif feature == 'five_prime_utr':
                    count_utr5 += 1
                elif feature == 'three_prime_utr':
                    count_utr3 += 1
                else:
                    count_other += 1
        
        print(f"\n输出统计:")
        print(f"  mRNA: {count_mrna}")
        print(f"  exon: {count_exon}")
        print(f"  CDS: {count_cds}")
        print(f"  5'UTR: {count_utr5}")
        print(f"  3'UTR: {count_utr3}")
        print(f"  其他: {count_other}")
        print(f"  翻转链方向: {flip_count}")
        print(f"\n输出文件: {self.output_gtf}")
    
    def run(self):
        """运行完整流程"""
        print("=" * 60)
        print("Step 1: 提取转录本信息")
        print("=" * 60)
        transcripts = self.get_transcript_info(self.input_gtf)
        print(f"总转录本: {len(transcripts)}")
        
        print("\n" + "=" * 60)
        print("Step 2: 过滤MSTRG和短转录本")
        print("=" * 60)
        keep_transcripts = self.filter_mstrg_and_short(transcripts)
        
        print("\n" + "=" * 60)
        print("Step 3: 运行featureCounts分析链方向")
        print("=" * 60)
        fwd_count, rev_count = self.run_featurecounts(keep_transcripts)
        
        print("\n" + "=" * 60)
        print("Step 4: 分析链方向")
        print("=" * 60)
        flip_dict = self.analyze_strand(fwd_count, rev_count)
        
        print("\n" + "=" * 60)
        print("Step 5: 应用过滤和链方向校正")
        print("=" * 60)
        self.apply_filter_and_correction(keep_transcripts, flip_dict)
        
        print("\n完成!")

if __name__ == '__main__':
    input_gtf = '/data/czh/reference_genome/cucumber/ChineseLong_v3.final.gtf'
    bam_fwd = '/data2/czh/TEL/cucumber/MeRIP_Seq_1/ssbamfiles/merged_BAM/cs_9300_Input_fwd.bam'
    bam_rev = '/data2/czh/TEL/cucumber/MeRIP_Seq_1/ssbamfiles/merged_BAM/cs_9300_Input_rev.bam'
    output_gtf = '/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.gtf'
    
    corrector = GTFFilterAndCorrector(input_gtf, bam_fwd, bam_rev, output_gtf, ratio_threshold=10)
    corrector.run()
