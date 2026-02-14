#!/usr/bin/env python3
"""
从final_annotation_v2.gtf开始，添加CDS和UTR
使用exons边界和参考CDS来计算UTR
"""
import re
from collections import defaultdict

def parse_ref_gtf(gtf_file):
    """解析参考GTF，提取转录本的CDS"""
    transcripts = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature = fields[2]
            start, end = int(fields[3]), int(fields[4])
            strand = fields[6]
            
            gene_id = None
            transcript_id = None
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    m = re.search(r'"([^"]+)"', attr)
                    if m:
                        gene_id = m.group(1)
                elif attr.startswith('transcript_id'):
                    m = re.search(r'"([^"]+)"', attr)
                    if m:
                        transcript_id = m.group(1)
            
            if feature == 'mRNA' and transcript_id:
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = {
                        'gene_id': gene_id,
                        'strand': strand,
                        'cds': []
                    }
            
            if feature == 'CDS' and transcript_id in transcripts:
                transcripts[transcript_id]['cds'].append((start, end))
    
    return transcripts

def parse_asm_gtf(gtf_file):
    """解析组装GTF"""
    transcripts = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature = fields[2]
            start, end = int(fields[3]), int(fields[4])
            strand = fields[6]
            
            transcript_id = None
            gene_id = None
            for attr in fields[8].split(';'):
                attr = attr.strip()
                if attr.startswith('transcript_id'):
                    m = re.search(r'"([^"]+)"', attr)
                    if m:
                        transcript_id = m.group(1)
                elif attr.startswith('gene_id'):
                    m = re.search(r'"([^"]+)"', attr)
                    if m:
                        gene_id = m.group(1)
            
            if transcript_id is None:
                continue
            
            if feature == 'mRNA':
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = {
                        'chrom': fields[0],
                        'source': fields[1],
                        'strand': strand,
                        'gene_id': gene_id,
                        'exons': []
                    }
            
            elif feature == 'exon' and transcript_id in transcripts:
                transcripts[transcript_id]['exons'].append((start, end))
    
    # 计算每个转录本的start和end（基于exons）
    for tr_id, tr in transcripts.items():
        if tr['exons']:
            tr['start'] = min(e[0] for e in tr['exons'])
            tr['end'] = max(e[1] for e in tr['exons'])
    
    return transcripts

def get_utr_from_cds(tr_start, tr_end, cds_list, strand):
    """根据CDS位置计算UTR"""
    if not cds_list:
        return None, None
    
    cds_sorted = sorted(cds_list)
    
    if strand == '+':
        first_cds = cds_sorted[0][0]
        last_cds = cds_sorted[-1][1]
        
        utr5 = (tr_start, first_cds - 1) if first_cds > tr_start else None
        utr3 = (last_cds + 1, tr_end) if last_cds < tr_end else None
    else:
        first_cds = cds_sorted[0][0]
        last_cds = cds_sorted[-1][1]
        
        utr3 = (tr_start, first_cds - 1) if first_cds > tr_start else None
        utr5 = (last_cds + 1, tr_end) if last_cds < tr_end else None
    
    return utr5, utr3

def write_final_gtf(ref_data, asm_data, output_file):
    """写入最终GTF"""
    utr5_count, utr3_count, cds_count = 0, 0, 0
    
    with open(output_file, 'w') as f:
        for tr_id in sorted(asm_data.keys()):
            tr = asm_data[tr_id]
            gene_id = tr['gene_id']
            strand = tr['strand']
            tr_start = tr.get('tr_start', tr['start'])
            tr_end = tr.get('tr_end', tr['end'])
            
            if not tr['exons']:
                continue
            
            # mRNA - 使用exons边界
            exons = sorted(tr['exons'])
            tr_start = min(e[0] for e in exons)
            tr_end = max(e[1] for e in exons)
            
            f.write(f"{tr['chrom']}\t{tr['source']}\tmRNA\t{tr_start}\t{tr_end}\t.\t{strand}\t.\ttranscript_id \"{tr_id}\"; gene_id \"{gene_id}\";\n")
            
            # exons
            for i, (start, end) in enumerate(exons, 1):
                f.write(f"{tr['chrom']}\t{tr['source']}\texon\t{start}\t{end}\t.\t{strand}\t.\ttranscript_id \"{tr_id}\"; gene_id \"{gene_id}\"; exon_number \"{i}\";\n")
            
            # CDS - 从参考获取
            cds_list = []
            if tr_id in ref_data and ref_data[tr_id]['cds']:
                cds_list = ref_data[tr_id]['cds']
            elif gene_id in ref_data:
                for tid, rd in ref_data.items():
                    if rd['gene_id'] == gene_id and rd['cds']:
                        cds_list = rd['cds']
                        break
            
            if cds_list:
                cds_sorted = sorted(cds_list)
                for i, (start, end) in enumerate(cds_sorted, 1):
                    f.write(f"{tr['chrom']}\t{tr['source']}\tCDS\t{start}\t{end}\t.\t{strand}\t0\ttranscript_id \"{tr_id}\"; gene_id \"{gene_id}\"; exon_number \"{i}\";\n")
                    cds_count += 1
                
                # UTR - 基于exons边界
                utr5, utr3 = get_utr_from_cds(tr_start, tr_end, cds_list, strand)
                
                if utr5:
                    s, e = utr5
                    f.write(f"{tr['chrom']}\t{tr['source']}\tfive_prime_utr\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id \"{tr_id}\"; gene_id \"{gene_id}\";\n")
                    utr5_count += 1
                
                if utr3:
                    s, e = utr3
                    f.write(f"{tr['chrom']}\t{tr['source']}\tthree_prime_utr\t{s}\t{e}\t.\t{strand}\t.\ttranscript_id \"{tr_id}\"; gene_id \"{gene_id}\";\n")
                    utr3_count += 1
    
    return utr5_count, utr3_count, cds_count

if __name__ == '__main__':
    ref_gtf = '/data/czh/reference_genome/cucumber/ChineseLong_v3.gtf'
    asm_gtf = '/data/czh/reference_genome/cucumber/final_annotation_v2.gtf'
    output_file = '/data/czh/reference_genome/cucumber/final_annotation_v3.gtf'
    
    print("Parsing reference GTF...")
    ref_data = parse_ref_gtf(ref_gtf)
    print(f"Found {len(ref_data)} reference transcripts")
    
    print("Parsing assembled GTF...")
    asm_data = parse_asm_gtf(asm_gtf)
    print(f"Found {len(asm_data)} assembled transcripts")
    
    print("Writing final GTF with CDS and UTR...")
    utr5, utr3, cds = write_final_gtf(ref_data, asm_data, output_file)
    print(f"Done! Output: {output_file}")
    print(f"CDS: {cds}, 5'UTR: {utr5}, 3'UTR: {utr3}")
