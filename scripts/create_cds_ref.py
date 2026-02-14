#!/usr/bin/env python3
"""
创建CDS-only参考GTF
从ChineseLong_v3.gtf中提取gene和CDS，去除mRNA/exon信息
避免StringTie过度依赖旧注释
"""
import re

def create_cds_ref(input_gtf, output_gtf):
    genes = {}
    cds_regions = {}
    
    with open(input_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            source = fields[1]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            frame = fields[7]
            attributes = fields[8]
            
            gene_id = None
            transcript_id = None
            
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    m = re.search(r'"([^"]+)"', attr)
                    if m:
                        gene_id = m.group(1)
                elif attr.startswith('transcript_id'):
                    m = re.search(r'"([^"]+)"', attr)
                    if m:
                        transcript_id = m.group(1)
            
            if feature == 'gene' and gene_id:
                genes[gene_id] = {'chrom': chrom, 'source': source,
                                  'start': start, 'end': end, 'strand': strand}
            elif feature == 'CDS' and transcript_id:
                if transcript_id not in cds_regions:
                    cds_regions[transcript_id] = []
                cds_regions[transcript_id].append({
                    'chrom': chrom, 'source': source, 'start': start,
                    'end': end, 'strand': strand, 'frame': frame, 'gene_id': gene_id
                })
    
    with open(output_gtf, 'w') as f:
        for gene_id, gene_info in genes.items():
            attr = f'gene_id "{gene_id}";'
            f.write(f"{gene_info['chrom']}\t{gene_info['source']}\tgene\t"
                   f"{gene_info['start']}\t{gene_info['end']}\t.\t"
                   f"{gene_info['strand']}\t.\t{attr}\n")
        
        for transcript_id, cds_list in cds_regions.items():
            gene_id = cds_list[0]['gene_id']
            for cds in cds_list:
                attr = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
                f.write(f"{cds['chrom']}\t{cds['source']}\tCDS\t"
                       f"{cds['start']}\t{cds['end']}\t.\t"
                       f"{cds['strand']}\t{cds['frame']}\t{attr}\n")
    
    print(f"基因数量: {len(genes)}, 转录本数量: {len(cds_regions)}")
    print(f"输出: {output_gtf}")

if __name__ == '__main__':
    create_cds_ref('/data/czh/reference_genome/cucumber/ChineseLong_v3.gtf',
                   '/data/czh/reference_genome/cucumber/ChineseLong_v3_CDS_only.gtf')
