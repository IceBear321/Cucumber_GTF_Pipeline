#!/usr/bin/env python3

import os
import subprocess
import tempfile
import re

GTF_FWD = "/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.fwd.gtf"
GTF_REV = "/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.rev.gtf"

PEAK_FWD = "exomePeak2_fwd/exomePeak2_fwd_whole_genome/peaks.csv"
PEAK_REV = "exomePeak2_rev/exomePeak2_rev_whole_genome/peaks.csv"
OUTPUT = "exomePeak2_annotated_peaks_cucumber_strand_corrected.tsv"

def csv2bed(csv_file, strand):
    peaks = []
    with open(csv_file, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 15:
                continue
            chr = parts[1].replace('"', '')
            start = parts[2].replace('"', '')
            end = parts[3].replace('"', '')
            log2fc = parts[12].replace('"', '')
            pval = parts[13].replace('"', '')
            fdr = parts[14].replace('"', '')
            peaks.append((chr, start, end, strand, log2fc, pval, fdr))
    return peaks

def extract_gtf_features(gtf_file):
    gene_regions = {}
    cds_start = {}
    cds_end = {}
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chr, source, feature, start, end, score, strand, frame, attr = parts
            start, end = int(start), int(end)
            
            gene_id = None
            m = re.search(r'gene_id "([^"]+)"', attr)
            if m:
                gene_id = m.group(1)
            
            if gene_id is None:
                continue
            
            key = (chr, strand, gene_id)
            
            if feature == 'mRNA':
                if key not in gene_regions:
                    gene_regions[key] = (start, end)
            elif feature == 'CDS':
                if strand == '+':
                    if key not in cds_start or start < cds_start[key]:
                        cds_start[key] = start
                    if key not in cds_end or end > cds_end[key]:
                        cds_end[key] = end
                else:
                    if key not in cds_start or start < cds_start[key]:
                        cds_start[key] = start
                    if key not in cds_end or end > cds_end[key]:
                        cds_end[key] = end
    
    return gene_regions, cds_start, cds_end

def write_bed_file(items, filename):
    with open(filename, 'w') as f:
        for item in items:
            if len(item) == 7:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\t{item[5]}\t{item[6]}\n")
            else:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\n")

def run_bedtools(cmd):
    subprocess.run(cmd, shell=True, check=True)

def annotate_peaks(peaks, gtf_file, out_strand):
    gene_regions, cds_start, cds_end = extract_gtf_features(gtf_file)
    
    priority = {'three_prime_utr': 1, 'stop_codon': 2, 'exon': 3, 'start_codon': 4, 'five_prime_utr': 5, 'intron': 6, 'intergenic': 7}
    results = {}
    
    tmp_dir = tempfile.mkdtemp()
    peak_file = os.path.join(tmp_dir, 'peaks.bed')
    with open(peak_file, 'w') as f:
        for p in peaks:
            f.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[3]}\t{p[4]}\t{p[5]}\t{p[6]}\n")
    
    remaining = peak_file
    
    features_to_check = [
        ('three_prime_utr', 'three_prime_utr'),
        ('stop_codon', 'stop_codon'),
        ('exon', 'exon'),
        ('start_codon', 'start_codon'),
        ('five_prime_utr', 'five_prime_utr')
    ]
    
    for feat_name, gtf_feature in features_to_check:
        feat_file = os.path.join(tmp_dir, f'{feat_name}.bed')
        
        feat_items = []
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                chr, source, feature, start, end, score, strand, frame, attr = parts
                if feature != gtf_feature:
                    continue
                start, end = int(start), int(end)
                gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    # Use name field (5th column) for gene_id
                    feat_items.append((chr, start, end, strand, '.', '.', gene_id))
        
        if feat_name in ['stop_codon', 'start_codon']:
            processed_items = []
            for key in (cds_start if feat_name == 'start_codon' else cds_end):
                chr, strand, gene_id = key
                if feat_name == 'stop_codon':
                    if strand == '+' and key in cds_end:
                        pos = cds_end[key]
                        processed_items.append((chr, pos - 10, pos + 10, strand, '.', '.', gene_id))
                    elif strand == '-' and key in cds_start:
                        pos = cds_start[key]
                        processed_items.append((chr, pos - 10, pos + 10, strand, '.', '.', gene_id))
                else:
                    if strand == '+' and key in cds_start:
                        pos = cds_start[key]
                        processed_items.append((chr, pos - 10, pos + 10, strand, '.', '.', gene_id))
                    elif strand == '-' and key in cds_end:
                        pos = cds_end[key]
                        processed_items.append((chr, pos - 10, pos + 10, strand, '.', '.', gene_id))
            feat_items = processed_items
        
        if not feat_items:
            continue
        
        write_bed_file(feat_items, feat_file)
        
        overlap_file = os.path.join(tmp_dir, f'overlap_{feat_name}.bed')
        run_bedtools(f"bedtools intersect -a {remaining} -b {feat_file} -wa -wb > {overlap_file}")
        
        with open(overlap_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # Peak has 7 cols, bed file has 7 cols, so gene_id is at column 14 (index 13)
                if len(parts) < 14:
                    continue
                chr, start, end = parts[0], int(parts[1]), int(parts[2])
                strand = parts[3]
                log2fc, pval, fdr = parts[4], parts[5], parts[6]
                
                # gene_id is in column 14 for bed files with 7 columns
                gene_id = parts[13]
                
                if not gene_id or gene_id == '' or gene_id == '.':
                    continue
                
                key = (chr, start, end, strand, log2fc, pval, fdr)
                if key not in results:
                    results[key] = {'gene': gene_id, 'feature': feat_name, 'priority': priority[feat_name]}
                else:
                    if gene_id not in results[key]['gene']:
                        results[key]['gene'] += ',' + gene_id
                    if priority[feat_name] < results[key]['priority']:
                        results[key]['feature'] = feat_name
                        results[key]['priority'] = priority[feat_name]
        
        new_remaining = os.path.join(tmp_dir, f'remaining_{feat_name}.bed')
        run_bedtools(f"bedtools intersect -a {remaining} -b {feat_file} -v > {new_remaining}")
        os.remove(remaining)
        remaining = new_remaining
    
    gene_file = os.path.join(tmp_dir, 'gene.bed')
    gene_items = [(chr, start, end, strand, '.', '.', gene_id) for (chr, strand, gene_id), (start, end) in gene_regions.items()]
    write_bed_file(gene_items, gene_file)
    
    intron_overlap = os.path.join(tmp_dir, 'intron_overlap.bed')
    run_bedtools(f"bedtools intersect -a {remaining} -b {gene_file} -wa -wb > {intron_overlap}")
    
    with open(intron_overlap, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # Peak has 7 cols, gene bed has 7 cols, so gene_id is at column 14 (index 13)
            if len(parts) < 14:
                continue
            chr, start, end = parts[0], int(parts[1]), int(parts[2])
            strand = parts[3]
            log2fc, pval, fdr = parts[4], parts[5], parts[6]
            gene_id = parts[13]
            
            key = (chr, start, end, strand, log2fc, pval, fdr)
            if key not in results:
                results[key] = {'gene': gene_id, 'feature': 'intron', 'priority': priority['intron']}
            else:
                if gene_id not in results[key]['gene']:
                    results[key]['gene'] += ',' + gene_id
    
    intergenic_file = os.path.join(tmp_dir, 'intergenic.bed')
    run_bedtools(f"bedtools intersect -a {remaining} -b {gene_file} -v > {intergenic_file}")
    
    with open(intergenic_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            chr, start, end = parts[0], int(parts[1]), int(parts[2])
            strand = parts[3]
            log2fc, pval, fdr = parts[4], parts[5], parts[6]
            
            key = (chr, start, end, strand, log2fc, pval, fdr)
            results[key] = {'gene': 'intergenic', 'feature': 'intergenic', 'priority': priority['intergenic']}
    
    import shutil
    shutil.rmtree(tmp_dir)
    
    return results

print("Processing forward strand...")
peaks_fwd = csv2bed(PEAK_FWD, '+')
results_fwd = annotate_peaks(peaks_fwd, GTF_FWD, '+')
print(f"  Found {len(results_fwd)} peaks")

print("Processing reverse strand...")
peaks_rev = csv2bed(PEAK_REV, '-')
results_rev = annotate_peaks(peaks_rev, GTF_REV, '-')
print(f"  Found {len(results_rev)} peaks")

all_results = {**results_fwd, **results_rev}

with open(OUTPUT, 'w') as f:
    f.write("chr\tpeak_start\tpeak_end\tstrand\tgeneid\tfeature\tlog2FC\tpvalue\tfdr\n")
    
    sorted_keys = sorted(all_results.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))
    prev_key = None
    for key in sorted_keys:
        chr, start, end, strand, log2fc, pval, fdr = key
        current_key_str = f"{chr}\t{start}\t{end}\t{strand}\t{log2fc}\t{pval}\t{fdr}"
        
        if prev_key is not None:
            prev_str = f"{prev_key[0]}\t{prev_key[1]}\t{prev_key[2]}\t{prev_key[3]}\t{prev_key[4]}\t{prev_key[5]}\t{prev_key[6]}"
            if current_key_str == prev_str:
                continue
        
        prev_key = key
        r = all_results[key]
        f.write(f"{chr}\t{start}\t{end}\t{strand}\t{r['gene']}\t{r['feature']}\t{log2fc}\t{pval}\t{fdr}\n")

print(f"Output saved to {OUTPUT}")
