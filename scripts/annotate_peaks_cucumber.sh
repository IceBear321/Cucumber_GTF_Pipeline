#!/usr/bin/env bash

GTF_FWD="/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.fwd.gtf"
GTF_REV="/data/czh/reference_genome/cucumber/ChineseLong_v3.final.strand_corrected.rev.gtf"

PEAK_FWD="exomePeak2_fwd/exomePeak2_fwd_whole_genome/peaks.csv"
PEAK_REV="exomePeak2_rev/exomePeak2_rev_whole_genome/peaks.csv"
OUTPUT="exomePeak2_annotated_peaks_cucumber_strand_corrected.tsv"

csv2bed() {
    local csv=$1
    local strand=$2
    awk -F',' -v s="$strand" 'NR>1 {
        gsub(/"/, "", $2)
        gsub(/"/, "", $3)
        gsub(/"/, "", $4)
        gsub(/"/, "", $13)
        gsub(/"/, "", $14)
        gsub(/"/, "", $15)
        print $2"\t"$3"\t"$4"\t"s"\t"$13"\t"$14"\t"$15
    }' "$csv"
}

annotate_single_strand() {
    local bed_file=$1
    local gtf_file=$2
    local out_strand=$3
    
    local tmp_exon=$(mktemp)
    local tmp_3utr=$(mktemp)
    local tmp_5utr=$(mktemp)
    local tmp_stop=$(mktemp)
    local tmp_start=$(mktemp)
    local tmp_gene=$(mktemp)
    local tmp_remain=$(mktemp)
    
    awk -F'\t' '$3 == "exon"' "$gtf_file" > "$tmp_exon"
    awk -F'\t' '$3 == "three_prime_utr"' "$gtf_file" > "$tmp_3utr"
    awk -F'\t' '$3 == "five_prime_utr"' "$gtf_file" > "$tmp_5utr"
    
    awk -F'\t' '$3 == "CDS"' "$gtf_file" | sort -k1,1 -k4,4n -k5,5n | awk -F'\t' '{
        gene=""; if ($0 ~ /gene_id "([^"]+)"/) {match($0, /gene_id "([^"]+)"/, m); gene=m[1]}
        key=$1"\t"$7"\t"gene
        if ($7 == "+") { if (!(key in cds_end) || $5 > cds_end[key]) cds_end[key] = $5 }
        else { if (!(key in cds_start) || $4 < cds_start[key]) cds_start[key] = $4 }
    } END{
        for (k in cds_end) { split(k,arr,"\t"); if (arr[2]=="+") print arr[1]"\t"(cds_end[k]-10)"\t"(cds_end[k]+10)"\t"arr[2]"\t"arr[3] }
        for (k in cds_start) { split(k,arr,"\t"); if (arr[2]=="-") print arr[1]"\t"(cds_start[k]-10)"\t"(cds_start[k]+10)"\t"arr[2]"\t"arr[3] }
    }' | sort -k1,1 -k2,2n > "$tmp_stop"
    
    awk -F'\t' '$3 == "CDS"' "$gtf_file" | sort -k1,1 -k4,4n -k5,5n | awk -F'\t' '{
        gene=""; if ($0 ~ /gene_id "([^"]+)"/) {match($0, /gene_id "([^"]+)"/, m); gene=m[1]}
        key=$1"\t"$7"\t"gene
        if ($7 == "+") { if (!(key in cds_start) || $4 < cds_start[key]) cds_start[key] = $4 }
        else { if (!(key in cds_end) || $5 > cds_end[key]) cds_end[key] = $5 }
    } END{
        for (k in cds_start) { split(k,arr,"\t"); if (arr[2]=="+") print arr[1]"\t"(cds_start[k]-10)"\t"(cds_start[k]+10)"\t"arr[2]"\t"arr[3] }
        for (k in cds_end) { split(k,arr,"\t"); if (arr[2]=="-") print arr[1]"\t"(cds_end[k]-10)"\t"(cds_end[k]+10)"\t"arr[2]"\t"arr[3] }
    }' | sort -k1,1 -k2,2n > "$tmp_start"
    
    awk -F'\t' '$3 == "mRNA"' "$gtf_file" | awk -F'\t' '{
        gene=""; if ($0 ~ /gene_id "([^"]+)"/) {match($0, /gene_id "([^"]+)"/, m); gene=m[1]}
        print $1"\t"$4"\t"$5"\t"$7"\t"gene
    }' | sort -k1,1 -k2,2n > "$tmp_gene"
    
    cp "$bed_file" "$tmp_remain"
    
    local -A genes features
    local -A priority
    local -A seen_peak
    
    priority["three_prime_utr"]=1
    priority["stop_codon"]=2
    priority["exon"]=3
    priority["start_codon"]=4
    priority["five_prime_utr"]=5
    priority["intron"]=6
    priority["intergenic"]=7
    
    for feat in three_prime_utr stop_codon exon start_codon five_prime_utr; do
        case $feat in
            three_prime_utr) local bed="$tmp_3utr" ;;
            stop_codon) local bed="$tmp_stop" ;;
            exon) local bed="$tmp_exon" ;;
            start_codon) local bed="$tmp_start" ;;
            five_prime_utr) local bed="$tmp_5utr" ;;
        esac
        
        local tmp_overlap=$(mktemp)
        bedtools intersect -a "$tmp_remain" -b "$bed" -wa -wb > "$tmp_overlap"
        
        while IFS=$'\t' read -r chr start end strand log2fc pval fdr rest; do
            key="${chr}:${start}:${end}:${strand}:${log2fc}:${pval}:${fdr}"
            gene_id=$(echo "$rest" | grep -oP 'gene_id "[^"]+' | head -1 | sed 's/gene_id "//')
            if [[ -z "$gene_id" ]]; then continue; fi
            
            if [[ -z "${seen_peak[$key]}" ]]; then
                seen_peak[$key]=1
                genes["$key"]=$gene_id
                features["$key"]=$feat
            fi
        done < <(awk -F'\t' -v feat="$feat" 'NF >= 16 {
            chr=$1; start=$2; end=$3; strand=$4; log2fc=$5; pval=$6; fdr=$7
            attr=""; for(i=16;i<=NF;i++) attr=attr (attr?" ":"") $i
            print chr,start,end,strand,log2fc,pval,fdr,attr
        }' "$tmp_overlap")
        
        local tmp_new=$(mktemp)
        bedtools intersect -a "$tmp_remain" -b "$bed" -v > "$tmp_new"
        mv "$tmp_new" "$tmp_remain"
        rm "$tmp_overlap"
    done
    
    local tmp_intron_overlap=$(mktemp)
    bedtools intersect -a "$tmp_remain" -b "$tmp_gene" -wa -wb > "$tmp_intron_overlap"
    
    while IFS=$'\t' read -r chr start end strand log2fc pval fdr gene_id; do
        key="${chr}:${start}:${end}:${strand}:${log2fc}:${pval}:${fdr}"
        if [[ -z "${seen_peak[$key]}" ]]; then
            seen_peak[$key]=1
            genes["$key"]=$gene_id
            features["$key"]="intron"
        fi
    done < <(awk -F'\t' 'NF >= 5 {print $1,$2,$3,$4,$5,$6,$7,$10}' "$tmp_intron_overlap")
    
    local tmp_new=$(mktemp)
    bedtools intersect -a "$tmp_remain" -b "$tmp_gene" -v | awk -v s="$out_strand" 'BEGIN{FS="\t"; OFS="\t"}
    {print $1,$2,$3,s,"intergenic","intergenic",$5,$6,$7}' > "$tmp_new"
    
    for key in "${!genes[@]}"; do
        IFS=':' read -r chr start end strand log2fc pval fdr <<< "$key"
        echo -e "${chr}\t${start}\t${end}\t${strand}\t${genes[$key]}\t${features[$key]}\t${log2fc}\t${pval}\t${fdr}"
    done | sort -k1,1 -k2,2n -k3,3n -k4,4
    
    cat "$tmp_new"
    
    rm "$tmp_exon" "$tmp_3utr" "$tmp_5utr" "$tmp_stop" "$tmp_start" "$tmp_gene" "$tmp_remain" "$tmp_intron_overlap" "$tmp_new"
}

echo -e "chr\tpeak_start\tpeak_end\tstrand\tgeneid\tfeature\tlog2FC\tpvalue\tfdr" > "$OUTPUT"

tmp_fwd=$(mktemp)
tmp_rev=$(mktemp)

csv2bed "$PEAK_FWD" "+" > "$tmp_fwd"
csv2bed "$PEAK_REV" "-" > "$tmp_rev"

annotate_single_strand "$tmp_fwd" "$GTF_FWD" "+" > "${tmp_fwd}.annotated"
annotate_single_strand "$tmp_rev" "$GTF_REV" "-" > "${tmp_rev}.annotated"

cat "${tmp_fwd}.annotated" "${tmp_rev}.annotated" \
| sort -k1,1 -k2,2n -k3,3n -k4,4 \
| awk 'BEGIN{FS="\t"; OFS="\t"}
{
    key = $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9
    if (key == prev) next
    prev = key
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9
}' \
>> "$OUTPUT"

rm "$tmp_fwd" "$tmp_rev" "${tmp_fwd}.annotated" "${tmp_rev}.annotated"

echo "Output saved to $OUTPUT"
