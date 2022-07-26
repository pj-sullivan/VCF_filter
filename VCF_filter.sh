##!/usr/bin/env bash

unsolved=0

while getopts "a:d:f:i:m:q:u" opt; do
    case $opt in
        a) min_altReads="$OPTARG";; # Minimum alternative reads for variant
        d) min_depth="$OPTARG";; # Minimum depth for variant
        f) allele_frequency="$OPTARG";; # Max Allele frequency for variant
        i) input_vcf="$OPTARG";; # Input VCF file
        m) metadata="$OPTARG";; # Input metadata file
        q) min_QUAL="$OPTARG";; # Minimum QUAL value for variants
        u) unsolved=1;; # Trigger unsolved filtering
    esac
done

if [ ! -z $metadata ]; then
    affected_column=$(head -1 $metadata | tr '\t' '\n' | grep -n "affected_status" | cut -f1 -d':')
    unsolved_column=$(head -1 $metadata | tr '\t' '\n' | grep -n "solve_state" | cut -f1 -d':')
    bcftools view -h $input_vcf | tail -1 | cut -f10- | tr '\t' '\n' > samples.txt

    cat $metadata | awk -F "\t" '{if ($'$affected_column' == "\"Unaffected\"")  {print $1}}' | tr -d '"' > unaffected.txt
    if (( $unsolved == 1 )); then
        cat $metadata | awk -F "\t" '{if ($'$affected_column' == "\"Affected\"" && ($'$unsolved_column' == "\"Unsolved\"" || $'$unsolved_column' == "\"Tier 2\""))  {print $1}}' | tr -d '"' > affected.txt
    else
        cat $metadata | awk -F "\t" '{if ($'$affected_column' == "\"Affected\"")  {print $1}}' | tr -d '"' > affected.txt
    fi

    unaffected_lines=($(grep -wn -f unaffected.txt samples.txt | cut -f1 -d':')) # unaffected
    affected_lines=($(grep -wn -f affected.txt samples.txt | cut -f1 -d':')) # affected

    query="AF<=$allele_frequency && (QUAL>=$min_QUAL || QUAL='.') && FILTER=='PASS' && MAX(FORMAT/DP[*])>=$min_depth && MAX(FORMAT/AD[*:1])>=$min_altReads && "
    for i in ${unaffected_lines[@]}; do
        index=$(($i-1))
        query+="FORMAT/GT[$index]!='1/1' && "
    done

    query+=" ("
    for i in ${affected_lines[@]}; do
        index=$(($i-1))
        query+="FORMAT/GT[$index]!='0/0' || "
    done

    final_query=$(echo "$query" | sed -e 's/\(.*\)|| /\1)/' -e 's/\(.*\) /\1/')
fi

bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i "$final_query" $input_vcf
