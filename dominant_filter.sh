##!/usr/bin/env bash

while getopts "f:i:m:" opt; do
    case $opt in
        f) dom_callset_AF="$OPTARG";; # Callset AF for dominant genes
        i) dominant_vcf="$OPTARG";; # Input dominant filtered VCF file
        m) metadata="$OPTARG";; # Input metadata file
    esac
done

# Collect information about what column affected and unsolved status is in from metadata
affected_column=$(head -1 $metadata | tr '\t' '\n' | grep -n "affected_status" | cut -f1 -d':')
unsolved_column=$(head -1 $metadata | tr '\t' '\n' | grep -n "solve_state" | cut -f1 -d':')
bcftools view -h $dominant_vcf | tail -1 | cut -f10- | tr '\t' '\n' > samples.txt

# Extract samples which are affected or unsolved and determine their index in the VCF
cat $metadata | awk -F "\t" '{if ($'$affected_column' == "\"Unaffected\"")  {print $1}}' | tr -d '"' > unaffected.txt
cat $metadata | awk -F "\t" '{if ($'$affected_column' == "\"Affected\"" && ($'$unsolved_column' == "\"Unsolved\"" || $'$unsolved_column' == "\"Tier 2\""))  {print $1}}' | tr -d '"' > affected.txt
unaffected_lines=($(grep -wn -f unaffected.txt samples.txt | cut -f1 -d':')) # unaffected
affected_lines=($(grep -wn -f affected.txt samples.txt | cut -f1 -d':')) # affected

# Create filter query with input variables
query="(AF<=$dom_callset_AF && "

# Filter out all variants which are heterozygous in unaffected individuals
for i in ${unaffected_lines[@]}; do
    index=$(($i-1))
    query+="FORMAT/GT[$index]!='0/1' && "
done

# Run filtering on input vcf
final_query=$(echo "$query" | sed -e 's/\(.*\)&& /\1)/' -e 's/\(.*\) /\1/') # Format query
bcftools filter --threads $(getconf _NPROCESSORS_ONLN) -i "$final_query" $dominant_vcf
