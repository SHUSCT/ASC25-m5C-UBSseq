#!/bin/bash
detected_file=passed_SRR23538290.tsv
true_file=GSE225614_HeLa-WT_sites.tsv
Precision=$(awk 'NR==FNR {a[$1,$2,$3]=1; next} ($1,$2,$3) in a' "$true_file" "$detected_file" | wc -l | awk -v total=$(wc -l < "$detected_file") '{printf "%.2f", ($1/total)*100}')
echo "Precision: $Precision"