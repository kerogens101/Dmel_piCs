#!/bin/bash
#SBATCH --mem=10gb

# copy the piC annotation files with their frequency in column 4th

## frequency-grouped enrichment analysis of INDELs in piCs

files=("DSPR_allMethods_multiIntersect_common_remergedV7_collapsed.repComb.250bp.bed" "DSPR_proTRACv3_merged_RPM_filterV4_curated_remap_collapsed_count.bed" "DSPR_restrictive_filterV5_curated_remap_collapsed_count.bed")

for file in "${files[@]}"; do
  for i in {1..8}; do
    awk 'int($4) == '$i'' "$file" >"${file%.*}_n$i.bed"
  done
done

for del_file in DSPR_allMethods_*DEL.bed; do
  parallel --progress --eta -j 3 './poverlap/poverlap.py poverlap -i ISO1_dm6_heterochromatin_coordinates.bed \
  -g ISO1_chrom.sizes -N 1000 --ncpus=6 -A {} -B DSPR_allMethods_multiIntersect_common_remergedV7_collapsed.repComb.250bp_n8.bed \
  >het_shuffled_results/{/.}.DEL.txt' ::: "$del_file"
done

for ins_file in DSPR_allMethods_*INS.bed; do
  parallel --progress --eta -j 3 './poverlap/poverlap.py poverlap -i ISO1_dm6_heterochromatin_coordinates.bed \
  -g ISO1_chrom.sizes -N 1000 --ncpus=6 -A {} -B DSPR_allMethods_multiIntersect_common_remergedV7_collapsed.repComb.250bp_n8.bed \
  >het_shuffled_results/{/.}.INS.txt' ::: "$ins_file"
done

# Create a header for the CSV file
echo "File,Observed,Metric,Simulated Mean Metric,Simulated_p" > output.csv

# Loop through all the JSON files in the current directory
for file in *.json; do
    # Extract values from the JSON file
    observed=$(jq -r ".['wc -l'].observed" "$file")
    metric=$(jq -r ".['wc -l'].metric" "$file")
    simulated_mean=$(jq -r ".['wc -l']['simulated mean metric']" "$file")
    simulated_p=$(jq -r ".['wc -l'].simulated_p" "$file")

    # Print the values in CSV format and append to the CSV file
    echo "$file,$observed,$metric,$simulated_mean,$simulated_p" >> output.csv
done

echo "CSV file 'output.csv' has been created with the concatenated values."

