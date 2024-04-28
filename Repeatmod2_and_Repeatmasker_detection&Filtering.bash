#!/bin/bash

RM_PATH="/workdir/repeatModeler2_and_dependencies/RM2"
CDHIT_PATH="/programs/cd-hit-4.8.1"

GENOMES=("A2" "A4" "A1" "A7" "B3" "B6" "ISO1" "OreR")

# Build RepeatModeler databases for multiple genomes
for GENOME in "${GENOMES[@]}"; do
    DB_NAME="${GENOME}_repmod2"
    GENOME_FA="${GENOME}_genome.fa"
    nohup "${RM_PATH}/BuildDatabase" -name "${DB_NAME}" "${GENOME_FA}" &
    nohup "${RM_PATH}/repeatModeler" -database "${DB_NAME}" -pa 6 -LTRStruct -LTRMaxSeqLen 10000 &
done

# Wait for background processes to finish
wait

# Run RepeatModeler on the built databases with specified parameters
for GENOME in "${GENOMES[@]}"; do
    DB_NAME="${GENOME}_repmod2"
    nohup "${RM_PATH}/repeatModeler" -database "${DB_NAME}" -pa 6 -LTRStruct -LTRMaxSeqLen 10000 &
done

# Wait for background processes to finish
wait

# Index the consensus sequences and perform processing
for file in *.fa; do
    samtools faidx "$file" > "$file.fai" 
done
for file in *.fai; do
    cut -f1 "$file" > "$file.txt"
done
sed -i 's/#/_/g; s/\//_/g; s/Unknown/Unknown Unknown /g' *.txt

# Create a summary file for each database
for file in *.fa.fai.txt; do
    printf "familyID class family size\n" | cat "$file" > "${file/.fa.fai.txt/_ID.txt}"
done
rm *.fai.txt

# Edit the families.fa files from RepeatModeler databases
# Remove specific types of repetitive elements
sed -i "/\b\(tRNA\|rRNA\|Satellite\|Simple_repeat\)\b/d" *.txt

# Save files with only IDs for each database
for GENOME in "${GENOMES[@]}"; do
    seqkit grep --pattern-file "${GENOME}_repmod2-families.fa.fai.txt" "${GENOME}_repmod2-families.fa" > "${GENOME}_repmod2_TEfam.fa"
done

# Process the edited files
for file in *rep.txt; do
    grep -n "*" "$file" > "${file/rep.txt/rep.txt}"
done
sed -i 's/, >/ /g; s/://g; s/...//g' *.txt
for file in *rep.txt; do
    cut -f1 "$file" > "${file/rep.txt/rep.txt}"
done

# Use SeqKit to filter the sequences based on pattern files
for GENOME in "${GENOMES[@]}"; do
    seqkit grep --pattern-file "${GENOME}_repmod2_TEfamrep.txtrep.txt" "${GENOME}_repmod2_TEfam_edited.fa" > "${GENOME}_repmod2_TEfam_rep.fa"
done

# Perform blastn for each genome against the repbase_seqs2020.fa to identify the Unknown families in case they are horizontally transferred
# Output format includes various information such as query sequence id, subject sequence id, alignment length, etc.
nohup blastn -subject repbase_seqs2020.fa -query A1_repmod2_TEfam_subfam_refTE_blastn_RMselectV5.fa \
-outfmt "6 qseqid sseqid qlen slen pident length mismatch qstart qend sstart send evalue bitscore" >A1_repmod2_TEfam_RMselectV5_repbase_blastn.txt &
# ... (similar commands for other genomes)

# Step 1: Filter out TEs based on blastn results
seqkit grep -v -f A1_flyBaseref_blastn_reduced_TEfam_ID.txt A1_repmod2_TEfam_subfam_removed.fa >A1_TEfam_subfam_blastn_reduced_TEfam_ID.fa &
# Similar commands for other samples (A2, A4, A5, A7, B3, B4, B6, ISO1, OreR)

# Step 2: Further filter based on another set of blastn results
seqkit grep -f A1_rm2_blastn_selected_ID.txt A1_TEfam_subfam_blastn_reduced_TEfam_ID.fa >A1_TEfam_subfam_blastn_RMselectV3.fa
# Similar commands for other samples

# Step 3: Rename headers in output files
sed -i 's/>/>OreR_/g' OreR_TEfam_subfam_blastn_RMselectV3.fa &
# Similar commands for other samples

# Step 4: Run RepeatMasker on genomes using a TE library
nohup /programs/RepeatMasker_4-1-0/RepeatMasker -lib DSPR_TEfam_RMSelectV2_cdhit.fa -no_is -nolow -pa 8 A1_genome.fa &
# Similar commands for other samples (A2, A4, A5, A7, B3, B4, B6, OreR, ISO1)

# Step 5: Run RepeatMasker on genomes using a clustered TE library
nohup /programs/RepeatMasker_4-1-0/RepeatMasker -lib DSPR_TEfam_RMSelectV2_cdhit_clstr.fa -no_is -nolow -pa 4 A1_genome.fa &
# Similar commands for other samples

# Step 6: Run RepeatMasker on genomes using a combined TE library
nohup /programs/RepeatMasker_4-1-0/RepeatMasker -lib combined_TEfam_RMSelect_cdhit_RMoutV2_clstr.fa -no_is -nolow -pa 4 A1_genome.fa &
# Similar commands for other samples

# Step 7: Parse RepeatMasker output to collect information
nohup perl parseRM.pl -i TE_sub-family/A1_genome.fa.out -l 50,1 &
# Similar commands for other samples

# Step 8: Collect new TE families based on divergence, length, and copy number selection in RStudio (done separately)
seqkit grep -f A1_genome_RM_CN_TEfam.txt A1_TEfam_subfam_blastn_reduced_TEfam_ID.fa >A1_TEfam_subfam_blastn_reduced_RM-select.fa &
# Similar commands for other samples (A2, A4, A5, A7, B3, B4, B6, ISO1, OreR)

# Step 9: Cluster and filter TE libraries using cd-hit
nohup /programs/cd-hit-4.8.1/cd-hit -i A2_TEfam_RMselectV3_repbase_compiled.fa -o A2_TEfam_RMselectV3_repbase_compiled_cdhit \
 -c 0.8 -aL 0.8 -aS 0.2 -A 80 -sc 1 -T 0 -d 1000 -g 1 -M 0 &
# Similar commands for other samples (A4, A7, B3, B4, B6, ISO1, OreR)

# Wait for background processes to finish
wait