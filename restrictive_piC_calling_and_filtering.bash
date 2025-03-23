#!/bin/sh

## this script consists of following stages: 
## mapping, filter by cutoff, analyze and filter by piRNA characteristics, and merge with low mappability scores

## Index your genomes and map respective piRNA reads to their genomes
## make sure the first two or three words in each genome file's name and correspdonding small RNA library name is the same.  
ls *genome.fasta | parallel -j 8 'bowtie-build {} {.} --threads 4' 
mkdir piC_calling/
## first add echo before -- 'bowtie -- and do a trial run. See if the individual commands make sense and change the cut and file extensions to match your inputs. Then remove the echo and run the changed command. 
ls *piRNA.fastq.gz | cut -d '_' -f1-2 | sort -u | parallel -j 4 'bowtie -n 1 -l 12 -a -m 1 -y -S {}_ovary_piRNA.fastq -x {}_genome -S piC_calling/{.}_genome.uniq.sam' 

## create sorted and indexed bam file for coverage
cd piC_calling/
ls *uniq.sam | parallel -j 5 'samtools view -bS -@ 10 {} >{.}.bam'
ls *uniq.bam | parallel -j 5 'samtools sort -@ 10 {} >{.}.sorted.bam'
ls *uniq.sorted.bam | parallel -j 5 'samtools index -@ 10 {}'

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
export PATH=/programs/deeptools-3.5.1/bin:$PATH
export PYTHONPATH=/programs/deeptools-3.5.1/lib64/python3.9/site-packages:/programs/deeptools-3.5.1/lib/python3.9/site-packages

## calculate piRNA coverage from bamCoverage and remove windows with less 2 CPM values. 
ls *sorted.bam | parallel -j 6 'bamCoverage -bs 500 -b {} -of bedgraph -p 2 -o {.}_500bp.bedGraph --normalizeUsing CPM' 
ls *500bp.CPM.bedGraph | parallel -j 10 "awk -v OFS='\t' '{ if (\$4>2) print \$0}' {} > {.}_RPM.bed"
for file in *5CPM.bed; do
    sort -k1,1 -k2,2n "$file" | bedtools merge -d 10 -c 4 -o sum > "${file%.bed}_merge10kb.bed"
done

## extract the piRNA-enriched domains from the respective genome files and remap the reads from the same fastq files used in line 11. 
ls *merge10kb.bed | cut -d '_' -f1-2 | sort -u | parallel -j 4 'bedtools getfasta -f ../{}_genome.fa -bed {}_ovary_piRNA_genome.uniq.5CPM_merge10kb.bed > {.}_ovary_piRNA_genome.uniq.5CPM_merge10kb_domain.fa' 
mkdir ../remap_and_filter
cp *domain.fa ../remap_and_filter
cd ../remap_and_filter
sed -i 's/:/_/g' *.fa
ls ovary_piRNA_genome.uniq.5CPM_merge10kb_domain.fa | parallel -j 5 'bowtie2-build -f {} {.}' 
ls *.fa | cut -d '_' -f1-4 | sort -u | parallel -j 5 'bowtie2 -N 1 -L 14 -k 1 -x {}_genome_uniq_500bp_5CPM.merge -q ../{}.fastq.gz -S {.}_remap.sam -p 5 --no-unal' 
ls *.sam | parallel -j 5 'samtools view -bS -@ 10 {} >{.}.bam'

# split this bam file by each merged domain and obtain characteristics of piRNAs from each domain to filter by 1U >60% of the reads and minimum 20 unique piRNAs per kb. 
mkdir bamtools_split
for file in *bam ; do bamtools split -in $file -reference ; done &
mv *REF*.bam bamtools_split
cd bamtools_split

#### from here onwards unless you use exactly my filename convention, you might want to carefully look at the input and output filenames to make sure YOUR correct input files will be used at each step
# fastqc runs
# mkdir fastqc_files
# for file in *.bam; do fastqc $file -t 30 -q --extract -f bam_mapped -o fastqc_files; done 
# cd fastqc_files
# rm *.zip *REF*.bam
# for f in $(find . -name '*data.txt'); do
    # mv "$f" "$(dirname "$f")/$(basename "$(dirname "$f")" | tr "/" "_")_$(basename "$f")"
# done
# find . -type f -name '*.data.txt' -exec cp {} ../ \;
# cd ../
# rm -r fastqc_files
# ls *fastqc.txt | parallel -j 20 'grep -A 1 "^#Base[[:blank:]]*G" {} | tail -n 1 | awk -v fname=$(basename {}) "{print \$0\"\t\"fname}" | cut -f4,6 | sed -E "s/\\.REF_|\\_fastqc|:|-/\\t/g; s/(rep[123])/\1\\t/g" | \
# cut -f1,2,4,5,6 > {.}.1U.txt'
# cat *1U.txt >fly_ovary_restrictive_piRNA_remapped_1U_percent.txt
# awk '$1 >= 60' fly_ovary_restrictive_piRNA_remapped_1U_percent.txt | awk '{print $3"\t"$4"\t"$5"\t"$2}' >fly_ovary_restrictive_piRNA_remapped_1U_percent_60p.bed

# deduplication (use dedupe.sh of bbmap if you want it to be faster and to specify the paremeters of deduplication)
# mkdir dedupe_files
# ls *.bam | parallel -j 15 'samtools rmdup -s {} - | samtools view -c -F 260 | awk -v fname=$(basename {}) "{print \$0\"\t\"fname}" > dedupe_files/{.}_uniq.txt' &
# cat dedupe_files/*uniq.txt >./fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled.txt
# sed -i 's/\.bam//g; s/:/\t/g; s/_piRNA_remap.REF_/\t/g; s/-/\t/g' fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled.txt
# awk -F'\t' -v OFS='\t' '{diff = $5 - $4; div = diff / 500; $6 = $1 / div} $6 >= 10 {print $3, $4, $5, $2}' fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled.txt \
# >fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb.bed

# generate the restrictive piC bed files by filtering the *merge10kb.bed entries using the dedupe counts and 1U percent reads. 
# bedtools intersect -a fly_ovary_restrictive_piRNA_remapped_1U_percent_60p.bed -b fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb.bed >fly_ovary_restrictive_piRNA_remapped_1U60p_10piRNAperKb.bed 

# copy the matching reference genome's mappability bedGraph file from the GEM mappability script 
# create a 10kb windows of the releveant genome and then map the mappability bedGraph values to the 10kb windows 
# bedtools makewindows -g genome.sizes -w 100000 >genome_10kb.bed
# bedtools map -a genome_10kb.bed -b genome_25nt_mappability.bedGraph >genome_25nt_mappability_10kb.bed
# awk '$4 =< 0.1' genome_25nt_mappability_10kb.bed | awk '{print $1"\t"$2"\t"$3"}' >genome_25nt_OnlyLowMap_10kb.bed
# bedtools closest -a genomic_25nt_OnlyLowMap_10kb.bed -b fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb.bed -d  | awk '$9 < 10 {print $1"\t"$2"\t"$3}' | sort -u  >fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb_LowMap10kbOnly.bed

# now combined those very close unmappable windows with the filtered restrictive piC files: 
# cat fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb_LowMap10kbOnly fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb.bed | bedtools sort -i \
- | bedtools merge -i - -d 11 > fly_ovary_restrictive_piRNA_remapped_dedupe_CT_compiled_10piRNAperKb_LowMap10kbOnly_merged.bed


