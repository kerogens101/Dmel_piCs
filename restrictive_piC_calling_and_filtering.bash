#!/bin/sh
## this script conssists of following stages: 
## mapping, filter by cutoff, analyze and filter by piRNA characteristics, and merge with low mappability scores
## Index your genomes and map respective piRNA reads to their genomes
## make sure the first two or thhree words in each genome file's name and correspdonding small RNA library name is the same.  
ls *genome.fasta | parallel -j 8 'bowtie2-build {} {.} --threads 4' 
mkdir piC_calling/
## first add echo before -- 'bowtie -- and do a trial run. See if the individual commands make sense and change the cut and file extensions to match your inputs. Then remove the echo and run the changed command. 
ls *piRNA.fastq.gz | cut -d '_' -f1-2 | sort -u | parallel -j 4 'bowtie2 -p 5 -N 1 -L 12 -k 1 -q {}_ovary_piRNA.fastq.gz -x {}_genome -S piC_calling/{.}_genome.uniq.sam --no-unal' 
## create sorted and indexed bam file for coverage
cd piC_calling/
ls *uniq.sam | parallel -j 5 'samtools view -bS -@ 10 {} >{.}.bam'
ls *uniq.bam | parallel -j 5 'samtools sort -@ 10 {} >{.}.sorted.bam'
ls *uniq.sorted.bam | parallel -j 5 'samtools index -@ 10 {}'
rm *uniq.sam *uniq.bam
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8
export PATH=/programs/deeptools-3.5.1/bin:$PATH
export PYTHONPATH=/programs/deeptools-3.5.1/lib64/python3.9/site-packages:/programs/deeptools-3.5.1/lib/python3.9/site-packages
## calculate piRNA coverage from bamCoverage and remove windows with less 5 CPM values. 
ls *sorted.bam | parallel -j 6 'bamCoverage -bs 500 -b {} -of bedgraph -p 2 -o {.}_500bp.bedGraph --normalizeUsing CPM' 
ls *500bp.CPM.bedGraph | parallel -j 10 "awk -v OFS='\t' '{ if (\$4>5) print \$0}' {} > {.}_5CPM.bed"
for file in *5CPM.bed; do
    sort -k1,1 -k2,2n "$file" | bedtools merge -d 10 -c 4 -o sum > "${file%.bed}_merge10kb.bed"
done
## extract the piRNA-enriched domains from the respective genome files 
ls *merge10kb.bed | cut -d '_' -f1-2 | sort -u | parallel -j 4 'bedtools getfasta -f ../{}_genome.fa -bed {}_ovary_piRNA_genome.uniq.5CPM_merge10kb.bed > {.}_ovary_piRNA_genome.uniq.5CPM_merge10kb_domain.fa' 
mkdir ../remap_and_filter
cp *domain.fa ../remap_and_filter
cd ../remap_and_filter
sed -i 's/:/_/g' *.fa
ls ovary_piRNA_genome.uniq.5CPM_merge10kb_domain.fa | parallel -j 5 'bowtie2-build -f {} {.}' 
ls *.fa | cut -d '_' -f1-4 | sort -u | parallel -j 5 'bowtie2 -N 1 -L 14 -k 1 -x {}_genome_uniq_500bp_5CPM.merge -q ../{}.fastq.gz -S {.}_remap.sam -p 5 --no-unal' 
ls *.sam | parallel -j 5 'samtools view -bS -@ 10 {} >{.}.bam'
## split this bam file by each merged domain and obtain characteristics of piRNAs from each domain to filter by 1U >60% of the reads and minimum 20 unique piRNAs per kb. 
mkdir bamtools_split
for file in *bam ; do bamtools split -in $file -reference ; done &
mv *REF*.bam bamtools_split
cd bamtools_split
## fastqc runs
mkdir fastqc_files
for file in *.bam; do fastqc $file -t 30 -q --extract -f bam_mapped -o fastqc_files; done 
cd fastqc_files
rm *.zip *REF*.bam
for f in $(find . -name '*data.txt'); do
    mv "$f" "$(dirname "$f")/$(basename "$(dirname "$f")" | tr "/" "_")_$(basename "$f")"
done
find . -type f -name '*.data.txt' -exec cp {} ../ \;
cd ../
rm -r fastqc_files
ls *fastqc.txt | parallel -j 20 'grep -A 1 "^#Base[[:blank:]]*G" {} | tail -n 1 | awk -v fname=$(basename {}) "{print \$0\"\t\"fname}" | cut -f4,6 | sed -E "s/\\.REF_|\\_fastqc|:|-/\\t/g; s/(rep[123])/\1\\t/g" | \
cut -f1,2,4,5,6 > {.}.1U.txt'
cat *1U.txt >Zfish_gonad_restrictive_piRNA_remapped_1U_percent.txt
awk '$1 >= 60' Zfish_gonad_restrictive_piRNA_remapped_1U_percent.txt | awk '{print $3"\t"$4"\t"$5"\t"$2}' > Zfish_gonad_restrictive_piRNA_remapped_1U_percent_60p.bed
## deduplication
mkdir dedupe_files
ls *.bam | parallel -j 15 'samtools rmdup -s {} - | samtools view -c -F 260 | awk -v fname=$(basename {}) "{print \$0\"\t\"fname}" > dedupe_files/{.}_uniq.txt' &
cat *dedupe_files/uniq.txt >gonad_restrictive_piRNA_remapped_dedupe_Cnt_compiled.txt
nohup sh -c 'for i in *.fastq; do /programs/bbmap-38.86/dedupe.sh in=$i out={$i}_dedupe.fastq threads=15 -Xmx60g ; done' >nohup_dedupe.txt &
wc -l *dedupe.fastq >unique_piRNA_domains_dedupe.txt
## now calculate
awk -F'\t' -v OFS='\t' '{diff = $5 - $4; div = diff / 500; $6 = $1 / div} $6 >= 10 {print $3, $4, $5, $2}' gonad_restrictive_piRNA_remapped_rmdup_count.txt >gonad_restrictive_piRNA_remapped_rmdup_count_10PRKb.bed