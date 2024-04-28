## essential steps are laid out for analyses of mapped long pacbio reads to iso-1 genome using minimap2 -x map-pb to detect INDELs and non-reference TE insertionsx 
## download long reads 


nohup sh -c 'for (( i =48; i <=56; i++ )); do fastq-dump --accession SRR78416$i --gzip ; done' &

## cuteSV
cuteSV --max_cluster_bias_INS 100 --diff_ratio_merging_INS	0.3 --max_cluster_bias_DEL	200 
--diff_ratio_merging_DEL 0.5 --min_read_len 1000 --min_support 5 --threads 10 & 

ls *.bam | parallel -j 5 'cuteSV --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 \
--min_support 5 --threads 10 -l 30 --batches 5000000 {} ISO1_genome_chrom.fa {.}.cuteSV.vcf ./' &

## svim
source /home/sps257/miniconda/bin/activate svim
conda activate svim_env

nohup svim alignment --min_mapq 40 --min_sv_size 30 --types DEL,INS,INV,DUP:TANDEM --heterozygous_threshold 0.05 \
--minimum_depth 11 A2_svim A2_pacbio_minimap2_ISO1_chrom.sort.bam ISO1_genome_chrom.fa >A2_svim/A2_svim.log &

bcftools view -i "QUAL >=20" variants.vcf >A2_svim_q20.vcf

## sniffles
/programs/sniffles-1.0.12/bin/sniffles-core-1.0.12/sniffles -m A2_pacbio_ngmlr.sort.bam -v A2_pacbio_ngmlr.sniffles.vcf -f 0.05 -l 30 -t 20 &

## merging step using SURVIVOR 

## SV genotyping
nohup cuteSV --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 \
--min_read_len 1000 --min_support 5 --threads 6 --min_size 30 --genotype A2_pacbio_minimap2_ISO1_chrom.sort.bam \
/fs/cbsuclarkfs1/storage/sps257/long_reads_wgs/ISO1_genome_chrom.fa A2_pacbio_minimap2_svim_gt.vcf A2_svim_gt/ -Ivcf Dmel_vcf_files_raw_svim_merged_150bpDist.vcf &

## after removal of TRA,BND SVs and SVs without AF
for file in *gt.bed; do bedtools merge -i $file -d 90 -c 5,6,7,8 -o distinct,mean,sum,mean >${file/.bed/collapsed.bed};done &
sed -i '/,/d' *collapsed.bed

# TLDR
 export PATH=/programs/exonerate-2.2.0-x86_64/bin:$PATH
 export PATH=/programs/minimap2-2.17:$PATH
 export PATH=/programs/htslib-1.14/bin:$PATH
 export PATH=/programs/samtools-1.14/bin:$PATH
 export PATH=/programs/mafft/bin:$PATH

nohup tldr/tldr -b A2_pacbio_minimap2_ISO1_chrom.sort.bam -r ISO1_genome_chrom.fa -e Dmel_repmod2_TEfam_library.fasta 
-p 12 -o A2_pacbio_minimap2_ISO1_tldr --color_consensus --detail_output --flanksize 200 --max_te_len 15000 
-m 2 >A2_pacbio_minimap2_ISO1_tldr.log &

## tldr selection
sed -i 's/|/\t/g' *.bed
cut -f2,3,4,5,6,7,13,14,15,19,20,21 A1_pacbio_minimap2_ISO1_tldr.table.txt >A1_pacbio_minimap2_ISO1_tldr.table.bed
MedianMapQ >20 and TEMatch>80



