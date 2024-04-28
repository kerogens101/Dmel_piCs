#!/bin/sh
## you'll require the proTRAC_version.pl script
## pay attention to how your information of sample/library/processing is stored in the filenames. Have atleast 
## first prepare the fastq files. Use the processed piRNA.fastq files from the small RNA processing script. 
## Now use it to collapse the duplicate reads get a unique set of piRNA sequences. 
ls *_piRNA.fastq.gz | cut -d '_' -f1-3 | sort -u | parallel -j 2 '/programs/bbmap-39.03/dedupe.sh in={} out={.}_dedupe.fastq.gz threads=10' &
## then index your genomes and map respective piRNA reads to their genomes
## make sure the first two or thhree words in each genome file's name and correspdonding small RNA library name is the same.  
ls *genome.fasta | parallel -j 10 'bowtie2-build {} {.} --threads 4' 
mkdir proTRAC_calling/
## first add echo before -- 'bowtie -- and do a trial run. See if the individual commands make sense and change the cut and file extensions to match your inputs. Then remove the echo and run the changed command. 
ls *.dedupe.fastq.gz | cut -d '_' -f1-2 | sort -u | parallel -j 4 'bowtie2 -p 5 -N 1 -L 12 -k 50 -q {}_ovary_piRNA_dedupe.fastq.gz -x {}_genome -S ./proTRAC_calling/{.}_genome.multi.sam' &
## once you get the SAM files, carry out proTRAC command. Adjust the parameters accordingly to your species, amount of reference genome piRNA clusters etc.
cd proTRAC_calling/
ls *.sam | cut -d '_' -f1-2 | sort -u | parallel -j 5 'perl ../proTRAC_2.4.4.pl -map {}_genome.multi.sam -g {}_genome.fasta -swsize 5000 -swincr 1000 -1Tor10A 0.5 \
-clsize 1000 -pimin 24 -pimax 32 -pdens 0.05 -clhitsn 10 -pti -nohtml -nofaspi' 
## once you have the proTRAC files, change the name of the clusters.gtf file to the name of the respective directory its in
for f in $(find . -name '*data.txt'); do mv "$f" "$(dirname "$f")/$(dirname "${f:2}" | tr "/" "_")_$(basename $f)" ; done 
## then copy the clusters.gtf files to current directory
find . -type f -name '*_clusters.gtf' -exec cp {} . \;
## remove the extra info from the filenames
for filename in proTRAC*clusters.gtf; do new_filename=$(echo "$filename" | sed 's/.sam.*_clusters/_clusters/');     mv "$filename" "$new_filename";     echo "File '$filename' renamed to '$new_filename'"; done
## In the gtf files, replace semicolon and colon with tab and then select columns to convert to bed file. 
sed -i 's/[:;]/\t/g' *clusters.gtf
for file in *clusters.gtf; do cut -f1,4,5,15 $file >${file/.gtf/_results.txt};done &
