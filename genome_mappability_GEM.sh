#!/bin/sh

# load the following version of GEM library (20130406-045632) and newest version of UCSC Kent utilities
export PATH=/programs/GEM/bin:$PATH
export PATH=/programs/kentUtils/bin:$PATH

ls *.fa | parallel -j 8 'gem-indexer -T 10 -i {}.fa -o {.}_gem.index'
ls *.gem | parallel -j 8 'gem-mappability -T 10 -I {} -l 25 -o {.}.25nt'
ls *.gem | awk -F "_gem.index.gem" '{print $1}' | sort -u | parallel -j 8 "gem-2-wig -i {}_gem.25nt.mappability -I {}_gem.index.gem -o {.}_25nt"
ls *.wig | awk -F "_25nt.wig" '{print $1}' | sort -u | parallel -j 8 "wigToBigWig {} {}.sizes {.}_25nt.bw"
ls *.bw | parallel -j 8 'bigWigToBedGraph {} {.}_25nt_mappability.bedGraph' 

