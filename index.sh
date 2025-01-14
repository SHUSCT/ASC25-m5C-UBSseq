#!/bin/bash -e
hisat-3n/hisat-3n-build -p 64 --base-change C,T ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./ref/index/hisat3n/Homo_sapiens.GRCh38.dna.primary_assembly
samtools/samtools faidx ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf
hisat-3n/hisat-3n-build -p 64 --base-change C,T ./ref/genome/Homo_sapiens.GRCh38.ncrna.fa ./ref/index/hisat3n/Homo_sapiens.GRCh38.ncrna
samtools/samtools faidx ./ref/genome/Homo_sapiens.GRCh38.ncrna.fa