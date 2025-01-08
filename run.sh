#!/bin/bash
hisat-3n/hisat-3n-build -p 64 --base-change C,T ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools/samtools faidx ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf
hisat-3n/hisat-3n-build -p 64 --base-change C,T ./ref/Homo_sapiens.GRCh38.ncrna.fa ./ref/Homo_sapiens.GRCh38.ncrna.fa
samtools/samtools faidx ./ref/Homo_sapiens.GRCh38.ncrna.fa

for DATA in "SRR23538290" "SRR23538291" "SRR23538292"
do
cutseq ./${DATA}/${DATA}.fastq -t 64 -A INLINE -m 20 --trim-polyA --ensure-inline-barcode \
    -o ./${DATA}/${DATA}.fastq_cut \
    -s ./${DATA}/${DATA}.fastq_tooshort \
    -u ./${DATA}/${DATA}.fastq_untrimmed
hisat-3n/hisat-3n --index ./ref/Homo_sapiens.GRCh38.ncrna.fa \
                  --summary-file ./${DATA}/map2ncrna.output.summary --new-summary -q \
                  -U ./${DATA}/${DATA}.fastq_cut \
                  -p 32  --base-change C,T --mp 8,2 --no-spliced-alignment --directional-mapping \
                  | samtools/samtools view -@ 32 -e '!flag.unmap' -O BAM -U ./${DATA}/${DATA}.ncrna.unmapped.bam \
                                                                           -o ./${DATA}/${DATA}.ncrna.mapped.bam
samtools/samtools fastq -@ 64 -O ./${DATA}/${DATA}.ncrna.unmapped.bam >./${DATA}/${DATA}.mRNA.fastq
hisat-3n/hisat-3n --index ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa -p 32 \
                  --summary-file ./${DATA}/map2genome.output.summary --new-summary -q \
                  -U ./${DATA}/${DATA}.mRNA.fastq --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 \
                  | samtools/samtools view -@ 32 -e '!flag.unmap' -O BAM -U ./${DATA}/${DATA}.mRNA.genome.unmapped.bam \
                                                                         -o ./${DATA}/${DATA}.mRNA.genome.mapped.bam
samtools/samtools sort -@ 64 --write-index  -O BAM -o ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam ./${DATA}/${DATA}.mRNA.genome.mapped.bam
samtools/samtools view -@ 64 -F 3980 -c ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam >./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam.tsv
java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=./${DATA} \
     -jar ./UMICollapse/umicollapse.jar bam -t 2 -T 32 --data naive --merge avgqual --two-pass \
     -i ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam \
     -o ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam \
     > ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.log
samtools/samtools view -@ 64 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" -O BAM -o ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.filtered.bam ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam
samtools/samtools index -@ 64 ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam.bai

samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam \
| hisat-3n/hisat-3n-table -p 64 -u --alignments - --ref ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_unfiltered_uniq.tsv.gz
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam \
| hisat-3n/hisat-3n-table -p 64 -m --alignments - --ref ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_unfiltered_multi.tsv.gz
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.filtered.bam \
| hisat-3n/hisat-3n-table -p 64 -u --alignments - --ref ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_filtered_uniq.tsv.gz
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.filtered.bam \
| hisat-3n/hisat-3n-table -p 64 -m --alignments - --ref ./ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_filtered_multi.tsv.gz

python ./bin/join_pileup.py -i ./${DATA}/${DATA}_unfiltered_uniq.tsv.gz \
                               ./${DATA}/${DATA}_unfiltered_multi.tsv.gz \
                               ./${DATA}/${DATA}_filtered_uniq.tsv.gz \
                               ./${DATA}/${DATA}_filtered_multi.tsv.gz \
                            -o ./${DATA}/${DATA}_genome.arrow
done

python ./bin/group_pileup.py -i ./SRR23538290/SRR23538290_genome.arrow ./SRR23538291/SRR23538291_genome.arrow ./SRR23538292/SRR23538292_genome.arrow -o WT.arrow
python ./bin/select_sites.py -i ./WT.arrow -o ./WT.prefilter.tsv
python ./bin/filter_sites.py -i  ./SRR23538290/SRR23538290_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538290/SRR23538290.bg.tsv -o ./SRR23538290/SRR23538290.filtered.tsv
python ./bin/filter_sites.py -i  ./SRR23538291/SRR23538291_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538291/SRR23538291.bg.tsv -o ./SRR23538291/SRR23538291.filtered.tsv
python ./bin/filter_sites.py -i  ./SRR23538292/SRR23538292_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538292/SRR23538292.bg.tsv -o ./SRR23538292/SRR23538292.filtered.tsv
