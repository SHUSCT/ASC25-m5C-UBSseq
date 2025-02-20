#!/bin/bash -e
eval "$(conda shell.bash hook)"
conda activate m5C
# hisat-3n/hisat-3n-build -p 64 --base-change C,T ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ./ref/index/hisat3n/Homo_sapiens.GRCh38.dna.primary_assembly
# samtools/samtools faidx ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai > ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.saf
# hisat-3n/hisat-3n-build -p 64 --base-change C,T ./ref/genome/Homo_sapiens.GRCh38.ncrna.fa ./ref/index/hisat3n/Homo_sapiens.GRCh38.ncrna
# samtools/samtools faidx ./ref/genome/Homo_sapiens.GRCh38.ncrna.fa

start_time=$(date +%s.%N)
parallel --linebuffer --tag '
DATA={}
step_start=$(date +%s.%N)
cutseq ./${DATA}/${DATA}.fastq -t 40 -A INLINE -m 20 --trim-polyA --ensure-inline-barcode \
    -o ./${DATA}/${DATA}.fastq_cut \
    -s ./${DATA}/${DATA}.fastq_tooshort \
    -u ./${DATA}/${DATA}.fastq_untrimmed --json-file ./${DATA}/${DATA}.trimming.json 2> /dev/null
step_end=$(date +%s.%N)
echo "cutseq took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
hisat-3n/hisat-3n --index ./ref/index/hisat3n/Homo_sapiens.GRCh38.ncrna \
                  --summary-file ./${DATA}/map2ncrna.output.summary --new-summary -q \
                  -U ./${DATA}/${DATA}.fastq_cut \
                  -p 36  --base-change C,T --mp 8,2 --no-spliced-alignment --directional-mapping \
                  | samtools/samtools view -@ 4 -e '!flag.unmap' -O BAM -U ./${DATA}/${DATA}.ncrna.unmapped.bam \
                                                                           -o ./${DATA}/${DATA}.ncrna.mapped.bam
step_end=$(date +%s.%N)
echo "hisat-3n rRNA, tRNA filtering took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
samtools/samtools fastq -@ 40 -O ./${DATA}/${DATA}.ncrna.unmapped.bam > ./${DATA}/${DATA}.mRNA.fastq
step_end=$(date +%s.%N)
echo "samtools fastq took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
hisat-3n/hisat-3n --index ./ref/index/hisat3n/Homo_sapiens.GRCh38.dna.primary_assembly -p 36 \
                  --summary-file ./${DATA}/map2genome.output.summary --new-summary -q \
                  -U ./${DATA}/${DATA}.mRNA.fastq --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 \
                  | samtools/samtools view -@ 4 -e '!flag.unmap' -O BAM -U ./${DATA}/${DATA}.mRNA.genome.unmapped.bam \
                                                                         -o ./${DATA}/${DATA}.mRNA.genome.mapped.bam
step_end=$(date +%s.%N)
echo "hisat-3n genome alignment took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
samtools/samtools sort -@ 40 --write-index  -O BAM -o ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam ./${DATA}/${DATA}.mRNA.genome.mapped.bam
step_end=$(date +%s.%N)
echo "samtools sort took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
samtools/samtools view -@ 40 -F 3980 -c ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam >./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam.tsv
step_end=$(date +%s.%N)
echo "samtools view took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=./${DATA} \
     -jar ./UMICollapse/umicollapse.jar bam -T 40 --data bktree --merge avgqual --two-pass \
     -i ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam \
     -o ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam \
     > ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.log
step_end=$(date +%s.%N)
echo "UMICollapse deduplication took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
samtools/samtools view -@ 40 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" -O BAM -o ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.filtered.bam ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam
step_end=$(date +%s.%N)
echo "samtools filter took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
samtools/samtools index -@ 40 ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam.bai
step_end=$(date +%s.%N)
echo "samtools index took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam \
| hisat-3n/hisat-3n-table -p 1 -u --alignments - --ref ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_unfiltered_uniq.tsv.gz &
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.bam \
| hisat-3n/hisat-3n-table -p 1 -m --alignments - --ref ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_unfiltered_multi.tsv.gz &
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.filtered.bam \
| hisat-3n/hisat-3n-table -p 1 -u --alignments - --ref ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_filtered_uniq.tsv.gz &
samtools/samtools view -e "rlen<100000" -h ./${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.filtered.bam \
| hisat-3n/hisat-3n-table -p 1 -m --alignments - --ref ./ref/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-name /dev/stdout --base-change C,T \
| cut -f 1,2,3,5,7 | gzip -c > ./${DATA}/${DATA}_filtered_multi.tsv.gz &
wait
step_end=$(date +%s.%N)
echo "hisat-3n-table site calling took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
python ./bin/join_pileup.py -i ./${DATA}/${DATA}_unfiltered_uniq.tsv.gz \
                               ./${DATA}/${DATA}_unfiltered_multi.tsv.gz \
                               ./${DATA}/${DATA}_filtered_uniq.tsv.gz \
                               ./${DATA}/${DATA}_filtered_multi.tsv.gz \
                            -o ./${DATA}/${DATA}_genome.arrow 2> /dev/null
step_end=$(date +%s.%N)
echo "join_pileup.py took $(echo "$step_end - $step_start" | bc) seconds."
' ::: SRR23538290 SRR23538291 SRR23538292

step_start=$(date +%s.%N)
python ./bin/group_pileup.py -i ./SRR23538290/SRR23538290_genome.arrow ./SRR23538291/SRR23538291_genome.arrow ./SRR23538292/SRR23538292_genome.arrow -o WT.arrow  2> /dev/null
step_end=$(date +%s.%N)
echo "group_pileup.py took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
python ./bin/select_sites.py -i ./WT.arrow -o ./WT.prefilter.tsv  2> /dev/null
step_end=$(date +%s.%N)
echo "select_sites.py took $(echo "$step_end - $step_start" | bc) seconds."

step_start=$(date +%s.%N)
python ./bin/filter_sites.py -i  ./SRR23538290/SRR23538290_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538290/SRR23538290.bg.tsv -o ./SRR23538290/SRR23538290.filtered.tsv 2> /dev/null
python ./bin/filter_sites.py -i  ./SRR23538291/SRR23538291_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538291/SRR23538291.bg.tsv -o ./SRR23538291/SRR23538291.filtered.tsv 2> /dev/null
python ./bin/filter_sites.py -i  ./SRR23538292/SRR23538292_genome.arrow -m ./WT.prefilter.tsv -b ./SRR23538292/SRR23538292.bg.tsv -o ./SRR23538292/SRR23538292.filtered.tsv 2> /dev/null
step_end=$(date +%s.%N)
echo "filter_sites.py took $(echo "$step_end - $step_start" | bc) seconds."

end_time=$(date +%s.%N)
echo "Total time: $(echo "$end_time - $start_time" | bc) seconds."

for DATA in SRR23538290 SRR23538291 SRR23538292
do
mv ${DATA}/${DATA}.filtered.tsv .
mv ${DATA}/${DATA}.trimming.json .
mv ${DATA}/${DATA}.mRNA.genome.mapped.sorted.bam.tsv .
mv ${DATA}/${DATA}.bg.tsv .
mv ${DATA}/map2genome.output.summary ${DATA}.map2genome.output.summary 
mv ${DATA}/map2ncrna.output.summary ${DATA}.map2ncrna.output.summary
mv ${DATA}/${DATA}.mRNA.genome.mapped.sorted.dedup.log .
done