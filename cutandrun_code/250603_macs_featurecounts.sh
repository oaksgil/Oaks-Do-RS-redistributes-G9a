conda activate snakemake
cd /home/briando/data/files/bd523_bd524/results
mkdir macs2
cd macs2


# Averaging
cd /home/briando/data/files/bd523_bd524/results/bigwigs
parallel -j 12 "sbatch -J map -N 1 --wrap \" \
bigwigAverage -bs 5 -p 10 --skipNAs -b clone32__{}__rep*_no_rmdup_norm.bw -o clone32__{}_avg_no_rmdup_norm.bw \" "  ::: \
HU__G9a HU__H3K27ac HU__H3K27me3 HU__H3K9ac HU__H3K9me2 HU__WIZ \
Unt__G9a Unt__H3K27ac Unt__H3K27me3 Unt__H3K9ac Unt__H3K9me2 Unt__WIZ \
HU__H3K9me3 Unt__H3K9me3

parallel -j 12 "sbatch -J map -N 1 --exclude=c[5-24] --wrap \" \
bigwigAverage -bs 5 -p 10 --skipNAs -b norm_to_mouse/clone32__{}__rep*_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/clone32__{}_avg_no_rmdup_norm_to_mouse.bw \" "  ::: \
HU__G9a HU__H3K27ac HU__H3K27me3 HU__H3K9ac HU__H3K9me2 HU__WIZ \
Unt__G9a Unt__H3K27ac Unt__H3K27me3 Unt__H3K9ac Unt__H3K9me2 Unt__WIZ \
HU__H3K9me3 Unt__H3K9me3

##ATAC
cd /net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs
parallel -j 12 "sbatch -J map -N 1 --exclude=c[5-24] --wrap \" \
bigwigAverage -bs 5 -p 10 --skipNAs -b {}_rep*.bw -o {}_avg.bw \" " ::: Unt_A Unt_F HU_A HU_F

cat *.narrowPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -i - > all_atac_peaks_merged.bed  
awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("atac_peak" NR,$1,$2,$3,"+")}}' all_atac_peaks_merged.bed > \
      all_atac_peaks_merged.saf

sbatch -J map -N 1 --wrap  "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
all_atac_peaks_merged.saf \
-o all_atac_peaks_merged_featurecounts.txt \
../bowtie2_bam/*.bam"

##clone A3
parallel -j 12 "sbatch -J map -N 1 --exclude=c[5-24] --wrap \" \
bigwigAverage -bs 5 -p 10 --skipNAs -b norm_to_mouse/cloneA3__{}__rep*_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/cloneA3__{}_avg_no_rmdup_norm_to_mouse.bw \" "  ::: \
Unt__H3K9ac HU__H3K9ac Unt__H3K9me2 HU__H3K9me2

parallel -j 12 "sbatch -J map -N 1 --exclude=c[5-24] --wrap \" \
bigwigAverage -bs 5 -p 10 --skipNAs -b cloneA3__{}__rep*_no_rmdup_norm.bw \
-o cloneA3__{}_avg_no_rmdup_norm.bw \" "  ::: \
Unt__H3K9ac HU__H3K9ac Unt__H3K9me2 HU__H3K9me2

sbatch -J map -N 1 --exclude=c[5-24] --wrap "\
bigwigAverage -bs 5 -p 10 --skipNAs -b \
norm_to_mouse/cloneA3__HU__H3K9me3__rep2_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__HU__H3K9me3__rep3_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__HU__H3K27ac__rep1_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/cloneA3__HU__H3K27me3_avg_no_rmdup_norm_to_mouse.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
norm_to_mouse/cloneA3__Unt__H3K27ac__rep1_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__Unt__H3K27ac__rep2_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__Unt__H3K27ac__rep3_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/cloneA3__Unt__H3K27me3_avg_no_rmdup_norm_to_mouse.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
norm_to_mouse/cloneA3__Unt__H3K9me3__rep1_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__Unt__H3K9me3__rep2_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__Unt__H3K27me3__rep1_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__Unt__H3K27me3__rep2_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/cloneA3__Unt__H3K9me3_avg_no_rmdup_norm_to_mouse.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
norm_to_mouse/cloneA3__Unt__H3K27me3__rep3_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__Unt__IgG__rep1_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/cloneA3__Unt__H3K27ac_avg_no_rmdup_norm_to_mouse.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
norm_to_mouse/cloneA3__HU__H3K27me3__rep1_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__HU__H3K27me3__rep2_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__HU__H3K27me3__rep3_no_rmdup_norm_to_mouse.bw \
norm_to_mouse/cloneA3__HU__IgG__rep1_no_rmdup_norm_to_mouse.bw \
-o norm_to_mouse/cloneA3__HU__H3K27ac_avg_no_rmdup_norm_to_mouse.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap "\
bigwigAverage -bs 5 -p 10 --skipNAs -b \
cloneA3__HU__H3K9me3__rep2_no_rmdup_norm.bw \
cloneA3__HU__H3K9me3__rep3_no_rmdup_norm.bw \
cloneA3__HU__H3K27ac__rep1_no_rmdup_norm.bw \
-o cloneA3__HU__H3K27me3_avg_no_rmdup_norm.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
cloneA3__Unt__H3K27ac__rep1_no_rmdup_norm.bw \
cloneA3__Unt__H3K27ac__rep2_no_rmdup_norm.bw \
cloneA3__Unt__H3K27ac__rep3_no_rmdup_norm.bw \
-o cloneA3__Unt__H3K27me3_avg_no_rmdup_norm.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
cloneA3__Unt__H3K9me3__rep1_no_rmdup_norm.bw \
cloneA3__Unt__H3K9me3__rep2_no_rmdup_norm.bw \
cloneA3__Unt__H3K27me3__rep1_no_rmdup_norm.bw \
cloneA3__Unt__H3K27me3__rep2_no_rmdup_norm.bw \
-o cloneA3__Unt__H3K9me3_avg_no_rmdup_norm.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
cloneA3__Unt__H3K27me3__rep3_no_rmdup_norm.bw \
cloneA3__Unt__IgG__rep1_no_rmdup_norm.bw \
-o cloneA3__Unt__H3K27ac_avg_no_rmdup_norm.bw"

sbatch -J map -N 1 --exclude=c[5-24] --wrap " \
bigwigAverage -bs 5 -p 10 --skipNAs -b \
cloneA3__HU__H3K27me3__rep1_no_rmdup_norm.bw \
cloneA3__HU__H3K27me3__rep2_no_rmdup_norm.bw \
cloneA3__HU__H3K27me3__rep3_no_rmdup_norm.bw \
cloneA3__HU__IgG__rep1_no_rmdup_norm.bw \
-o cloneA3__HU__H3K27ac_avg_no_rmdup_norm.bw"


# MACS2 narrow peaks for G9a, WIZ, and H3K27ac

parallel -j 2 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c /net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bowtie2_bam/IgG_0h_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs \
-n {/.} 2> {/.}.log \"" ::: \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bowtie2_bam/BCL11A_24h_aligned_sorted_noblacklist.bam \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bowtie2_bam/BCL11A_0h_aligned_sorted_noblacklist.bam

parallel -j 7 --dryrun "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__HU__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__HU__G9a*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__HU__WIZ*aligned_sorted_noblacklist.bam

parallel -j 7 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__Unt__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__Unt__G9a*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__Unt__WIZ*aligned_sorted_noblacklist.bam 

parallel -j 7 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__HU__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__HU__H3K27ac*aligned_sorted_noblacklist.bam 

parallel -j 7 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__Unt__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__Unt__H3K27ac*aligned_sorted_noblacklist.bam

# broad peaks for H3K9ac
parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__HU__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__HU__H3K9ac*aligned_sorted_noblacklist.bam

parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__Unt__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__Unt__H3K9ac*aligned_sorted_noblacklist.bam

# broad peaks for H3K9me2 and H3K9me3
parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__HU__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__HU__H3K9me2*aligned_sorted_noblacklist.bam

parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__Unt__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__Unt__H3K9me2*aligned_sorted_noblacklist.bam

parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__HU__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__HU__H3K9me3*aligned_sorted_noblacklist.bam

parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__Unt__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__Unt__H3K9me3*aligned_sorted_noblacklist.bam

# broad peaks for H3K27me3
parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__HU__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__HU__H3K27me3*aligned_sorted_noblacklist.bam

parallel -j 3 "sbatch -J map -N 1 --wrap \"macs2 callpeak -t {} \
-c ../bowtie2_bam/clone32__Unt__IgG__rep1_aligned_sorted_noblacklist.bam \
-f BAMPE -g hs --broad \
-n {/.} 2> {/.}.log \"" ::: ../bowtie2_bam/clone32__Unt__H3K27me3*aligned_sorted_noblacklist.bam

# take the union of all peaks and then merges them
cat *G9a*.narrowPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -i - > g9a_peaks_merged.bed
cat *WIZ*.narrowPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -i - > wiz_peaks_merged.bed
cat *H3K9ac*.broadPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -i - > h3k9ac_peaks_merged.bed
cat *H3K9me2*.broadPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -d 5000 -i - > h3k9me2_peaks_merged.bed
cat *H3K9me3*.broadPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -d 5000 -i - > h3k9me3_peaks_merged.bed
cat *H3K27me3*.broadPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -d 5000 -i - > h3k27me3_peaks_merged.bed

# make SAF file
awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("g9a_peak" NR,$1,$2,$3,"+")}}' g9a_peaks_merged.bed > \
      g9a_peaks_merged.saf
awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("wiz_peak" NR,$1,$2,$3,"+")}}' wiz_peaks_merged.bed > \
      wiz_peaks_merged.saf
awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("h3k9ac_peak" NR,$1,$2,$3,"+")}}' h3k9ac_peaks_merged.bed > \
      h3k9ac_peaks_merged.saf

awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("h3k9me2_peak" NR,$1,$2,$3,"+")}}' h3k9me2_peaks_merged.bed > \
      h3k9me2_peaks_merged.saf

awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("h3k9me3_peak" NR,$1,$2,$3,"+")}}' h3k9me3_peaks_merged.bed > \
      h3k9me3_peaks_merged.saf

awk 'BEGIN {{OFS="\t"; print("GeneID", "Chr", "Start", "End", "Strand")}}; \
      {{print("h3k27me3_peak" NR,$1,$2,$3,"+")}}' h3k27me3_peaks_merged.bed > \
      h3k27me3_peaks_merged.saf

# featureCounts for everything over these peaks

sbatch -J map -N 1 --wrap  "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
g9a_peaks_merged.saf \
-o g9a_featurecounts.txt \
../bowtie2_bam/clone32__*G9a*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*WIZ*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me2*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me3*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27me3*aligned_sorted_noblacklist.bam"

sbatch -J map -N 1 --wrap "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
wiz_peaks_merged.saf \
-o wiz_featurecounts.txt \
../bowtie2_bam/clone32__*G9a*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*WIZ*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me2*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me3*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27me3*aligned_sorted_noblacklist.bam"

sbatch -J map -N 1 --wrap "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
h3k9ac_peaks_merged.saf \
-o h3k9ac_featurecounts.txt \
../bowtie2_bam/clone32__*G9a*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*WIZ*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me2*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me3*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27me3*aligned_sorted_noblacklist.bam"

sbatch -J map -N 1 --wrap "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
h3k9me2_peaks_merged.saf \
-o h3k9me2_featurecounts.txt \
../bowtie2_bam/clone32__*H3K9me2*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me3*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27me3*aligned_sorted_noblacklist.bam"

sbatch -J map -N 1 --wrap "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
h3k9me3_peaks_merged.saf \
-o h3k9me3_featurecounts.txt \
../bowtie2_bam/clone32__*H3K9me2*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me3*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27me3*aligned_sorted_noblacklist.bam"

sbatch -J map -N 1 --wrap "featureCounts -T 10 -p -O -M --fraction -F SAF -a \
h3k27me3_peaks_merged.saf \
-o h3k27me3_featurecounts.txt \
../bowtie2_bam/clone32__*H3K9me2*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K9me3*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27ac*aligned_sorted_noblacklist.bam \
../bowtie2_bam/clone32__*H3K27me3*aligned_sorted_noblacklist.bam"


# deepTools to generate plots
# computeMatrix
sbatch -J map -N 1 -n8 --wrap "computeMatrix reference-point -S \
../bigwigs/clone32__HU__G9a__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__G9a__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__G9a__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__G9a__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__G9a__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__G9a__rep3_no_rmdup_norm.bw \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed \
--referencePoint center -a 2000 -b 2000 -bs 50 -p 8 -o 250917_g9a_cr.gz"

sbatch -J map -N 1 -n8 --wrap "computeMatrix reference-point -S \
../bigwigs/clone32__HU__WIZ__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__WIZ__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__WIZ__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__WIZ__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__WIZ__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__WIZ__rep3_no_rmdup_norm.bw \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed \
--referencePoint center -a 2000 -b 2000 -bs 50 -p 8 -o 250917_wiz_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix reference-point -S \
../bigwigs/clone32__HU__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me2__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me2__rep3_no_rmdup_norm.bw \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed h3k9ac_peaks_merged.bed \
--referencePoint center --maxThreshold 5 -a 5000 -b 5000 -bs 50 -p 8 -o 250917_h3k9me2_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix reference-point -S \
../bigwigs/cloneA3__HU__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/cloneA3__HU__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/cloneA3__HU__H3K9me2__rep3_no_rmdup_norm.bw \
../bigwigs/cloneA3__Unt__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/cloneA3__Unt__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/cloneA3__Unt__H3K9me2__rep3_no_rmdup_norm.bw \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed h3k9ac_peaks_merged.bed \
--referencePoint center --maxThreshold 5 -a 5000 -b 5000 -bs 50 -p 8 -o 250917_cloneA3_h3k9me2_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix reference-point -S \
../bigwigs/clone32__HU__H3K9me3__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me3__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me3__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me3__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me3__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me3__rep3_no_rmdup_norm.bw \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed h3k9ac_peaks_merged.bed \
--referencePoint center --maxThreshold 5 -a 5000 -b 5000 -bs 50 -p 8 -o 250917_h3k9me3_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix reference-point -S \
../bigwigs/clone32__HU__H3K9ac__rep1_no_rmdup_norm.bw  \
../bigwigs/clone32__HU__H3K9ac__rep2_no_rmdup_norm.bw  \
../bigwigs/clone32__HU__H3K9ac__rep3_no_rmdup_norm.bw  \
../bigwigs/clone32__Unt__H3K9ac__rep1_no_rmdup_norm.bw  \
../bigwigs/clone32__Unt__H3K9ac__rep2_no_rmdup_norm.bw  \
../bigwigs/clone32__Unt__H3K9ac__rep3_no_rmdup_norm.bw  \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed h3k9ac_peaks_merged.bed \
--referencePoint center -a 2000 -b 2000 -p 8 -bs 50 -o 250917_h3k9ac_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix reference-point -S \
../bigwigs/clone32__HU__H3K27ac__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K27ac__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K27ac__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K27ac__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K27ac__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K27ac__rep3_no_rmdup_norm.bw \
-R wiz_peaks_merged.bed g9a_peaks_merged.bed h3k9ac_peaks_merged.bed \
--referencePoint center -a 2000 -b 2000 -p 8 -bs 50 -o 250917_h3k27ac_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix scale-regions -S \
../bigwigs/cloneA3__HU__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/cloneA3__HU__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/cloneA3__HU__H3K9me2__rep3_no_rmdup_norm.bw \
../bigwigs/cloneA3__Unt__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/cloneA3__Unt__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/cloneA3__Unt__H3K9me2__rep3_no_rmdup_norm.bw \
-R h3k9me2_peaks_merged.bed h3k9me3_peaks_merged.bed \
-a 5000 -b 5000 -m 10000 -p 8 -bs 50 -o 250917_cloneA3_h3k9me2_repressed_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix scale-regions -S \
../bigwigs/clone32__HU__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me2__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me2__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me2__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me2__rep3_no_rmdup_norm.bw \
-R h3k9me2_peaks_merged.bed h3k9me3_peaks_merged.bed \
-a 5000 -b 5000 -m 10000 -p 8 -bs 50 -o 250917_h3k9me2_repressed_cr.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix scale-regions -S \
../bigwigs/clone32__HU__H3K9me3__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me3__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__HU__H3K9me3__rep3_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me3__rep1_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me3__rep2_no_rmdup_norm.bw \
../bigwigs/clone32__Unt__H3K9me3__rep3_no_rmdup_norm.bw \
-R h3k9me2_peaks_merged.bed h3k9me3_peaks_merged.bed \
-a 5000 -b 5000 -m 10000 --maxThreshold 5 -p 8 -bs 50 -o 250917_h3k9me3_repressed_cr.gz"

# RPA at H3K9me3 peaks
bedtools merge -d 10000 -i h3k9me3_peaks_merged.bed > h3k9me3_peaks_merged_10kb.bed
sbatch -J map -N 1 -n 8 --wrap  "computeMatrix scale-regions -S \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/RPA_0h_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/RPA_6h_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/RPA_12h_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/RPA_24h_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/RPA_48h_norm.bw \
-R h3k9me3_peaks_merged_100kb.bed \
-a 20000 -b 20000 -m 40000 -bs 1000 -p 8 -o 250811_rpa_at_k9me3_10kb.gz \
--outFileNameMatrix 250811_rpa_at_k9me3_10kb.tab"

# computeRegions for regions that change H3K9me2/3 and H3K27me3
sbatch -J map -N 1 -n 8 --wrap  "computeMatrix scale-regions -S \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/unt_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/Unt_A_avg.bw \
-R 250809_k9me2_down.bed 250809_k9me2_unch.bed 250809_k9me2_up.bed \
--missingDataAsZero --skipZeros -a 5000 -b 5000 -m 10000 -p 8 -bs 50 -o 250918_k9me2_chg.gz"

sbatch -J map -N 1 -n 8 --wrap  "computeMatrix scale-regions -S \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/unt_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/Unt_A_avg.bw \
-R 250809_k27me3_down.bed 250809_k27me3_unch.bed 250809_k27me3_up.bed \
--missingDataAsZero --skipZeros -a 5000 -b 5000 -m 10000 -bs 50 -p 8 -o 250918_k27me3_chg.gz"

sbatch -J map -N 1 -n 15 --wrap  "computeMatrix reference-point -S \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/hu50_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/unt_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/hu50_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/hu50_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/Unt_A_avg.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/HU_A_avg.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__G9a_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__G9a_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__WIZ_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__WIZ_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/BCL11A_0h_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/BCL11A_24h_no_rmdup_norm.bw \
-R 250809_k9me2_down_atac.bed 250809_k9me2_unch_atac.bed 250809_k9me2_up_atac.bed \
--referencePoint center --missingDataAsZero --skipZeros -a 2000 -b 2000 -bs 50 -p 15 -o 250918_k9me2_chg_withhu_atacpeaks.gz"


sbatch -J map -N 1 -n 15 --wrap  "computeMatrix reference-point -S \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/hu50_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/unt_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/hu50_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/hu50_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/Unt_A_avg.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/HU_A_avg.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__G9a_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__G9a_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__WIZ_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__WIZ_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/BCL11A_0h_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/BCL11A_24h_no_rmdup_norm.bw \
-R 250809_k27me3_down_atac.bed 250809_k27me3_unch_atac.bed 250809_k27me3_up_atac.bed \
--referencePoint center --missingDataAsZero --skipZeros -a 2000 -b 2000 -p 15 -bs 50 -o 250918_k27me3_chg_withhu_atacpeaks.gz"

sbatch -J map -N 1 -n 15 --wrap  "computeMatrix reference-point -S \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/hu50_H3K4me1_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/unt_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/201220Van/results/bigwigs/hu50_H3K4me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/unt_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210219Van/results/bigwigs/hu50_H3K27me3_A_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K27ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9me2_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9me3_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__H3K9ac_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/Unt_A_avg.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210305Van/bd186/code/results/norm_bigwigs/HU_A_avg.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__G9a_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__G9a_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__Unt__WIZ_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/bd523_bd524/results/bigwigs/clone32__HU__WIZ_avg_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/BCL11A_0h_no_rmdup_norm.bw \
/net/bmc-pub17/mirror/lab/vanderheiden/briando/files/210723Van/results/bigwigs/BCL11A_24h_no_rmdup_norm.bw \
-R wiz_peaks_without_g9a.bed wiz_peaks_with_g9a.bed g9a_peaks_without_wiz.bed \
--referencePoint center --missingDataAsZero --skipZeros -a 2000 -b 2000 -bs 50 -p 15 -o 250918_wiz_g9a_cobinding.gz"

sbatch -J map -N 1 --wrap "plotProfile -m 250918_wiz_g9a_cobinding.gz -o 250918_wiz_g9a_cobinding.pdf \
--averageType mean --colors red black blue \
--samplesLabel \"H3K4me1\" \"H3K4me1 HU\" \"H3K4me3\" \"H3K4me3 HU\" \
\"H3K27me3\" \"H3K27me3 HU\" \
\"H3K27ac\" \"H3K27ac HU\" \"H3K9me2\" \"H3K9me2 HU\" \
\"H3K9me3\" \"H3K9me3 HU\" \"H3K9ac\" \"H3K9ac HU\" \
\"ATAC\" \"ATAC HU\" \"G9a\" \"G9a HU\" \
\"WIZ\" \"WIZ HU\" \"BCL11a\" \"BCL11a HU\" \
--regionsLabel \"WIZ peaks w/o G9a\" \"WIZ-G9a cobinding\" \"G9a peaks without WIZ\" \
--yMin 0 --yMax 0.25 0.25 0.6 0.6 0.05 0.05 \
0.8 0.8 0.2 0.2 0.1 0.1 0.4 0.4 140 140 0.8 0.8 1.5 1.5 0.6 0.6 \
--refPointLabel Center --numPlotsPerRow 4"

sbatch -J map -N 1 --wrap "plotHeatmap -m 250918_wiz_g9a_cobinding.gz -o 250918_wiz_g9a_cobinding_heatmap.pdf \
--averageType mean --colorMap Greys Greys Purples Purples \
Blues Blues Greens Greens Oranges Oranges Reds Reds Purples Purples \
Blues Blues Greens Greens Oranges Oranges Reds Reds \
--samplesLabel \"H3K4me1\" \"H3K4me1 HU\" \"H3K4me3\" \"H3K4me3 HU\" \
\"H3K27me3\" \"H3K27me3 HU\" \
\"H3K27ac\" \"H3K27ac HU\" \"H3K9me2\" \"H3K9me2 HU\" \
\"H3K9me3\" \"H3K9me3 HU\" \"H3K9ac\" \"H3K9ac HU\" \
\"ATAC\" \"ATAC HU\" \"G9a\" \"G9a HU\" \
\"WIZ\" \"WIZ HU\" \"BCL11a\" \"BCL11a HU\" \
--regionsLabel \"WIZ peaks w/o G9a\" \"WIZ-G9a cobinding\" \"G9a peaks without WIZ\" \
--yMin 0 --yMax 0.25 0.25 0.6 0.6 0.05 0.05 \
0.8 0.8 0.1 0.1 0.1 0.1 0.4 0.4 140 140 0.8 0.8 1.5 1.5 0.6 0.6 \
--zMax 0.25 0.25 0.6 0.6 0.05 0.05 \
0.8 0.8 0.1 0.1 0.1 0.1 0.4 0.4 140 140 0.8 0.8 1.5 1.5 0.6 0.6 \
--refPointLabel Center"

sbatch -J map -N 1 --wrap "plotProfile -m 250918_k9me2_chg_withhu_atacpeaks.gz -o 250918_k9me2_chg_withhu_atacpeaks.pdf \
--averageType mean --colors red black blue \
--samplesLabel \"H3K4me1\" \"H3K4me1 HU\" \
\"H3K4me3\" \"H3K4me3 HU\" \"H3K27me3\" \"H3K27me3 HU\" \
\"H3K27ac\" \"H3K27ac HU\" \"H3K9me2\" \"H3K9me2 HU\" \
\"H3K9me3\" \"H3K9me3 HU\" \"H3K9ac\" \"H3K9ac HU\" \
\"ATAC\" \"ATAC HU\" \"G9a\" \"G9a HU\" \
\"WIZ\" \"WIZ HU\" \"BCL11a\" \"BCL11a HU\" \
--regionsLabel \"Down\" \"Unchanged\" \"Up\" \
--yMin 0 --yMax 0.25 0.25 0.1 0.1 0.1 0.1 \
0.5 0.5 0.5 0.5 0.3 0.3 0.2 0.2 25 25 0.3 0.3 0.3 0.3 0.2 0.2 \
--refPointLabel Center --numPlotsPerRow 4"


sbatch -J map -N 1 --wrap "plotHeatmap -m 250918_k9me2_chg_withhu_atacpeaks.gz -o 250918_k9me2_chg_withhu_atacpeaks_heatmap50.pdf \
--averageType mean --colorMap Greys Greys Purples Purples \
Blues Blues Greens Greens Oranges Oranges Reds Reds Purples Purples \
Blues Blues Greens Greens Oranges Oranges Reds Reds \
--samplesLabel \"H3K4me1\" \"H3K4me1 HU\" \
\"H3K4me3\" \"H3K4me3 HU\" \"H3K27me3\" \"H3K27me3 HU\" \
\"H3K27ac\" \"H3K27ac HU\" \"H3K9me2\" \"H3K9me2 HU\" \
\"H3K9me3\" \"H3K9me3 HU\" \"H3K9ac\" \"H3K9ac HU\" \
\"ATAC\" \"ATAC HU\" \"G9a\" \"G9a HU\" \
\"WIZ\" \"WIZ HU\" \"BCL11a\" \"BCL11a HU\" \
--regionsLabel \"Down\" \"Unchanged\" \"Up\" \
--yMin 0 --yMax 0.25 0.25 0.1 0.1 0.1 0.1 \
0.5 0.5 0.5 0.5 0.3 0.3 0.2 0.2 25 25 0.3 0.3 0.3 0.3 0.2 0.2 \
--zMin 0 --zMax 0.25 0.25 0.1 0.1 0.1 0.1 \
0.5 0.5 0.5 0.5 0.3 0.3 0.2 0.2 25 25 0.3 0.3 0.3 0.3 0.2 0.2 \
--refPointLabel 0 --heatmapHeight 50"


sbatch -J map -N 1 --wrap "plotProfile -m 250918_k27me3_chg_withhu_atacpeaks.gz -o 250918_k27me3_chg_withhu_atacpeaks.pdf \
--averageType mean --colors red black blue \
--samplesLabel \"H3K4me1\" \"H3K4me1 HU\" \
\"H3K4me3\" \"H3K4me3 HU\" \"H3K27me3\" \"H3K27me3 HU\" \
\"H3K27ac\" \"H3K27ac HU\" \"H3K9me2\" \"H3K9me2 HU\" \
\"H3K9me3\" \"H3K9me3 HU\" \"H3K9ac\" \"H3K9ac HU\" \
\"ATAC\" \"ATAC HU\" \"G9a\" \"G9a HU\" \
\"WIZ\" \"WIZ HU\" \"BCL11a\" \"BCL11a HU\" \
--regionsLabel \"Down\" \"Unchanged\" \"Up\" \
--yMin 0 --yMax 0.1 0.1 0.1 0.1 0.15 0.15 \
0.1 0.1 0.1 0.1 0.05 0.05 0.15 0.15 20 20 \
0.2 0.2 0.3 0.30.1 0.1 \
--refPointLabel Center --numPlotsPerRow 4"


sbatch -J map -N 1 --wrap "plotProfile -m 250918_k9me2_chg.gz -o 250809_k9me2_chg.pdf \
--averageType mean --plotType se --colors red black blue \
--samplesLabel \"H3K4me1\" \"H3K4me3\" \"H3K27me3\" \"H3K27ac\" \"H3K9me2\" \"H3K9me3\" \"H3K9ac\" \"ATAC\" \
--regionsLabel \"Down\" \"Unchanged\" \"Up\" \
--yMin 0 --yMax 0.05 0.05 0.1 0.15 0.5 0.3 0.2 3 --startLabel Start --endLabel End --numPlotsPerRow 4"

sbatch -J map -N 1 --wrap "plotProfile -m 250809_k27me3_chg.gz -o 250809_k27me3_chg.pdf \
--averageType mean --plotType se --colors red black blue \
--samplesLabel \"H3K4me1\" \"H3K4me3\" \"H3K27me3\" \"H3K27ac\" \"H3K9me2\" \"H3K9me3\" \"H3K9ac\" \"ATAC\" \
--regionsLabel \"Down\" \"Unchanged\" \"Up\" \
--yMin 0 --yMax 0.05 0.05 0.1 0.05 0.1 0.1 0.2 5 --startLabel Start --endLabel End --numPlotsPerRow 4"

sbatch -J map -N 1 --wrap "plotProfile -m 250811_rpa_at_k9me3_10kb.gz -o 250811_rpa_at_k9me3_10kb.pdf \
--averageType mean --colors red orange green blue black \
--samplesLabel \"0h\" \"6h\" \"12h\" \"24h\" \"48h\" \
--yMin 0.01 --yMax 0.02 --startLabel Start --endLabel End --perGroup"

# can also do computeRegions for ATAC peaks within this set of regions

cd /home/briando/data/files/210305Van/bd186/code/results/macs2

# take the union of all peaks and then merges them
cat Unt*.narrowPeak | sort -k1,1 -k2,2n | cut -f1-3 | bedtools merge -i - > unt_atac_peaks_merged.bed  # 130973 peaks
cp /home/briando/data/files/210305Van/bd186/code/results/macs2/unt_atac_peaks_merged.bed /home/briando/data/files/bd523_bd524/results/macs2


# now intersect everything with it
parallel -j 12 "bedtools intersect -a unt_atac_peaks_merged.bed -b {} -wa > {.}_atac.bed" ::: 250809*.bed

# Plot profiles
sbatch -J map -N 1 --wrap "plotProfile -m 250917_g9a_cr.gz -o 250917_g9a_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--yMin 0 --yMax 0.3 0.8 0.2 --perGroup \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \
--refPointLabel Center"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_wiz_cr.gz -o 250917_wiz_cr1.pdf \
--averageType mean --plotType se --colors red red red black black black \
--yMin 0 --yMax 1 --perGroup \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \
--refPointLabel Center"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_h3k27ac_cr.gz -o 250917_h3k27ac_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \"H3K9ac peaks\" \
--refPointLabel Center --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_h3k9ac_cr.gz -o 250917_h3k9ac_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \"H3K9ac peaks\" \
--refPointLabel Center --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_h3k9me2_cr.gz -o 250917_h3k9me2_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \"H3K9ac peaks\" \
--refPointLabel Center --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_h3k9me3_cr.gz -o 250917_h3k9me3_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \"H3K9ac peaks\" \
--refPointLabel Center --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_cloneA3_h3k9me2_cr.gz -o 250917_cloneA3_h3k9me2_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"WIZ peaks\" \"G9a peaks\" \"H3K9ac peaks\" \
--refPointLabel Center --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_cloneA3_h3k9me2_repressed_cr.gz -o 250917_cloneA3_h3k9me2_repressed_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"H3K9me2 peaks\" \"H3K9me3 peaks\" \
--startLabel Start --endLabel End --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_h3k9me2_repressed_cr.gz -o 250917_h3k9me2_repressed_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"H3K9me2 peaks\" \"H3K9me3 peaks\" \
--startLabel Start --endLabel End --perGroup"

sbatch -J map -N 1 --wrap "plotProfile -m 250917_h3k9me3_repressed_cr.gz -o 250917_h3k9me3_repressed_cr.pdf \
--averageType mean --plotType se --colors red red red black black black \
--samplesLabel \"HU 1\" \"HU 2\" \"HU 3\" \"Unt 1\" \"Unt 2\" \"Unt 3\" \
--regionsLabel \"H3K9me2 peaks\" \"H3K9me3 peaks\" \
--startLabel Start --endLabel End --perGroup"

# Not all H3K9me2 peaks change... how do the ones near WIZ/G9a or associated with H3K9me3 or active vs repressed regions differ?
bedtools closest -a h3k9me2_peaks_merged.bed -b wiz_peaks_merged.bed -d -t first > h3k9me2_peaks_dist_wiz.bed
bedtools closest -a h3k9me2_peaks_merged.bed -b g9a_peaks_merged.bed -d -t first > h3k9me2_peaks_dist_g9a.bed
bedtools closest -a h3k9me2_peaks_merged.bed -b h3k9me3_peaks_merged.bed -d -t first > h3k9me2_peaks_dist_h3k9me3.bed

bedtools closest -a h3k9me3_peaks_merged.bed -b wiz_peaks_merged.bed -d -t first > h3k9me3_peaks_dist_wiz.bed
bedtools closest -a h3k9me3_peaks_merged.bed -b g9a_peaks_merged.bed -d -t first > h3k9me3_peaks_dist_g9a.bed
bedtools closest -a h3k9me3_peaks_merged.bed -b h3k9me2_peaks_merged.bed -d -t first > h3k9me3_peaks_dist_h3k9me2.bed

# HOMER
conda activate intel_env
perl /Users/briando/miniconda3/envs/intel_env/share/homer-4.10-0/.//configureHomer.pl -install hg38
cd ~/Dropbox\ \(MIT\)/MVH_Lab/Replication\ stress/Data/BD524/results

findMotifsGenome.pl 250809_g9a_all.bed hg38 homer_g9a -p 4
nohup findMotifsGenome.pl 250809_wiz_all.bed hg38 homer_wiz -p 4 &
nohup findMotifsGenome.pl 250809_k9me2_down.bed  hg38 homer_k9me2_down -bg 250809_k9me2_all.bed -p 4 &
nohup findMotifsGenome.pl 250809_k9me2_up.bed   hg38 homer_k9me2_up -bg 250809_k9me2_all.bed -p 4 &
nohup findMotifsGenome.pl 250809_k27me3_down.bed  hg38 homer_k27me3_down  -bg 250809_k27me3_all.bed -p 4 &
nohup findMotifsGenome.pl 250809_k27me3_up.bed  hg38  homer_k27me3_up  -bg 250809_k27me3_all.bed -p 4 &

nohup findMotifsGenome.pl 250809_g9a_all_atac.bed hg38 homer_g9a_atac -p 4 &
nohup findMotifsGenome.pl 250809_wiz_all_atac.bed hg38 homer_wiz_atac -p 4 &
nohup findMotifsGenome.pl 250809_k9me2_down_atac.bed  hg38 homer_k9me2_down_atac -bg 250809_k9me2_all_atac.bed -p 4 &
nohup findMotifsGenome.pl 250809_k9me2_down_atac.bed  hg38 homer_k9me2_down_atac_vs_genome -bg unt_atac_peaks_merged.bed -p 4 &
nohup findMotifsGenome.pl 250809_k9me2_up_atac.bed   hg38 homer_k9me2_up_atac -bg 250809_k9me2_all_atac.bed -p 4 &
nohup findMotifsGenome.pl 250809_k27me3_down_atac.bed  hg38 homer_k27me3_down_atac  -bg 250809_k27me3_all_atac.bed -p 4 &
nohup findMotifsGenome.pl 250809_k27me3_down_atac.bed  hg38 homer_k27me3_down_atac_vs_genome -bg unt_atac_peaks_merged.bed -p 4 &
nohup findMotifsGenome.pl 250809_k27me3_up_atac.bed  hg38  homer_k27me3_up_atac  -bg 250809_k27me3_all_atac.bed -p 4 &

nohup findMotifsGenome.pl 250809_k9me2_down_atac.bed  hg38  homer_k9me2down_vs_k27me3down_atac -bg 250809_k27me3_down_atac.bed -p 4 &
nohup findMotifsGenome.pl 250809_k27me3_down_atac.bed  hg38  homer_k27me3down_vs_k9me2down_atac -bg  250809_k9me2_down_atac.bed -p 4 &

nohup findMotifsGenome.pl BCL11A_0h_aligned_sorted_noblacklist_peaks.narrowPeak hg38 homer_bcl11a -p 4


#next steps
#Look at effect of H3K9me2 and H3K9me3 proximity to H3K27ac -- "active" vs repressed regions
#poised genes? -- can overlay H3K4me1/3

#do deepTools with mouse comparisons?
#Better H3K9me2 and H3K9me3 peak finding
#Which h3K9me2/3 sites change most and why?
#Motif analysis for WIZ and G9a -- which ones dont change? which ones are most sensitive?
#Think about H3K27me3 also
