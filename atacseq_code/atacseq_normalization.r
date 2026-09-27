# atacseq_normalization.r
# Usage: Rscript atacseq_normalization.r <featureCounts output> <MDSplot.pdf> <norm_factors.txt>

library(edgeR)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
	print("Usage: Rscript atacseq_normalization.r <featureCounts output> <MDSplot.pdf> <norm_factors.txt>")
	quit(save="no", 1)
}

# Read in table
df = read.table(args[1], header = TRUE, row.names=1)
df = df[,6:dim(df)[2]]
colnames(df) = gsub("results.bowtie2_bam.", "", 
                    gsub("_aligned_reads_sorted_nodup_noblacklist_noM.bam", "", colnames(df)))

# edgeR to calculate norm factors and get peaks
d0 <- DGEList(df)
d0 <- calcNormFactors(d0, method="TMM")
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

# plot 
pdf(args[2], width=8, height=8)
plotMDS(d)
dev.off()

# norm factors
d0$samples$div.factors = d0$samples$lib.size * d0$samples$norm.factors / 1e7
write.table(d0$samples, args[3], quote=F, sep='\t')

quit(save="no", 0)