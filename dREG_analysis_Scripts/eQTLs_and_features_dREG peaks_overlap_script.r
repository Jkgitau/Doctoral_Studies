library(regioneR)

# Load genome sizes
32 read.table("rice.genome", header=FALSE, sep="\t",
                          col.names=c("chrom", "size"))
rice.genome <- toGRanges(rice.genome)


# Load eQTLs and features (dREG peaks)
eqtls <- toGRanges("eQTLs.bed")
features <- toGRanges("features.bed")


# Run the Correlation between the two files
pt <- overlapPermTest(A=eqtls,
                      B=features,
                      ntimes=1000,
                      genome=rice.genome,             # rice genome
                      alternative="greater",
                      evaluate.function=numOverlaps,
                      randomize.function=randomizeRegions,
                      allow.overlaps=FALSE,           # avoid overlaps
                      per.chromosome=TRUE,            # keep within chromosome
                      universe=intergenic)            # restrict to intergenic

# Plot results
print(pt)
plot(pt)
