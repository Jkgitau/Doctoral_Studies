
# Working directory: "C:/Users/hp/Desktop/marina/quick" - This is my HP computer

# A script to arrange chromosome naming convention and order by co-ordinates --------


### Read in the dREG peaks overlapping pore-C (before processing)

Overlapping_all <- read.table('dREG_signals_overlapping_intergenic_Pore_contacts.tsv', header = TRUE, sep = '\t')


# Define the correct chromosome order
chrom_order <- paste0("chr", 1:12)

# Convert chr column to a factor with the correct order
all_dREG$chr <- factor(all_dREG$chr, levels = chrom_order, ordered = TRUE)

# Writye the table tp file
write.table(
  Overlapping,
  file = 'OverlappingdREGs_wit_Pore-c_sorted_by_chrm_and_coordinated_CM.tsv'
)


# Now extract all the overlapping dREGs in the all dREG files, and get theor proximal or distal status


# Read in a file with distal dREG peaks
distal_dREGs <- read.table('distal_dREGs.txt', header = TRUE)

# Read in the file with ALL dREG peaks overlapping with non-genic part of gene-no_gene
Overlapping <- read.table('OverlappingdREGs_wit_Pore-c_sorted_by_chrm_and_coordinated_CM.tsv', header = TRUE)


# Ensure DE_intergenic_proximity_status_calc is the last column
merged_df <- merged_df[, c(names(Overlapping), "DE_intergenic_proximity_status_calc")]

#subset dREG peaks that are also distal
common <- semi_join(
  Overlapping_dREG,
  distal_dREGs_final,
  by = c("Distal_dREG", "Start", "End")
)




