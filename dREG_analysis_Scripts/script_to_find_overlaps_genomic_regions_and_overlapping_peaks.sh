
#!/bin/bash


*****************************************************************************************


# Part 1: Finding genic - genic pore-C contacts with overlapping dREG peaks (Control and drought) 

# Input files:
# 1104 final_common_peaks_CM.bed
# 12875 loops_gene_gene.bed


## This section is for analysis of gene-no_gene loops
12875 loops_gene_gene.bed # The number of contacts between gene-gene location 


# Step 1: Extract the relevant columns from the loop files shared by Marina
cut -f1-3 loops_gene_gene.bed > loops_anchor1.bed -Anchor one
cut -f1,7,8 loops_gene_gene.bed > loops_anchor2.bed - Anchor two

# Make sure the number of contacts is equal to total contacts (in the original file),i.e, in anchor 1, anchor 2 and the main file, as # the only thing changing is the arrangements - contacts separated into different files based on their anchors (1,2)


## rename the chr column in the second anchor file so that it's consistent with start and stop naming convention:
sed '1s/^chrbin1/chrbin2/' loops_anchor2.bed > loops_anchor2_mod.bed # change the file name to loops_anchor2.bed


# Rename chromosome names to be consistent with the peak file shred by Marina
sed -e 's/^CM020633.1\t/chr1\t/' \
    -e 's/^CM020634.1\t/chr2\t/' \
    -e 's/^CM020635.1\t/chr3\t/' \
    -e 's/^CM020636.1\t/chr4\t/' \
    -e 's/^CM020637.1\t/chr5\t/' \
    -e 's/^CM020638.1\t/chr6\t/' \
    -e 's/^CM020639.1\t/chr7\t/' \
    -e 's/^CM020640.1\t/chr8\t/' \
    -e 's/^CM020641.1\t/chr9\t/' \
    -e 's/^CM020642.1\t/chr10\t/' \
    -e 's/^CM020643.1\t/chr11\t/' \
    -e 's/^CM020644.1\t/chr12\t/' \
    coordinates.txt > coordinates_CM.tsv

# Rename chromosome names to be consistent with the peak file shred by Marina
sed -e 's/^chr1\t/CM020633.1\t/' \
    -e 's/^chr2\t/CM020634.1\t/' \
    -e 's/^chr3\t/CM020635.1\t/' \
    -e 's/^chr4\t/CM020636.1\t/' \
    -e 's/^chr5\t/CM020637.1\t/' \
    -e 's/^chr6\t/CM020638.1\t/' \
    -e 's/^chr7\t/CM020639.1\t/' \
    -e 's/^chr8\t/CM020640.1\t/' \
    -e 's/^chr9\t/CM020641.1\t/' \
    -e 's/^chr10\t/CM020642.1\t/' \
    -e 's/^chr11\t/CM020643.1\t/' \
    -e 's/^chr12\t/CM020644.1\t/' \
    all_dREG_CM.tsv > all_dREG_CM_ordered.tsv


# Remove the first row from the two files (to be consistent with files having peaks in control and drought)

tail -n +2 loops_anchor1_CM.bed > loops_anchor1_CM.noheader.bed # This is now good to go
tail -n +2 loops_anchor2_CM.bed > loops_anchor2_CM.noheader.bed


# Combined the two anchor files, and sort them accordingly 
cat loops_anchor1_CM.noheader.bed loops_anchor2_CM.noheader.bed | sort -k1,1 -k2,2n > genic_genic_CM_regions.bed


# genic_genic_CM_regions.bed - This should be the sum of anchor 1 and 2 - wc-l to confirm


# Find the overlaps between the final genic file and final peak files
bedtools intersect -a final_common_peaks_CM.bed -b genic_genic_CM_regions.bed -wo | wc -l

# 408 contacts - genic_genic_contacts_with_overlapping_peaks_a_b_option.bed
#  408 contacts - genic_genic_contacts_with_overlapping_peaks-wo_option.bed


# Check the bedtools intersect for explanation of parameters used: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html




***************************************************************************************************************



# Part 2: Finding genic - intergenic pore-C contacts with overlapping dREG peaks (Control and drought) 
  
# Input files:

# 1104 final_common_peaks_CM.bed
# 14199 loops_gene_nogene.bed



# Step 1: Extract the relevant columns from the loop files shared by Marina
cut -f1-3 loops_gene_nogene.bed > loops_anchor1.bed #-Anchor one
cut -f1,7,8 loops_gene_nogene.bed > loops_anchor2.bed #- Anchor two

# Make sure the number of contacts is equal to total contacts (in the original file),i.e, in anchor 1, anchor 2 and the main file, as # the only thing changing is the arrangements - contacts separated into different files based on their anchors (1,2)


## rename the chr column in the second anchor file so that it's consistent with start and stop naming convention:
sed '1s/^chrbin1/chrbin2/' loops_anchor2.bed > loops_anchor2_mod.bed # change the file name to loops_anchor2.bed



# Rename chromosome names to be consistent with the peak file shred by Marina
sed -e 's/^chr.1\t/CM020633.1\t/' \
    -e 's/^chr.2\t/CM020634.1\t/' \
    -e 's/^chr.3\t/CM020635.1\t/' \
    -e 's/^chr.4\t/CM020636.1\t/' \
    -e 's/^chr.5\t/CM020637.1\t/' \
    -e 's/^chr.6\t/CM020638.1\t/' \
    -e 's/^chr.7\t/CM020639.1\t/' \
    -e 's/^chr.8\t/CM020640.1\t/' \
    -e 's/^chr.9\t/CM020641.1\t/' \
    -e 's/^chr.10\t/CM020642.1\t/' \
    -e 's/^chr.11\t/CM020643.1\t/' \
    -e 's/^chr.12\t/CM020644.1\t/' \
    loops_anchor1.bed > loops_anchor1_CM.bed

# Rename chromosome names to be consistent with the peak file shred by Marina
sed -e 's/^chr.1\t/CM020633.1\t/' \
    -e 's/^chr.2\t/CM020634.1\t/' \
    -e 's/^chr.3\t/CM020635.1\t/' \
    -e 's/^chr.4\t/CM020636.1\t/' \
    -e 's/^chr.5\t/CM020637.1\t/' \
    -e 's/^chr.6\t/CM020638.1\t/' \
    -e 's/^chr.7\t/CM020639.1\t/' \
    -e 's/^chr.8\t/CM020640.1\t/' \
    -e 's/^chr.9\t/CM020641.1\t/' \
    -e 's/^chr.10\t/CM020642.1\t/' \
    -e 's/^chr.11\t/CM020643.1\t/' \
    -e 's/^chr.12\t/CM020644.1\t/' \
    loops_anchor2.bed > loops_anchor2_CM.bed


# Remove the first row from the two files (to be consistent with files having peaks in control and drought)

tail -n +2 loops_anchor1_CM.bed > loops_anchor1_CM.noheader.bed # This is now good to go
tail -n +2 loops_anchor2_CM.bed > loops_anchor2_CM.noheader.bed


# Combined the two anchor files, and sort them accordingly 
cat loops_anchor1_CM.noheader.bed loops_anchor2_CM.noheader.bed | sort -k1,1 -k2,2n > genic_intergenic_CM_regions.bed


# genic_intergenic_CM_regions.bed - This should be the sum of anchor 1 and 2 - wc-l to confirm


# Find the overlaps between the final genic file and final peak files
bedtools intersect -a final_common_peaks_CM.bed -b genic_intergenic_CM_regions.bed -wo | wc -l
bedtools intersect -a final_common_peaks_CM.bed -b genic_intergenic_CM_regions.bed | wc -l

# 256 genic_intergenic_contacts_with_overlapping_peaks_a_b_option.bed
# 256 genic_intergenic_contacts_with_overlapping_peaks_wo_option.bed






***********************************************************************************************************************************
# Part 3: Finding intergenic - intergenic pore-C contacts with overlapping dREG peaks (Control and drought) 
  
# Input files:

# 1104 final_common_peaks_CM.bed
# 6708 loops_nogene_nogene.bed



# Step 1: Extract the relevant columns from the loop files shared by Marina
cut -f1-3 loops_nogene_nogene.bed > loops_anchor1.bed #-Anchor one
cut -f1,7,8 loops_nogene_nogene.bed > loops_anchor2.bed #- Anchor two

# Make sure the number of contacts is equal to total contacts (in the original file),i.e, in anchor 1, anchor 2 and the main file, as # the only thing changing is the arrangements - contacts separated into different files based on their anchors (1,2)


## rename the chr column in the second anchor file so that it's consistent with start and stop naming convention:
sed '1s/^chrbin1/chrbin2/' loops_anchor2.bed > loops_anchor2_mod.bed # change the file name to loops_anchor2.bed



# Rename chromosome names to be consistent with the peak file shred by Marina
sed -e 's/^chr.1\t/CM020633.1\t/' \
    -e 's/^chr.2\t/CM020634.1\t/' \
    -e 's/^chr.3\t/CM020635.1\t/' \
    -e 's/^chr.4\t/CM020636.1\t/' \
    -e 's/^chr.5\t/CM020637.1\t/' \
    -e 's/^chr.6\t/CM020638.1\t/' \
    -e 's/^chr.7\t/CM020639.1\t/' \
    -e 's/^chr.8\t/CM020640.1\t/' \
    -e 's/^chr.9\t/CM020641.1\t/' \
    -e 's/^chr.10\t/CM020642.1\t/' \
    -e 's/^chr.11\t/CM020643.1\t/' \
    -e 's/^chr.12\t/CM020644.1\t/' \
    loops_anchor1.bed > loops_anchor1_CM.bed

# Rename chromosome names to be consistent with the peak file shred by Marina
sed -e 's/^chr.1\t/CM020633.1\t/' \
    -e 's/^chr.2\t/CM020634.1\t/' \
    -e 's/^chr.3\t/CM020635.1\t/' \
    -e 's/^chr.4\t/CM020636.1\t/' \
    -e 's/^chr.5\t/CM020637.1\t/' \
    -e 's/^chr.6\t/CM020638.1\t/' \
    -e 's/^chr.7\t/CM020639.1\t/' \
    -e 's/^chr.8\t/CM020640.1\t/' \
    -e 's/^chr.9\t/CM020641.1\t/' \
    -e 's/^chr.10\t/CM020642.1\t/' \
    -e 's/^chr.11\t/CM020643.1\t/' \
    -e 's/^chr.12\t/CM020644.1\t/' \
    loops_anchor2.bed > loops_anchor2_CM.bed


# Remove the first row from the two files (to be consistent with files having peaks in control and drought)

tail -n +2 loops_anchor1_CM.bed > loops_anchor1_CM.noheader.bed # This is now good to go
tail -n +2 loops_anchor2_CM.bed > loops_anchor2_CM.noheader.bed


# Combined the two anchor files, and sort them accordingly 
cat loops_anchor1_CM.noheader.bed loops_anchor2_CM.noheader.bed | sort -k1,1 -k2,2n > intergenic_intergenic_CM_regions.bed


# intergenic_intergenic_CM_regions.bed - This should be the sum of anchor 1 and 2 - wc-l to confirm


# Find the overlaps between the final genic file and final peak files
bedtools intersect -a final_common_peaks_CM.bed -b intergenic_intergenic_CM_regions.bed -wo | wc -l
bedtools intersect -a final_common_peaks_CM.bed -b intergenic_intergenic_CM_regions.bed | wc -l

# 68 genic_intergenic_contacts_with_overlapping_peaks_a_b_option.bed
# 68 genic_intergenic_contacts_with_overlapping_peaks_wo_option.bed



