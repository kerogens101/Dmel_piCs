library(rstudioapi)

# Function to read tree data and perform labeling
read_and_label_tree <- function(tree_file, label_file) {
  tree <- read.tree(tree_file)
  treelabs <- read.table(label_file, header = FALSE, sep = "\t")
  names(treelabs) <- c("group", "label")
  tree_labelled <- full_join(as_tibble(tree), treelabs, by = c('label'))
  tree_data <- as.treedata(tree_labelled)
  return(list(tree = tree_data, label = treelabs))
}

# Set working directory
cur_dir <- dirname(getSourceEditorContext()$path)
setwd(cur_dir)

# Process ISO1_RMv8blood_tree
blood_tree <- read_and_label_tree("ISO1_genome_RMv8_classif_flank50bp_blood_nonpiC-piC_compiled_mafft_FastTree.nwk",
                                  "ISO1_genome_RMv8_classif_flank50bp_blood_nonpiC-piC_compiled_mafft_FastTreeNames.txt")
ggtree(blood_tree$tree) + geom_tippoint(aes(color = group))

# Process ISO1_RMv8blastopia_tree
blastopia_tree <- read_and_label_tree("ISO1_RMv8_classif_flank50bp_Blastopia-like_nonpiC-piC_compiled_mafft_FastTreeV1.nwk",
                                      "ISO1_RMv8Blastopia_tree_data.txt")
ggtree(blastopia_tree$tree) + geom_tippoint(aes(color = group))

# Process ISO1_RMv8DM_412_tree
dm_412_tree <- read_and_label_tree("/local/workdir/sps257/dspr_piRNA/TEfam_CN_trees/DM_412_tree_files/cdhit_collapse/DSPR_RMv8_Ty1Copia_classif_flank50bp_DM412_NONpiC-piC_cdhit0.64kb_compiled_mafEINSi_FastTree.nwk",
                                   "/local/workdir/sps257/dspr_piRNA/TEfam_CN_trees/DM_412_tree_files/cdhit_collapse/DSPR_RMv8_Ty1Copia_classif_flank50bp_DM412_NONpiC-piC_cdhit0.64kb_compiled_mafEINSi_FastTree_names.txt")
ggtree(dm_412_tree$tree) + geom_tippoint(aes(color = group))

# Plotting the first ggplot
ggplot(data = ISO1_genome_RMv8_classif_flank50bp_blood_nonpiC_piC_compiled_mafft_FastTree.nwk.txt,
       aes(x = V1, y = V3, fill = V1)) +
  geom_violin(position = dodge) +
  geom_boxplot(position = dodge, width = 0.2, alpha = 0.0) +
  scale_y_continuous(limits = c(0, 0.12), breaks = c(0, 0.05, 0.1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_rect(), legend.position = "none",
        axis.title = element_text(size = 13, color = "black", face = "bold")) +
  xlab("blood") + ylab("terminal branch length")

# Plotting the second ggplot
ggplot(data = ISO1_RMv8_Tc1Mariner_classif_unflank50bp_cdhit_repComb_nonpiC_piC_compiled_mafft_FastTree.nwk.txt,
       aes(x = V1, y = V3, fill = V1)) +
  geom_violin(position = dodge, alpha = 0.5, lwd = 0.6, width = 0.9) +
  geom_boxplot(position = dodge, width = 0.2, alpha = 0.0) +
  scale_y_continuous(limits = c(0, 0.21), breaks = c(0, 0.1, 0.2)) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_rect(), legend.position = "none",
        axis.title = element_text(size = 13, color = "black", face = "bold")) +
  xlab("Ty1/Copia") + ylab("terminal branch length") + scale_x_discrete(limits = c("piC", "nonpiC"))
