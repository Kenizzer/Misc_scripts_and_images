# Making nested relative abundance barplots with Fantaxtic
# Code by Ceyda Kural and Joel Swift


# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(phyloseq); packageVersion('phyloseq')
library(fantaxtic); packageVersion("fantaxtic")


#### Load datasets for use throughout ####
# Non-CLR transformed here since we are making relative abundance barplots
drt.bact.late <- readRDS('phyloseq_b_asv_clean_inoculated_late.RDS')

# Make our nested dataframe for plotting
# Here you can set the top and bottom levels of the nesting along with the #s of each
# So, we are taking the top 6 phyla and from each the top 4 classes.
# e.g. Proteobacteria | Alpha, Delta, Gamma, Beta
# For cases where only 1 class for a phyla exists nested_merged_label argument will capture this
top_nested <- nested_top_taxa(drt.bact.late,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Class",
                              n_top_taxa = 6, 
                              n_nested_taxa = 4, 
                              nested_merged_label = "<tax>")

# plot with a facet by Drought (D) or Well-Watered (W) plants 
P1 <- plot_nested_bar(top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Class",
                nested_merged_label = "<tax>",
                legend_title = "Phylum and Class") +
  labs(y = "Relative Abuance") +
  facet_wrap(~Drought.or.Watered, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        legend.key.size = unit(10, "points")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')

ggsave("Example_plot1.svg", P1, height = 6, width = 10)

# What if we want to sort the samples by the abundance of a particular taxonomic rank?
# All we need to do is create a vector that contains the sample names in the order 
# we wish to plot them.

# sort by a particular phyla
sample_order.vec <- psmelt(top_nested$ps_obj) %>%
  data.frame() %>%
  
  # Calculate relative abundances
  group_by(Sample) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  
  # Sort by taxon of interest
  filter(Phylum == "Actinobacteria") %>%
  group_by(Sample) %>%
  summarise(Abundance = sum(Abundance)) %>%
  arrange(Abundance) %>% 
  
  # Extract the sample order
  pull(Sample) %>%
  as.character()

# Add the sample_order.vec to the sample_order argument in plot_nested_bar
P2 <- plot_nested_bar(top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Class",
                nested_merged_label = "<tax>",
                legend_title = "Phylum and Class",
                sample_order = rev(sample_order.vec)) +
  labs(y = "Relative Abuance") +
  facet_wrap(~Drought.or.Watered, scales = "free_x") +
  theme(plot.title = element_text(hjust = 0.5,  size = 8,  face = "bold"),
        legend.key.size = unit(10, "points")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')
ggsave("Example_plot2.svg", P2, height = 6, width = 10)


# Want to make it a single barplot for each facet?
# Give sample_order = character(0)
P3 <- plot_nested_bar(top_nested$ps_obj,
                top_level = "Phylum",
                nested_level = "Class",
                nested_merged_label = "<tax>",
                legend_title = "Phylum and Class",
                sample_order = character(0)) +
  labs(y = "Relative Abuance") +
  facet_wrap(~Drought.or.Watered) +
  theme(plot.title = element_text(hjust = 0.5, size = 8,  face = "bold"),
        legend.key.size = unit(10, "points")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'right')
ggsave("Example_plot3.svg", P3, height = 6, width = 10)