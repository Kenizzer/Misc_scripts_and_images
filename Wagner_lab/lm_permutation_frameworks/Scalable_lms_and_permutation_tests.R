# The code here provides scalable frameworks to do linear modeling and permutation
# testing across any taxonomic rank, all the way to individual ASVs! The frameworks
# make heavy use of phyloseq and its functions. 

# The code is still pretty rough but my hope is someone in the future finds it 
# useful... even if that someone is me :)

# Date: 11-06-2022
# Code by: Joel Swift

# Packages w/ version numbers.
library(tidyverse); packageVersion('tidyverse')
library(phyloseq); packageVersion('phyloseq')
library(lme4); packageVersion('lme4')
library(lmerTest); packageVersion('lmerTest')
library(emmeans); packageVersion('emmeans')
library(ggpubr); packageVersion('ggpubr')
# Theme set 
theme_set(theme_pubr())

#### Load dataset for use throughout ####
# This was made by taking the 101 taxa from a phyloseq data frame and
# running data.frame(matrix(rnorm(19897, 100, 25), nrow = 101)) to generate a 
# new otu matrix.
Phyloseq_example <- readRDS('Phyloseq_example.RDS')
# I replaced 5 taxa abundances to show strong differences across the treatment
Phyloseq_example@otu_table[1,] <- c(76.999995, 7.066666, 42.14286, 62.33334, 14.95238, 10.488888, 22.85715, 6.857142, 19.466666, 11.35, 37.5, 21.380952, 8.714286, 18, 52.5, 72, 53.108115, 24.648648, 1.351352, 30.2, 71.70732, 5.272728, 17.314286, 80.131575, 72.162165, 10.578948, 48, 93.375, 88.536585, 68.108115, 95.921055, 100.65789, 14.714286, 45.375, 21.526316, 33.571425, 27.75, 62.25, 24.5, 18.243902, 26.133334, 32.714286, 23.55, 59.66667, 79.736835, 21.75, 15.25, 24.486486, 12.594594, 31.216215, 137.837835, 43.333335, 21.4, 9.5, 17.333334, 20.918918, 115.13514, 90.7317, 22.7, 19.8, 3.333334, 20.81081, 96.375, 104.21052, 16.108108, 6.756756, 75.88236, 18.1, 22.45, 25.027028, 21.263158, 22.04878, 21.351352, 165, 20.25, 30.888888, 39.375, 91.62162, 65.625, 10.108108, 10.432432, 15.028572, 26.685714, 18.864864, 91.973685, 22.85, 15.35, 92.625, 78.55263, 26.378378, 81.891885, 22.263158, 30.162162, 80.121945, 4.2, 15.65, 56.571435, 59.18919, 83.25, 20.97561, 70.000005, 14.60526, 45.35715, 31.72973, 19.422222, 81, 45, 76.875, 84.51219, 9.513514, 74.594595, 18.142858, 22.8, 29.405406, 83.25, 15.6, 26.378378, 118.125, 16.628572, 21.36363, 23.9, 70.13514, 50.66667, 107.43243, 77.142855, 76.216215, 29.297298, 29.894736, 98.51352, 24.3, 81.891885, 17.47619, 90.81081, 24.108108, 56.70732, 27.804878, 15.85, 27.675676, 24.052632, 49.125, 70.65789, 74.571435, 42.428565, 45.78948, 101.75676, 25.72973, 28.864864, 59.33334, 22.15, 42.499995, 62.43243, 64.459455, 45.66666, 58.815795, 22.216216, 24.756756, 14.162162, 18.809524, 27.2, 99.729735, 76.216215, 73.783785, 65.357145, 60.33333, 69.33333, 16.761904, 25.190476, 47.368425, 95.625, 49.125, 63.181815, 18.228572, 15.714286, 94.125, 79.875, 71.78571, 9.75, 23.333334, 30.75, 48.24324, 28.324324, 78.75, 85.13514, 18.216216, 15.466666, 29.72973, 104.594595, 30.39474, 10.190476, 71.78571, 19.095238, 96, 25.6, 32.864864, 16.619048, 13.380952, 17.428572)
Phyloseq_example@otu_table[4,] <- c(61.599996, 5.6533328, 33.714288, 49.866672, 11.961904, 8.3911104, 18.28572, 5.4857136, 15.5733328, 9.08, 30, 17.1047616, 6.9714288, 14.4, 42, 57.6, 42.486492, 19.7189184, 1.0810816, 24.16, 57.365856, 4.2181824, 13.8514288, 64.10526, 57.729732, 8.4631584, 38.4, 74.7, 70.829268, 54.486492, 76.736844, 80.526312, 11.7714288, 36.3, 17.2210528, 26.85714, 22.2, 49.8, 19.6, 14.5951216, 20.9066672, 26.1714288, 18.84, 47.733336, 63.789468, 17.4, 12.2, 19.5891888, 10.0756752, 24.972972, 110.270268, 34.666668, 17.12, 7.6, 13.8666672, 16.7351344, 92.108112, 72.58536, 18.16, 15.84, 2.6666672, 16.648648, 77.1, 83.368416, 12.8864864, 5.4054048, 60.705888, 14.48, 17.96, 20.0216224, 17.0105264, 17.639024, 17.0810816, 132, 16.2, 24.7111104, 31.5, 73.297296, 52.5, 8.0864864, 8.3459456, 12.0228576, 21.3485712, 15.0918912, 73.578948, 18.28, 12.28, 74.1, 62.842104, 21.1027024, 65.513508, 17.8105264, 24.1297296, 64.097556, 3.36, 12.52, 45.257148, 47.351352, 66.6, 16.780488, 56.000004, 11.684208, 36.28572, 25.383784, 15.5377776, 64.8, 36, 61.5, 67.609752, 7.6108112, 59.675676, 14.5142864, 18.24, 23.5243248, 66.6, 12.48, 21.1027024, 94.5, 13.3028576, 17.090904, 19.12, 56.108112, 40.533336, 85.945944, 61.714284, 60.972972, 23.4378384, 23.9157888, 78.810816, 19.44, 65.513508, 13.980952, 72.648648, 19.2864864, 45.365856, 22.2439024, 12.68, 22.1405408, 19.2421056, 39.3, 56.526312, 59.657148, 33.942852, 36.631584, 81.405408, 20.583784, 23.0918912, 47.466672, 17.72, 33.999996, 49.945944, 51.567564, 36.533328, 47.052636, 17.7729728, 19.8054048, 11.3297296, 15.0476192, 21.76, 79.783788, 60.972972, 59.027028, 52.285716, 48.266664, 55.466664, 13.4095232, 20.1523808, 37.89474, 76.5, 39.3, 50.545452, 14.5828576, 12.5714288, 75.3, 63.9, 57.428568, 7.8, 18.6666672, 24.6, 38.594592, 22.6594592, 63, 68.108112, 14.5729728, 12.3733328, 23.783784, 83.675676, 24.315792, 8.1523808, 57.428568, 15.2761904, 76.8, 20.48, 26.2918912, 13.2952384, 10.7047616, 13.9428576) 
Phyloseq_example@otu_table[12,] <- c(107.799993, 9.8933324, 59.000004, 87.266676, 20.933332, 14.6844432, 32.00001, 9.5999988, 27.2533324, 15.89, 52.5, 29.9333328, 12.2000004, 25.2, 73.5, 100.8, 74.351361, 34.5081072, 1.8918928, 42.28, 100.390248, 7.3818192, 24.2400004, 112.184205, 101.027031, 14.8105272, 67.2, 130.725, 123.951219, 95.351361, 134.289477, 140.921046, 20.6000004, 63.525, 30.1368424, 46.999995, 38.85, 87.15, 34.3, 25.5414628, 36.5866676, 45.8000004, 32.97, 83.533338, 111.631569, 30.45, 21.35, 34.2810804, 17.6324316, 43.702701, 192.972969, 60.666669, 29.96, 13.3, 24.2666676, 29.2864852, 161.189196, 127.02438, 31.78, 27.72, 4.6666676, 29.135134, 134.925, 145.894728, 22.5513512, 9.4594584, 106.235304, 25.34, 31.43, 35.0378392, 29.7684212, 30.868292, 29.8918928, 231, 28.35, 43.2444432, 55.125, 128.270268, 91.875, 14.1513512, 14.6054048, 21.0400008, 37.3599996, 26.4108096, 128.763159, 31.99, 21.49, 129.675, 109.973682, 36.9297292, 114.648639, 31.1684212, 42.2270268, 112.170723, 5.88, 21.91, 79.200009, 82.864866, 116.55, 29.365854, 98.000007, 20.447364, 63.50001, 44.421622, 27.1911108, 113.4, 63, 107.625, 118.317066, 13.3189196, 104.432433, 25.4000012, 31.92, 41.1675684, 116.55, 21.84, 36.9297292, 165.375, 23.2800008, 29.909082, 33.46, 98.189196, 70.933338, 150.405402, 107.999997, 106.702701, 41.0162172, 41.8526304, 137.918928, 34.02, 114.648639, 24.466666, 127.135134, 33.7513512, 79.390248, 38.9268292, 22.19, 38.7459464, 33.6736848, 68.775, 98.921046, 104.400009, 59.399991, 64.105272, 142.459464, 36.021622, 40.4108096, 83.066676, 31.01, 59.499993, 87.405402, 90.243237, 63.933324, 82.342113, 31.1027024, 34.6594584, 19.8270268, 26.3333336, 38.08, 139.621629, 106.702701, 103.297299, 91.500003, 84.466662, 97.066662, 23.4666656, 35.2666664, 66.315795, 133.875, 68.775, 88.454541, 25.5200008, 22.0000004, 131.775, 111.825, 100.499994, 13.65, 32.6666676, 43.05, 67.540536, 39.6540536, 110.25, 119.189196, 25.5027024, 21.6533324, 41.621622, 146.432433, 42.552636, 14.2666664, 100.499994, 26.7333332, 134.4, 35.84, 46.0108096, 23.2666672, 18.7333328, 24.4000008)
Phyloseq_example@otu_table[64,] <- c(123.199992, 11.3066656, 67.428576, 99.733344, 23.923808, 16.7822208, 36.57144, 10.9714272, 31.1466656, 18.16, 60, 34.2095232, 13.9428576, 28.8, 84, 115.2, 84.972984, 39.4378368, 2.1621632, 48.32, 114.731712, 8.4363648, 27.7028576, 128.21052, 115.459464, 16.9263168, 76.8, 149.4, 141.658536, 108.972984, 153.473688, 161.052624, 23.5428576, 72.6, 34.4421056, 53.71428, 44.4, 99.6, 39.2, 29.1902432, 41.8133344, 52.3428576, 37.68, 95.466672, 127.578936, 34.8, 24.4, 39.1783776, 20.1513504, 49.945944, 220.540536, 69.333336, 34.24, 15.2, 27.7333344, 33.4702688, 184.216224, 145.17072, 36.32, 31.68, 5.3333344, 33.297296, 154.2, 166.736832, 25.7729728, 10.8108096, 121.411776, 28.96, 35.92, 40.0432448, 34.0210528, 35.278048, 34.1621632, 264, 32.4, 49.4222208, 63, 146.594592, 105, 16.1729728, 16.6918912, 24.0457152, 42.6971424, 30.1837824, 147.157896, 36.56, 24.56, 148.2, 125.684208, 42.2054048, 131.027016, 35.6210528, 48.2594592, 128.195112, 6.72, 25.04, 90.514296, 94.702704, 133.2, 33.560976, 112.000008, 23.368416, 72.57144, 50.767568, 31.0755552, 129.6, 72, 123, 135.219504, 15.2216224, 119.351352, 29.0285728, 36.48, 47.0486496, 133.2, 24.96, 42.2054048, 189, 26.6057152, 34.181808, 38.24, 112.216224, 81.066672, 171.891888, 123.428568, 121.945944, 46.8756768, 47.8315776, 157.621632, 38.88, 131.027016, 27.961904, 145.297296, 38.5729728, 90.731712, 44.4878048, 25.36, 44.2810816, 38.4842112, 78.6, 113.052624, 119.314296, 67.885704, 73.263168, 162.810816, 41.167568, 46.1837824, 94.933344, 35.44, 67.999992, 99.891888, 103.135128, 73.066656, 94.105272, 35.5459456, 39.6108096, 22.6594592, 30.0952384, 43.52, 159.567576, 121.945944, 118.054056, 104.571432, 96.533328, 110.933328, 26.8190464, 40.3047616, 75.78948, 153, 78.6, 101.090904, 29.1657152, 25.1428576, 150.6, 127.8, 114.857136, 15.6, 37.3333344, 49.2, 77.189184, 45.3189184, 126, 136.216224, 29.1459456, 24.7466656, 47.567568, 167.351352, 48.631584, 16.3047616, 114.857136, 30.5523808, 153.6, 40.96, 52.5837824, 26.5904768, 21.4095232, 27.885715)
Phyloseq_example@otu_table[22,] <- c(192.4999875, 17.666665, 105.35715, 155.83335, 37.38095, 26.22222, 57.142875, 17.142855, 48.666665, 28.375, 93.75, 53.45238, 21.785715, 45, 131.25, 180, 132.7702875, 61.62162, 3.37838, 75.5, 179.2683, 13.18182, 43.285715, 200.3289375, 180.4054125, 26.44737, 120, 233.4375, 221.3414625, 170.2702875, 239.8026375, 251.644725, 36.785715, 113.4375, 53.81579, 83.9285625, 69.375, 155.625, 61.25, 45.609755, 65.333335, 81.785715, 58.875, 149.166675, 199.3420875, 54.375, 38.125, 61.216215, 31.486485, 78.0405375, 344.5945875, 108.3333375, 53.5, 23.75, 43.333335, 52.297295, 287.83785, 226.82925, 56.75, 49.5, 8.333335, 52.027025, 240.9375, 260.5263, 40.27027, 16.89189, 189.7059, 45.25, 56.125, 62.56757, 53.157895, 55.12195, 53.37838, 412.5, 50.625, 77.22222, 98.4375, 229.05405, 164.0625, 25.27027, 26.08108, 37.57143, 66.714285, 47.16216, 229.9342125, 57.125, 38.375, 231.5625, 196.381575, 65.945945, 204.7297125, 55.657895, 75.405405, 200.3048625, 10.5, 39.125, 141.4285875, 147.972975, 208.125, 52.439025, 175.0000125, 36.51315, 113.392875, 79.324325, 48.555555, 202.5, 112.5, 192.1875, 211.280475, 23.783785, 186.4864875, 45.357145, 57, 73.513515, 208.125, 39, 65.945945, 295.3125, 41.57143, 53.409075, 59.75, 175.33785, 126.666675, 268.581075, 192.8571375, 190.5405375, 73.243245, 74.73684, 246.2838, 60.75, 204.7297125, 43.690475, 227.027025, 60.27027, 141.7683, 69.512195, 39.625, 69.18919, 60.13158, 122.8125, 176.644725, 186.4285875, 106.0714125, 114.4737, 254.3919, 64.324325, 72.16216, 148.33335, 55.375, 106.2499875, 156.081075, 161.1486375, 114.16665, 147.0394875, 55.54054, 61.89189, 35.405405, 47.02381, 68, 249.3243375, 190.5405375, 184.4594625, 163.3928625, 150.833325, 173.333325, 41.90476, 62.97619, 118.4210625, 239.0625, 122.8125, 157.9545375, 45.57143, 39.285715, 235.3125, 199.6875, 179.464275, 24.375, 58.333335, 76.875, 120.6081, 70.81081, 196.875, 212.83785, 45.54054, 38.666665, 74.324325, 261.4864875, 75.98685, 25.47619, 179.464275, 47.738095, 240, 64, 82.16216, 41.54762, 33.45238, 43.57143)
# Extract the sample data from phyloseq, for whenever we to plot or model.
Phyloseq_sample_example <- data.frame(sample_data(Phyloseq_example))


#===========================================================#
#### Linear modeling example, scalable to each rank/taxa ####
#===========================================================#


##### 1) Run a model w/o the factor you want to test #####
# First we can isolate to the factor we are interested by pulling the residuals 
# from a model with the other factors controlled for. Here we are looking at 
# plant shoot mass rate (i.e., the grams of shoot mass / number of days of 
# growth). We are interested in whether a particular ASV abundance is associated
# positively or negatively with shoot mass rate and whether this is specific to a 
# particular treatment (i.e., drought or well-watered plants).

# THIS WILL BE SPECIFC TO WHATEVER YOU STUDY DESIGN IS (Modify as required/desired)
# Mixed model with genotype + soilInoculum + useable reads (log/Zscore)
mod <- lmerTest::lmer(sqrt(ShootMassRate) ~ Genotype + SoilInoculum + logObs.z +
                        (1|Block), data = Phyloseq_sample_example)
# Add residuals of the model for SMR to sample data and phyloseq obj
Phyloseq_sample_example$ShootMassRateResid <- resid(mod)
sample_data(Phyloseq_example)$ShootMassRateResid <- resid(mod)


##### 2) Prepare phyloseq dataframes by taxonomic rank #####
# We want to to test whether the abundance at a particular rank associates with
# Shoot mass rate, so we need to make subset our phyloseq obj by taxonomic rank
# and use psmelt to melt the phyloseq object into a dataframe.
Phyloseq_example_phy <- psmelt(tax_glom(Phyloseq_example, "Phylum")) # 5 taxa
Phyloseq_example_cla <- psmelt(tax_glom(Phyloseq_example, "Class")) # 10 taxa
Phyloseq_example_ord <- psmelt(tax_glom(Phyloseq_example, "Order")) # 18 taxa
Phyloseq_example_fam <- psmelt(tax_glom(Phyloseq_example, "Family")) # 40 taxa
Phyloseq_example_gen <- psmelt(tax_glom(Phyloseq_example, "Genus")) # 51 taxa
# We now have all of the dataframes per taxonomic rank in long form.


##### 3) Linear models run on each taxa by taxonomic rank #####
# Here I will break down one of the one liners that runs these models
Phyloseq_example_aovs_phy <- Phyloseq_example_phy %>%
  # make a tibble for each taxa by the given rank (phylum here)
  nest_by(Phylum) %>% 
  # create a tibble to hold anova results from a linear model by rank
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>% 
  # Use tidy to turn the anovas into tibbles and summarize to condense the 
  # tibbles into a single data frame grouping by phylum
  summarize(broom::tidy(mod))

# Class 
Phyloseq_example_aovs_cla <- Phyloseq_example_cla %>% nest_by(Class) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
# Order
Phyloseq_example_aovs_ord <- Phyloseq_example_ord %>% nest_by(Order) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
# Family
Phyloseq_example_aovs_fam <- Phyloseq_example_fam %>% nest_by(Family) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))
#Genus
Phyloseq_example_aovs_gen <- Phyloseq_example_gen %>% nest_by(Genus) %>%
  mutate(mod = list(anova(lm(ShootMassRateResid ~ Abundance*Drought.or.Watered, data=data)))) %>%
  summarize(broom::tidy(mod))


##### 4) Correct p-values and count the significant taxa by term #####
# Lets correct p-values and count the number of significant taxa within a rank
# by term. Our terms we care about are the ASV abundance and the drought treatment
# Phylum
sum(p.adjust(Phyloseq_example_aovs_phy[Phyloseq_example_aovs_phy$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(Phyloseq_example_aovs_phy[Phyloseq_example_aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)
### We can check which were significant with a couple lines
temp <- which(p.adjust(Phyloseq_example_aovs_phy[Phyloseq_example_aovs_phy$term == "Abundance", ]$p.value, method = "BH") < 0.05)
Phyloseq_example_aovs_phy[Phyloseq_example_aovs_phy$term == "Abundance", ][c(temp), ]

# Class
sum(p.adjust(Phyloseq_example_aovs_cla[Phyloseq_example_aovs_cla$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(Phyloseq_example_aovs_cla[Phyloseq_example_aovs_cla$term == "Abundance", ]$p.value, method = "BH") < 0.05)

# Order
sum(p.adjust(Phyloseq_example_aovs_ord[Phyloseq_example_aovs_ord$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(Phyloseq_example_aovs_ord[Phyloseq_example_aovs_ord$term == "Abundance", ]$p.value, method = "BH") < 0.05)

# Family
sum(p.adjust(Phyloseq_example_aovs_fam[Phyloseq_example_aovs_fam$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(Phyloseq_example_aovs_fam[Phyloseq_example_aovs_fam$term == "Abundance", ]$p.value, method = "BH") < 0.05)

# Genus
sum(p.adjust(Phyloseq_example_aovs_gen[Phyloseq_example_aovs_gen$term == "Abundance:Drought.or.Watered", ]$p.value, method = "BH") < 0.05)
sum(p.adjust(Phyloseq_example_aovs_gen[Phyloseq_example_aovs_gen$term == "Abundance", ]$p.value, method = "BH") < 0.05)


## Plotting function to improve
# only works for phylum at the moment.
plot_asvs <- function(DF, growth_measure = NULL, Phylum, Wet.or.Dry){
  if(is.null(growth_measure)){stop("Please specify a growth measure")
  } else if(growth_measure == "ShootMassRate"){
    x <- ggplot(DF[DF$Phylum  == Phylum & DF$Drought.or.Watered == Wet.or.Dry,],
                aes( x = Abundance, y = ShootMassRate)) +
      geom_point() + geom_smooth(method = 'lm') + ylab("Shoot Mass Rate (g/day)") + xlab("CLR Transformed Abundance") +
      #ggtitle(paste(ASV_number, sep = " "))
      ggtitle(paste("Phylum:", Phylum, sep = " "))
    return(x)

  }
}
plot_asvs(Phyloseq_example_phy, growth_measure = "ShootMassRate", Phylum ="Firmicutes", Wet.or.Dry = "D")
plot_asvs(Phyloseq_example_phy, growth_measure = "ShootMassRate", Phylum ="Other", Wet.or.Dry = "D")
plot_asvs(Phyloseq_example_phy, growth_measure = "ShootMassRate", Phylum ="Proteobacteria", Wet.or.Dry = "D")
plot_asvs(Phyloseq_example_phy, growth_measure = "ShootMassRate", Phylum ="Actinobacteria", Wet.or.Dry = "D")


##### Ideas for improvement #####
# + Function run the p-value corrections and pull the significant results
#     - Would need the rank_df
# + Function for quick plotting of the significant results 
#     - Would need the rank_df, Taxa #, and many a grouping factor.
###


#=======================================================#
#### Permutation example, scalable to each rank/taxa ####
#=======================================================#
# Warning this is still pretty experimental at the moment :)

# Lets cut the dataset down to only 10 taxa to save time
Phyloseq_example <- readRDS('Phyloseq_example.RDS')
temp <- (fantaxtic::top_taxa(Phyloseq_example, n = 10)) # Fantaxtic is an awesome package, check it out!
Phyloseq_example <- temp$ps_obj


##### 1) Create our observed value list #####
Phyloseq_example_long <- psmelt(Phyloseq_example)
# Borrowing some of the code from above to 
Phyloseq_example_aovs <- Phyloseq_example_long %>% nest_by(OTU, Drought.or.Watered) %>%
  mutate(mod = list(anova(lmer(sqrt(ShootMassRate) ~ Abundance + logObs.z + (1|Block), data=data)))) %>%
  summarize(broom::tidy(mod))

# Create a dataframe to store the observed values along with the OTU number
drought_obs.vec <- data.frame(F_stat=Phyloseq_example_aovs[Phyloseq_example_aovs$term == "Abundance" &
                                                             Phyloseq_example_aovs$Drought.or.Watered == "D",]$statistic,
                              OTU = Phyloseq_example_aovs[Phyloseq_example_aovs$term == "Abundance" &
                                                             Phyloseq_example_aovs$Drought.or.Watered == "D",]$OTU)


##### 2) Create data frame to permute #####
# dataframes for permutation
dfs_for_perm <- Phyloseq_example_long %>% nest_by(OTU, Drought.or.Watered)
# only drought plants
dfs_for_perm <- dfs_for_perm[dfs_for_perm$Drought.or.Watered == "D",]


##### 3) Set up parallelization using "foreach" and "doParallel" #####
# can be done without parallelization... but this way is around ~3x as fast
# and with higher core CPUs would be even faster.
library(foreach); packageVersion('foreach')
library(doParallel); packageVersion('doParallel')

# Set seed and number of simulations
set.seed(0195739)
n.sims <- 999 # number of simulations to run

#setup parallel backend to use all CPU cores - 1 (So it doesn't overload your computer)
cores <- detectCores() # Does this detect hyper threaded cores? 
cl <- makeCluster(cores[1]-1) # make a set of parallel R instances 
registerDoParallel(cl) # start other R instances, THESE MUST BE STOPPED LATER ON!

# Pre-reqs lists and variables
ASV <- 1:length(dfs_for_perm$data) # Vector 1 to 11 (indexing any column should work here)
permuted_F_list_all <- c() # List to store all the permutations, List of lists with 612x1000
permuted_F_list <- c() # List per ASV, temporary object will be 1x1000


##### 4) Run the parallel permutation loop #####
parallel_ASV_permuted_Fs <- foreach(ASV=1:length(dfs_for_perm$data)) %dopar% {
  library(tidyverse) # Need to load as daughter R instances won't have it
  library(lmerTest)  # Need to load as daughter R instances won't have it
  if (ASV==1) {(permuted_F_list_all <- c())} # Master List
  for (simulation in 1:n.sims){
    if (simulation==1) {(permuted_F_list <- c())} # Simulation list per ASV
    #print(simulation)
    data.permuted  <- mutate(dfs_for_perm$data[[ASV]], OTU = sample(Abundance, size = length(Abundance), replace = FALSE))
    model.permuted <- lmer(sqrt(ShootMassRate) ~ OTU + logObs.z + (1|Block), data = data.permuted)
    anova.permuted <- rownames_to_column(as.data.frame(anova(model.permuted)), var = "Term")
    F.permuted <- filter(anova.permuted, Term == "OTU")$`F value`
    permuted_F_list <- append(permuted_F_list, values=F.permuted) # Can I name by ASV here?
  }
  print(paste("Finished permuting ASV ", ASV))
  permuted_F_list_all[[ASV]] <- permuted_F_list
}
stopCluster(cl)

# Save the permuted values for quick reloading
saveRDS(parallel_ASV_permuted_Fs, "ASV_permutations_SMR.rds")
parallel_ASV_permuted_Fs <- readRDS("ASV_permutations_SMR.rds")


# some example histograms and stats
drought_obs.vec

# bASV_1
hist(parallel_ASV_permuted_Fs[[1]], breaks = 30)
abline(v = drought_obs.vec[drought_obs.vec$OTU == "bASV_1",]$F_stat, col="red", lwd=3, lty=2)
sum(parallel_ASV_permuted_Fs[[1]] >= drought_obs.vec[drought_obs.vec$OTU == "bASV_1",]$F_stat, na.rm = TRUE)/(n.sims+1)
# non-significant

# bASV_2
hist(parallel_ASV_permuted_Fs[[5]], breaks = 30)
abline(v = drought_obs.vec[drought_obs.vec$OTU == "bASV_2",]$F_stat, col="red", lwd=3, lty=2)
sum(parallel_ASV_permuted_Fs[[5]] >= drought_obs.vec[drought_obs.vec$OTU == "bASV_2",]$F_stat, na.rm = TRUE)/(n.sims+1)
# non-significant

# bASV_55
hist(parallel_ASV_permuted_Fs[[8]], breaks = 30)
abline(v = drought_obs.vec[drought_obs.vec$OTU == "bASV_55",]$F_stat, col="red", lwd=3, lty=2)
sum(parallel_ASV_permuted_Fs[[8]] >= drought_obs.vec[drought_obs.vec$OTU == "bASV_55",]$F_stat, na.rm = TRUE)/(n.sims+1)
# non-significant

# bASV_95
hist(parallel_ASV_permuted_Fs[[11]], breaks = 30)
abline(v = drought_obs.vec[drought_obs.vec$OTU == "bASV_95",]$F_stat, col="red", lwd=3, lty=2)
sum(parallel_ASV_permuted_Fs[[11]] >= drought_obs.vec[drought_obs.vec$OTU == "bASV_95",]$F_stat, na.rm = TRUE)/(n.sims+1)
# non-significant


##### Ideas for improvement #####
# + The initial passing of the ASV=1:length(dfs_for_perm$data) to foreach is weird
#   It would be better to name each of these lists with the OTU/ASV number...
#   This seems easy enough to work out.
# + Make a function to plot histograms
#    - Would need to take the observed vector, permuted vector, ASV number
###