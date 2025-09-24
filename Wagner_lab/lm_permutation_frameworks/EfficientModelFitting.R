####### Description #######
"Our goal is to fit the same model to a large number of traits.
This code uses lists and a custom function to do so efficiently.
The example data is a spreadsheet of plant traits that were measured by
Dr. Nichole Ginnan. In this demo we will answer the question, 
'which traits differ between maize and eastern gamagrass?'"

####### Load packages #######
library(tidyverse)
library(lme4)
library(lmerTest)

####### Load and clean data #######
data <- read.delim('example_plant_data.csv',sep=',')
summary(data) # Make sure all variables are the correct type
# Fix the ones that aren't:
data$PlantID <- as.character(data$PlantID)
data$Species <- as.factor(data$Species)
data$Block <- as.factor(data$Block)
data$Inoculum_type <- as.factor(data$Inoculum_type)
summary(data) # Looks better now

####### Make a template model and output for a single trait #######
model0.mixed <- lmer(LAT.Average.area.of.metaxylem ~ Species + (1|Block),data=data) # example mixed-effects model

results0.mixed <- as.data.frame(anova(model0.mixed)) # extract the results

# Make a "template" for the output by replacing all the values with 0s
blank.results.mixed <- results0.mixed # copy it over
blank.results.mixed[,c('Sum Sq', 'Mean Sq', 'NumDF', 'DenDF', 'F value', 'Pr(>F)')] <- 0

####### Write a function to do this process for any trait #######
myModel = function(trait, myData){
  #trait is trait name (must be part of colnames(df))
  #myData is the dataframe we are analyzing
  form = as.formula(paste(trait, '~ Species + (1|Block)')) # this is where you specify what model to use
  assign('form',form, envi = globalenv())
  mod = tryCatch(     
    lmer(form, data = myData), # here is where you actually fit the model to the data
    error = function(cond) return(NULL)) # if the model fails to fit for some reason, mod will have a value of NULL
  #get the output
  if (!is.null(mod)) {results = as.data.frame(anova(mod))}
  else {results = blank.results.mixed} # if mod is NULL, the function will return the blank template
  return(results)
}

####### Test out the function #######
myModel('LAT.Average.area.of.metaxylem',data)

####### Apply the function to multiple traits #######
# Make a list of traits you want to model:
colnames(data) # look at the options
rootTraits <- colnames(data)[c(5,50:68)] # specify by subsetting
shootTraits <- list('Shoot_mass_g','LICOR_E','LICOR_A','LICOR_water_use_eff','iWUE','LICOR_Ci','LICOR_gsw_stomata') # specify by listing

# Apply the function to all traits in the list, store results in a new list called rootResults or shootResults
rootResults <- lapply(rootTraits, FUN = myModel, myData = data)
names(rootResults) <- rootTraits # name each result after the trait it describes
shootResults <- lapply(shootTraits, FUN = myModel, myData = data)
names(shootResults) <- shootTraits # name each result after the trait it describes

# View all results:
rootResults
shootResults
# But you're not done yet! Your type I error rate is inflated!
# You also need to check whether the assumptions are met by inspecting the residuals!

####### *** adjust P-values to correct for multiple comparisons *** #######
# Make a new function:
list.padj <- function(results_list,padjMethod='fdr') { # will use FDR by default
  # Extract the Pr(>F) values
  pr_values <- sapply(results_list, function(x) x[['Pr(>F)']])
  
  # Apply p.adjust to the extracted p-values
  adjusted_pvalues <- p.adjust(pr_values,method=padjMethod)
  
  # Add the adjusted p-values back to the lists
  updated_results <- lapply(seq_along(results_list), function(i) {
    results_list[[i]]$`Pr(>F)_adjusted` <- adjusted_pvalues[i]
    return(results_list[[i]])
  })
  
  return(updated_results)
}

# Apply this function to your list of results to add the corrected p-values
rootResults <- list.padj(rootResults)
shootResults <- list.padj(shootResults)

####### *** Check that assumptions are met for each trait *** #######
# coming soon