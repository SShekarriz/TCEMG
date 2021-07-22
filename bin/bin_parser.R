
#!/usr/bin/env Rscript
# SShekarriz Jun03,2021
# Making a Bin info file to be imported into anvio
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[1] = "HMT141_bin_qc.txt"
  args[2] = "HMT141.bac120.summary.tsv"
  args[3] = "bins_across_samples"
  args[4] = "HMT141_bins_sum_cov.txt"
}

#################
library(tidyverse)
#################

checkM_i=args[1] #check output table
taxa_i=args[2] #GTDB taxonomy table
#coverage tables per patients
path=args[3] #path to where all the coverage files are
output=args[4] #output summary file

#######################################################

#Here is a script to merge all related bin Info into one
checkM<-read.csv(checkM_i, sep = "\t")
taxa<-read.csv(taxa_i, sep = "\t")
#importing GC taxonomy
taxa%>% 
  separate(classification, c("Kingdom", "Phylum", "Class",
                           "Order", "Family", "Genus",
                           "Species"), ";.__") %>%
  select(user_genome, Kingdom:Species) %>%
  mutate(Bin.Id=user_genome)-> taxa
#merge them into one
checkM %>%
  select("Bin.Id","Marker.lineage",
         "Completeness","Contamination","Strain.heterogeneity") %>%
  left_join(taxa) %>%
  mutate(Bin.Id = gsub("bin.", "bin_", Bin.Id))-> bins

#coverage tables per patients
patt=".txt"
data.frame(type = paste(dir(path, 
                            pattern = patt))) %>%
# read each file from the directory (current path) and then unnest
mutate(file_contents = map(type, ~ read_tsv(
  file.path(path, .),
            col_names = T))) %>%
         unnest() %>%
  rename(Bin.Id=bins) %>%
  left_join(bins) %>%
  write.table(output, sep = "\t", row.names = F, quote = F)




