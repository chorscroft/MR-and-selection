# Adds in the chromosome and position data to the METAL output

require(dplyr)
require(jsonlite)

# Get file path to project file 
config <- read_json("~/config.json")
mr_sel_path <- config$mr_sel_path

# Read in data files
setwd(paste0(mr_sel_path,"/Meta_analysis_no_children_SD"))
meta<-read.table("METAANALYSIS1.TBL",header = TRUE)
ukb_b_1209<-read.table("ukb-b-1209.tsv",header=TRUE)
ukb_b_2227<-read.table("ukb-b-2227.tsv",header=TRUE)

# only keep columns needed
ukb_b_1209<-select(ukb_b_1209,variant_id,chromosome,base_pair_location)
ukb_b_2227<-select(ukb_b_2227,variant_id,chromosome,base_pair_location)

# join all snp locations into one dataframe
snps<-unique(full_join(ukb_b_1209,ukb_b_2227))

# join snp locations to meta file
new_meta<-left_join(meta,snps,by=c("MarkerName"="variant_id"))
write.table(new_meta,"METAANALYSIS1_NEW.TBL",quote = FALSE,row.names = FALSE)
