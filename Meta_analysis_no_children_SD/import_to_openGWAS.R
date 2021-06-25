require(GwasDataImport)
require(jsonlite)

# Get file path to project file 
config <- read_json("~/config.json")
mr_sel_path <- config$mr_sel_path

# Define file name
filename=paste0(mr_sel_path,"/Meta_analysis_no_children_SD/METAANALYSIS1_NEW.TBL")

# Authenticate
ieugwasr::get_access_token()

# pick ID
#id <- "ieu-b-4747"

# Initialise
x <- Dataset$new(filename=filename)#, igd_id=id)

# Specify columns
x$determine_columns(list(chr_col=16, snp_col=1, pos_col=17, oa_col=3, ea_col=2, eaf_col=4, beta_col=8, se_col=9, pval_col=10))

# Process dataset
x$format_dataset()

# Input metadata
# sample size = 209,872 + 250,782
x$collect_metadata(list(trait="Number of children", build="HG19/GRCh37", category="Categorical Ordered", subcategory="NA", group_name="public", population="European", sex="Males and Females",unit="SD",ontology="NA",sample_size="460654",author="Clare Horscroft",year="2021",consortium="MRC-IEU"))

# Upload metadata
x$api_metadata_upload()

# Upload summary data
x$api_gwasdata_upload()

# View report
x$api_report()

# Release dataset
x$api_gwas_release()

# Delete wd
x$delete_wd()
