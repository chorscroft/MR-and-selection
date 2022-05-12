library(rms2)
library(TwoSampleMR)

no_fathered <- rms2$new("ukb-b-2227")
no_fathered$extract_gwashits()
no_fathered$gwashits
r_g_y_bsen<-data.frame(rsid=no_fathered$gwashits$rsid,r2=get_r_from_bsen(no_fathered$gwashits$beta,no_fathered$gwashits$se,no_fathered$gwashits$n)^2)
traits_filter<-list()
no_fathered_mr_filtered<-list()
r_x_y_bsen<-list()
for (hit in 1:nrow(no_fathered$gwashits)){
  # Scan OpenGWAS for associations with each GWAS hit
  no_fathered$scan_rsid(no_fathered$gwashits$rsid[hit])
  if (nrow(no_fathered$rsid_scan[[no_fathered$gwashits$rsid[hit]]])>0){
    # Perform colocalisation for each of the candidate associations
    no_fathered$coloc_scan(no_fathered$gwashits$rsid[hit])
    # Perform MR for all the candidate traits
    no_fathered$mr(no_fathered$gwashits$rsid[hit], exclude_rsid_region=TRUE)
    # Filter on colocalisation, traits synonymous with the outcome, eqtls
    traits_filter[[no_fathered$gwashits$rsid[hit]]]<-data.frame(id=no_fathered$rsid_scan[[no_fathered$gwashits$rsid[hit]]]$id,keep=no_fathered$coloc_result[[no_fathered$gwashits$rsid[hit]]]$PP.H4.abf>0.8 &
      no_fathered$rsid_scan[[no_fathered$gwashits$rsid[hit]]]$id!="ieu-b-4760" &
      no_fathered$rsid_scan[[no_fathered$gwashits$rsid[hit]]]$id!="ukb-a-304" &
      no_fathered$rsid_scan[[no_fathered$gwashits$rsid[hit]]]$id!="ukb-b-1209" &
      substr(no_fathered$rsid_scan[[no_fathered$gwashits$rsid[hit]]]$id,1,4)!="eqtl"
      )
    # MR analysis of traits that pass the filter and have a pval < 0.05
    no_fathered_mr_filtered[[no_fathered$gwashits$rsid[hit]]]<-no_fathered$mr_scan[[no_fathered$gwashits$rsid[hit]]][no_fathered$mr_scan[[no_fathered$gwashits$rsid[hit]]]$id.exposure %in% traits_filter[[no_fathered$gwashits$rsid[hit]]]$id[traits_filter[[no_fathered$gwashits$rsid[hit]]]$keep] & no_fathered$mr_scan[[no_fathered$gwashits$rsid[hit]]]$pval<0.05,]
    
    # calculate r squared value
    r_x_y_bsen[[no_fathered$gwashits$rsid[hit]]]<-data.frame(id=no_fathered_mr_filtered[[no_fathered$gwashits$rsid[hit]]]$id.exposure,outcome=no_fathered_mr_filtered[[no_fathered$gwashits$rsid[hit]]]$exposure,r2=get_r_from_bsen(no_fathered_mr_filtered[[no_fathered$gwashits$rsid[hit]]]$b,no_fathered_mr_filtered[[no_fathered$gwashits$rsid[hit]]]$se,no_fathered_mr_filtered[[no_fathered$gwashits$rsid[hit]]]$nsnp)^2)
  }
}

# Do MVMR on all traits that passed the filter
exposure_ids<-NULL
for (i in 1:length(no_fathered_mr_filtered)){
  exposure_ids<-c(exposure_ids,no_fathered_mr_filtered[[i]]$id.exposure)
}
exposure_ids<-unique(exposure_ids)

d <- mv_extract_exposures(c(exposure_ids))
o <- extract_outcome_data(d$SNP, "ukb-b-2227")
d <- mv_harmonise_data(d, o)
mvmr1<-mv_multiple(d)

# lasso
lasso<-mv_lasso_feature_selection(d)

# mvmr with lassoed exposures
d <- mv_extract_exposures(c(lasso$exposure))
o <- extract_outcome_data(d$SNP, "ukb-b-2227")
d <- mv_harmonise_data(d, o)
mvmr2<-mv_multiple(d)

# Factors that are the most important
imp_traits<-mvmr2$result$exposure[mvmr2$result$pval<0.05]

# sankey plot
r2<-list()
traits<-NULL
rsid<-NULL
base_r2<-NULL
value_r2<-NULL
for (hit in 1:nrow(no_fathered$gwashits)){
  r2[[no_fathered$gwashits$rsid[hit]]]<-r_x_y_bsen[[no_fathered$gwashits$rsid[hit]]][r_x_y_bsen[[no_fathered$gwashits$rsid[hit]]]$outcome %in% imp_traits,]
  traits<-c(traits,substr(r2[[no_fathered$gwashits$rsid[hit]]]$outcome,1,regexpr("\\|",r2[[no_fathered$gwashits$rsid[hit]]]$outcome)-2),"Unknown")
  rsid<-c(rsid,rep(no_fathered$gwashits$rsid[hit],length(r2[[no_fathered$gwashits$rsid[hit]]]$outcome)+1))
  base_r2<-c(base_r2,r2[[no_fathered$gwashits$rsid[hit]]]$r2,1-sum(r2[[no_fathered$gwashits$rsid[hit]]]$r2))
  value_r2<-c(value_r2,c(r2[[no_fathered$gwashits$rsid[hit]]]$r2,1-sum(r2[[no_fathered$gwashits$rsid[hit]]]$r2))*r_g_y_bsen$r2[r_g_y_bsen$rsid==no_fathered$gwashits$rsid[hit]])
}

sankey_df<-data.frame(N1=traits,N2=rsid,r2=base_r2,Value=value_r2)

library(riverplot)
library(RColorBrewer)
nodes <- data.frame(ID=unique(c(sankey_df$N1, sankey_df$N2)),
                    x=rep(c(1,2), times=c(length(unique(sankey_df$N1)),length(unique(sankey_df$N2)))))#, y=c(length(unique(sankey_df$N1)):1,length(unique(sankey_df$N2)):1))
palette = paste0(brewer.pal(12, "Set3"), "60")
styles = lapply(1:nrow(nodes), function(n) {
  list(col = ifelse(n<=12,palette[n],palette[n-12]), lty = 0, textcol = "black")
})
names(styles) = nodes$ID
r <- makeRiver(nodes=nodes, edges=select(sankey_df,c("N1","N2","Value")), styles=styles)
d_sty <- list(srt=0, textcex=0.5) # default style
plot(r, plot_area=1, nodewidth=10, default_style=d_sty)


# Save output
jpeg("~/Documents/rms2_analysis/sankey/number_of_children_fathered_sankey.jpeg",width = 20,height = 12,units = "cm",quality = 100,res=300)
plot(r, plot_area=1, nodewidth=10, default_style=d_sty)
dev.off()
library(openxlsx)
excelsheets<-list(mvmr1=mvmr1$result,lasso=lasso,mvmr2=mvmr2$result,r_g_y=r_g_y_bsen)
for (i in 1:length(no_fathered_mr_filtered)){
  excelsheets[[paste0("mr_",names(no_fathered_mr_filtered)[i])]]<-no_fathered_mr_filtered[[i]]
}
for (i in 1:length(r_x_y_bsen)){
  excelsheets[[paste0("r_x_y_",names(r_x_y_bsen)[i])]]<-r_x_y_bsen[[i]]
}
write.xlsx(excelsheets, "~/Documents/rms2_analysis/sankey/number_of_children_fathered_sankey.xlsx", asTable = FALSE, overwrite = TRUE)



























r2<-r_x_y_bsen$rs62250711[r_x_y_bsen$rs62250711$outcome %in% imp_traits,]

# sankey plot
traits<-substr(r2$outcome,1,regexpr("\\|",r2$outcome)-2)
sankey_df<-data.frame(N1=c(traits,rep("Unknown",3)),N2=rep(no_fathered$gwashits$rsid,times=c(6,1,1)),r2=c(r2$r2,1-sum(r2$r2),1,1),Value=c(r2$r2,1-sum(r2$r2),1,1)*rep(r_g_y_bsen$r2,times=c(6,1,1)))

library(riverplot)
nodes <- data.frame(ID=unique(c(sankey_df$N1, sankey_df$N2)),
                    x=rep(c(1,2), times=c(6,3)), y=c(6:1,3:1))
r <- makeRiver(nodes=nodes, edges=select(sankey_df,c("N1","N2","Value")))
d_sty <- list(srt=0, textcex=0.5) # default style
plot(r, plot_area=1, nodewidth=10, default_style=d_sty)

# Save output
jpeg("~/Documents/rms2_analysis/sankey/number_of_children_fathered_sankey.jpeg",width = 12,height = 12,units = "cm",quality = 100,res=300)
plot(r, plot_area=1, nodewidth=10, default_style=d_sty)
dev.off()
library(openxlsx)
excelsheets<-list(mvmr1=mvmr1$result,lasso=lasso,mvmr2=mvmr2$result,r_g_y=r_g_y_bsen)
for (i in 1:length(no_fathered_mr_filtered)){
  excelsheets[[paste0("mr_",names(no_fathered_mr_filtered)[i])]]<-no_fathered_mr_filtered[[i]]
}
for (i in 1:length(r_x_y_bsen)){
  excelsheets[[paste0("r_x_y_",names(r_x_y_bsen)[i])]]<-r_x_y_bsen[[i]]
}
write.xlsx(excelsheets, "~/Documents/rms2_analysis/sankey/number_of_children_fathered_sankey.xlsx", asTable = FALSE, overwrite = TRUE)

     