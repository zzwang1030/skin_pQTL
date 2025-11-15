library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=16, face = "bold")
  )
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

## proteome process ####
# 1/2 minimum impution, VSN normalizaiton, average of different samples of the same age
qq2 <-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv",header  = T)
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- tidyr::separate_rows(qq3, PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn); fit <- vsn2(mm); mm1 <- predict(fit, mm) 
mm2 <- as.data.table(mm1); mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); names(mm3)<-c("pro_id","sample","intensity")

## protein level vs age: rho ####
hh3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Sample_info.csv", header = T)
jj<-merge(mm3, hh3, by.x="sample", by.y="No"); jj<-jj[order(pro_id, age)]
jj[, intensity2:=mean(intensity), by=.(pro_id, age)] # 同一年龄不同样本取平均
jj1<-unique(jj[,.(pro_id, intensity2, age)])
jj1[, c("rho_age", "rho_p_age") := .(cor.test(intensity2, age, method = "spearman")$estimate, cor.test(intensity2, age, method = "spearman")$p.value), by = .(pro_id)]
fwrite(jj1, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", col.names = T, row.names = F, sep=",", quote = F)
jj2<-unique(jj1[,.(pro_id, rho_age, rho_p_age)])
jj2<-merge(jj2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
nrow(jj2); nrow(jj2[rho_p_age<0.05]); nrow(jj2[rho_p_age<0.05 & rho_age>0]); nrow(jj2[rho_p_age<0.05 & rho_age<0]) # 5956; 891; 461; 430

## protein level vs age: heatmap ####
library(ComplexHeatmap); library(circlize)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
jj_mat <- dcast(jj1, pro_id ~ age, value.var = "intensity2")
mat <- as.matrix(jj_mat[, -1]); rownames(mat) <- jj_mat$pro_id
mat_scaled <- t(scale(t(mat)))
age_order <- sort(as.numeric(colnames(mat_scaled))); mat_scaled <- mat_scaled[, as.character(age_order)]

jj_pos<-jj2[rho_p_age<0.05 & rho_age>0][order(rho_p_age, -rho_age)]
jj_neg<-jj2[rho_p_age<0.05 & rho_age<0][order(rho_p_age, rho_age)]
jj_sig<-rbind(jj_pos, jj_neg)
tmp <- mat_scaled[match(jj_pos$pro_id, rownames(mat_scaled)), ] 
tmp <- mat_scaled[match(jj_sig$pro_id, rownames(mat_scaled)), ]

age_labels <- rep("", ncol(tmp)); idx_to_label <- seq(1, ncol(tmp), by = 5)
age_labels[idx_to_label] <-  colnames(tmp)[idx_to_label]; age_labels
ht<-Heatmap(tmp, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, 
            column_title = "Age", column_title_side = "bottom", row_title = "Protein level",
            column_names_rot = 45, column_names_gp = gpar(fontface = "bold"), column_labels = age_labels, heatmap_legend_param = list(direction = "horizontal"))
draw(ht, heatmap_legend_side = "top")

## protein level vs age: case plot ####
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
jj2<-unique(jj1[,.(pro_id, rho_age, rho_p_age)])
jj_pos<-jj2[rho_p_age<0.05 & rho_age>0][order(rho_p_age, -rho_age)]
jj_neg<-jj2[rho_p_age<0.05 & rho_age<0][order(rho_p_age, rho_age)]
jj2<-merge(jj2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")

case_test<-c(jj_pos$pro_id[1:5], jj_neg$pro_id[1:5]) # TIMP3, APCS, HTRA1, NT5E, EMILIN2;; DDX46, COL1A1, COL1A2, CPXM1, ECM2 = P35625, P02743, Q92743, P21589, Q9BXX0;; Q7L014, P02452, P08123, Q96SM3, O94769
for (i in 1:length(case_test)) {
  p<-ggplot(jj1[pro_id == case_test[i]], aes(x = age, y = intensity2)) + geom_point(color = "grey30", size = 1.5, alpha = 1) + 
    geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth = 1.2) + labs(x = "Age", y = "log2(Protein abundance)") + mytheme
  print(p)
}

## protein level vs age: GO for age-dependent proteins ####
## Protein–protein interaction (PPI) analysis by metascape online
library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
jj2<-unique(jj1[,.(pro_id, rho_age, rho_p_age)])
jj_pos<-jj2[rho_p_age<0.05 & rho_age>0][order(rho_p_age, -rho_age)]
jj_neg<-jj2[rho_p_age<0.05 & rho_age<0][order(rho_p_age, rho_age)]

# positive-ce protein; similar for negative-cor proteins
enrich_go1 <- enrichGO(jj_pos$pro_id, universe= jj2$pro_id,
                        OrgDb = org.Hs.eg.db, keyType = "UNIPROT", ont = "BP",
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
barplot(enrich_go1, showCategory = 20, x = "GeneRatio", color = "p.adjust")
# Semantic simplification
gg1<-as.data.table(enrich_go1@result); gg1<-gg1[order(p.adjust)]; gg1<-gg1[p.adjust<0.05]
library(rrvgo); library(GOSemSim); library(org.Hs.eg.db)
go_ids <- gg1$ID; scores <- setNames(-log10(gg1$`p.adjust`), gg1$ID)
simMatrix <- calculateSimMatrix(go_ids, orgdb = "org.Hs.eg.db", ont = "BP", method = "Rel")
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = "org.Hs.eg.db")
reducedTerms<-as.data.table(reducedTerms); reducedTerms<-reducedTerms[order(cluster, -score)]
tmp<-reducedTerms
gg11<-gg1[ID%in%tmp$go]
gg11$negLogP <- -log10(gg11$p.adjust)
gg11$Description <- factor(gg11$Description, levels = rev(gg11$Description[order(gg11$negLogP, decreasing = TRUE)]))
ggplot(gg11, aes(x = negLogP, y = Description, fill = Count)) + geom_bar(stat = "identity") +
  scale_fill_gradient(name = "Gene Count", low = "mistyrose", high = "brown4") +
  labs(x = "-log10(p.adj)", y = NULL) + mytheme + 
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.title = element_text(size = 14), legend.text = element_text(size = 14))

## protein variance vs age ####
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
jj1[, age_group := cut(age, breaks = c(13, 34.5, 52, 69.5, 95), labels = c("age1", "age2", "age3", "age4"), include.lowest = TRUE, right = TRUE)]
table(jj1[pro_id=="A0A075B6H7"]$age_group) # age1 age2 age3 age4  = 18   18   17   18 
jj1_cv <- jj1[, .(SD = sd(intensity2, na.rm = TRUE), CV = sd(intensity2, na.rm = TRUE) / mean(intensity2, na.rm = TRUE)), by = .(pro_id, age_group)]

tmp<-jj1_cv[age_group=="age1"|age_group=="age4"]
tmp$age_group<-factor(tmp$age_group,levels = c("age1","age4"))
wilcox.test(jj1_cv[age_group=="age1"]$SD, jj1_cv[age_group=="age4"]$SD, paired = T)$`p.value` # 1.221947e-257
ggplot(tmp, aes(x=age_group, y=SD))+
  geom_line(aes(group = pro_id), color = "grey60", linewidth = 0.2, alpha = 0.2) +
  geom_jitter(aes(color = age_group, group = age_group), width = 0.25, size = 0.8, alpha = 0.2) +
  geom_boxplot(aes(group = age_group), width = 0.4,  outlier.shape = NA, fill = NA, color = "black", lwd = 1.1) +
  scale_color_manual(values = c("indianred", "steelblue"))+ # scale_fill_manual(values = alpha(c("#E04643", "#686F76"), 0.8))+
  labs(x = ' ', y = 'SD of protein level') + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(labels = c("age1" = "Young", "age4" = "Old"))+
  mytheme + theme(legend.position="none")

## Abundance correlation of RNA  vs protein ####
rr<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # a total table contain raw cnt, log2Normalized_cnt and FPKM value; The complete file is too large, so only a demo was provided here.
# rr<-fread("~/Desktop/省皮/project/pQTL/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # complete data
sel_col<-grep("gene_ID|.+SRNA.+fpkm",names(rr),value=T); rr[,..sel_col]->rr_rpkm; names(rr_rpkm)<-gsub("SRNA_fpkm","",names(rr_rpkm))
sel_col<-grep("gene_ID|.+SRNA$",names(rr),value=T); rr[,..sel_col]->rr_cnt; names(rr_cnt)<-gsub("SRNA","",names(rr_cnt))
names(rr_rpkm)[1]<-"gene_id";names(rr_cnt)[1]<-"gene_id"
identical(names(rr_rpkm),names(rr_cnt)) # TRUE
rr_rpkm<-rr_rpkm[grep('ENSG',gene_id)]; rr_cnt<-rr_cnt[grep('ENSG',gene_id)] # remove predicted genes
identical(rr_rpkm$gene_id, rr_cnt$gene_id); nrow(rr_cnt) # TRUE; 60564
non_zero_sample_cols <- colSums(rr_rpkm[, -1, with = FALSE]) != 0
rr_rpkm <- rr_rpkm[, c(TRUE, non_zero_sample_cols), with = FALSE] # Remove the sample column containing all zeros
rr_rpkm <- rr_rpkm[rowSums(rr_rpkm[, -1, with = FALSE]) != 0]; dim(rr_rpkm) # Remove gene rows with all values of 0; 52080 gene; 208-1 smp
# defining expressed gene:  RPKM >0.1 in at least 10 smps; summed cnt > 5 for each gene. Adapted from ref: GTEx 29022597
rr_cnt[, `:=`(mean_expr = rowMeans(.SD), sum_expr = rowSums(.SD)), .SDcols = !'gene_id'] # Calculate the mean and sum of cnt for each gene (excluding the gene_id column).
melt(rr_rpkm, id.vars = "gene_id")->rr_rpkm2; names(rr_rpkm2)<-c("gene_id","sample","value")
rr_rpkm2[, rpkm01_gene:=sum(value>0.1), by=.(gene_id)]
rr_rpkm3 <- rr_rpkm2[rpkm01_gene >= 10]; length(unique(rr_rpkm3$gene_id)) # 33359；RPKM >0.1 in at least 10 smps
rr_rpkm3 <- rr_rpkm3[gene_id %in% rr_cnt[sum_expr >= 6, gene_id],]; length(unique(rr_rpkm3$gene_id)); length(unique(rr_rpkm3$sample)) # summed cnt > 5 for each gene; 33359,207

## rho for overlapped of samples and genes
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
qq2 <-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv",header  = T)
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; names(qq3)<-gsub("SP", "", names(qq3))
qq3<-as.data.table(qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";")) # 5956
qq4<-merge(qq3, name_id_reviewd[,1:2], by.x="PG.ProteinAccessions", by.y="pro_id") # 6659；有些 pro_id对应多个 ENSG_id
qq5<-melt.data.table(qq4[,-1],id.vars = "gene_id"); names(qq5)<-c("gene_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(gene_id)]; qq5[,intensity:=ifelse(is.na(intensity),min/2,intensity)] # 填补 NA

jj<-merge(qq5[,1:3], rr_rpkm3[,1:3], by=c("gene_id","sample")); names(jj)[3:4]<-c("protein","rna")
length(unique(qq5$gene_id)); length(unique(qq5$sample)) # 6653; 248
length(unique(rr_rpkm3$gene_id)); length(unique(rr_rpkm3$sample)) # 33359; 207
length(unique(jj$gene_id)); length(unique(jj$sample)) # 5860; 199
jj[,c("protein_log","rna_log"):=.(log(protein), log(rna+0.01))]
jj[,c("protein_log_median","protein_log_mean","rna_log_median","rna_log_mean"):= .(median(protein_log), mean(protein_log), median(rna_log), mean(rna_log)), by=.(gene_id)]
jj[, c("rho", "rho_p") := .(cor.test(protein_log, rna_log, method = "spearman")$estimate, 
                            cor.test(protein_log, rna_log, method = "spearman")$p.value), by = .(gene_id)]
jj1<-unique(jj[,.(gene_id, rho, rho_p, protein_log_median, protein_log_mean, rna_log_median,rna_log_mean)])
cor.test(jj1[protein_log_mean>0 & rna_log_mean>0]$protein_log_mean, jj1[protein_log_mean>0 & rna_log_mean>0]$rna_log_mean, method = "spearman") # 0.323
fwrite(jj1[protein_log_mean>0 & rna_log_mean>0], file="~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/expressed_RNA_protein_rho.csv", sep=",", col.names = T, row.names = F,quote=F)

jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/expressed_RNA_protein_rho.csv")
library(viridis); library(ggpointdensity); library(ggpubr)
nrow(jj1[protein_log_mean>0 & rna_log_mean>0]) # 5127
ggplot(data = jj1[protein_log_mean>0 & rna_log_mean>0], mapping = aes(x = protein_log_mean, y = rna_log_mean)) +
  geom_pointdensity() + scale_color_viridis() +
  labs(x = 'Protein abundance', y = 'RNA expression', color = NULL)+
  geom_smooth(method='lm', se=T,formula= y~x,col="red", linetype="dashed")+  
  coord_cartesian(xlim = c(0, 11), ylim = c(0, 7)) + mytheme 

## RNA vs protein: by age group ####
hh3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Sample_info.csv", header = T)

qq2 <-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv",header  = T)
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA
mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn); fit <- vsn2(mm); mm1 <- predict(fit, mm) 
mm2 <- as.data.table(mm1); mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); names(mm3)<-c("pro_id","sample","intensity")

rr<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # a total table contain raw cnt, log2Normalized_cnt and FPKM value
# rr<-fread("~/Desktop/省皮/project/pQTL/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # complete data
sel_col<-grep("gene_ID|.+SRNA.+fpkm",names(rr),value=T); rr_rpkm<-rr[,..sel_col]; names(rr_rpkm)<-gsub("SRNA_fpkm","",names(rr_rpkm))
names(rr_rpkm)[1]<-"gene_id"; rr_rpkm<-rr_rpkm[grep('ENSG',gene_id)]
rr_rpkm2<-melt(rr_rpkm, id.vars = "gene_id"); names(rr_rpkm2)<-c("gene_id","sample","value") 
rr_rpkm4<-merge(rr_rpkm2, unique(name_id_reviewd[,1:2]), by="gene_id")

jj<-merge(mm3, rr_rpkm4[, .(pro_id, sample, value)], by=c("pro_id", "sample")); names(jj)[3:4]<-c("pro_raw","rna_raw")
jj<-jj[order(pro_id, sample, -pro_raw, -rna_raw)]; jj[,id:=paste(pro_id, sample, sep="_")]; jj<-jj[!duplicated(id)]
length(unique(mm3$pro_id)); length(unique(rr_rpkm4$pro_id)); length(unique(jj$pro_id)) # 5956; 19275; 5929
length(unique(mm3$sample)); length(unique(rr_rpkm4$sample)); length(unique(jj$sample)) # 248; 207; 199
jj[,c("pro_log","rna_log"):=.(pro_raw, log2(rna_raw+0.001))] # The protein has already undergone VSN log2 transformation; performing log first, then mean, is more robust and reduces the impact of extremes
jj[,c("pro_log_mean", "pro_log_median", "rna_log_mean", "rna_log_median"):=.(mean(pro_log), median(pro_log), mean(rna_log), median(rna_log)), by=.(pro_id)]

jj<-merge(jj, hh3, by.x="sample", by.y="No"); jj<-jj[order(pro_id, age)]
jj[, age_group := cut(age, breaks = c(13, 38, 52, 62, 95), labels = c("age1", "age2", "age3", "age4"), include.lowest = TRUE, right = TRUE)]
table(jj[pro_id=="A0FGR8"]$age_group) # age1 age2 age3 age4 = 59   51   44   45  
table(jj[pro_id=="A0FGR8"]$gender) # female   male = 72    127 
jj[, c("rho_age", "rho_p_age", "R_age", "R_p_age") := .(cor.test(pro_log, rna_log, method = "spearman")$estimate, cor.test(pro_log, rna_log, method = "spearman")$p.value,
                                                        cor.test(pro_log, rna_log, method = "pearson")$estimate,  cor.test(pro_log, rna_log, method = "pearson")$p.value), by = .(pro_id, age_group)]
jj[, c("rho_sex", "rho_p_sex", "R_sex", "R_p_sex") := .(cor.test(pro_log, rna_log, method = "spearman")$estimate, cor.test(pro_log, rna_log, method = "spearman")$p.value,
                                                        cor.test(pro_log, rna_log, method = "pearson")$estimate,  cor.test(pro_log, rna_log, method = "pearson")$p.value), by = .(pro_id, gender)]

tmp<-jj[pro_log_mean>0 & rna_log_mean>0]; length(unique(tmp$pro_id)) # 5115
jj_age<-unique(jj[pro_id%in%tmp$pro_id, .(pro_id, age_group, rho_age, rho_p_age, R_age, R_p_age)]); jj_age<-jj_age[!is.na(rho_age)]

# count barplot
nrow(jj_age[age_group=="age1"]); nrow(jj_age[age_group=="age1" & rho_p_age<0.05]); nrow(jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age>0]); nrow(jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age<0]) # 5115; 1026; 940; 86 /1026=91.62 
nrow(jj_age[age_group=="age4"]); nrow(jj_age[age_group=="age4" & rho_p_age<0.05]); nrow(jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age>0]); nrow(jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age<0]) # 5115; 598; 526; 72 /598=87.96
chisq.test(matrix(c(5115, 940, 5115, 526), byrow = T, ncol = 2))$`p.value` # 6.069992e-24
mymat<-data.table(group = rep(c("Young", "Old"), each = 2), status = rep(c("Positive", "Negative"), times = 2),count = c(940, 86, 526, 72))
mymat$group<-factor(mymat$group, levels = c("Young", "Old")); mymat$status<-factor(mymat$status, levels = c("Negative", "Positive"))
ggplot(mymat, aes(x = group, y = count, fill = status)) + 
  geom_bar(stat = "identity", width = 0.7) + labs(x = NULL, y = "Number of genes") + 
  scale_fill_manual(values = c("Positive" = "indianred", "Negative" = "steelblue"), labels = c("Negative cor", "Positive cor")) +
  mytheme+ theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank())

# decoupling GO
tmp1<-jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age>0] # positive in young
tmp2<-jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age>0] # positive in old
old_loss<-setdiff(tmp1$pro_id, tmp2$pro_id); old_gain<-setdiff(tmp2$pro_id, tmp1$pro_id)
nrow(jj_age[age_group=="age1"]); nrow(tmp1); nrow(tmp2); length(old_loss); length(old_gain) # 5115; 940; 526; 687; 273

df <- data.frame(age_group = rep(c("Young", "Old"), each = 2), status = rep(c("Positive cor", "Positive loss in old"), 2), count = c(940, 0, 253, 687))
df$age_group <- factor(df$age_group, levels = c("Young", "Old"))
ggplot(df, aes(x = age_group, y = count, fill = status)) + geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual( values = c("Positive cor" = "indianred", "Positive loss in old" = "grey70"), name = NULL, labels = c("Couping in young", "Decoupling in Old")) +
  labs(x = NULL, y = "Number of Genes") + mytheme +
  theme(legend.position = "top",legend.direction = "horizontal", legend.text = element_text(size = 16))

library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
enrich_go3 <- enrichGO( old_loss, OrgDb = org.Hs.eg.db, keyType = "UNIPROT", universe= unique(jj_age$pro_id), 
                        ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
tmp3<-as.data.table(enrich_go3@result); tmp3<-tmp3[order(p.adjust)]; tmp3<-tmp3[p.adjust<0.05] # Manual concise
df_plot <- tmp3[order(p.adjust)]; df_plot$negLogP <- -log10(df_plot$p.adjust)
df_plot$Description <- factor(df_plot$Description, levels = rev(df_plot$Description[order(df_plot$negLogP, decreasing = TRUE)]))
ggplot(df_plot, aes(x = negLogP, y = Description, fill = Count)) + geom_bar(stat = "identity") +
  scale_fill_gradient(name = "Gene Count", low = "mistyrose", high = "brown4") +
  labs(x = "-log10(p.adj)", y = NULL) + mytheme + 
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.title = element_text(size = 14), legend.text = element_text(size = 14))

# rho compare ECDF
wilcox.test(jj_age[age_group=="age1"]$rho_age, jj_age[age_group=="age4"]$rho_age, paired = TRUE) # 3.964501e-20,用这个
ggplot(jj_age[age_group %in% c("age1", "age4")], aes(x = rho_age, color = age_group)) + 
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_color_manual(values = c("age1" = "indianred", "age4" = "steelblue"), labels = c("age1" = "Young", "age4" = "Old")) +
  labs(x = "Spearman rho", y = "Cumulative Distribution Fraction", color = NULL) + coord_cartesian(xlim = c(-0.3, 0.6)) +
  mytheme+theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.key.width = unit(2.5, "lines"), legend.text=element_text(size=18, face = "bold"))

## protein predict age model ####
library(glmnet); library(iml)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
thresholds <- 10^seq(0, -13, by = -1)
results <- data.frame(threshold = thresholds, n_proteins = NA, MAE = NA, R = NA)

for (i in seq_along(thresholds)) {
  thr <- thresholds[i]
  expr_mat <- dcast(jj1[rho_p_age < thr], age ~ pro_id, value.var = "intensity2_zs"); expr_mat <- na.omit(expr_mat)
  if (nrow(expr_mat) < 5 || ncol(expr_mat) <= 2) next  # Skip if sample size is too small or protein content is too low
  
  y <- expr_mat$age; X <- as.matrix(expr_mat[, !c("age")])
  n <- nrow(X); pred_age <- rep(NA, n)
  for (j in 1:n) { # Modeling + Cross-validation: Choosing lambda
    train_idx <- setdiff(1:n, j)
    X_train <- X[train_idx, ]; y_train <- y[train_idx]
    X_test <- X[j, , drop = FALSE]
    
    cvfit <- cv.glmnet(X_train, y_train, alpha = 0.5, nfolds = 5)
    best_lambda <- cvfit$lambda.min
    model <- glmnet(X_train, y_train, alpha = 0.5, lambda = best_lambda)
    pred_age[j] <- predict(model, newx = X_test) }
  
  mae <- mean(abs(pred_age - y)); R <- cor(y, pred_age)  # Evaluation indicators
  results$n_proteins[i] <- ncol(X); results$MAE[i] <- mae; results$R[i] <- R 
}

# Plotting: The left y-axis is MAE, and the right y-axis is R.
library(ggplot2); library(scales); library(patchwork)  
results<-as.data.table(results)
results$threshold_label <- format(results$threshold, scientific = TRUE) 
results_long <- melt(results, id.vars = c("threshold", "threshold_label", "n_proteins"), measure.vars = c("MAE", "R"), variable.name = "Metric", value.name = "Value")
results_long[, Value_scaled := ifelse(Metric == "R", (Value - 0.5) * 10 + 5, Value)]

p1 <- ggplot(results_long, aes(x = threshold_label, y = Value_scaled, group = Metric, color = Metric)) + 
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_y_continuous(name = "MAE", sec.axis = sec_axis(~ (. - 5) / 10 + 0.5, name = "Pearson R")) +
  scale_x_discrete(limits = results$threshold_label) +
  labs(x = NULL, color = NULL, title = "Model Performance") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "top")

p2 <- ggplot(results, aes(x = threshold_label, y = log10(n_proteins)))+ geom_col(fill = "grey70") + 
  scale_x_discrete(limits = results$threshold_label) +
  labs(x = "Rho threshold", y = "log10(Count)", title = "Included protein count") +
  theme_minimal(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 / p2 + plot_layout(heights = c(3, 1)) 


## Hyperparameter search; Merging samples of the same age
library(mlr3verse); library(mlr3tuning); library(paradox); library(mlr3learners)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
expr_mat <- dcast(jj1[rho_p_age<1e-9], age ~ pro_id, value.var = "intensity2_zs")  
expr_mat <- na.omit(expr_mat)  

# 0. Modeling part of LOOCV and hyperparameter tuning
y <- expr_mat$age; X <- as.data.frame(expr_mat[, !c("age")]) # Split X and y
n <- nrow(X); pred_age <- rep(NA, n)

for (i in 1:n) { # ~15 min
  # 1. LOOCV partitioning
  train_idx <- setdiff(1:n, i); X_train <- X[train_idx, ]; y_train <- y[train_idx]; X_test <- X[i, , drop = FALSE]
  # 2. build Task
  train_data <- data.table::data.table(X_train); train_data$age <- y_train
  task <- TaskRegr$new("age_prediction", backend = train_data, target = "age")
  # 3. define learner (elastic net via glmnet)
  learner <- lrn("regr.glmnet", predict_type = "response")
  # 4. Search space
  param_space <- ps(
    alpha = p_dbl(lower = 0, upper = 1),
    lambda = p_dbl(lower = -4, upper = 1, trafo = function(x) 10^x)) # The actual value is 1e-4~10. Log-scale search with lambda sampling is more uniform and the search is more efficient.
  # 5. Tuning
  instance <- TuningInstanceBatchSingleCrit$new(task = task, learner = learner, resampling = rsmp("cv", folds = 5), 
                                                measure = msr("regr.mae"), search_space = param_space, terminator = trm("evals", n_evals = 50),  # 尝试 30 个组合
  )
  tuner <- tnr("random_search"); tuner$optimize(instance)
  # 6. Optimal parameter modeling
  learner$param_set$values <- instance$result_learner_param_vals
  learner$train(task)
  #7. Prediction
  pred_age[i] <- learner$predict_newdata(X_test)$response
}
mae <- mean(abs(pred_age - y)); mae # 7.63
cor.test(y, pred_age, method = "pearson") # R=0.882, P=3.49e-24

# 1. Train the final model using all the data.
# Full Data Construction Task
final_data <- data.table::data.table(X); final_data$age <- y
task_final <- TaskRegr$new("age_prediction", backend = final_data, target = "age")
# learner and search space
learner_final <- lrn("regr.glmnet", predict_type = "response")
param_space <- ps(
  alpha = p_dbl(lower = 0, upper = 1),
  lambda = p_dbl(lower = -4, upper = 1, trafo = function(x) 10^x)  # log-scale lambda
)
# Parameter Tuning Examples
instance_final <- TuningInstanceBatchSingleCrit$new(task = task_final, learner = learner_final, resampling = rsmp("cv", folds = 5), 
                                                    measure = msr("regr.mae"), search_space = param_space, terminator = trm("evals", n_evals = 30))
# Tuning and setting optimal parameters
tuner <- tnr("random_search")
tuner$optimize(instance_final)
learner_final$param_set$values <- instance_final$result_learner_param_vals
# Training the final model
learner_final$train(task_final)
# model performace
pred_all <- learner_final$predict(task_final) 
pred_age_all <- pred_all$response
mae_all <- mean(abs(pred_age_all - y)) # MAE= 6.86
cor.test(pred_age_all, y, method = "pearson") # R=0.91; P=1.50e-27

tmp<-as.data.table(cbind(y, pred_age_all)); names(tmp)<-c("actual", "predict")
ggplot(tmp, aes(x = actual, y = predict)) + geom_point(color = "grey30", size = 1.5, alpha = 1) + geom_smooth(method = "lm", se = TRUE, color = "indianred",  size = 1.2) + 
  labs(x = "Actual age", y = "Predicted age") + coord_cartesian(xlim=c(13,100), ylim = c(13, 100))+ mytheme

# 2. Identifying important proteins using the iml package
library(iml)
# Construct the Predictor object of iml
predictor <- Predictor$new(model = learner_final, data = X, y = y) 
# use permutation importance
imp <- FeatureImp$new(predictor, loss = "mae")  # or "rmse"
imp_dt_top <- head(imp$results[order(-imp$results$importance), ], 20); imp_dt_top<-as.data.table(imp_dt_top)
imp_dt_top<-merge(imp_dt_top, unique(name_id_reviewd[,c(1,4)]), by.x="feature", by.y="pro_id"); imp_dt_top<-imp_dt_top[order(-importance)]
# 3. visiluze top features
ggplot(imp_dt_top, aes(x = reorder(gene_name, importance), y = importance)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95), linewidth = 1.1, width = 0.2, color = "steelblue") +
  coord_flip() + labs(x = "Gene", y = "Permutation Importance") +
  mytheme + theme(axis.text.y = element_text(face = "italic"),legend.position = "none")

## time-resolved DEPs ####
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Protein_ageMean_cor.csv", header = T)
jj2<-jj1[rho_p_age<0.05] # rho_p_age<1 for all protiens

# Step 0: Preset parameters
window_size <- 20; step_size <- 5
age_min <- min(jj2$age); age_max <- max(jj2$age); age_windows <- seq(age_min, age_max - window_size, by = step_size)
# Step 1: protein list (all proteins + total number of proteins in each tissue).
all_proteins <- unique(jj2$pro_id); n_detected_proteins <- length(all_proteins)
time_resolved_DEPs2 <- list() # initialization
# Step 2: Sliding window loop
for (start_age in age_windows) {
  window_label <- paste0("Age_", start_age + window_size/2)
  end_age <- start_age + window_size
  dat_win <- jj2[age >= start_age & age <= end_age] 
  if (length(unique(dat_win$age)) < 5) next 
  
  dep_result <- dat_win[, {
    cor_test <- cor.test(intensity2, age, method = "spearman")
    list(rho = cor_test$estimate, p = cor_test$p.value)}, by = pro_id]
  
  dep_sig <- dep_result[p < 0.05, pro_id] 
  time_resolved_DEPs2[[window_label]] <- unique(dep_sig)
}
# Step 3: Constructing the cumulative DEP curve
cumulative_set <- c(); cumulative_counts2 <- data.frame()

for (i in seq_along(time_resolved_DEPs2)) {
  window_name <-  names(time_resolved_DEPs2)[i]
  current_prots <- time_resolved_DEPs2[[i]]
  
  cumulative_set <- union(cumulative_set, current_prots)
  count <- length(cumulative_set)
  
  cumulative_counts2 <- rbind(cumulative_counts2, data.frame(
    window = window_name,
    midpoint_age = as.numeric(gsub("Age_", "", window_name)),
    cumulative_unique_DEPs = count,
    cumulative_percent = count / n_detected_proteins * 100))
}

# Step 4: plot
library(ggplot2); library(patchwork)

setDT(cumulative_counts2)
cumulative_counts2$count_DEPs<-as.numeric(sapply(time_resolved_DEPs2, length))
cumulative_counts2[, delta_DEPs := c(0, diff(cumulative_unique_DEPs))] # cumulative_unique_DEPs[1]
cumulative_counts2[, delta_cumperc := c(0, diff(cumulative_percent))] # cumulative_percent[1]

scale_factor <- max(cumulative_counts2$delta_DEPs) / max(cumulative_counts2$count_DEPs) # Calculate scaling factor
cumulative_counts2[midpoint_age==min(midpoint_age)]$delta_DEPs<-NA 
ggplot(cumulative_counts2, aes(x = midpoint_age)) +
  geom_line(aes(y = delta_DEPs), color = "firebrick", size = 1, linetype = "solid") +
  geom_point(aes(y = delta_DEPs), color = "firebrick", size = 2) + 
  geom_col(aes(y = count_DEPs * scale_factor), fill = "gray50", alpha = 0.6) +  
  scale_x_continuous(breaks = cumulative_counts2$midpoint_age) +
  scale_y_continuous(name = "Newly emerged ACPs per window", sec.axis = sec_axis(~ . / scale_factor, name = "Total ACPs per window")) +
  labs(x = "Age") + theme_minimal(base_size = 18) +
  theme(text = element_text(face = "bold"),axis.title.y.left = element_text(color = "firebrick"), axis.text.y.left = element_text(color = "firebrick"), 
        axis.title.y.right = element_text(color = "gray50"), axis.text.y.right = element_text(color = "gray50"))


