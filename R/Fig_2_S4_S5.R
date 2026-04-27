library(data.table);library(ggplot2);
library("tidyr"); library(viridis);library(ggpointdensity);library(ggpubr)

std <- function(x) sd(x,na.rm = T)/sqrt(length(x[!is.na(x)]))
mytheme <- theme_classic(base_size = 22) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5), 
  legend.text=element_text(size=18,face = "bold")
)

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

#### protein with aging ####
hh3<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableSX_sample_information.xlsx",sheet=1); 
hh3<-as.data.table(hh3)

qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
mm<-as.matrix(qq3); 
rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm) 
mm2 <- as.data.table(mm1); 
mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); 
names(mm3)<-c("pro_id","sample","intensity")

jj<-merge(mm3, hh3, by.x="sample", by.y="sample_id"); 
jj<-jj[order(pro_id, age)]
jj[, intensity2:=mean(intensity), by=.(pro_id, age)] 
jj1<-unique(jj[,.(pro_id, intensity2, age)])
jj1[, c("rho_age", "rho_p_age") := .(cor.test(intensity2, age, method = "spearman")$estimate, 
                                     cor.test(intensity2, age, method = "spearman")$p.value), by = .(pro_id)]

## protein vs age: LOESS-fitted trajectories: Fig2A ####
jj1[, intensity2_zs := scale(intensity2), by = .(pro_id)]
jj2<-unique(jj1[,.(pro_id, rho_age, rho_p_age)])

sig_apc<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableSX_protein_age_correlation_adjustedP.csv", header = T)
jj_pos<-jj2[pro_id%in%sig_apc[rho_age>0]$pro_id][order(rho_p_age, -rho_age)] 
jj_neg<-jj2[pro_id%in%sig_apc[rho_age<0]$pro_id][order(rho_p_age, rho_age)] 

tmp1<-jj1[pro_id%in%jj_pos$pro_id,] 
tmp_avg <- tmp1[, .(mean_z = mean(intensity2_zs, na.rm = TRUE)), by = age] 
ggplot() +
  geom_smooth(data = tmp1, 
              aes(x = age, y = intensity2_zs, group = pro_id), 
              method = "loess", se = FALSE, color = "grey60", linewidth = 0.1, alpha = 0.01) +
  geom_smooth(data = tmp_avg, 
              aes(x = age, y = mean_z), 
              method = "lm", se = T, color = "indianred", linewidth = 1.2, fill = "indianred", alpha = 0.2)+ 
  labs(x = "Age", y = "Scaled protein level") + scale_x_continuous(breaks = seq(0, 100, by = 15))+
  mytheme

tmp1<-jj1[pro_id%in%jj_neg$pro_id,] 
tmp_avg <- tmp1[, .(mean_z = mean(intensity2_zs, na.rm = TRUE)), by = age] 
ggplot() +
  geom_smooth(data = tmp1, 
              aes(x = age, y = intensity2_zs, group = pro_id), 
              method = "loess", se = FALSE, color = "grey60", linewidth = 0.1, alpha = 0.01) +
  geom_smooth(data = tmp_avg, 
              aes(x = age, y = mean_z), 
              method = "lm", se = T, color = "indianred", linewidth = 1.2, fill = "indianred", alpha = 0.2)+ 
  labs(x = "Age", y = "Scaled protein level") + scale_x_continuous(breaks = seq(0, 100, by = 15))+
  mytheme

## heatmap protein vs age: Fig. 2B ####
library(ComplexHeatmap); 
library(circlize)
jj_mat <- dcast(jj1, pro_id ~ age, value.var = "intensity2")
mat <- as.matrix(jj_mat[, -1]); 
rownames(mat) <- jj_mat$pro_id
mat_scaled <- t(scale(t(mat))) 
age_order <- sort(as.numeric(colnames(mat_scaled))); 
mat_scaled <- mat_scaled[, as.character(age_order)]; 
colnames(mat_scaled) 

sig_apc<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableSX_protein_age_correlation_adjustedP.csv", header = T)
jj_pos<-jj2[pro_id%in%sig_apc[rho_age>0]$pro_id][order(rho_p_age, -rho_age)] 
jj_neg<-jj2[pro_id%in%sig_apc[rho_age<0]$pro_id][order(rho_p_age, rho_age)]
jj_sig<-rbind(jj_pos, jj_neg)
tmp <- mat_scaled[match(jj_pos$pro_id, rownames(mat_scaled)), ]
tmp <- mat_scaled[match(jj_sig$pro_id, rownames(mat_scaled)), ]

age_labels <- rep("", ncol(tmp)); 
idx_to_label <- seq(1, ncol(tmp), by = 5)
age_labels[idx_to_label] <-  colnames(tmp)[idx_to_label]; 
ht<-Heatmap(tmp, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, 
            column_title = "Age", column_title_side = "bottom", row_title = "Protein level",
            column_names_rot = 45, column_names_gp = gpar(fontface = "bold"), column_labels = age_labels, heatmap_legend_param = list(direction = "vertical"))
draw(ht, heatmap_legend_side = "right")

## protein level vs age: case plot: Fig2C & FigS4A ####
ggplot(jj1[pro_id == "P02452"], aes(x = age, y = intensity2)) + 
  geom_point(color = "grey30", size = 1.5, alpha = 1) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth = 1.2) + 
  labs(x = "Age", y = "log2(Protein abundance)") + 
  mytheme

ggplot(jj1[pro_id == "P02743"], aes(x = age, y = intensity2)) + 
  geom_point(color = "grey30", size = 1.5, alpha = 1) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth = 1.2) + 
  labs(x = "Age", y = "log2(Protein abundance)") + 
  mytheme

ggplot(jj1[pro_id == "Q7L014"], aes(x = age, y = intensity2)) + 
  geom_point(color = "grey30", size = 1.5, alpha = 1) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth = 1.2) + 
  labs(x = "Age", y = "log2(Protein abundance)") + 
  mytheme

## GSEA: Fig. 2D ####
library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr)

jj2<-unique(jj1[,.(pro_id, rho_age, rho_p_age)])
jj3<-merge(jj2, unique(name_id_reviewd[,1:2]), by="pro_id") 
jj3<-jj3[order(gene_id, rho_p_age, rho_age)]; 
jj3<-jj3[!duplicated(gene_id)];
set.seed(123); 
jj3[, rho_jitter:= rho_age + rnorm(nrow(jj3), mean = 0, sd = 1e-6)]; 
jj3<-jj3[!duplicated(rho_jitter)]
jj3<-jj3[order(rho_jitter, decreasing = T)]; 
gene_list <- jj3$rho_jitter; 
names(gene_list) <- jj3$gene_id  

msigdb_human <- msigdbr( species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gene_sets <- msigdb_human[,c("gs_name","ensembl_gene")]; names(gene_sets)<-c("term", "gene")
gene_sets$term <- tolower(gsub("^GOBP_", "", gene_sets$term));  length(unique(gene_sets$term)) # 不带 GOBP_ 前缀; 7583

gsea_result <- GSEA(geneList = gene_list, 
                    exponent = 1, 
                    minGSSize = 10, 
                    maxGSSize = 500, 
                    eps = 0,
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH", 
                    TERM2GENE = gene_sets, 
                    by = "fgsea")
p0 <- dotplot(gsea_result, 
              showCategory = 15, 
              font.size = 16, 
              x="NES", 
              split=".sign", 
              orderBy="NES", 
              color="p.adjust", 
              label_format=50)
df <- p0$data; 
df$logFDR <- -log10(df$p.adjust)
ggplot(df, aes(x = NES, y = gsub("_"," ",Description))) + 
  geom_point(aes(size = setSize, color = logFDR)) +
  facet_grid(.sign ~ ., scales = "free_y", space = "free_y")+
  scale_color_gradient(low = "#56B1F7", high = "#CA0020", breaks = c(2, 4, 6, 8), labels = c("1e-2", "1e-4", "1e-6", "1e-8")) +
  theme_bw(base_size = 18) + 
  theme(text = element_text(color = "black",face = "bold"))+
  labs(y = NULL, color = "p.adjust", x = "NES")

## protein SD vs age: Fig. S4E ####
jj1[, age_group := cut(age, breaks = c(13, 36, 51, 61, 95), labels = c("age1", "age2", "age3", "age4"), include.lowest = TRUE, right = TRUE)]
jj1_cv <- jj1[, .(SD = sd(intensity2, na.rm = TRUE), 
                  CV = sd(intensity2, na.rm = TRUE) / mean(intensity2, na.rm = TRUE)), 
              by = .(pro_id, age_group)]
tmp<-jj1_cv[age_group=="age1"|age_group=="age4"]
tmp$age_group<-factor(tmp$age_group,levels = c("age1","age4"))

wilcox.test(jj1_cv[age_group=="age1"]$SD, jj1_cv[age_group=="age4"]$SD, paired = T)$`p.value` 
ggplot(tmp, aes(x=age_group, y=SD))+
  geom_line(aes(group = pro_id), color = "grey60", linewidth = 0.2, alpha = 0.2) +
  geom_jitter(aes(color = age_group, group = age_group), width = 0.25, size = 0.8, alpha = 0.2) +
  geom_boxplot(aes(group = age_group), width = 0.4,  outlier.shape = NA, fill = NA, color = "black", lwd = 1.1) +
  scale_color_manual(values = c("indianred", "steelblue"))+ 
  labs(x = ' ', y = 'SD of protein level') + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_x_discrete(labels = c("age1" = "Young", "age4" = "Old"))+
  mytheme + 
  theme(legend.position="none")

#### RNA-protein decoulping with aging ####
## RNA vs protein: overall correlation: FigS4B ####
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
mm<-as.matrix(qq3); 
rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm)
mm2 <- as.data.table(mm1); 
mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); 
names(mm3)<-c("pro_id","sample","intensity")

rr_rpkm <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/gene_sample_FPKM.csv")
rr_rpkm2<-melt(rr_rpkm, id.vars = "gene_id"); 
names(rr_rpkm2)<-c("gene_id","sample","value") 
rr_rpkm4<-merge(rr_rpkm2, unique(name_id_reviewd[,1:2]), by="gene_id") 

jj<-merge(mm3, rr_rpkm4[, .(pro_id, sample, value)], by=c("pro_id", "sample")); 
names(jj)[3:4]<-c("pro_raw","rna_raw")
jj<-jj[order(pro_id, sample, -pro_raw, -rna_raw)]; 
jj[,id:=paste(pro_id, sample, sep="_")]; 
jj<-jj[!duplicated(id)]
jj[,c("pro_log","rna_log"):=.(pro_raw, log2(rna_raw+0.001))] 
jj[,c("pro_log_mean", "pro_log_median", "rna_log_mean", "rna_log_median"):=.(
  mean(pro_log), median(pro_log), mean(rna_log), median(rna_log)), 
  by=.(pro_id)]

jj1<-unique(jj[, .(pro_id, pro_log_mean, pro_log_median, rna_log_mean, rna_log_median)])
tmp<-jj1[pro_log_mean>0 & rna_log_mean>0]; 
cor.test(tmp$pro_log_mean, tmp$rna_log_mean, method = "spearman") 

library(viridis); library(ggpointdensity); library(ggpubr)
ggplot(data = tmp, mapping = aes(x = pro_log_mean, y = rna_log_mean)) +
  geom_pointdensity() + 
  scale_color_viridis() +
  labs(x = 'Protein abundance', y = 'RNA expression', color = NULL)+
  geom_smooth(method='lm', se=T,formula= y~x,col="red", linetype="dashed")+  
  mytheme + 
  theme(legend.text = element_text(face = "bold", size = 16))

## RNA-protein cor by age group ####
hh3<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableSX_sample_information.xlsx",sheet=1); 
hh3<-as.data.table(hh3)

jj<-merge(jj, hh3[, c(1:3,8)], by.x="sample", by.y="sample_id"); 
jj<-jj[order(pro_id, age, sample)]
jj[, age_group := cut(age, breaks = c(13, 36, 51, 61, 95), labels = c("age1", "age2", "age3", "age4"), include.lowest = TRUE, right = TRUE)]
table(jj[pro_id=="P05091"]$age_group)   
jj[, c("rho_age", "rho_p_age", "R_age", "R_p_age") := .(cor.test(pro_log, rna_log, method = "spearman")$estimate, 
                                                        cor.test(pro_log, rna_log, method = "spearman")$p.value,
                                                        cor.test(pro_log, rna_log, method = "pearson")$estimate,  
                                                        cor.test(pro_log, rna_log, method = "pearson")$p.value), 
   by = .(pro_id, age_group)]
tmp<-jj[pro_log_mean>0 & rna_log_mean>0];
jj_age<-unique(jj[pro_id%in%tmp$pro_id, .(pro_id, age_group, rho_age, rho_p_age, R_age, R_p_age)]); 
jj_age<-jj_age[!is.na(rho_age)]

## count barplot: Fig2E ####
nrow(jj_age[age_group=="age1"]); 
nrow(jj_age[age_group=="age1" & rho_p_age<0.05]); 
nrow(jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age>0]); 
nrow(jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age<0]) 
nrow(jj_age[age_group=="age4"]); 
nrow(jj_age[age_group=="age4" & rho_p_age<0.05]); 
nrow(jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age>0]); 
nrow(jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age<0]) 
chisq.test(matrix(c(5115, 940, 5115, 526), byrow = T, ncol = 2))$`p.value` # 6.069992e-24

mymat<-data.table(group = rep(c("Young", "Old"), each = 2), 
                  status = rep(c("Positive", "Negative"), times = 2),
                  count = c(940, 86, 526, 72))
mymat$group<-factor(mymat$group, levels = c("Young", "Old")); 
mymat$status<-factor(mymat$status, levels = c("Negative", "Positive"))

ggplot(mymat, aes(x = group, y = count, fill = status)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  labs(x = NULL, y = "Number of genes") + 
  scale_fill_manual(values = c("Positive" = "indianred", "Negative" = "steelblue"), labels = c("Negative cor", "Positive cor")) +
  mytheme + 
  theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank())

## decoupling GO: Fig. S4C & S4D ####
tmp1<-jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age>0] 
tmp2<-jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age>0] 
old_loss<-setdiff(tmp1$pro_id, tmp2$pro_id); 

df <- data.frame(age_group = rep(c("Young", "Old"), each = 2), 
                 status = rep(c("Positive cor", "Positive loss in old"), 2), 
                 count = c(940, 0, 253, 687))
df$age_group <- factor(df$age_group, levels = c("Young", "Old"))

ggplot(df, aes(x = age_group, y = count, fill = status)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual( values = c("Positive cor" = "indianred", "Positive loss in old" = "grey70"), 
                     name = NULL, 
                     labels = c("Couping in young", "Decoupling in Old")) +
  labs(x = NULL, y = "Number of Genes") + 
  mytheme +
  theme(legend.position = "top", 
        legend.direction = "horizontal", 
        legend.text = element_text(size = 16))

library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
enrich_go3 <- enrichGO( old_loss, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "UNIPROT", 
                        universe= unique(jj_age$pro_id), 
                        ont = "ALL", 
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
barplot(enrich_go3, showCategory = 20, x = "GeneRatio", color = "p.adjust")

## rho compare ECDF: Fig. 2F ####
wilcox.test(jj_age[age_group=="age1"]$rho_age, jj_age[age_group=="age4"]$rho_age, paired = TRUE)$`p.value`
ggplot(jj_age[age_group %in% c("age1", "age4")], aes(x = rho_age, color = age_group)) + 
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_color_manual(values = c("age1" = "indianred", "age4" = "steelblue"), 
                     labels = c("age1" = "Young", "age4" = "Old")) +
  labs(x = "Spearman rho", y = "Cumulative Distribution Fraction", color = NULL) + 
  coord_cartesian(xlim = c(-0.3, 0.6)) +
  mytheme + 
  theme(legend.position = c(0.95, 0.05), 
        legend.justification = c("right", "bottom"), 
        legend.key.width = unit(2.5, "lines"), 
        legend.text=element_text(size=18, face = "bold"))

#### age-related correlation changes ####
jj_age<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/RNA_protein_correlation_by_age_group.csv", header = T)
res <- dcast(jj_age[age_group %in% c("age1", "age4")], pro_id ~ age_group, value.var = "rho_age")
res<-merge(res, name_id_reviewd[,.(pro_id, gene_name)], by="pro_id", all.x=T)
res[, delta_r := abs(age4 - age1)]

## select top n proteins: Fig. S5A ####
res <- res[order(-delta_r)]; 
top_ids <- res$pro_id
k_list <- c(5, 8, 10, 13, 15, 20) 
nn2<-t(mm1); 
nn2<-as.data.frame(nn2); 
nn2$sample<-rownames(nn2)

results <- lapply(k_list, function(k) {
  ids <- top_ids[1:k]
  nn3 <- nn2[, c("sample", ids)]  
  nn3<-merge(nn3, hh3[,1:4], by.x="sample", by.y="sample_id", all.x=T)
  formula <- as.formula(paste("age ~", paste(ids, collapse = " + "), "+ sex + sampling_site"))
  fit <- lm(formula, data = nn3)
  list(k = k, 
       R2 = summary(fit)$r.squared, 
       adjR2 = summary(fit)$adj.r.squared, 
       AIC = AIC(fit), 
       BIC = BIC(fit))
}
) 
tmp<-do.call(rbind, lapply(results, as.data.frame)) 

# plot
scale_factor <- (max(tmp$adjR2) - min(tmp$adjR2)) / (max(tmp$AIC) - min(tmp$AIC)) 
ggplot(tmp, aes(x = k)) +
  geom_line(aes(y = adjR2, color = "Adjusted R²"), linewidth = 1.2) +
  geom_point(aes(y = adjR2, color = "Adjusted R²"), size = 2) + 
  geom_line(aes(y = (AIC - min(AIC)) * scale_factor + min(adjR2), color = "AIC"), linetype = "dashed", linewidth = 1.2) +
  geom_point(aes(y = (AIC - min(AIC)) * scale_factor + min(adjR2), color = "AIC"),size = 2) +
  scale_y_continuous(name = "Adjusted R²",
                     sec.axis = sec_axis( ~ (. - min(tmp$adjR2)) / scale_factor + min(tmp$AIC), name = "AIC")) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Adjusted R²" = "steelblue", "AIC" = "indianred")) +
  labs(x = "Number of proteins (k)", color = "") +
  mytheme + 
  theme(axis.title.y.left = element_text(color = "steelblue"), 
        axis.title.y.right = element_text(color = "indianred"))

## multi-protein regression model using these top 10 proteins: Age ∼ Protein1 + Protein2 + ⋯ + Proteink + covariates: Fig. S5B & S5C ####  
nn2<-t(mm1); 
nn2<-as.data.frame(nn2); 
nn2$sample<-rownames(nn2); 
nn2<-as.data.table(nn2)
res <- res[order(-delta_r)]; 
top_ids <- res$pro_id[1:10]; 
nn3 <- nn2[, c("sample", top_ids)]
nn3<-merge(nn3, hh3[,1:4], by.x="sample", by.y="sample_id", all.x=T) 

formula <- as.formula(paste("age ~", paste(top_ids, collapse = " + "), "+ sex + sampling_site"))
fit <- lm(formula, data = nn3)
summary(fit)

# Standardized regression coefficients
library(lm.beta); 
fit_beta <- lm.beta(fit)
dt <- as.data.table(summary(fit_beta)$coefficients, keep.rownames = "Protein")

# plot Fig. S5B
dt1<-merge(dt, name_id_reviewd[,.(pro_id, gene_name)], by.x="Protein", by.y="pro_id")
setorder(dt1, Standardized)
dt1[, gene_name := factor(gene_name, levels = gene_name)]
ggplot(dt1, aes(x = Standardized, y = gene_name)) +
  geom_col(aes(fill = Standardized > 0), width = 0.7) +
  geom_vline(xintercept = 0, color = "black") + 
  geom_text(aes(x = 0, label = gene_name, hjust = ifelse(Standardized > 0, 1.1, -0.1)), 
            size = 6, fontface = "italic") + 
  scale_fill_manual(values = c("TRUE" = "indianred", "FALSE" = "steelblue")) +
  labs(x = "Standardized coefficient (β)", y = NULL) +
  mytheme + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank()) 

# partial R²
library(rsq)
rsq <- rsq.partial(fit)
dt_rsq <- data.table(Protein = rsq$variable, partial_rsq = rsq$partial.rsq) # 10 proteins + sex + sampling_site
dt_rsq1<-merge(dt_rsq, name_id_reviewd[,.(pro_id, gene_name)], by.x="Protein", by.y="pro_id")
setorder(dt_rsq1, partial_rsq)
dt_rsq1[, gene_name := factor(gene_name, levels = gene_name)]

# leave-one-out ΔR^2
full_r2 <- summary(fit)$r.squared
delta_r2 <- sapply(top_ids, function(p) {
  vars <- setdiff(top_ids, p)
  formula_sub <- as.formula( paste("age ~", paste(vars, collapse = " + "), "+ sex + sampling_site"))
  fit_sub <- lm(formula_sub, data = nn3)
  full_r2 - summary(fit_sub)$r.squared })
dt_delta <- as.data.table(delta_r2, keep.rownames = "Protein")
dt_delta1<-merge(dt_delta, name_id_reviewd[,.(pro_id, gene_name)], by.x="Protein", by.y="pro_id")
setorder(dt_delta1, delta_r2)
dt_delta1[, gene_name := factor(gene_name, levels = gene_name)]

# plot Fig. S5C
dt_plot <- merge(dt_rsq1, dt_delta1, by = "gene_name")
setorder(dt_plot, -delta_r2)
dt_plot[, gene_name := factor(gene_name, levels = rev(gene_name))]
dt_long <- melt(dt_plot, id.vars = "gene_name", measure.vars = c("partial_rsq", "delta_r2"), variable.name = "metric", value.name = "value")
ggplot(dt_long, aes(x = gene_name, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("partial_rsq" = "steelblue", "delta_r2" = "indianred"),
                    labels = c("partial_rsq" = expression(Partial~R^2), "delta_r2"   = expression(Delta~R^2)), name = NULL) +
  labs(x = NULL, y = expression(R^2)) +
  mytheme + theme(axis.text.y = element_text(face = "italic"), legend.position = "top")

## multi-protein regression model using these top 10 proteins: Age ∼ Residual + Residual + ⋯ + Residual + covariates: Fig. S5D ####  
jj<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/RNA_protein_correlation_by_age_group_full.csv", header  = T)
tmp<-jj[pro_log_mean>0 & rna_log_mean>0]; 
jj2<-jj[pro_id%in%tmp$pro_id]
jj2[, resid := {
  fit <- lm(pro_log ~ rna_log + sex + sampling_site)
  residuals(fit)
  }, 
  by = pro_id] 

resid_wide <- dcast( jj2, sample + age + sex + sampling_site ~ pro_id, value.var = "resid")
jj_age<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/RNA_protein_correlation_by_age_group.csv", header = T)
res <- dcast(jj_age[age_group %in% c("age1", "age4")], pro_id ~ age_group, value.var = "rho_age")
res<-merge(res, name_id_reviewd[,.(pro_id, gene_name)], by="pro_id", all.x=T)
res[, delta_r := abs(age4 - age1)]; 
res <- res[order(-delta_r)]; 
top_ids <- res$pro_id[1:10] # |Δr| top 10
resid_wide_sub <- resid_wide[, c("sample", "age", "sex", "sampling_site", top_ids), with = FALSE]  

formula_resid <- as.formula(paste("age ~", paste(top_ids, collapse = " + "), "+ sex + sampling_site"))
fit_resid <- lm(formula_resid, data = resid_wide_sub)
summary(fit_resid)

# Standardized regression coefficients
library(lm.beta); 
fit_resid_beta <- lm.beta(fit_resid)
summary(fit_resid_beta) 
dt <- as.data.table(summary(fit_resid_beta)$coefficients, keep.rownames = "Protein")

# plot
dt1<-merge(dt, name_id_reviewd[,.(pro_id, gene_name)], by.x="Protein", by.y="pro_id")
setorder(dt1, Standardized)
dt1[, gene_name := factor(gene_name, levels = gene_name)]
ggplot(dt1, aes(x = Standardized, y = gene_name)) +
  geom_col(aes(fill = Standardized > 0), width = 0.7) +
  geom_vline(xintercept = 0, color = "black") + 
  geom_text(aes(x = 0, label = gene_name, hjust = ifelse(Standardized > 0, 1.1, -0.1)), 
            size = 6, fontface = "italic") + 
  scale_fill_manual(values = c("TRUE" = "indianred", "FALSE" = "steelblue")) +
  labs(x = "Standardized coefficient (β)", y = NULL) +
  mytheme + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank())

#### protemoic clock: Fig. 2G ####
library(mlr3verse); 
library(mlr3tuning); 
library(paradox); 
library(mlr3learners)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/Protein_ageMean_cor.csv", header = T)
expr_mat <- dcast(jj1[rho_p_age < 1e-9], age ~ pro_id, value.var = "intensity2_zs")
expr_mat <- na.omit(expr_mat)  

y <- expr_mat$age; 
X <- as.data.frame(expr_mat[, !c("age")]) 
final_data <- data.table::data.table(X); 
final_data$age <- y
task_final <- TaskRegr$new("age_prediction", backend = final_data, target = "age")
# learner and search space
learner_final <- lrn("regr.glmnet", predict_type = "response")
param_space <- ps(
  alpha = p_dbl(lower = 0, upper = 1),
  lambda = p_dbl(lower = -4, upper = 1, trafo = function(x) 10^x)  # log-scale lambda
)
# Parameter Tuning Examples
instance_final <- TuningInstanceBatchSingleCrit$new(task = task_final, 
                                                    learner = learner_final, 
                                                    resampling = rsmp("cv", folds = 5), 
                                                    measure = msr("regr.mae"), 
                                                    search_space = param_space, 
                                                    terminator = trm("evals", n_evals = 30))
# Tune parameters and set optimal parameters
tuner <- tnr("random_search")
tuner$optimize(instance_final)
learner_final$param_set$values <- instance_final$result_learner_param_vals
# Training the final model
learner_final$train(task_final)
# model performance
pred_all <- learner_final$predict(task_final) 
pred_age_all <- pred_all$response
mae_all <- mean(abs(pred_age_all - y)) 
cor.test(pred_age_all, y, method = "pearson")

tmp<-as.data.table(cbind(y, pred_age_all));
names(tmp)<-c("actual", "predict")
ggplot(tmp, aes(x = actual, y = predict)) + 
  geom_point(color = "grey30", size = 1.5, alpha = 1) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  size = 1.2) + 
  labs(x = "Actual age", y = "Predicted age") + 
  coord_cartesian(xlim=c(13,100), ylim = c(13, 100)) + 
  mytheme

## Identifying important proteins: Fig. 2H ####
library(iml)
predictor <- Predictor$new(model = learner_final, data = X, y = y) 
imp <- FeatureImp$new(predictor, loss = "mae")
imp_dt_top <- head(imp$results[order(-imp$results$importance), ], 20); 
imp_dt_top<-as.data.table(imp_dt_top)
imp_dt_top<-merge(imp_dt_top, unique(name_id_reviewd[,c(1,4)]), by.x="feature", by.y="pro_id"); 
imp_dt_top<-imp_dt_top[order(-importance)]
ggplot(imp_dt_top, aes(x = reorder(gene_name, importance), y = importance)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95), linewidth = 1.1, width = 0.2, color = "steelblue") +
  coord_flip() + 
  labs(x = "Gene", y = "Permutation Importance") +
  mytheme + 
  theme(axis.text.y = element_text(face = "italic"),legend.position = "none")

#### independent validation of proteomic clock: Fig. S5F ####
# Extract the coefficients and intercept of the final model
instance_final$result_learner_param_vals
learner_final$model
coef_mat <- coef(learner_final$model)
coef_df <- data.frame(feature = rownames(as.matrix(coef_mat)),coef = as.numeric(as.matrix(coef_mat)))
coef_df_nonzero <- subset(coef_df, coef != 0)

dia<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/DIA_蛋白质鉴定列表.xlsx"); 
dia<-as.data.table(dia)
dia<-dia[,c(2, 6:15)]; 
names(dia)[1]<-"pro_id"
dia2 <- as.data.table(dia %>% tidyr::separate_rows(pro_id, sep = ";")); 
dia2 <- dia2[!duplicated(pro_id)]
dia3<-melt.data.table(dia2, id.vars = "pro_id"); 
names(dia3)<-c("pro_id","sample","value")
dia4 <- dia3[, .SD[sum(!is.na(value)) >= 5], by = .(pro_id)]
dia4[,min:=min(value, na.rm = T),by=.(pro_id)]; 
dia4<-dia4[!is.infinite(min)]
if (1) {dia4[, intensity := ifelse(is.na(value), min/2, value)]} 

mm<-dcast(dia4,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); 
rownames(mm)<-mm[,1]; 
mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm) 
mm2 <- as.data.table(mm1); 
mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); 
names(mm3)<-c("pro_id","sample","intensity")

hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA 数据验证/sample_info.xlsx", sheet = 2) 
hh2<-as.data.table(hh2); 
hh2<-hh2[,c(1,4)]; 
names(hh2)<-c("sample","age") 
jj<-merge(mm3, hh2, by="sample"); 
jj<-jj[order(pro_id, age)]
jj[, intensity2:=mean(intensity), by=.(pro_id, age)] 
jj1<-unique(jj[,.(pro_id, intensity2, age)])

jj1[, intensity2_zs := scale(intensity2), by = .(pro_id)]
jj2 <- dcast(jj1, age ~ pro_id, value.var = "intensity2_zs")
feature_coef <- coef_df_nonzero[coef_df_nonzero$feature != "(Intercept)", ]
intercept <- coef_df_nonzero$coef[coef_df_nonzero$feature == "(Intercept)"]
cols <- feature_coef$feature
ext_data <- jj2[, c("age", ..cols)]
all(colnames(ext_data[,-1]) == feature_coef$feature) 

pred_age_ext <- intercept + as.matrix(ext_data[, -1]) %*% feature_coef$coef
ext_data$predicted_age <- pred_age_ext
ggplot(ext_data, aes(x = age, y = predicted_age)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "indianred", size = 1, linetype = "dashed",linewidth = 1, alpha = 0.2) + 
  coord_cartesian(xlim = c(30, 80), ylim = c(30, 80)) +
  mytheme + 
  labs(x = "Actual age", y = "Predicted age")
