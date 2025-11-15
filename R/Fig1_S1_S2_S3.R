#### base prepare ####
library(data.table);library(ggplot2);
library(tidyr); library(viridis);library(ggpointdensity);library(ggpubr)

std <- function(x) sd(x,na.rm = T)/sqrt(length(x[!is.na(x)]))
mytheme <- theme_classic(base_size = 22) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=18,face = "bold"))

#### DIA QC: hela ####
## Hela protein overlap ##
aa<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/M-GSGC0290947_13hela_Report_pro.tsv",header = T)
names(aa)[-1:-5]<-gsub("DIA","",gsub("\\[\\d+\\] ","",gsub("_S.+.PG.Quantity","",names(aa)[-1:-5])))
aa[,PG.ProteinAccessions:=gsub(";.+","",PG.ProteinAccessions)]

library("UpSetR")
listInput<-list(aa[!is.na(hela_1)]$PG.ProteinAccessions, aa[!is.na(hela_2)]$PG.ProteinAccessions, aa[!is.na(hela_3)]$PG.ProteinAccessions, 
                aa[!is.na(hela_4)]$PG.ProteinAccessions, aa[!is.na(hela_5)]$PG.ProteinAccessions, aa[!is.na(hela_6)]$PG.ProteinAccessions, 
                aa[!is.na(hela_7)]$PG.ProteinAccessions, aa[!is.na(hela_8)]$PG.ProteinAccessions, aa[!is.na(hela_9)]$PG.ProteinAccessions, 
                aa[!is.na(hela_10)]$PG.ProteinAccessions, aa[!is.na(hela_11)]$PG.ProteinAccessions, aa[!is.na(hela_12)]$PG.ProteinAccessions, aa[!is.na(hela_12)]$PG.ProteinAccessions)
names(listInput)<-paste0("Hela_", seq(1,13,1))
upset(fromList(listInput), order.by = "freq",sets = rev(names(listInput)),keep.order=T,mainbar.y.label = "Intersection of identified proteins", nintersects = 10, 
      sets.x.label = "Identified protein count", point.size = 2, line.size = 1.2, text.scale = c(2, 2, 2, 2, 1.5, 3))

# diag
library(GGally) 
my_lower_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.3, size = 0.1) +
    geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
    theme_minimal()
}

my_diag_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(aes(y = ..density..), binwidth = 0.5, fill = "#4B0082", color = "white", ...) +
    geom_density(color = "red", size = 1) + theme_minimal()
}

tmp<-aa[,-1:-5]
myorder<-paste0("hela_", seq(1,13,1)); setcolorder(tmp, myorder)
tmp[, names(tmp) := lapply(.SD, log10)]
ggpairs(tmp,
        upper = list(continuous = wrap("cor", size = 5, colour = "black", fontface = "bold")),
        lower = list(continuous =wrap(my_lower_lm)),
        diag = list(continuous = wrap(my_diag_hist))) +
  mytheme + theme(axis.text = element_text(size = 12))

#### PRM vs DIA ####
library(data.table); library(ggplot2); library(ggpubr)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=22),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  axis.line = element_line(linewidth = 1.2),axis.ticks=element_line(linewidth=1.2),
  #axis.text.x = element_text(angle =  45,vjust = 0.9, hjust=1)
)

dia<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/DIAvsPRM_DIA蛋白质鉴定列表.xlsx"); dia<-as.data.table(dia)
dia<-dia[,c(2, 6:15)]; names(dia)[1]<-"pro_id"
dia2 <- as.data.table(dia %>% tidyr::separate_rows(pro_id, sep = ";")); dia2 <- dia2[!duplicated(pro_id)]
dia3<-melt.data.table(dia2, id.vars = "pro_id"); names(dia3)<-c("pro_id","sample","value")
dia3[,min:=min(value, na.rm = T),by=.(pro_id)]; dia3<-dia3[!is.infinite(min)]
if (1) {dia3[, value2 := ifelse(is.na(value), min/2, value)]}

prm<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/DIAvsPRM_PRM蛋白质定量结果.xlsx"); prm<-as.data.table(prm)
prm<-prm[,c(1, 15:24)]; names(prm)[1]<-"pro_id"
prm2 <- as.data.table(prm %>% tidyr::separate_rows(pro_id, sep = ";")); prm2 <- prm2[!duplicated(pro_id)]
prm3<-melt.data.table(prm2, id.vars = "pro_id"); names(prm3)<-c("pro_id","sample","value")
prm3[,min:=min(value, na.rm = T),by=.(pro_id)]; prm3<-prm3[!is.infinite(min)]
if (1) {prm3[, value2 := ifelse(is.na(value), min/2, value)]}
prm3[, sample := gsub(".*\\] (\\d+)_.*", "NCN\\1", sample)]; prm3[, sample := sprintf("NCN%03d", as.numeric(gsub("NCN", "", sample)))]

hh<-merge(prm3[,c(1,2,5)], dia3[,c(1,2,5)], all.x=T, by=c("pro_id","sample")); names(hh)[3:4]<-c("prm", "dia")
ggplot(hh, mapping = aes(x = log2(prm), y = log2(dia))) +
  geom_point(color="black", size=3, shape=20, alpha=1)+
  labs(x = 'log2(Protein abundance) by PRM', y = 'log2(Protein abundance) by DIA', color = NULL)+
  geom_smooth(aes(group = sample), method='lm', se=T, formula= y~x, col="red", linetype="dashed")+
  stat_cor(method = "spearman", size = 5, cor.coef.name="rho", label.x.npc="left", label.y.npc="top", label.sep = "\n")+
  facet_wrap(~sample, nrow=2, scales = "free", labeller=labeller(sample = c("NCN004" = "Sample 1","NCN005" = "Sample 2","NCN006" = "Sample 3","NCN007" = "Sample 4","NCN008" = "Sample 5",
                                                                            "NCN009" = "Sample 6","NCN010" = "Sample 7","NCN011" = "Sample 8","NCN012" = "Sample 9","NCN013" = "Sample 10")))+
  mytheme

#### Quality control of full-DIA proteomes ####
## remove sample and protein based on NA number ####
# for protein, identified at least 10 smps; for smp, < 80% median
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
ff <- as.data.table(ff %>% tidyr::separate_rows(protein, sep = ";"))
length(unique(ff$protein)); ncol(ff)-1 # 5956 individual proteins; 248 samples

## overall pheatmap and GO of filtered proteomes ####
qq2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv"); qq2<-qq2[,-2]
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]
qq3<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq3<-as.data.table(qq3)

melt.data.table(qq3,id.vars = "PG.ProteinAccessions")->qq4; names(qq4)<-c("pro_id","sample","intensity")
qq5<-qq4; qq5[grep("^S",sample,invert = T)]->qq5
qq5[,min:=min(intensity,na.rm = T),by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity),min/2,intensity)] # minimum/2填补 NA
qq5[,intensity_log:=ifelse(is.na(intensity),NA,log10(intensity))]
qq5[,intensity_log_median:=median(intensity_log),by=.(pro_id)]
tmp<-dcast(qq5, formula = sample~pro_id, value.var = "intensity_log")

qq6<-unique(qq5[,.(pro_id,intensity_log_median)]); qq6<-qq6[order(intensity_log_median,decreasing = T)]
setcolorder(tmp, c("sample",qq6$pro_id))
pheatmap::pheatmap(tmp[,-1],cluster_rows=1,cluster_cols=0,col=colorRampPalette(c("blue","white","red"))(50), fontsize=20,
                   treeheight_row = 0,treeheight_col = 0,show_rownames = 0,show_colnames = 0)

# top 10% most-abundant proteins GO
library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
tmp3 <- as.data.table(qq6 %>% tidyr::separate_rows(pro_id, sep = ";")); tmp3 <- tmp3[order(intensity_log_median,decreasing = T)];
tmp3[,bin:=ceiling(rank(intensity_log_median)*100/length(intensity_log_median))] 
enrich_go <- enrichGO( gene = tmp3[bin>=90]$pro_id,
                       OrgDb = org.Hs.eg.db, keyType = "UNIPROT", universe= tmp3$pro_id, ont = "ALL",
                       pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
dotplot(enrich_go, x='GeneRatio',showCategory = 10, font.size = 18) 
barplot(enrich_go, showCategory = 20, x = "GeneRatio")

## Keratin, collegen ####
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
KRT_gene<-name_id_reviewd[grep("^KRT",name_id_reviewd$gene_name)]
COL_gene<-name_id_reviewd[grep("^COL",name_id_reviewd$gene_name)]
qq6<-unique(qq5[pro_id!="",c(1,6)])
qq6[,rank_median:=rank(-intensity_log_median)]

qq6[,class:=ifelse(pro_id%in%KRT_gene$pro_id,"KRT", ifelse(pro_id%in%COL_gene$pro_id,"COL","Others"))]
qq6$class<-factor(qq6$class,levels = c("KRT","COL","Others"))
ggplot(qq6, aes(x=rank_median, y=intensity_log_median)) + 
  geom_point(data=qq6[class=="Others"],colour="gray",size=2,shape=20,alpha=0.1)+
  geom_point(data=qq6[class=="KRT"],colour="red",size=3,shape=20,alpha=0.8)+
  geom_point(data=qq6[class=="COL"],colour="blue",size=3,shape=20,alpha=0.8)+
  labs(x = 'Protein abundance rank', y = 'Protein abundance', color = NULL)+
  mytheme

ggplot(qq6, aes(x = class, y = intensity_log_median, color = class))+
  geom_jitter(size=2,alpha = 0.4,position = position_jitter(width = 0.25))+
  geom_boxplot(outlier.shape = NA, fill = '#00000000',width = 0.7, lwd=1.1,position = position_dodge(width = 0.85))+
  scale_color_manual(values = c("red","blue","gray"))+
  labs(x = '', y = 'Protein abundance (log2)', color = NULL)+
  scale_x_discrete(labels=c("KRT" = "Keratins", "COL" = "Collagen","Others" = "Others"))+
  mytheme + theme(legend.position="none",axis.text.x = element_text(face="bold", size=14,angle = 45,vjust = 0.95, hjust=1))
wilcox.test(qq6[class=="KRT"]$intensity_log_median,qq6[class=="Others"]$intensity_log_median)
wilcox.test(qq6[class=="COL"]$intensity_log_median,qq6[class=="Others"]$intensity_log_median)

#  known-correlated KRT in our data: KRT5 vs KRT14  KRT1 vs KRT10 from PMID: 36285538
tmp1<-qq5[pro_id=="P13647"]; tmp2<-qq5[pro_id=="P02533"] # OR tmp1<-qq5[pro_id=="P04264"]; tmp2<-qq5[pro_id=="P13645"]
tmp3<-merge(tmp1[,1:3], tmp2[,1:3], by="sample")
ggplot(tmp3, aes(x = log2(intensity.x), y = log2(intensity.y))) + 
  geom_point(color = "black", size = 1.5, alpha = 0.8) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + labs(x = NULL, y = NULL) + mytheme
cor.test(tmp3$intensity.x, tmp3$intensity.y, method = "pearson") 

#### Our skin proteome vs. PMID40713952 multi-tissue aging proteomes ####
## self data processing: 1/2 minimum impution, VSN normalizaiton, log2 transform ####
qq2 <-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv",header  = T)
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
library(tidyr); qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn); fit <- vsn2(mm); mm1 <- predict(fit, mm)
mm2<-as.data.table(mm1); mm2$pro_id<-rownames(mm1)
mm2_long<-melt(mm2, id.vars = "pro_id"); names(mm2_long)<-c("pro_id","sample","intensity_log")
mm2_long[, intensity:=2^intensity_log]
mm3<-mm2_long[,.(intensity_median=median(intensity,na.rm = T), intensity_mean=mean(intensity,na.rm = T)), by=.(pro_id)] # 同一 tissue 的不同 rep 之间 aggregate
mm3[,c("intensity_median_log", "intensity_mean_log"):=.(ifelse(is.na(intensity_median), NA, log10(intensity_median)), ifelse(is.na(intensity_mean), NA, log10(intensity_mean)))]

## overlap and cor with skin in PMID: 40713952 multi-tissue aging proteomes ####
aa<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/PMID40713952_aging_proteome.xlsx", sheet = 4, skip = 2); aa<-as.data.table(aa) 
mycounts <- colSums(aa[, 4:ncol(aa)] > 0) 
# Skin    Muscle     Heart  Pancreas     Liver     Aorta    Spleen     Lymph Intestine   Adrenal   Adipose      Lung 
# 6198      6690      8797      9734      9864      8929     10491     10642     10399     10420      9023     10824 
aa<-aa[,-2:-3]; aa[aa == 0] <- NA
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
mm4<-merge(mm3, unique(name_id_reviewd[,c(1,4)]), by="pro_id") 
jj<-merge(aa, mm4[,-2:-3], by.x="Protein", by.y="gene_name", all = T);
jj[,pro_id:=NULL]; names(jj)[c(ncol(jj)-1, ncol(jj))]<-c("Skin_self_median","Skin_self_mean")
jj<-jj[order(Protein, -Skin_self_median)]; jj<-jj[!duplicated(Protein)]; names(jj)[1]<-"gene_name"


cols <- setdiff(names(jj), "gene_name"); n <- length(cols)
overlap_mat <- matrix(NA_integer_, nrow = n, ncol = n, dimnames = list(cols, cols))
for (i in seq_along(cols)) {
  for (j in seq_along(cols)) {
    overlap_mat[i, j] <- sum(!is.na(jj[[cols[i]]]) & !is.na(jj[[cols[j]]]))
  }
}

overlap_mat2<-as.data.table(overlap_mat); overlap_mat2[, tissue:=rownames(overlap_mat)]; overlap_mat2<-overlap_mat2[order(-Skin_self_median)]
overlap_mat2<-overlap_mat2[!tissue%in%c("Skin_self_median", "Skin_self_mean", "reference")]; overlap_mat2<-overlap_mat2[order(-Skin_self_median)]

## paired spearman cor ####
library(corrplot)
mycor<-cor(jj[,-1], method = "spearman", use = "pairwise.complete.obs") # spearman 0.73, pearson 0.78
corrplot(mycor, method = "circle",type = "upper",addCoef.col="black", cl.cex=1.5, number.digits=2, number.cex =0.5,col.lim=c(min(mycor),1), is.cor=F, order = "original", tl.col = "black", tl.srt = 45, tl.cex=0.8, mar = c(0,0,0,0))
mycor2<-as.data.table(mycor); mycor2[, tissue:=rownames(mycor)]
mycor2<-mycor2[!tissue%in%c("Skin_self_median", "Skin_self_mean", "reference")]

mycor3<-merge(mycor2[, .(Skin_self_mean, tissue)], overlap_mat2[, .(Skin_self_mean, tissue)], by="tissue") # mean更好
names(mycor3)[2:3]<-c("rho","overlap_num"); mycor3<-mycor3[order(-rho)]
mycor3$tissue<-factor(mycor3$tissue, levels = rev(mycor3$tissue))
## 画图
ggplot(mycor3, aes(x = tissue, y = rho)) + geom_col(fill = "steelblue", width = 0.7) + 
  labs(x = NULL, y = "Spearman rho") + coord_cartesian(ylim = c(0.3, 0.8))+
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 18))

library(viridis); library(ggpointdensity)
ggplot(data = jj, mapping = aes(x = Skin_self_mean, y = Skin)) +
  geom_pointdensity(size = 2, alpha = 0.5) + scale_color_viridis(option = "D") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", linewidth = 1) + 
  labs(x = 'Protein abundance in this study', y = 'Skin protein abundance in Ding et al.', color = NULL)+
  # coord_cartesian(xlim = c(0, 6), ylim = c(0, 6))+
  mytheme

cor.test(jj$Skin_self_mean, jj$Skin, method = "spearman"); cor.test(jj$Skin_self_mean, jj$Skin, method = "pearson")  # 0.7252399; 0.7787051

#### co-expression tSNE for protein pairs ####
library(data.table); library(treeClust); library(WGCNA)
set.seed(42)
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
qq2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv"); qq2<-qq2[,-2]
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]
qq3 <- as.data.table(qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"))

prohd<-copy(qq3)
melt.data.table(prohd,id.vars = "PG.ProteinAccessions")->gg1
names(gg1)<-c("pro_id","sample","intensity")
gg1[,median_inten:=median(intensity,na.rm = T),by=.(pro_id)] # "converted the protein intensities into log2 fold-changes over the median intensity measured for each protein across all cell lines"
gg1[,nor_inten:=log2(intensity/median_inten)] 
gg1[,nonNA_num:=sum(!is.na(intensity)),by=.(pro_id)]
dcast(gg1[nonNA_num>124,],pro_id ~ sample, value.var = "nor_inten")-> prohd_ratios_min95 # 4868 proteins detected in > half samples
prohd_ratios_min95 <- data.frame(prohd_ratios_min95, row.names = "pro_id") 

tc_distances <- treeClust.dist( prohd_ratios_min95, d.num = 2,verbose = TRUE,rcontrol = rpart.control(cp = 0.01),control = treeClust.control(serule = 0) ) # Obtain treeClust distances using default settings
tc_sim_symm <- 1-as.matrix(tc_distances) # Turn the distance matrix into a similarity matrix
adj_mat <- sigmoidAdjacencyFunction(tc_sim_symm, mu = 0.8, alpha = 20) # Calculate the adjacency matrix using the sigmoid function using default settings
tom_sim <- TOMsimilarity( adj_mat, TOMDenom = "mean" ) # Get the Topological Overlap Matrix
colnames(tom_sim) <- colnames(tc_sim_symm)
row.names(tom_sim) <- colnames(tc_sim_symm)

mydt<-c("tc_sim_symm","adj_mat","tom_sim")
mydt2<-c("tc_sim_dt","adj_mat_dt","tc_tom_dt")
for (i in 1:length(mydt)) { # Turn similarity matrices into long-format, remove duplicates and merge into final table
  x<-get(mydt[i])
  y <- as.data.table( reshape2::melt(x) )   
  y <- y[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_value = value ) ]
  y <- y[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl.self-comparisons)
  assign(mydt2[i],y)
}
names(tc_sim_dt)[3]<-"tc_sim";names(adj_mat_dt)[3]<-"tc_adj";names(tc_tom_dt)[3]<-"tc_tom"
tc_dt <- merge( tc_sim_dt, adj_mat_dt, by = c("Protein_1", "Protein_2"))
tc_dt <- merge( tc_dt, tc_tom_dt,by = c("Protein_1", "Protein_2"))
tc_dt_final <- tc_dt[, .(Protein_1, Protein_2, coregulation_score = tc_tom)]
fwrite(tc_dt_final, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/pro_coregulation_scores.csv") # large data are not provided, run above code every time to obtain full data
fwrite(tc_dt_final[1:100000,], "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/pro_coregulation_scores_demo.csv")

# tSNE_map
library(Rtsne); library(gridExtra)
DT <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/pro_coregulation_scores.csv", header  = T)  
DT[, coreg_distance := (1 - log2(coregulation_score)) ]  # Turn co-regulation score back into a distance metric (distance smaller，protein more correlated) and log2-transform for better tSNE performance
DTm <- dcast( data = rbind( DT[, .(Protein_1, Protein_2, coreg_distance)], DT[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coreg_distance)]),
              Protein_1 ~ Protein_2 , value.var = "coreg_distance")  # Turn the melted pairwise table back into a matrix
DTm <- as.data.frame(DTm) ; rownames(DTm) <- DTm$Protein_1 ; DTm$Protein_1 <- NULL
DTm <- as.dist( as.matrix( DTm ))         # Turn into numeric matrix then dist object
protein_IDs <- attr(DTm, "Labels")        # Extract protein IDs from dist object

set.seed(123)
SNE <- Rtsne(DTm, is_distance = TRUE, theta = 0.0, perplexity = 50, max_iter = 1000, verbose = TRUE) # long time
SNE <- as.data.table( SNE$Y );SNE[, ID := protein_IDs ]
SNE2 <- as.data.table(SNE %>% tidyr::separate_rows(ID, sep = ";"))

# Organelle annotation on tSNE map
my_annotations <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/tSNE_map_annotation_PMID31690884.csv") # protein annotation file from PMID31690884
my_annotations[ Organelle == "" , Organelle := NA ]; my_annotations[ Manual_annotation == "" , Manual_annotation := NA ] # Replace empty strings with NAs
my_annotations <- as.data.table(my_annotations %>% tidyr::separate_rows(ID, sep = ";"))
SNE3 <- merge(SNE2, my_annotations, by = "ID", all.x = TRUE )

ggplot(SNE3, aes(x=V1, y=V2))+
  geom_point(data = SNE3[ Organelle == "Cytoplasm"],     shape = 16, size= 1.2, alpha=1, colour = "green")+
  geom_point(data = SNE3[ Organelle == "Mitochondrion"], shape = 16, size= 1.2, alpha=0.8, colour = "red")+
  geom_point(data = SNE3[ Organelle == "ER"],            shape = 16, size= 1.2, alpha=1, colour = "magenta")+
  geom_point(data = SNE3[ Organelle == "Nucleus"],       shape = 16, size= 1.2, alpha=1, colour = "blue")+
  geom_point(data = SNE3[ Organelle == "Nucleolus"],     shape = 16, size= 1.2, alpha=0.8, colour = "purple")+
  geom_point(data = SNE3[ Organelle == "Secreted"],      shape = 16, size= 1.2, alpha=0.8, colour = "cyan")+
  geom_point(data = SNE3[ Organelle == "Ribosome"],      shape = 16, size= 1.2, alpha=0.8, colour = "orange")+
  geom_point(data = SNE3[ is.na(Organelle)],             shape = 16, size= 0.8, alpha=0.2, colour = "gray")+
  labs(x = 'tSNE dimension1', y = 'tSNE dimension2') + mytheme +theme (axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())

#### variation of populational proteome ####
## PCA ####
# proteome  process: 1/2 minimum impution, VSN normalizaiton, log2 transform
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
qq2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv"); qq2<-qq2[,-2]
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5, formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn); fit <- vsn2(mm) ; mm1 <- predict(fit, mm) 
pca_result <- prcomp(t(mm1), center = TRUE, scale. = TRUE) 
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$sample <- rownames(pca_df); pca_df$sample<-gsub("SP", "", pca_df$sample)

# sample information
hh3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Sample_info.csv", header = T)
hh3$age<-as.numeric(hh3$age)
hh3$gender <- factor(hh3$gender, levels=c("female","male"))
hh3$Sun_exposure<-factor(hh3$Sun_exposure, levels=c("Sun_exposed", "Sun_unexposed"))

# PCA plot
pca_df <- merge(pca_df, hh3, by.x = "sample", by.y="No")
pca_df$age_group <- cut(pca_df$age,  breaks = 5,  labels = c("Group1","Group2","Group3","Group4","Group5"))
ggplot(pca_df, aes(x = PC1, y = PC2, shape = gender, color = Sun_exposure, size = age_group)) + 
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Sun_exposed" = "indianred", "Sun_unexposed" = "steelblue")) +
  scale_size_manual(values = seq(2, 6, 1)) +  # 对应 5 级大小
  coord_cartesian(xlim=c(-50, 80), ylim=c(-50,90))+
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
       color = "Sun Exposure", shape = "Gender", size = "Age") +
  guides(shape = guide_legend(override.aes = list(size = 3)), color = guide_legend(override.aes = list(size = 3))) + 
  mytheme + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14))


## Variation partition by ANOVA Eta: variance ~ gender + age + sampling_site ####
library(data.table); library(limma); library(car); set.seed(123)

## lm for loop, use type II sum of squares ANOVA
meta<-as.data.frame(hh3); rownames(meta)<-meta$No
common_samples <- intersect(colnames(mm1), rownames(meta)); length(common_samples) # 248
mm2 <- mm1[, common_samples, drop = FALSE]; meta2 <- meta[common_samples, , drop = FALSE]; 
all(rownames(meta2) == colnames(mm2))  # TRUE

res2 <- lapply(1:nrow(mm2), function(i) {
  fit <- lm(mm2[i, ] ~ age + gender + Sun_exposure, data = meta2) # Protein-by-protein fitting
  ss  <- Anova(fit, type = 2) # Type II ANOVA
  SS_err <- ss["Residuals", "Sum Sq"] # Sum of Squares of Residuals
  pe2 <- function(term) {  # calculate Partial Eta² by the sum of squares
    if (!term %in% rownames(ss)) return(NA_real_)
    SS_eff <- ss[term, "Sum Sq"]
    SS_eff / (SS_eff + SS_err) }
  
  coefs <- coef(fit)
  data.frame(pro_id = rownames(mm2)[i],
             beta_age  = coefs["age"], beta_sexM = coefs["gendermale"], beta_Sun = coefs["Sun_exposureSun_unexposed"],
             p_age  = ss["age", "Pr(>F)"], p_sex  = ss["gender", "Pr(>F)"], p_Sun = ss["Sun_exposure", "Pr(>F)"],
             eta2_age  = pe2("age"), eta2_sex  = pe2("gender"), eta2_Sun = pe2("Sun_exposure"))
}
)

res2 <- do.call(rbind, res2); res2<-as.data.table(res2)
res2[, q_age:=p.adjust(p_age, method = "BH")]; res2[, q_sex:=p.adjust(p_sex, method = "BH")]; res2[, q_Sun:=p.adjust(p_Sun, method = "BH")]
nrow(res2[q_age<0.05]); nrow(res2[q_sex<0.05]); nrow(res2[q_Sun<0.05]) # 174; 8; 8
res2<-merge(res2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
summary(res2$eta2_age) # 0.0000000 0.0007747 0.0032590 0.0082210 0.0097746 0.4229756 
summary(res2$eta2_sex) # 0.0000000 0.0005166 0.0022327 0.0049919 0.0064123 0.2046370
summary(res2$eta2_Sun) # 0.0000000 0.0004835 0.0020969 0.0047260 0.0060475 0.1232338

## plot boxplot
res2_long_eta <- melt(res2, id.vars = "pro_id", measure.vars = c("eta2_age", "eta2_sex", "eta2_Sun"), variable.name = "factor", value.name   = "eta2")
res2_long_eta[, factor := fcase(factor == "eta2_age", "Age", factor == "eta2_sex", "Gender", factor == "eta2_Sun", "Sun exposure")]
res2_long_eta$factor <- factor(res2_long_eta$factor, levels = c("Sun exposure", "Gender", "Age"))
library(scales); res2_long_eta[, eta3:=ifelse(eta2>0.05, 0.05, eta2)]
ggplot(res2_long_eta, aes(x = factor, y = 100*eta2, fill = factor)) +
  geom_violin(width = 1, trim = TRUE, alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.size = 0.7) +
  labs(x = NULL, y = "Percentage of variance (%)") + 
  mytheme+ theme(legend.position = "none")

## case plot
mm3<-as.data.table(mm2); mm3$pro_id<-rownames(mm2)
mm3_long<-melt.data.table(mm3, id.vars = "pro_id"); names(mm3_long)<-c("pro_id", "sample", "intensity")
jj<-merge(mm3_long, hh3, by.x = "sample", by.y="No")
# sun exposed increased 增加表达：SAA1，P0DJI8，PMID: 26900010
wilcox.test(jj[pro_id=="P0DJI8" & Sun_exposure=="Sun_exposed"]$intensity, jj[pro_id=="P0DJI8" & Sun_exposure=="Sun_unexposed"]$intensity)
ggplot(jj[pro_id=="P0DJI8"], aes(x=Sun_exposure, y=intensity))+
  geom_jitter(aes(color = Sun_exposure), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = Sun_exposure, group=factor(Sun_exposure)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("Sun_exposed"="Unexposed", "Sun_unexposed"="Exposed"))+
  scale_color_manual(values = c("Sun_exposed"="steelblue", "Sun_unexposed"="indianred"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none")

# associated with gender: RPS4Y1, P22090, PMID: 33670450,14583743; PZP, P20742, PMID: 32913072, 40360480
wilcox.test(jj[pro_id=="P22090" & gender=="female"]$intensity, jj[pro_id=="P22090" & gender=="male"]$intensity)  
ggplot(jj[pro_id=="P22090"], aes(x=gender, y=intensity))+
  geom_jitter(aes(color = gender), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = gender, group=factor(gender)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("female"="Female", "male"="Male"))+
  scale_color_manual(values = c("female"="indianred", "male"="steelblue"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none")

wilcox.test(jj[pro_id=="P20742" & gender=="female"]$intensity, jj[pro_id=="P20742" & gender=="male"]$intensity) 
ggplot(jj[pro_id=="P20742"], aes(x=gender, y=intensity))+
  geom_jitter(aes(color = gender), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = gender, group=factor(gender)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("female"="Female", "male"="Male"))+
  scale_color_manual(values = c("female"="indianred", "male"="steelblue"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none",axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1))
