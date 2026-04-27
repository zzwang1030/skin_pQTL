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

#### process clinical information ####
hh<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS1_sample_information.xlsx", sheet=1, n_max = 264, col_names = T); hh<-as.data.table(hh)
nrow(hh[proteomics==1]);
# 248
summary(hh[proteomics==1]$age); 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 13.00   38.00   52.00   50.08   62.00   95.00 
table(hh[proteomics==1]$sex); 
# female   male = 96(38.7%), 152(61.3%); 248
table(hh[proteomics==1]$sampling_site)
# lower_limb      trunk upper_limb  = 123         97         28 

## plot: Fig. S1A ####
hh3<-hh[proteomics==1]
ggplot(hh3, aes(x=age)) + geom_histogram(bins=20, fill="steelblue",color='grey90', linewidth=0.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(x ="Age", y="Numbers of individuals") + mytheme

## plot: Fig. S1B ####
mypie<-as.data.table(table(hh3$sex)); names(mypie)<-c("category","value")
ggplot(mypie, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1, color = "white") + coord_polar(theta = "y") + 
  scale_fill_manual(values = c("male" = "steelblue","female" = "indianred")) +
  theme_void() + theme(legend.position = "none")

#### DIA QC: hela ####
aa<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/M-GSGC0290947_13hela_Report_pro.tsv",header = T)
names(aa)[-1:-5]<-gsub("DIA","",gsub("\\[\\d+\\] ","",gsub("_S.+.PG.Quantity","",names(aa)[-1:-5])))
aa[,PG.ProteinAccessions:=gsub(";.+","",PG.ProteinAccessions)]

## Hela protein overlap: Fig. S1E ####
library("UpSetR")
listInput<-list(aa[!is.na(hela_1)]$PG.ProteinAccessions, aa[!is.na(hela_2)]$PG.ProteinAccessions, aa[!is.na(hela_3)]$PG.ProteinAccessions, 
                aa[!is.na(hela_4)]$PG.ProteinAccessions, aa[!is.na(hela_5)]$PG.ProteinAccessions, aa[!is.na(hela_6)]$PG.ProteinAccessions, 
                aa[!is.na(hela_7)]$PG.ProteinAccessions, aa[!is.na(hela_8)]$PG.ProteinAccessions, aa[!is.na(hela_9)]$PG.ProteinAccessions, 
                aa[!is.na(hela_10)]$PG.ProteinAccessions, aa[!is.na(hela_11)]$PG.ProteinAccessions, aa[!is.na(hela_12)]$PG.ProteinAccessions, aa[!is.na(hela_12)]$PG.ProteinAccessions)
names(listInput)<-paste0("Hela_", seq(1,13,1))
upset(fromList(listInput), order.by = "freq",sets = rev(names(listInput)),keep.order=T,mainbar.y.label = "Intersection of identified proteins", nintersects = 10, 
      sets.x.label = "Identified protein count", point.size = 2, line.size = 1.2, text.scale = c(2, 2, 2, 2, 1.5, 3))

## correlation of protein abundace: Fig. S1D ####
library(corrplot)
mycor<-cor(aa[,-1:-5], method = "pearson",use = "pairwise.complete.obs")
library(GGally) 
my_lower_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.3, size = 0.1) +
    geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
    theme_minimal()
  }

my_diag_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(aes(y = ..density..),
                   binwidth = 0.5,
                   fill = "#4B0082", 
                   color = "white", ...) +
    geom_density(color = "red", size = 1) +
    theme_minimal()
  }

tmp<-aa[,-1:-5]
myorder<-paste0("hela_", seq(1,13,1)); setcolorder(tmp, myorder) 
tmp[, names(tmp) := lapply(.SD, log10)]
ggpairs(tmp,
        upper = list(continuous = wrap("cor", size = 5, colour = "black", fontface = "bold")),  
        lower = list(continuous =wrap(my_lower_lm)),  
        diag = list(continuous = wrap(my_diag_hist)))+  
  mytheme + 
  theme(axis.text = element_text(size = 12))


#### PRM vs DIA ####
## quantification correlation: Fig. S2B ####
hh<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS11_PRM_vs_DIA_correlation.csv", header = T)
ggplot(hh, mapping = aes(x = log2(PRM), y = log2(DIA))) +
  geom_point(color="black", size=3, shape=20, alpha=1)+
  labs(x = 'log2(Protein abundance) by PRM', y = 'log2(Protein abundance) by DIA', color = NULL)+
  geom_smooth(aes(group = Sample_ID), method='lm', se=T, formula= y~x, col="red", linetype="dashed")+
  stat_cor(method = "spearman", size = 5, cor.coef.name="rho", label.x.npc="left", label.y.npc="top", label.sep = "\n")+
  facet_wrap(~Sample_ID, nrow=2, scales = "free")+
  mytheme

#### chraracterization of skin proteomes ####
## pheatmap: Fig. 1B ####
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
qq5<-melt.data.table(qq3,id.vars = "pro_id"); 
names(qq5)<-c("pro_id","sample","intensity")
qq5[,intensity_log:=ifelse(is.na(intensity),NA,log10(intensity))]
qq5[,intensity_log_median:=median(intensity_log),by=.(pro_id)]
summary(qq5$intensity_log) # -1.552   1.318   1.830   1.754   2.247   5.957  
tmp<-dcast(qq5, formula = sample~pro_id, value.var = "intensity_log")

qq6<-unique(qq5[,.(pro_id,intensity_log_median)]); 
qq6<-qq6[order(intensity_log_median,decreasing = T)]
setcolorder(tmp, c("sample",qq6$pro_id))
pheatmap::pheatmap(tmp[,-1], 
                   cluster_rows=1, 
                   cluster_cols=0, 
                   col=colorRampPalette(c("blue","white","red"))(50), 
                   fontsize=20,
                   treeheight_row = 0, 
                   treeheight_col = 0, 
                   show_rownames = 0, 
                   show_colnames = 0)

## top 10%  GO: Fig. S3A ####
library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
tmp3 <- as.data.table(qq6 %>% tidyr::separate_rows(pro_id, sep = ";")); 
tmp3 <- tmp3[order(intensity_log_median,decreasing = T)];
tmp3[,bin:=ceiling(rank(intensity_log_median)*100/length(intensity_log_median))] 
enrich_go <- enrichGO(gene = tmp3[bin>=90]$pro_id, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "UNIPROT", 
                      universe= tmp3$pro_id,
                      ont = "ALL", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1, 
                      qvalueCutoff = 1, 
                      readable = TRUE)
dotplot(enrich_go, x='GeneRatio',showCategory = 10, font.size = 18) 
aa<-as.data.table(enrich_go); 

## Keratin, collegen: Fig. 1D ####
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
qq5<-melt.data.table(qq3,id.vars = "pro_id"); 
qq6<-merge(qq5, unique(name_id_reviewd[,.(pro_id, gene_name)]), by="pro_id", all.x = T)
names(qq6)<-c("pro_id", "sample", "intensity", "gene_id")
qq6[, intensity_gene:=sum(intensity[!is.na(intensity)]), by=.(sample, gene_id)] 
qq6<-unique(qq6[,c(-1,-3)]); 
qq6[,c("mean_intensity_gene", "std_intensity_gene", "median_intensity_gene"):= 
      .(mean(intensity_gene), 
        std(intensity_gene), 
        median(intensity_gene)),
    by=.(gene_id)]
qq6[,c("mean_intensity_gene_log", "std_intensity_gene_log", "median_intensity_gene_log"):= 
      .(log2(mean_intensity_gene), 
        log2(std_intensity_gene), 
        log(median_intensity_gene)),
    by=.(gene_id)]
qq7<-unique(qq6[gene_id!="",c(-1,-3)])
qq7[,c("rank_mean" ,"rank_median"):=
      .(rank(-mean_intensity_gene_log), 
        rank(-median_intensity_gene_log))]
qq7[,class:=ifelse(grepl("KRT",gene_id), "KRT", 
                   ifelse(grepl("COL",gene_id), "COL", "Others"))]
qq7$class<-factor(qq7$class,levels = c("KRT","COL","Others"))

ggplot(qq7, aes(x=rank_mean, y=mean_intensity_gene_log)) + 
  geom_point(data=qq7[class=="Others"],colour="gray",size=2,shape=20,alpha=0.1)+
  geom_point(data=qq7[class=="KRT"],colour="red",size=3,shape=20,alpha=0.8)+
  geom_point(data=qq7[class=="COL"],colour="blue",size=3,shape=20,alpha=0.8)+
  labs(x = 'Protein abundance rank', y = 'Protein abundance', color = NULL)+
  mytheme

ggplot(qq7, aes(x = class, y = mean_intensity_gene_log,color = class))+
  geom_jitter(size=2,alpha = 0.4,position = position_jitter(width = 0.25))+
  geom_boxplot(outlier.shape = NA, fill = '#00000000',width = 0.7, 
               lwd=1.1,position = position_dodge(width = 0.85))+
  scale_color_manual(values = c("red","blue","gray"))+
  labs(x = '', y = 'Protein abundance (log2)', color = NULL)+
  scale_x_discrete(labels=c("KRT" = "Keratins", "COL" = "Collagen","Others" = "Others"))+
  mytheme + theme(legend.position="none", 
                  axis.text.x = element_text(face="bold", size=14,angle = 45,vjust = 0.95, hjust=1))

wilcox.test(qq7[class=="KRT"]$mean_intensity_gene_log,qq7[class=="Others"]$mean_intensity_gene_log)
wilcox.test(qq7[class=="COL"]$mean_intensity_gene_log,qq7[class=="Others"]$mean_intensity_gene_log)

## KRT5 vs KRT14: Fig. 1E ####
tmp1<-qq6[gene_id=="KRT5"]; tmp2<-qq6[gene_id=="KRT14"]
tmp3<-merge(tmp1[,1:3], tmp2[,1:3], by="sample")
ggplot(tmp3, aes(x = log2(intensity_gene.x), y = log2(intensity_gene.y))) + 
  geom_point(color = "black", size = 1.5, alpha = 0.8) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + labs(x = NULL, y = NULL) + 
  mytheme
cor.test(tmp3$intensity_gene.x, tmp3$intensity_gene.y, method = "pearson") 

## KRT1 vs KRT10: Fig. S3B ####
tmp1<-qq6[gene_id=="KRT1"]; tmp2<-qq6[gene_id=="KRT10"]
tmp3<-merge(tmp1[,1:3], tmp2[,1:3], by="sample")
ggplot(tmp3, aes(x = log2(intensity_gene.x), y = log2(intensity_gene.y))) + 
  geom_point(color = "black", size = 1.5, alpha = 0.8) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + labs(x = NULL, y = NULL) + 
  mytheme
cor.test(tmp3$intensity_gene.x, tmp3$intensity_gene.y, method = "pearson") 

## PKP1 vs DESP: Fig. S3C ####
tmp1<-qq6[gene_id=="PKP1"]; tmp2<-qq6[gene_id=="DSP"]
tmp3<-merge(tmp1[,1:3], tmp2[,1:3], by="sample")
ggplot(tmp3, aes(x = log2(intensity_gene.x), y = log2(intensity_gene.y))) + 
  geom_point(color = "black", size = 1.5, alpha = 0.8) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + labs(x = NULL, y = NULL) + 
  mytheme
cor.test(tmp3$intensity_gene.x, tmp3$intensity_gene.y, method = "pearson") 

## nonNA for each protein: Fig. S2D ####
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_NOimputed.csv",header = T)
qq4<-melt.data.table(qq3,id.vars = "PG.ProteinAccessions"); 
names(qq4)<-c("pro_id","sample","intensity")
qq4[, nonNA_num_pro:=sum(!is.na(intensity)), by=.(pro_id)]; 
tmp1<-unique(qq4[,.(pro_id, nonNA_num_pro)])
summary(tmp1$nonNA_num_pro)  # 10.0   167.0   241.0   200.9   248.0   248.0
ggplot(tmp1, aes(x = nonNA_num_pro)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(x = "Number of samples in which the protein was quantified", y = "Protein count") +
  mytheme

## nonNA for each samples: Fig. S2E ####
qq4[, nonNA_num_smp:=sum(!is.na(intensity)), by=.(sample)]; 
tmp2<-unique(qq4[,.(sample, nonNA_num_smp)])
summary(tmp2$nonNA_num_smp)
tmp2<-tmp2[order(nonNA_num_smp,decreasing = T)]; 
tmp2[,sample_id:=1:.N]
ggplot(tmp2, aes(x = "", y = nonNA_num_smp)) + 
  geom_boxplot(width = 0.3, fill = 'grey90', color = "black", lwd=1.1 , fatten = 2) +  
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = "steelblue") +  
  labs(x = "", y = "Number of quantified proteins per sample") +
  mytheme

#### skin proteome comapare with other study PMID: 40713952 ####
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
mm<-as.matrix(qq3); 
rownames(mm)<-mm[,1]; 
mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm)
mm1 <- predict(fit, mm)
mm2<-as.data.table(mm1); 
mm2$pro_id<-rownames(mm1)
mm2_long<-melt(mm2, id.vars = "pro_id"); 
names(mm2_long)<-c("pro_id","sample","intensity_log")
mm2_long[, intensity:=2^intensity_log] 
mm3<-mm2_long[,.(intensity_median = median(intensity,na.rm = T), 
                 intensity_mean = mean(intensity,na.rm = T)), 
              by=.(pro_id)]
mm3[,c("intensity_median_log", "intensity_mean_log"):=.(ifelse(is.na(intensity_median), NA, log10(intensity_median)), ifelse(is.na(intensity_mean), NA, log10(intensity_mean)))]
mm4<-merge(mm3, unique(name_id_reviewd[,c(1,4)]), by="pro_id") 

aa<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/mmc1_from_PMID40713952.xlsx", sheet = 4, skip = 2); 
aa<-as.data.table(aa) 
aa<-aa[,-2:-3]; 
aa[aa == 0] <- NA

jj<-merge(aa, mm4[,-2:-3], by.x="Protein", by.y="gene_name", all = T);
jj[,pro_id:=NULL]; 
names(jj)[c(ncol(jj)-1, ncol(jj))]<-c("Skin_self_median","Skin_self_mean")
jj<-jj[order(Protein, -Skin_self_median)]; 
jj<-jj[!duplicated(Protein)]; 
names(jj)[1]<-"gene_name"

## paired-tissue spearman cor: Fig. 1C-upper ####
library(corrplot)
mycor<-cor(jj[,-1], method = "spearman", use = "pairwise.complete.obs") # spearman 0.73, pearson 0.78
mycor2<-as.data.table(mycor); 
mycor2[, tissue:=rownames(mycor)]
mycor2<-mycor2[!tissue%in%c("Skin_self_median", "Skin_self_mean", "reference")]
mycor2<-mycor2[order(-Skin_self_mean)]
mycor2$tissue<-factor(mycor2$tissue, levels = rev(mycor2$tissue))

ggplot(mycor2, aes(x = tissue, y = Skin_self_mean)) + 
  geom_col(fill = "steelblue", width = 0.7) + 
  labs(x = NULL, y = "Spearman rho") +
  coord_cartesian(ylim = c(0.3, 0.8))+
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 18))

## our skin vs punlished skin: Fig. 1C-lower ####
library(viridis); library(ggpointdensity)
ggplot(data = jj, mapping = aes(x = Skin_self_mean, y = Skin)) +
  geom_pointdensity(size = 2, alpha = 0.5) + scale_color_viridis(option = "D") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", linewidth = 1) + 
  labs(x = 'Protein abundance in this study', y = 'Skin protein abundance in Ding et al.', color = NULL)+
  mytheme

#### co-regulation tSNE_map:  Fig. 1F ####
qq2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_NOimputed.csv",header = T)
prohd <- copy(qq2)
prohd <- data.frame(prohd, row.names = 1)
feature_count <- apply(prohd, 1, function(x) {sum(!is.na(x))})
prohd_ratios_min95 <- t(prohd[feature_count >= 50,]) 

library(treeClust); library(WGCNA)
set.seed(42)
prohd <- copy(qq2)
gg1<-melt.data.table(prohd, id.vars = "PG.ProteinAccessions")
names(gg1) <- c("pro_id","sample","intensity")
gg1[, median_inten := median(intensity,na.rm = T), by=.(pro_id)] 
gg1[,nor_inten := log2(intensity/median_inten)] 
gg1[,nonNA_num := sum(!is.na(intensity)),by=.(pro_id)]
prohd_ratios_min95 <- dcast(gg1[nonNA_num>124,], pro_id ~ sample, value.var = "nor_inten") 
prohd_ratios_min95 <- data.frame(prohd_ratios_min95, row.names = "pro_id") 
# Obtain treeClust distances using default settings
tc_distances <- treeClust.dist(prohd_ratios_min95, 
                               d.num = 2, 
                               verbose = TRUE, 
                               rcontrol = rpart.control(cp = 0.01), 
                               control = treeClust.control(serule = 0)) 
tc_sim_symm <- 1-as.matrix(tc_distances) # Turn the distance matrix into a similarity matrix
adj_mat <- sigmoidAdjacencyFunction(tc_sim_symm, mu = 0.8, alpha = 20) # Calculate the adjacency matrix using the sigmoid function using default settings
tom_sim <- TOMsimilarity( adj_mat, TOMDenom = "mean" ) # Get the Topological Overlap Matrix
colnames(tom_sim) <- colnames(tc_sim_symm)
row.names(tom_sim) <- colnames(tc_sim_symm)

mydt <- c("tc_sim_symm","adj_mat","tom_sim")
mydt2 <- c("tc_sim_dt","adj_mat_dt","tc_tom_dt")
for (i in 1:length(mydt)) { # Turn similarity matrices into long-format, remove duplicates and merge into final table
  x<-get(mydt[i])
  y <- as.data.table(x, keep.rownames = "Var1")
  y <- melt(y, id.vars = "Var1", variable.name = "Var2", value.name = "tc_value")
  y <- y[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_value = tc_value ) ]
  y <- y[ Protein_1 > Protein_2 ]
  assign(mydt2[i],y)
}

names(tc_sim_dt)[3]<-"tc_sim"; 
names(adj_mat_dt)[3]<-"tc_adj"; 
names(tc_tom_dt)[3]<-"tc_tom"
tc_dt <- merge( tc_sim_dt, adj_mat_dt, by = c("Protein_1", "Protein_2"))
tc_dt <- merge( tc_dt, tc_tom_dt,by = c("Protein_1", "Protein_2"))
tc_dt_final <- tc_dt[, .(Protein_1, Protein_2, coregulation_score = tc_tom)]

library(Rtsne); library(gridExtra)
DT <- copy(tc_dt_final)
DT[, coreg_distance := (1 - log2(coregulation_score)) ]  # Turn co-regulation score back into a distance metric and log2-transform for better tSNE performance
tmp<-merge(DT, unique(name_id_reviewd[,c(1,4)]), by.x="Protein_1", by.y="pro_id", all.x=T, allow.cartesian = TRUE)

DTm <- dcast( data = rbind( DT[, .(Protein_1, Protein_2, coreg_distance)],
                            DT[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coreg_distance)]),
              Protein_1 ~ Protein_2 , value.var = "coreg_distance")  # Turn the melted pairwise table back into a matrix

DTm <- as.data.frame(DTm); 
rownames(DTm) <- DTm$Protein_1; 
DTm$Protein_1 <- NULL
DTm <- as.dist(as.matrix( DTm ))         # Turn into numeric matrix then dist object
protein_IDs <- attr(DTm, "Labels")        # Extract protein IDs from dist object

set.seed(123)
SNE <- Rtsne(DTm, is_distance = TRUE, theta = 0.0, perplexity = 50, max_iter = 1000, verbose = TRUE) 
SNE <- as.data.table( SNE$Y ); SNE[, ID := protein_IDs ]
SNE2 <- as.data.table(SNE %>% tidyr::separate_rows(ID, sep = ";"))
ggplot(SNE2, aes(x = V1, y = V2))+
  geom_point(shape = 16, size= 1, alpha=0.5)+labs(x = 'tSNE dimension1', y = 'tSNE dimension2')+
  geom_point(data = SNE2[ ID == "P15924"], shape = 3, size= 4, alpha=1, colour = "red")+ # DESP
  geom_point(data = SNE2[ ID == "Q13835"], shape = 1, size= 4, alpha=1, colour = "blue")+ # PKP1
  mytheme+theme (axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())

# Organelle annotation on tSNE map
my_annotations <- fread("~/Desktop/省皮/project/pQTL/NBT2020 co-regulation map/tSNE_map_annotation.csv") # Read in annotation file
my_annotations[ Organelle == "" , Organelle := NA ]; 
my_annotations[ Manual_annotation == "" , Manual_annotation := NA ] # Replace empty strings with NAs
my_annotations <- as.data.table(my_annotations %>% tidyr::separate_rows(ID, sep = ";"))
SNE3 <- merge(SNE2, my_annotations, by = "ID", all.x = TRUE )

ggplot(SNE3, aes(x=V1, y=V2))+
  geom_point(data = SNE3[ is.na(Organelle)],             shape = 16, size= 0.8, alpha=0.2, colour = "gray")+
  geom_point(data = SNE3[ Organelle == "Cytoplasm"],     shape = 16, size= 1.2, alpha=1, colour = "green")+
  geom_point(data = SNE3[ Organelle == "Mitochondrion"], shape = 16, size= 1.2, alpha=0.8, colour = "red")+
  geom_point(data = SNE3[ Organelle == "ER"],            shape = 16, size= 1.2, alpha=1, colour = "magenta")+
  geom_point(data = SNE3[ Organelle == "Nucleus"],       shape = 16, size= 1.2, alpha=1, colour = "cyan")+
  geom_point(data = SNE3[ Organelle == "Nucleolus"],     shape = 16, size= 1.2, alpha=0.8, colour = "purple")+
  geom_point(data = SNE3[ Organelle == "Secreted"],      shape = 16, size= 1.2, alpha=0.8, colour = "blue")+
  geom_point(data = SNE3[ Organelle == "Ribosome"],      shape = 16, size= 1.2, alpha=0.8, colour = "orange")+
  geom_point(data = SNE2[ ID == "P15924"], shape = 3, size= 4, alpha=1, colour = "red")+ # DESP
  geom_point(data = SNE2[ ID == "Q13835"], shape = 1, size= 4, alpha=1, colour = "blue")+ # PKP1
  labs(x = 'tSNE dimension1', y = 'tSNE dimension2') + 
  mytheme +theme (axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())

#### inter-individual variance ####
## PCA: Fig. S3D ####
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
mm<-copy(qq3)
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm)

pca_result <- prcomp(t(mm1), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$sample <- rownames(pca_df); pca_df$sample<-gsub("SP", "", pca_df$sample)

pca_df <- merge(pca_df, hh3, by.x = "sample", by.y="No")
pca_df$age_group <- cut(pca_df$age,  breaks = 5,  labels = c("Group1","Group2","Group3","Group4","Group5"))
ggplot(pca_df, aes(x = PC1, y = PC2, shape = gender, color = Sun_exposure, size = age_group)) + 
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Sun_exposed" = "indianred", "Sun_unexposed" = "steelblue")) +
  scale_size_manual(values = seq(2, 6, 1)) +  
  coord_cartesian(xlim=c(-50, 80), ylim=c(-50,90))+
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
       color = "Sun Exposure", shape = "Gender", size = "Age") +
  guides(shape = guide_legend(override.aes = list(size = 3)), color = guide_legend(override.aes = list(size = 3))) + 
  mytheme + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14))

## PCA_site: Fig. S3E ####
hh3<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableSX_sample_information.xlsx",sheet=1); 
hh3<-as.data.table(hh3)
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
mm<-copy(qq3)
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm)

meta<-as.data.frame(hh3); rownames(meta)<-meta$sample_id
common_samples <- intersect(colnames(mm1), rownames(meta)); length(common_samples) # 248
mm2 <- mm1[, common_samples, drop = FALSE]; meta2 <- meta[common_samples, , drop = FALSE]
all(rownames(meta2) == colnames(mm2))  # TRUE

library(limma)
expr_res <- removeBatchEffect(mm2, covariates = model.matrix(~ age + sex, data = meta2)[, -1]) 
pca_result <- prcomp(t(expr_res), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$sample <- rownames(pca_df); pca_df$sample<-gsub("SP", "", pca_df$sample)

pca_df <- merge(pca_df, hh3, by.x = "sample", by.y="sample_id", all.x = TRUE)
ggplot(pca_df, aes(x = PC1, y = PC2, color = sampling_site)) + 
  geom_point(alpha = 1, size=1.5) +
  scale_color_manual(values = c("lower_limb" = "indianred", "upper_limb" = "steelblue", "trunk" = "darkolivegreen3")) +
  #coord_cartesian(xlim=c(-50, 80), ylim=c(-50,90))+
  labs(x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),  color = "Sampling site") +
  mytheme + theme(legend.text = element_text(size = 16), legend.title = element_text(size = 16))

## variance partition: Fig. 1G ####
library(data.table); library(limma); library(car); set.seed(123)

hh3<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableSX_sample_information.xlsx",sheet=1); 
hh3<-as.data.table(hh3)
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)
mm<-copy(qq3)
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm)

meta<-as.data.frame(hh3); rownames(meta)<-meta$No
common_samples <- intersect(colnames(mm1), rownames(meta)); length(common_samples) # 248
mm2 <- mm1[, common_samples, drop = FALSE]; meta2 <- meta[common_samples, , drop = FALSE]; 
all(rownames(meta2) == colnames(mm2))  # TRUE

res2 <- lapply(1:nrow(mm2), function(i) {
  fit <- lm(mm2[i, ] ~ age + gender + Sun_exposure, data = meta2) 
  ss  <- Anova(fit, type = 2) 
  SS_err <- ss["Residuals", "Sum Sq"] 
  pe2 <- function(term) { 
    if (!term %in% rownames(ss)) return(NA_real_)
    SS_eff <- ss[term, "Sum Sq"]
    SS_eff / (SS_eff + SS_err) }
  
  coefs <- coef(fit)
  data.frame(pro_id = rownames(mm2)[i],
             beta_age  = coefs["age"], beta_sexM = coefs["gendermale"], beta_Sun = coefs["Sun_exposureSun_unexposed"],
             p_age  = ss["age", "Pr(>F)"], p_sex  = ss["gender", "Pr(>F)"], p_Sun = ss["Sun_exposure", "Pr(>F)"],
             eta2_age  = pe2("age"), eta2_sex  = pe2("gender"), eta2_Sun = pe2("Sun_exposure"))
})

res2 <- do.call(rbind, res2); res2<-as.data.table(res2)
res2[, q_age:=p.adjust(p_age, method = "BH")]; res2[, q_sex:=p.adjust(p_sex, method = "BH")]; res2[, q_Sun:=p.adjust(p_Sun, method = "BH")]
nrow(res2[q_age<0.05]); nrow(res2[q_sex<0.05]); nrow(res2[q_Sun<0.05]) 
res2<-merge(res2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")

res2_long_eta <- melt(res2, id.vars = "pro_id", measure.vars = c("eta2_age", "eta2_sex", "eta2_Sun"), variable.name = "factor", value.name   = "eta2")
res2_long_eta[, factor := fcase(factor == "eta2_age", "Age", factor == "eta2_sex", "Gender", factor == "eta2_Sun", "Sun exposure")]
res2_long_eta$factor <- factor(res2_long_eta$factor, levels = c("Sun exposure", "Gender", "Age"))
library(scales); res2_long_eta[, eta3:=ifelse(eta2>0.05, 0.05, eta2)]
ggplot(res2_long_eta, aes(x = factor, y = 100*eta2, fill = factor)) +
  geom_violin(width = 1, trim = TRUE, alpha = 0.3) +
  geom_boxplot(width = 0.1, outlier.size = 0.7) +
  labs(x = NULL, y = "Percentage of variance (%)") + 
  # coord_cartesian(ylim = c(0, 2))+
  mytheme+ theme(legend.position = "none")

## SAA case plot
mm3<-as.data.table(mm2); mm3$pro_id<-rownames(mm2)
mm3_long<-melt.data.table(mm3, id.vars = "pro_id"); names(mm3_long)<-c("pro_id", "sample", "intensity")
jj<-merge(mm3_long, hh3[,c(1,2,3,11)], by.x = "sample", by.y="No")
wilcox.test(jj[pro_id=="P0DJI8" & Sun_exposure=="Sun_exposed"]$intensity, jj[pro_id=="P0DJI8" & Sun_exposure=="Sun_unexposed"]$intensity)  # p=2.349e-08
ggplot(jj[pro_id=="P0DJI8"], aes(x=Sun_exposure, y=intensity))+
  geom_jitter(aes(color = Sun_exposure), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = Sun_exposure, group=factor(Sun_exposure)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("Sun_exposed"="Unexposed", "Sun_unexposed"="Exposed"))+
  scale_color_manual(values = c("Sun_exposed"="steelblue", "Sun_unexposed"="indianred"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none")

#### cell composition by BisqueRNA: Fig. S3F ####
library(Biobase); library(BisqueRNA); library(Matrix)
## expression matrix、meta and gene list from h5ad from https://spatial-skin-atlas.cellgeni.sanger.ac.uk/
## expression matrix and meta
if(1){ # python
  conda install -c conda-forge scanpy anndata scipy pandas numpy
  python
  import scanpy as sc
  import scipy.io as sio
  adata = sc.read_h5ad("/sh2/home/sunyuanqiang/projects/pQTL/scRNA/PMID38165934/bcc_and_normal-CG_portal_fat.h5ad")
  adata.X = adata.raw.X
  adata.write_h5ad("/sh2/home/sunyuanqiang/projects/pQTL/scRNA/PMID38165934/subset_raw.h5ad")
  sio.mmwrite("/sh2/home/sunyuanqiang/projects/pQTL/scRNA/PMID38165934/expr.mtx", adata.X)
  adata.obs.to_csv("meta.csv")
  adata.obs.to_csv("/sh2/home/sunyuanqiang/projects/pQTL/scRNA/PMID38165934/meta.csv")
  adata.var.to_csv("/sh2/home/sunyuanqiang/projects/pQTL/scRNA/PMID38165934/genes.csv")
} 

expr <- readMM("/Volumes/ElementsSE_syq/scRNA/PMID38165934/expr.mtx")
expr <- t(expr); dim(expr) # gene × cell = 32983 155401
# meta data
meta <- read.csv("/Volumes/ElementsSE_syq/scRNA/PMID38165934/meta.csv",row.names = 1)
dim(meta) # 155401 cells x  9 annotations
# gene list
genes <- read.csv("/Volumes/ElementsSE_syq/scRNA/PMID38165934/genes.csv", header = FALSE)
head(genes) # 32983 genes
# Data alignment
rownames(expr) <- genes$V1
colnames(expr) <- rownames(meta)
stopifnot(nrow(expr) == length(genes$V1)); 
stopifnot(ncol(expr) == nrow(meta))
all(colnames(expr) == rownames(meta)) # TRUE

## only use samples from normal body，20562 cells
meta1<-copy(meta); 
meta1$barcode<-rownames(meta1); 
meta1<-as.data.table(meta1) 
meta1<-meta1[X02_group=="body",]; 
dim(meta1)
meta1[, mysample:= ifelse(grepl("-", barcode), sub(".*-", "", barcode), 
                          sub("^(SC\\d+control)_.*\\.(\\d+_\\d+)$", "\\1_\\2", barcode))]
table(meta1$X03_location) # abdomen     arm = 14033    6529 

meta1$cellType <- meta1$X04_celltypes
meta1[X04_celltypes %in% c("B cells", "Plasma Cells")]$cellType <- "B cells"
meta1[X04_celltypes %in% c("Basal keratinocytes", "Suprabasal keratinocytes")]$cellType <- "Keratinocytes"
meta1[X04_celltypes %in% c("ILC_NK", "NK cells")]$cellType <- "NK cells"
meta1[X04_celltypes %in% c("DC", "Macrophages", "Monocytes")]$cellType <- "Myeloid cells"
meta1[X04_celltypes == "Neuronal_Schwann Cells"]$cellType <- "Schwann Cells "
meta1 <- meta1[!X04_celltypes %in% c("Chondrocytes", "Skeletal muscle cells", "Mast cells"), ]
table(meta1$cellType)

## build scRNA input
sc_meta <- data.frame(cellType = meta1$cellType, 
                      SubjectName = meta1$mysample, 
                      row.names = meta1$barcode, 
                      stringsAsFactors = FALSE)
expr1 <- expr[, rownames(sc_meta)]; 
keep_genes <- Matrix::rowSums(expr1) > 0
expr1 <- expr1[keep_genes, ]; 
head(rownames(sc_meta)); head(colnames(expr1)); 
all(colnames(expr1) == rownames(sc_meta)) # TRUE

expr1_dense <- as.matrix(expr1)
sc.eset <- ExpressionSet(assayData = expr1_dense, 
                         phenoData = AnnotatedDataFrame(sc_meta))

## build bulk RNA input
rr<-fread("~/Desktop/省皮/project/pQTL/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") 
sel_col<-grep("gene_ID|.+SRNA$",names(rr),value=T); 
rr[,..sel_col]->rr_cnt; 
names(rr_cnt)<-gsub("SRNA","",names(rr_cnt))
names(rr_cnt)[1]<-"gene_id"; 
rr_cnt<-rr_cnt[grep('ENSG',gene_id)]; 
rr_cnt[, c("n_detect", "total_count", "mean_exp") := .(rowSums(.SD > 0), rowSums(.SD), rowMeans(.SD)), .SDcols = -1]
rr_cnt2 <- rr_cnt[n_detect >= 3 & total_count > 5] 

tss<-fread("~/Desktop/reference/ensembl_allGene_pos_biomart.txt")
tss1<-unique(tss[,c(1,14,19)]); 
names(tss1)<-c("gene_id", "gene_name", "gene_type")
rr_cnt3<-merge(rr_cnt2, tss1[,.(gene_id, gene_name)], by="gene_id", all.x=T)
rr_cnt3[, gene_name:=ifelse(is.na(gene_name), gene_id, gene_name)] 
rr_cnt3[, gene_name:=ifelse(gene_name=="", gene_id, gene_name)] 
sum(duplicated(rr_cnt3$gene_name)) 
rr_cnt3<-rr_cnt3[order(gene_name, -total_count)]
rr_cnt4<-rr_cnt3[!duplicated(gene_name)]

bulk_expr_mat <- as.matrix(rr_cnt4[, 2:208])
rownames(bulk_expr_mat) <- rr_cnt4$gene_name
class(bulk_expr_mat); 
typeof(bulk_expr_mat) 
dim(bulk_expr_mat) # 45345 genes × 207 samples
bulk.eset <- ExpressionSet(assayData = bulk_expr_mat, 
                           phenoData = AnnotatedDataFrame(data.frame(row.names = colnames(bulk_expr_mat))))

## check bulk RNA and scRNA correspondance
length(intersect(rownames(sc.eset), rownames(bulk.eset))) # 15926; >10000 good
table(pData(sc.eset)$cellType) # without NA, without extremely rare cell types（<20 cells）
exprs(sc.eset)[1:5,1:5]; exprs(bulk.eset)[1:5,1:5]

## perform BisqueRNA 
res<-BisqueRNA::ReferenceBasedDecomposition(bulk.eset = bulk.eset, 
                                            sc.eset = sc.eset, 
                                            markers = NULL, 
                                            use.overlap = FALSE)
tmp<-res$bulk.props

## Deconvolution Results Evaluation: The results are reliable and no further adjustments are needed.
colSums(tmp)
mycor<-cor(t(tmp)) 
apply(tmp, 1, function(x) mean(x > 0)) 
cell_summary <- data.frame(cell_type = rownames(tmp), 
                           mean = apply(tmp, 1, mean), 
                           median = apply(tmp, 1, median),
                           min = apply(tmp, 1, min), 
                           max = apply(tmp, 1, max)) 
## plot
library(ggnewscale)
aa<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/BisqueRNA_cell_proportion.txt", header = T)
aa_long <- melt(aa, id.vars = "CellType", variable.name = "Sample", value.name = "Proportion")

hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS1_sample_information.xlsx",  sheet=1); 
hh2<-as.data.table(hh2);
names(hh2)[1]<-"sample"
aa_long2<-merge(aa_long, hh2[,.(sample, sampling_site)], by.x="Sample", by.y="sample", all.x=T); 
aa_long2[is.na(sampling_site)]
aa_long2<-aa_long2[order(sampling_site, Sample),] 
aa_long2$Sample<-factor(aa_long2$Sample, levels = unique(aa_long2$Sample))

median_df <- aa_long2[, .(median_prop = median(Proportion)), by = CellType]; 
median_df <- median_df[order(-median_prop)]
aa_long2$CellType <- factor(aa_long2$CellType, levels = median_df$CellType) 

ggplot(aa_long2, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity",width = 1, colour = "white", linewidth = 0.05) +
  labs(fill = "Cell type") + 
  theme_bw() +  
  ylab("Cell proportion") + 
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 14, face = "bold"))


