setwd("~/Desktop/省皮/project/pQTL/")
#### base prepare ####
library(data.table);library(ggplot2);
library("tidyr"); library(viridis);library(ggpointdensity);library(ggpubr)

std <- function(x) sd(x,na.rm = T)/sqrt(length(x[!is.na(x)]))
mytheme <- theme_classic(base_size = 22) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=18,face = "bold") #设置图例字体的大小及字体
  #,axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1,color="black"),
  #,axis.line = element_line(size = 1.2),axis.ticks=element_line(size=1.2),#坐标轴及 ticks 的粗度
)

name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

sample_id<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息.xlsx",col_names  = T); sample_id<-as.data.table(sample_id)
sample_id<-sample_id[!is.na(sample_id$皮肤蛋白样本名称),c(15,17)]; names(sample_id)<-c("names1","names2")
sample_id[,names3:=gsub("\\[\\d+\\] ","",names2)]

sample_info<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息.xlsx",sheet = 2,col_names = T); sample_info<-as.data.table(sample_info)

#### Clinical info & sampling sites ####
## process clinical information ####
hh<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_修.xlsx", sheet=2,col_names = T); hh<-as.data.table(hh)
hh1<-hh[!is.na(skin_protein_sample_No),][grep("^S",No,invert = T),c(1,3:11)]
hh1[, sample_pos:=gsub(" ", "_", sample_pos)]; hh1[, sample_pos:=gsub("/", "_", sample_pos)]; hh1[, Body_Region:=gsub(" ", "_", Body_Region)] # 去掉空格，省掉麻烦
Sun_exposed_pos<-c("Forearm", "Upper_Arm", "Shoulder", "Neck", "Neck_and_Upper_Back", "Lower_Leg_Calf", "Knee", "Kneecap", "Elbow", "Back")
Sun_Unexposed_pos<-c("Lower_Back_Lumbar", "Thigh", "Abdomen", "Chest", "Back")
hh1[, Sun_exposure := ifelse(sample_pos %in% Sun_exposed_pos, "Sun_exposed", ifelse(sample_pos %in% Sun_Unexposed_pos,"Sun_unexposed", NA_character_))]
table(hh1$sample_pos); table(hh1$Body_Region); table(hh1$Anatomical_Category); table(hh1$Sun_exposure) # 11; 3; 8; 2个 levels

hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_省立正常对照_外送20200623.xls", sheet = 1) #里面含有采样时间信息
hh2<-as.data.table(hh2); names(hh2)[1:2]<-c("sample","date"); hh2<-hh2[,1:2]
hh2[,date2:=as.Date(date, format = "%Y-%m-%d")]; hh2[,date_duration:=as.numeric(max(date2, na.rm = TRUE)-date2)] # 距离最新采样时间的间隔时间
hh3<-merge(hh1, hh2[,c(1,4)], by.x="No", by.y="sample", all.x=T)
hh3$age<-as.numeric(hh3$age); hh3$date_duration<-as.numeric(hh3$date_duration)
hh3$gender <- factor(hh3$gender, levels=c("female","male")); hh3$Sun_exposure<-factor(hh3$Sun_exposure, levels=c("Sun_exposed", "Sun_unexposed")) # 以 female 和 Sun_exposed作为基线
hh3$sample_pos<-factor(hh3$sample_pos); hh3$Body_Region<-factor(hh3$Body_Region); hh3$Anatomical_Category<-factor(hh3$Anatomical_Category); 

hh3[, region:=gsub("_.+","",region)]; hh3[, region:=paste("China",region, sep="_")]; nrow(hh3) # 261
summary(hh3$age); table(hh3$gender); table(hh3$Sun_exposure)
# 13.00   38.00   52.00   50.31   63.00   95.00; female   male = 103(39.5%), 158(60.5%); Sun_exposed Sun_unexposed = 158           103
write.csv(hh3[,c(1:5, 11)],"~/Desktop/省皮/project/pQTL/manuscript/TableS1_261smp_info.csv",quote = F,row.names = F)

ggplot(hh3, aes(x=age)) + geom_histogram(bins=20, fill="steelblue",color='grey90', linewidth=0.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) +
  labs(x ="Age", y="Numbers of individuals") + mytheme

mypie<-as.data.table(table(hh3$gender)); names(mypie)<-c("category","value")
ggplot(mypie, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1, color = "white") + coord_polar(theta = "y") + 
  scale_fill_manual(values = c("male" = "steelblue","female" = "indianred")) +
  theme_void() + theme(legend.position = "none")

#### DIA QC: hela ####
sample_id<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息.xlsx", col_names  = T)
## Hela protein number ##
aa<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/Hela质控数据/M-GSGC0290947_13hela_Report_pro.tsv",header = T)
names(aa)[-1:-5]<-gsub("DIA","",gsub("\\[\\d+\\] ","",gsub("_S.+.PG.Quantity","",names(aa)[-1:-5])))
aa[,PG.ProteinAccessions:=gsub(";.+","",PG.ProteinAccessions)]
tmp<-aa; tmp[is.na(tmp)]<-1
pheatmap::pheatmap(log10(tmp[,-1:-5]),cluster_rows=1,cluster_cols=1,col=colorRampPalette(c("white","red"))(50),
                   treeheight_row = 0,treeheight_col = 0,show_rownames = 0,show_colnames = 0)

melt.data.table(aa[,-2:-5],id.vars = "PG.ProteinAccessions")->aa1
names(aa1)<-c("pro_id","sample","intensity")
aa1[,nonNA_num:=sum(intensity!="NaN"),by=.(sample)]
unique(aa1[order(nonNA_num,decreasing = F),c(2,4)])->aa2
ggplot(aa2, aes(x = sample, y = nonNA_num))+
  geom_bar(position="dodge", stat="identity")+
  coord_cartesian(ylim = c(0, 7000)) +
  # coord_flip()+ #将垂直的条形图变为竖直的条形图，此时注意横坐标的顺序，用 rev 反向一下
  labs(x = 'Hela sample', y = 'Number of identified proteins', color = NULL)+
  mytheme + theme(axis.text.x = element_text(angle = 45,vjust = 0.9, hjust=1))

## Hela protein overlap ##
library("UpSetR")
listInput<-list(aa[!is.na(hela_1)]$PG.ProteinAccessions, aa[!is.na(hela_2)]$PG.ProteinAccessions, aa[!is.na(hela_3)]$PG.ProteinAccessions, 
                aa[!is.na(hela_4)]$PG.ProteinAccessions, aa[!is.na(hela_5)]$PG.ProteinAccessions, aa[!is.na(hela_6)]$PG.ProteinAccessions, 
                aa[!is.na(hela_7)]$PG.ProteinAccessions, aa[!is.na(hela_8)]$PG.ProteinAccessions, aa[!is.na(hela_9)]$PG.ProteinAccessions, 
                aa[!is.na(hela_10)]$PG.ProteinAccessions, aa[!is.na(hela_11)]$PG.ProteinAccessions, aa[!is.na(hela_12)]$PG.ProteinAccessions, aa[!is.na(hela_12)]$PG.ProteinAccessions)
names(listInput)<-paste0("Hela_", seq(1,13,1))
upset(fromList(listInput), order.by = "freq",sets = rev(names(listInput)),keep.order=T,mainbar.y.label = "Intersection of identified proteins", nintersects = 10, 
      sets.x.label = "Identified protein count", point.size = 2, line.size = 1.2, text.scale = c(2, 2, 2, 2, 1.5, 3))

## correlation of protein abundace ##
library(corrplot)
mycor<-cor(aa[,-1:-5], method = "pearson",use = "pairwise.complete.obs")
corrplot(mycor, method = "circle",type = "upper", col.lim=c(0.95,1), is.cor=F, order = "hclust", tl.col = "black", tl.srt = 45)

library(GGally) 
# 自定义 lower 面板函数：散点图 + lm 拟合线
my_lower_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.3, size = 0.1) +
    geom_smooth(method = "lm", color = "red", se = FALSE, size = 0.5) +
    theme_minimal()
}

# 自定义 diag 面板函数：靛蓝直方图 + 密度曲线
my_diag_hist <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_histogram(aes(y = ..density..), binwidth = 0.5,
                   fill = "#4B0082", color = "white", ...) +
    geom_density(color = "red", size = 1) +
    theme_minimal()
}

tmp<-aa[,-1:-5]
myorder<-paste0("hela_", seq(1,13,1)); setcolorder(tmp, myorder) #调整顺序
tmp[, names(tmp) := lapply(.SD, log10)] # 转变为 log10
ggpairs(tmp,
        upper = list(continuous = wrap("cor", size = 5, colour = "black", fontface = "bold")),  # 显示相关系数
        lower = list(continuous =wrap(my_lower_lm)),  # 下三角为散点图 + 回归线
        diag = list(continuous = wrap(my_diag_hist)))+  # 对角线为直方图
  mytheme+theme(axis.text = element_text(size = 12))

#### Quality control of DIA proteomes ####
## remove sample and protein based on NA number ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro.xls",header  = T) 
names(qq2)[6:ncol(qq2)]<-gsub("\\[\\d+\\] ","",names(qq2)[6:ncol(qq2)])
names(qq2)[6:ncol(qq2)]<-sample_id$names1[match(names(qq2)[6:ncol(qq2)],sample_id$names3)] # rename sample names
fwrite(qq2[,-2:-4], "~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv", col.names = T, row.names = F, sep=",", quote = F)

qq3<-as.data.frame(qq2[, -2:-5]); rownames(qq3)<-qq2$PG.ProteinAccessions; qq3<-qq3[,-1]
qq3[] <- lapply(qq3, function(x) { # 0/1值替换原值
  if (is.numeric(x)) {
    as.numeric(!is.na(x))
  } else { x }
})

# for smp, < 80% median
sample_counts <- colSums(qq4 == 1); median_count <- median(sample_counts) # 4787
qq5 <- qq4[, sample_counts >= 0.8 * median_count] # 248 out of 249
# for protein, identified at least 10 smps;
sum(rowSums(qq5) >= 10) # 5850; sum(rowSums(qq5) > ceiling(ncol(qq5)*0.5)), 4867
protein_counts <- rowSums(qq5 == 1); keep_proteins <- protein_counts >= 10
qq6 <- qq5[keep_proteins, ] # 5850 out of 5873
tmp<-copy(qq6); tmp$protein<-rownames(tmp)
fwrite(tmp[,c(249, 1:248)], "~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv", col.names = T, row.names = F, sep=",", quote = F)

ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
length(unique(ff$protein)); length(unique(ff[grep(";",protein)]$protein)); length(unlist(strsplit(ff$protein, ";"))) # 5850; 74 groups; 5956 indicidual proteins

## kept peptide number followed filtered protein ####
sample_info<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息.xlsx",sheet = 2,col_names = T); sample_info<-as.data.table(sample_info)
sample_info[!is.na(sample_info$skin_protein_sample_name),c(15,17)]->sample_id # differ with blood
names(sample_id)<-c("names1","names2"); sample_id[,names3:=gsub("\\[\\d+\\] ","",names2)]; sample_id[,names4:=gsub(".PG.Quantity","",names3)]

aa<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pep.xls",header  = T)
sel_col<-grep("PG.ProteinGroups|PEP.StrippedSequence|TotalQuantity",names(aa),value = T); aa <- aa[,..sel_col]
names(aa)[3:ncol(aa)]<-gsub("\\[\\d+\\] ","",names(aa)[3:ncol(aa)]); names(aa)[3:ncol(aa)]<-gsub(".EG.TotalQuantity \\(Settings\\)", "", names(aa)[3:ncol(aa)])
names(aa)[3:ncol(aa)]<-sample_id$names1[match(names(aa)[3:ncol(aa)],sample_id$names4)] # rename sample names
fwrite(aa, "~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pep_rename.csv", col.names = T, row.names = F, sep=",", quote = F)

qq6<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(aa)[3:ncol(aa)] %in% qq6$protein; selected_cols <- c(TRUE, TRUE, cols_to_keep)
rows_to_keep <- aa$PG.ProteinGroups %in% qq6$protein
aa2 <- aa[rows_to_keep, ..selected_cols]
dim(aa2); length(unique(aa2$PG.ProteinGroups)) # 40753, 250; 5859
length(unique(aa2$PEP.StrippedSequence)) # 34041

## pheatmap and GO of filtered proteomes ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]

melt.data.table(qq3,id.vars = "PG.ProteinAccessions")->qq4; names(qq4)<-c("pro_id","sample","intensity")
qq5<-qq4; qq5[grep("^S",sample,invert = T)]->qq5
qq5[,min:=min(intensity,na.rm = T),by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity),min/2,intensity)] # minimum/2填补 NA
qq5[,intensity_log:=ifelse(is.na(intensity),NA,log10(intensity))]
qq5[,intensity_log_median:=median(intensity_log),by=.(pro_id)]
hist(qq5$intensity_log); summary(qq5$intensity_log) # -1.552   1.318   1.832   1.758   2.252   5.957 
tmp<-dcast(qq5, formula = sample~pro_id, value.var = "intensity_log")

# top 10% most-abundant proteins GO：
qq6<-unique(qq5[,.(pro_id,intensity_log_median)]); qq6<-qq6[order(intensity_log_median,decreasing = T)]
library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
tmp3 <- as.data.table(qq6 %>% tidyr::separate_rows(pro_id, sep = ";")); tmp3 <- tmp3[order(intensity_log_median,decreasing = T)];
tmp3[,bin:=ceiling(rank(intensity_log_median)*100/length(intensity_log_median))] 
enrich_go <- enrichGO( gene = tmp3[bin>=90]$pro_id, # top 10% most-abundant proteins
                       OrgDb = org.Hs.eg.db, keyType = "UNIPROT", universe= tmp3$pro_id, # background genes
                       ont = "ALL",
                       pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
dotplot(enrich_go, x='GeneRatio',showCategory = 10, font.size = 18) 
barplot(enrich_go, showCategory = 20, x = "GeneRatio")
## Keratin, collegen ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro.xls",header  = T) 
names(qq2)[6:ncol(qq2)]<-gsub("\\[\\d+\\] ","",names(qq2)[6:ncol(qq2)])
names(qq2)[6:ncol(qq2)]<-sample_id$names1[match(names(qq2)[6:ncol(qq2)],sample_id$names3)] # rename sample names
melt.data.table(qq2[,c(-1,-3:-5)],id.vars = "PG.Genes")->qq4; names(qq4)<-c("gene_id","sample","intensity")
qq4<-qq4[grep("^S",sample,invert = T)]; qq5<- qq4 %>% tidyr::separate_rows(gene_id, sep = ";"); qq5<-as.data.table(qq5); qq5<-qq5[gene_id!=""]
qq5[, intensity_gene:=sum(intensity[!is.na(intensity)]), by=.(sample,gene_id)] # sum intensity by gene in each sample, 因为有些基因对应多个蛋白
unique(qq5[,-3])->qq5; qq5[intensity_gene==0]$intensity_gene<-NA

qq5[,c("mean_intensity_gene","std_intensity_gene","median_intensity_gene"):= .(mean(intensity_gene[!is.na(intensity_gene)]), std(intensity_gene[!is.na(intensity_gene)]), median(intensity_gene[!is.na(intensity_gene)])),by=.(gene_id)]
qq5[,c("mean_intensity_gene_log","std_intensity_gene_log","median_intensity_gene_log"):= .(log2(mean_intensity_gene),log2(std_intensity_gene),log(median_intensity_gene)),by=.(gene_id)]
unique(qq5[gene_id!="",c(-2,-3)])->qq6
qq6[,c("rank_mean","rank_median"):=.(rank(-mean_intensity_gene_log),rank(-median_intensity_gene_log))]

qq6[,class:=ifelse(grepl("KRT",gene_id),"KRT", ifelse(grepl("COL",gene_id),"COL","Others"))]
qq6$class<-factor(qq6$class,levels = c("KRT","COL","Others"))
ggplot(qq6, aes(x=rank_mean, y=mean_intensity_gene_log)) + 
  geom_point(data=qq6[class=="Others"],colour="gray",size=2,shape=20,alpha=0.1)+
  geom_point(data=qq6[class=="KRT"],colour="red",size=3,shape=20,alpha=0.8)+
  geom_point(data=qq6[class=="COL"],colour="blue",size=3,shape=20,alpha=0.8)+
  #geom_point(data=qq6[class=="BloodProtein"],colour="black",size=3,shape=20)+
  labs(x = 'Protein abundance rank', y = 'Protein abundance', color = NULL)+
  mytheme

ggplot(qq6, aes(x = class, y = mean_intensity_gene_log,color = class))+
  geom_jitter(size=2,alpha = 0.4,position = position_jitter(width = 0.25))+
  geom_boxplot(outlier.shape = NA, fill = '#00000000',width = 0.7, lwd=1.1,position = position_dodge(width = 0.85))+
  scale_color_manual(values = c("red","blue","gray"))+
  labs(x = '', y = 'Protein abundance (log2)', color = NULL)+
  #scale_color_discrete(name = "Protein class", labels = c("Keratins(40)", "Collagen(30)", "Others(5956)"))+
  scale_x_discrete(labels=c("KRT" = "Keratins", "COL" = "Collagen","Others" = "Others"))+
  mytheme + theme(legend.position="none",axis.text.x = element_text(face="bold", size=14,angle = 45,vjust = 0.95, hjust=1))
wilcox.test(qq6[class=="KRT"]$mean_intensity_gene_log,qq6[class=="Others"]$mean_intensity_gene_log)
wilcox.test(qq6[class=="COL"]$mean_intensity_gene_log,qq6[class=="Others"]$mean_intensity_gene_log)

# correlated KRT or collagen: KRT5 vs KRT14  KRT1 vs KRT10 from PMID: 36285538
tmp1<-qq5[gene_id=="KRT5"]; tmp2<-qq5[gene_id=="KRT14"]
tmp1<-qq5[gene_id=="KRT1"]; tmp2<-qq5[gene_id=="KRT10"]
tmp3<-merge(tmp1[,1:3], tmp2[,1:3], by="sample")
ggplot(tmp3, aes(x = log2(intensity_gene.x), y = log2(intensity_gene.y))) + 
  geom_point(color = "black", size = 1.5, alpha = 0.8) + 
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + labs(x = NULL, y = NULL) + mytheme
cor.test(tmp3$intensity_gene.x, tmp3$intensity_gene.y, method = "pearson") # KRT5 vs KRT14 R=0.985, p<2.2e-16; KRT1 vs KRT10 R=0.966, p<2.2e-16

## nonNA for each protein and sample ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]
melt.data.table(qq3,id.vars = "PG.ProteinAccessions")->qq4; names(qq4)<-c("pro_id","sample","intensity")

# for protein
qq4[, nonNA_num_pro:=sum(!is.na(intensity)), by=.(pro_id)]; tmp1<-unique(qq4[,.(pro_id, nonNA_num_pro)])
nrow(tmp1[nonNA_num_pro >= 124,])/nrow(tmp1) # 4877/5850=83.37%
summary(tmp1$nonNA_num_pro)  # 10.0   167.0   241.0   200.9   248.0   248.0
ggplot(tmp1, aes(x = nonNA_num_pro)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(x = "Number of samples in which the protein was quantified", y = "Protein count") +
  mytheme
# for sample
qq4[, nonNA_num_smp:=sum(!is.na(intensity)), by=.(sample)]; tmp2<-unique(qq4[,.(sample, nonNA_num_smp)])
summary(tmp2$nonNA_num_smp) # 3858    4519    4788    4740    4978    5395
tmp2<-tmp2[order(nonNA_num_smp,decreasing = T)]; tmp2[,sample_id:=1:.N]
# barplot
ggplot(tmp2, aes(x = sample_id, y = nonNA_num_smp)) + 
  geom_col(width = 1, fill = "steelblue", color = "white", linewidth = 0.2) + 
  labs(x = "Sample IDs", y = "Number of quantified proteins") +
  mytheme
# 箱线图
ggplot(tmp2, aes(x = "", y = nonNA_num_smp)) + 
  geom_boxplot(width = 0.3, fill = 'grey90', color = "black", lwd=1.1 , fatten = 2) +  
  geom_jitter(width = 0.2, size = 3, alpha = 0.7, color = "steelblue") +  
  labs(x = "", y = "Number of quantified proteins per sample") +
  mytheme

#### PRM vs DIA ####

## quantification correlation ####
library(data.table); library(ggplot2); library(ggpubr)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=22),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  axis.line = element_line(linewidth = 1.2),axis.ticks=element_line(linewidth=1.2),
  #axis.text.x = element_text(angle =  45,vjust = 0.9, hjust=1)
)

name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

dia<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA 数据验证/DIA_蛋白质鉴定列表.xlsx"); dia<-as.data.table(dia)
dia<-dia[,c(2, 6:15)]; names(dia)[1]<-"pro_id"
dia2 <- as.data.table(dia %>% tidyr::separate_rows(pro_id, sep = ";")); dia2 <- dia2[!duplicated(pro_id)]
dia3<-melt.data.table(dia2, id.vars = "pro_id"); names(dia3)<-c("pro_id","sample","value")
dia3[,min:=min(value, na.rm = T),by=.(pro_id)]; dia3<-dia3[!is.infinite(min)]
if (1) {dia3[, value2 := ifelse(is.na(value), min/2, value)]}

prm<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA 数据验证/M-GSGC0447052_PRM实验报告/项目结果展示/附件3_4D-PRM蛋白质定量结果.xlsx"); prm<-as.data.table(prm)
prm<-prm[,c(1, 15:24)]; names(prm)[1]<-"pro_id"
prm2 <- as.data.table(prm %>% tidyr::separate_rows(pro_id, sep = ";")); prm2 <- prm2[!duplicated(pro_id)]
prm3<-melt.data.table(prm2, id.vars = "pro_id"); names(prm3)<-c("pro_id","sample","value")
prm3[,min:=min(value, na.rm = T),by=.(pro_id)]; prm3<-prm3[!is.infinite(min)]
if (1) {prm3[, value2 := ifelse(is.na(value), min/2, value)]}
prm3[, sample := gsub(".*\\] (\\d+)_.*", "NCN\\1", sample)]; prm3[, sample := sprintf("NCN%03d", as.numeric(gsub("NCN", "", sample)))]

hh<-merge(prm3[,c(1,2,5)], dia3[,c(1,2,5)], all.x=T, by=c("pro_id","sample")); names(hh)[3:4]<-c("prm", "dia")
table(hh$sample) # n=93 for all samples

ggplot(hh, mapping = aes(x = log2(prm), y = log2(dia))) +
  geom_point(color="black", size=3, shape=20, alpha=1)+
  labs(x = 'log2(Protein abundance) by PRM', y = 'log2(Protein abundance) by DIA', color = NULL)+
  geom_smooth(aes(group = sample), method='lm', se=T, formula= y~x, col="red", linetype="dashed")+
  stat_cor(method = "spearman", size = 5, cor.coef.name="rho", label.x.npc="left", label.y.npc="top", label.sep = "\n")+
  facet_wrap(~sample, nrow=2, scales = "free", labeller=labeller(sample = c("NCN004" = "Sample 1","NCN005" = "Sample 2","NCN006" = "Sample 3","NCN007" = "Sample 4","NCN008" = "Sample 5",
                                                                            "NCN009" = "Sample 6","NCN010" = "Sample 7","NCN011" = "Sample 8","NCN012" = "Sample 9","NCN013" = "Sample 10")))+
  mytheme

hh1<-merge(hh, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
jj<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA 数据验证/十块正常组织样本信息.xlsx",sheet=2)
hh2<-merge(hh1, jj[,c(1,3,4)], by.x="sample", by.y="no")
fwrite(hh2[,.(sample, Sex, Age, gene_name, prm, dia)],"~/Desktop/省皮/project/pQTL/manuscript/Table_SX_PRM_vs_DIA_correlation.csv", col.names = T, row.names = F, sep=",", quote = F)

#### skin proteome comapare with other study ####
## self data processing: 1/2 minimum impution, VSN normalizaiton, log2 transform ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
library(tidyr); qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm) 
mm2<-as.data.table(mm1); mm2$pro_id<-rownames(mm1)
mm2_long<-melt(mm2, id.vars = "pro_id"); names(mm2_long)<-c("pro_id","sample","intensity_log")
mm2_long[, intensity:=2^intensity_log] # 将 log2 值恢复原样，以进行样本间的 aggregate
mm3<-mm2_long[,.(intensity_median=median(intensity,na.rm = T), intensity_mean=mean(intensity,na.rm = T)), by=.(pro_id)] # 同一 tissue 的不同 rep 之间 aggregate
mm3[,c("intensity_median_log", "intensity_mean_log"):=.(ifelse(is.na(intensity_median), NA, log10(intensity_median)), ifelse(is.na(intensity_mean), NA, log10(intensity_mean)))]

## overlap and cor with skin in PMID: 40713952 多组织衰老蛋白质组: DIA ####
aa<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/Cell_proteome_aging/mmc1.xlsx", sheet = 4, skip = 2); aa<-as.data.table(aa) 
mycounts <- colSums(aa[, 4:ncol(aa)] > 0) 
# Skin    Muscle     Heart  Pancreas     Liver     Aorta    Spleen     Lymph Intestine   Adrenal   Adipose      Lung 
# 6198      6690      8797      9734      9864      8929     10491     10642     10399     10420      9023     10824 
aa<-aa[,-2:-3]; aa[aa == 0] <- NA
mm4<-merge(mm3, unique(name_id_reviewd[,c(1,4)]), by="pro_id") 
jj<-merge(aa, mm4[,-2:-3], by.x="Protein", by.y="gene_name", all = T);
jj[,pro_id:=NULL]; names(jj)[c(ncol(jj)-1, ncol(jj))]<-c("Skin_self_median","Skin_self_mean")
jj<-jj[order(Protein, -Skin_self_median)]; jj<-jj[!duplicated(Protein)]; names(jj)[1]<-"gene_name"

## 计算两两间 spearman cor
library(corrplot)
mycor<-cor(jj[,-1], method = "spearman", use = "pairwise.complete.obs") # spearman 0.73, pearson 0.78
# tmp<-cor(jj[,-1], method = "pearson", use = "pairwise.complete.obs")
# corrplot(tmp, method = "circle",type = "upper", addCoef.col="black", cl.cex=1.5, number.digits=2, number.cex =0.5,col.lim=c(min(tmp),1), is.cor=F, order = "original", tl.col = "black", tl.srt = 45, tl.cex=0.8, mar = c(0,0,0,0))
corrplot(mycor, method = "circle",type = "upper",addCoef.col="black", cl.cex=1.5, number.digits=2, number.cex =0.5,col.lim=c(min(mycor),1), is.cor=F, order = "original", tl.col = "black", tl.srt = 45, tl.cex=0.8, mar = c(0,0,0,0))
mycor2<-as.data.table(mycor); mycor2[, tissue:=rownames(mycor)]
mycor2<-mycor2[!tissue%in%c("Skin_self_median", "Skin_self_mean", "reference")]

mycor3<-merge(mycor2[, .(Skin_self_mean, tissue)], overlap_mat2[, .(Skin_self_mean, tissue)], by="tissue") # mean更好
names(mycor3)[2:3]<-c("rho","overlap_num"); mycor3<-mycor3[order(-rho)]
mycor3$tissue<-factor(mycor3$tissue, levels = rev(mycor3$tissue))

ggplot(mycor3, aes(x = tissue, y = rho)) + # 横向 
  geom_col(fill = "steelblue", width = 0.7) + # 不画 overlap的数目了，要不然 overlap 的太多了，没有 skin-specific了
  labs(x = NULL, y = "Spearman rho") +
  coord_cartesian(ylim = c(0.3, 0.8))+
  mytheme + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 18))

library(viridis); library(ggpointdensity)
ggplot(data = jj, mapping = aes(x = Skin_self_mean, y = Skin)) +
  geom_pointdensity(size = 2, alpha = 0.5) + scale_color_viridis(option = "D") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", linewidth = 1) + 
  labs(x = 'Protein abundance in this study', y = 'Skin protein abundance in Ding et al.', color = NULL)+
  # coord_cartesian(xlim = c(0, 6), ylim = c(0, 6))+
  mytheme

cor.test(jj$Skin_self_mean, jj$Skin, method = "spearman"); cor.test(jj$Skin_self_mean, jj$Skin, method = "pearson")  # 0.7252399; 0.7787051


#### variation of populational proteome ####
## PCA ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) # 建立 vsn 模型并归一化，不能提前 log2 转化
meanSdPlot(fit) # 标准差(sd)在不同均值区间(rank mean)近乎恒定，说明拟合很好
mm1 <- predict(fit, mm) # 已经是log2 转化了

# 运行 PCA
pca_result <- prcomp(t(mm1), center = TRUE, scale. = TRUE) # 列是样本，需转置并 scale
summary(pca_result) # PC1+PC2=14.22%
pca_df <- as.data.frame(pca_result$x[, 1:2])
summary(pca_df$PC1); summary(pca_df$PC2) # -46.986 -17.998  -2.527   0.000  12.163  86.310 # -43.4719 -10.9111  -0.9038   0.0000  10.8494  60.5146
pca_df$sample <- rownames(pca_df); pca_df$sample<-gsub("SP", "", pca_df$sample)

# 处理样本信息
hh<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_修.xlsx", sheet=2,col_names = T); hh<-as.data.table(hh)
hh1<-hh[!is.na(skin_protein_sample_No),][grep("^S",No,invert = T),c(1,3:11)]
hh1[, sample_pos:=gsub(" ", "_", sample_pos)]; hh1[, sample_pos:=gsub("/", "_", sample_pos)]; hh1[, Body_Region:=gsub(" ", "_", Body_Region)] # 去掉空格，省掉麻烦
Sun_exposed_pos<-c("Forearm", "Upper_Arm", "Shoulder", "Neck", "Neck_and_Upper_Back", "Lower_Leg_Calf", "Knee", "Kneecap", "Elbow", "Back")
Sun_Unexposed_pos<-c("Lower_Back_Lumbar", "Thigh", "Abdomen", "Chest", "Back")
hh1[, Sun_exposure := ifelse(sample_pos %in% Sun_exposed_pos, "Sun_exposed", ifelse(sample_pos %in% Sun_Unexposed_pos,"Sun_unexposed", NA_character_))]
table(hh1$sample_pos); table(hh1$Body_Region); table(hh1$Anatomical_Category); table(hh1$Sun_exposure) # 11; 3; 8; 2个 levels

hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_省立正常对照_外送20200623.xls", sheet = 1) #里面含有采样时间信息
hh2<-as.data.table(hh2); names(hh2)[1:2]<-c("sample","date"); hh2<-hh2[,1:2]
hh2[,date2:=as.Date(date, format = "%Y-%m-%d")]; hh2[,date_duration:=as.numeric(max(date2, na.rm = TRUE)-date2)] # 距离最新采样时间的间隔时间
hh3<-merge(hh1, hh2[,c(1,4)], by.x="No", by.y="sample", all.x=T)
hh3$age<-as.numeric(hh3$age); hh3$date_duration<-as.numeric(hh3$date_duration)
hh3$gender <- factor(hh3$gender, levels=c("female","male")); hh3$Sun_exposure<-factor(hh3$Sun_exposure, levels=c("Sun_exposed", "Sun_unexposed")) # 以 female 和 Sun_exposed作为基线
hh3$sample_pos<-factor(hh3$sample_pos); hh3$Body_Region<-factor(hh3$Body_Region); hh3$Anatomical_Category<-factor(hh3$Anatomical_Category); 

# 绘图
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

## 变异贡献度测量 by ANOVA Eta: variance ~ gender + age + sampling_site ####
library(data.table); library(limma); library(car); set.seed(123)
hh<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_修.xlsx", sheet=2,col_names = T); hh<-as.data.table(hh)
hh1<-hh[!is.na(skin_protein_sample_No),][grep("^S",No,invert = T),c(1,3:11)]
hh1[, sample_pos:=gsub(" ", "_", sample_pos)]; hh1[, sample_pos:=gsub("/", "_", sample_pos)]; hh1[, Body_Region:=gsub(" ", "_", Body_Region)] # 去掉空格，省掉麻烦
Sun_exposed_pos<-c("Forearm", "Upper_Arm", "Shoulder", "Neck", "Neck_and_Upper_Back", "Lower_Leg_Calf", "Knee", "Kneecap", "Elbow", "Back")
Sun_Unexposed_pos<-c("Lower_Back_Lumbar", "Thigh", "Abdomen", "Chest", "Back")
hh1[, Sun_exposure := ifelse(sample_pos %in% Sun_exposed_pos, "Sun_exposed", ifelse(sample_pos %in% Sun_Unexposed_pos,"Sun_unexposed", NA_character_))]
table(hh1$sample_pos); table(hh1$Body_Region); table(hh1$Anatomical_Category); table(hh1$Sun_exposure) # 11; 3; 8; 2个 levels

hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_省立正常对照_外送20200623.xls", sheet = 1) #里面含有采样时间信息
hh2<-as.data.table(hh2); names(hh2)[1:2]<-c("sample","date"); hh2<-hh2[,1:2]
hh2[,date2:=as.Date(date, format = "%Y-%m-%d")]; hh2[,date_duration:=as.numeric(max(date2, na.rm = TRUE)-date2)] # 距离最新采样时间的间隔时间
hh3<-merge(hh1, hh2[,c(1,4)], by.x="No", by.y="sample", all.x=T)
hh3$age<-as.numeric(hh3$age); hh3$date_duration<-as.numeric(hh3$date_duration)
hh3$gender <- factor(hh3$gender, levels=c("female","male")); hh3$Sun_exposure<-factor(hh3$Sun_exposure, levels=c("Sun_exposed", "Sun_unexposed")) # 以 female 和 Sun_exposed作为基线
hh3$sample_pos<-factor(hh3$sample_pos); hh3$Body_Region<-factor(hh3$Body_Region); hh3$Anatomical_Category<-factor(hh3$Anatomical_Category); 

## proteome process: 1/2 minimum impution, VSN normalizaiton, log2 transform
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) # 建立 vsn 模型并归一化，不能提前 log2 转化
meanSdPlot(fit) # 标准差(sd)在不同均值区间(rank mean)近乎恒定，说明拟合很好
mm1 <- predict(fit, mm) # 已经是log2 转化了

## lm for循环, 使用type II sum of squares ANOVA
meta<-as.data.frame(hh3); rownames(meta)<-meta$No
common_samples <- intersect(colnames(mm1), rownames(meta)); length(common_samples) # 248
mm2 <- mm1[, common_samples, drop = FALSE]; meta2 <- meta[common_samples, , drop = FALSE]; 
all(rownames(meta2) == colnames(mm2))  # TRUE

res2 <- lapply(1:nrow(mm2), function(i) {
  fit <- lm(mm2[i, ] ~ age + gender + Sun_exposure, data = meta2) # 逐个蛋白拟合
  ss  <- Anova(fit, type = 2) # Type II ANOVA
  SS_err <- ss["Residuals", "Sum Sq"] # 残差平方和
  pe2 <- function(term) {  # 按平方和计算 Partial Eta²
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
nrow(res2[q_age<0.05]); nrow(res2[q_sex<0.05]); nrow(res2[q_Sun<0.05]) # 174; 8; 8
res2<-merge(res2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")

## 作图
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

## case plot
mm3<-as.data.table(mm2); mm3$pro_id<-rownames(mm2)
mm3_long<-melt.data.table(mm3, id.vars = "pro_id"); names(mm3_long)<-c("pro_id", "sample", "intensity")
jj<-merge(mm3_long, hh3[,c(1,2,3,11)], by.x = "sample", by.y="No")
# 阳光暴露增加表达：SAA1，P0DJI8，PMID: 26900010；CRP，P02741， PMID: 40595199
wilcox.test(jj[pro_id=="P0DJI8" & Sun_exposure=="Sun_exposed"]$intensity, jj[pro_id=="P0DJI8" & Sun_exposure=="Sun_unexposed"]$intensity)  # p=2.349e-08
ggplot(jj[pro_id=="P0DJI8"], aes(x=Sun_exposure, y=intensity))+
  geom_jitter(aes(color = Sun_exposure), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = Sun_exposure, group=factor(Sun_exposure)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("Sun_exposed"="Unexposed", "Sun_unexposed"="Exposed"))+
  scale_color_manual(values = c("Sun_exposed"="steelblue", "Sun_unexposed"="indianred"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none")

# 与 gender 紧密相关: RPS4Y1, P22090, PMID: 33670450,PMID: 14583743; PZP, P20742, PMID: 32913072, PMID: 40360480
wilcox.test(jj[pro_id=="P22090" & gender=="female"]$intensity, jj[pro_id=="P22090" & gender=="male"]$intensity)  # p= 6.693e-14
ggplot(jj[pro_id=="P22090"], aes(x=gender, y=intensity))+
  geom_jitter(aes(color = gender), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = gender, group=factor(gender)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("female"="Female", "male"="Male"))+
  scale_color_manual(values = c("female"="indianred", "male"="steelblue"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none")

wilcox.test(jj[pro_id=="P20742" & gender=="female"]$intensity, jj[pro_id=="P20742" & gender=="male"]$intensity)  # p= 2.426e-13
ggplot(jj[pro_id=="P20742"], aes(x=gender, y=intensity))+
  geom_jitter(aes(color = gender), size=2.5, alpha = 0.6, position = position_jitter(width = 0.25))+
  geom_boxplot(aes(color = gender, group=factor(gender)), outlier.shape = NA, fill = '#00000000', width = 0.7, fatten = 2, lwd=1.5, position = position_dodge(width = 0.85))+
  scale_x_discrete(labels = c("female"="Female", "male"="Male"))+
  scale_color_manual(values = c("female"="indianred", "male"="steelblue"))+
  labs(x = NULL, y = 'Normalized protein abundance', color = NULL)+
  mytheme+theme(legend.position="none",axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1))

#### co-expression/variation analysis ####
## co-regulation tSNE_map ####
## co-variation,ref : https://github.com/Rappsilber-Laboratory/ProteomeHD; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6901355/#R37
library(data.table); library(treeClust); library(WGCNA)
set.seed(42)
qq2[,-2:-5]->prohd
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
  y <- as.data.table( melt( x )) 
  y <- y[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_value = value ) ]
  y <- y[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl.self-comparisons)
  assign(mydt2[i],y)
}
names(tc_sim_dt)[3]<-"tc_sim";names(adj_mat_dt)[3]<-"tc_adj";names(tc_tom_dt)[3]<-"tc_tom"
tc_dt <- merge( tc_sim_dt, adj_mat_dt, by = c("Protein_1", "Protein_2"))
tc_dt <- merge( tc_dt, tc_tom_dt,by = c("Protein_1", "Protein_2"))

fwrite(tc_dt, "~/Desktop/pro_treeClust_similarities.csv") #不能直接写入移动硬盘。mv过去
tc_dt_final <- tc_dt[, .(Protein_1, Protein_2, coregulation_score = tc_tom)]
fwrite(tc_dt_final, "/Volumes/ElementsSE_syq/pQTL/pro_coregulation_scores.csv")

library(Rtsne); library(gridExtra)
DT <- fread("/Volumes/ElementsSE_syq/pQTL/pro_coregulation_scores.csv")  # 11846278
DT[, coreg_distance := (1 - log2(coregulation_score)) ]  # Turn co-regulation score back into a distance metric (distance越小，蛋白越相关) and log2-transform for better tSNE performance
tmp<-merge(DT, unique(name_id_reviewd[,c(1,4)]), by.x="Protein_1", by.y="pro_id", all.x=T)

DTm <- dcast( data = rbind( DT[, .(Protein_1, Protein_2, coreg_distance)],
                            DT[, .(Protein_1 = Protein_2, Protein_2 = Protein_1, coreg_distance)]),
              Protein_1 ~ Protein_2 , value.var = "coreg_distance")  # Turn the melted pairwise table back into a matrix
DTm <- as.data.frame(DTm) ; rownames(DTm) <- DTm$Protein_1 ; DTm$Protein_1 <- NULL
DTm <- as.dist(as.matrix( DTm ))         # Turn into numeric matrix then dist object
protein_IDs <- attr(DTm, "Labels")        # Extract protein IDs from dist object

set.seed(123)
SNE <- Rtsne(DTm, is_distance = TRUE, theta = 0.0, perplexity = 50, max_iter = 1000, verbose = TRUE) # 时间较长,10 min
SNE <- as.data.table( SNE$Y );SNE[, ID := protein_IDs ]
SNE2 <- as.data.table(SNE %>% tidyr::separate_rows(ID, sep = ";"))
ggplot(SNE2, aes(x = V1, y = V2))+
  geom_point(shape = 16, size= 1, alpha=0.5)+labs(x = 'tSNE dimension1', y = 'tSNE dimension2')+
  geom_point(data = SNE2[ ID == "P15924"], shape = 3, size= 4, alpha=1, colour = "red")+ # DESP
  geom_point(data = SNE2[ ID == "Q13835"], shape = 1, size= 4, alpha=1, colour = "blue")+ # PKP1
  mytheme+theme (axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())

# Organelle annotation on tSNE map
my_annotations <- fread("~/Desktop/省皮/project/pQTL/NBT2020 co-regulation map/tSNE_map_annotation.csv") # Read in annotation file
my_annotations[ Organelle == "" , Organelle := NA ]; my_annotations[ Manual_annotation == "" , Manual_annotation := NA ] # Replace empty strings with NAs
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
