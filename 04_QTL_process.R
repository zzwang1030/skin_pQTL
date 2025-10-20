setwd("~/Desktop/省皮/project/pQTL/QTL_calling_formal/")
######## common data ########
library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 24) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
  #, legend.text=element_text(size=9,face = "italic") #设置图例字体的大小及字体https://ggplot2.tidyverse.org/articles/ggplot2-specs.html
  #,axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1,color="black"),
  #,axis.line = element_line(size = 1.2),axis.ticks=element_line(size=1.2),#坐标轴及 ticks 的粗度
)
## protein uniprot_id2gene_name
name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

name_id_reviewd2<-unique(name_id_reviewd[,c(1,4)]); name_id_reviewd2<-name_id_reviewd2[!duplicated(pro_id)] # 去重，有些 pro_id对应多个基因 name，容易对不上数目
length(unique(name_id_reviewd2$pro_id)); length(name_id_reviewd2$pro_id) # 42433 = 42433

## transcriptional start sites info
tss<-fread("~/Desktop/reference/ensembl_uniprot_id_pos_biomart.txt",header = T) # download from ensembl biomart and manually curated/add ~50 proteins by syq
tss[,strand:=ifelse(strand==1,"+","-")]
tss2<-tss[,.(chr=chromosome,start=transcription_start_site-1,end=transcription_start_site,id=uniprot_id,gene_id=gene_id,strand=strand)]
tss2<-tss2[order(-id,chr,start),]
tss2[,sel_tss:=ifelse(length(unique(chr))>1,start[1],ifelse(strand=="+",min(start),max(start))),by=.(id)]
tss3<-tss2[start==sel_tss,.SD,by=.(id)]; tss3<-tss3[!duplicated(id),]; tss3<-tss3[id!="",.(chr,start,end,id,gene_id,strand)]

if (0) { # manually curated/add proteins in ensembl_uniprot_id_pos_biomart.txt
  merge(name_id_reviewd,tss3,by.x="pro_id",by.y="id",all.x = T)->tmp
  tmp[grep("-",pro_id,invert = T),]->tmp 
  tmp[is.na(chr)]->tmp1
}
 

sel_smp<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/tmp_sample",header = F)
sel_chr<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/sel_chr",header = F)

######################### skin protein ######################### 
######## common data ########
## donor basic info
sample_info <- readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息.xlsx",sheet = 2,col_names = T)
sample_info <- as.data.table(sample_info)
sample_id <- sample_info[!is.na(sample_info$skin_protein_sample_name),c(15,17)]# differ with blood
names(sample_id)<-c("names1","names2")
sample_id[,names3:=gsub("\\[\\d+\\] ","",names2)]

######## Preapre protein matrix ########
## QC of skin protein matrix ####
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro.xls",header  = T)
names(qq2)[6:ncol(qq2)]<-gsub("\\[\\d+\\] ","",names(qq2)[6:ncol(qq2)])
names(qq2)[6:ncol(qq2)]<-sample_id$names1[match(names(qq2)[6:ncol(qq2)],sample_id$names3)] # rename sample names
names(qq2)[1]<-"pro_id"
qq2 <- as.data.table(qq2 %>% tidyr::separate_rows(pro_id, sep = ";")); nrow(qq2) # 5979
qq2<-qq2[order(pro_id,PG.Qvalue),];qq2<-qq2[!duplicated(pro_id)] # remove duplicated proteins based on Qvalue, 5979
qq2[,-2:-5]->qq2;names(qq2)<-gsub("SP","",names(qq2)) # remove suffix of sample names to match genotype sample
my_sel <- names(qq2)%in%c("pro_id",sel_smp$V1)
qq3 <- qq2[,..my_sel]; ncol(qq3)-1 # select pheno samples overlapped with geno samoles (186 = 249 overlap 198)
qq3<-melt(qq3, id.vars = "pro_id"); names(qq3)<-c("pro_id","sample","value")
qq3[, NA_num_pro:=sum(is.na(value)), by=.(pro_id)]; qq3[,NA_num_smp:=sum(is.na(value)),by=.(sample)] # NA counts for each protein and sample
qq3<-qq3[NA_num_pro < ceiling((ncol(qq2)-1)*0.5)] # remove proteins with NA in more than half of samples。这一步其实没过滤掉样本。
unique(qq3[,.(sample,NA_num_smp)])->tmp; tmp[,pro_num:=nrow(qq2) - NA_num_smp]
qq3<-qq3[NA_num_smp < nrow(qq2)-0.8*median(tmp$pro_num)] # remove samples with protein identifications below 80% of the median value of all samples. ref: PMID_36797296
length(unique(qq3$pro_id)); length(unique(qq3$sample)) # remaining 5391/5979 proteins and 186 smp

## Normalization of skin protein matrix  ####
qq3[,min:=min(value,na.rm = T), by=.(pro_id)]
if (1) {qq3[,value:=ifelse(is.na(value),min/2,value)]}
qq3[,value_log:=ifelse(is.na(value),NA,log2(value))] # log transformed
qq3[,value_log_scale:=as.numeric(scale(value_log,center = T,scale = T)), by=.(pro_id)] # scale for each protein across sample
qq3[,value_log_scale_rank:=rank(value_log_scale,na.last="keep",ties.method = "average"), by=.(sample)] # rank within each sample across  proteins
qq3[,value_log_scale_rank_qnorm:= qnorm((value_log_scale_rank-0.5)/sum(!is.na(value_log_scale_rank))), by=.(sample)] # quantile normalized (rank inverse normalized) within each smp
qq4<-dcast(qq3,pro_id~sample,value.var = "value_log_scale_rank_qnorm") 

## add TSS info to protein matrix ####
tmp<-merge(qq4, tss3[,-5:-6], by.x = "pro_id", by.y="id"); nrow(tmp); ncol(tmp)-4 # 5381 proteins and 186 smp
my_sel2<-c(ncol(tmp)-2,ncol(tmp)-1,ncol(tmp),1,2:(ncol(tmp)-3))
qq5<-tmp[,..my_sel2] # adjust column order to bed format
qq5[,chr:=paste0("chr",chr)]; qq5<-qq5[order(chr,start)] # 后续index 需要排序
qq6<-copy(qq5); qq6<-qq6[chr%in%sel_chr$V1]; dim(qq6)# 选择常染色体上的基因。5192 proteins and 186 smp
names(qq6)[1]<-"#Chr"; fwrite(qq6, "~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt", sep = "\t",col.names = T)
# 后续在 linux中 进行 expression phenotype的 PCA: ~/Desktop/省皮/project/pQTL/QTL_calling_formal/pQTL_formal_code.docx

## skin protein covariats processing ########
library(findPC)
pca<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca",header = T)
names(pca)[1]<-"PC_id"
pca[,PC_id:=paste0("geno_",gsub("hg38.+svd_","",PC_id))]
pca1<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.forPCA.gz.pca",header = T)
names(pca1)[1]<-"PC_id"
pca1[,PC_id:=paste0("pheno_",gsub("skin.+svd_","",PC_id))]

# select K PCs for use
pve<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca_stats",header = F,nrows = 3)
pve<-transpose(pve,keep.names = "PC",make.names = "V1");pve[,PC:=1:.N];pve[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve$PC,pve$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve$prop_var,number = c(30,60,100),method = 'all',aggregate = 'median',figure = T)
sel_genoPCs<-5 

pve1<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.forPCA.gz.pca_stats",header = F,nrows = 3)
pve1<-transpose(pve1,keep.names = "PC",make.names = "V1"); pve1[,PC:=1:.N]; pve1[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve1$PC,pve1$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve1$prop_var,number = c(40,60,100,150),method = 'all',aggregate = 'median',figure = T)
sel_phenoPCs<-11 

transpose(pca,keep.names = "sample", make.names = "PC_id")->tmp
transpose(pca1,keep.names = "sample", make.names = "PC_id")->tmp1
merge(tmp1[,1:(sel_phenoPCs+1)],tmp[,1:(sel_genoPCs+1)],by="sample")->tmp2 
merge(tmp2,sample_info[,.(No,gender,age,Other_disease,sample_pos)],by.x="sample",by.y="No",all.x = T)->tmp3 # sample_pos should be included for skin. All geno_pheno overlapped sample should be retained. The missing values should be encoded as NA.
# convert catagory variables to dummy variables
tmp3[,gender:=ifelse(gender=="male",1,0)]
Other_disease<-tmp3$Other_disease; Other_disease <- model.matrix(~Other_disease); Other_disease<-dplyr::as_tibble(Other_disease)
sample_pos<-tmp3$sample_pos; sample_pos <- model.matrix(~sample_pos); sample_pos<-dplyr::as_tibble(sample_pos)

tmp3<-cbind(tmp3,Other_disease[,-1], sample_pos[,-1])
tmp3[,Other_disease:=NULL]; tmp3[,sample_pos:=NULL]; str(tmp3)

transpose(tmp3,keep.names = "PC_id", make.names = "sample")->tmp4
tmp4[is.na(tmp4)]<-"NA" #  The missing values should be encoded as NA 
fwrite(tmp4,"~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.12pheno_5geno.PCs",col.names = T,quote = F,sep = "\t")
# 后续在 linux 里进行 QTL calling: ~/Desktop/省皮/project/pQTL/QTL_calling_formal/pQTL_formal_code.docx

#### pQTL charactering ####
## 计算显著的 SNP-protein pair 数目 ####
# all pairs太大，在服务器跑
library(data.table)
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header=T); nrow(nn) # 15,669,687
length(unique(nn$gene_id)); length(unique(nn$variant_id)) # 5152 proteins; 3461609 snps。其实有5192-5152 个蛋白因为 variation=0，跑的时候被 fastQTL 过滤掉了
pp<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
length(unique(pp$gene_id)); length(unique(pp$variant_id)) # 5152 proteins; 5061 snps，说明有的 SNP 是多个蛋白的 lead pQTL(但不一定显著)
names(pp)[-1]<-paste0(names(pp)[-1],"_perm"); nrow(pp[qval_perm<=0.05]) # 34 Sig proteins 

hh<-merge(nn,pp[,.(gene_id, num_var_perm, pval_perm_perm, qval_perm, pval_nominal_threshold_perm)],by="gene_id")
hh1<-hh[pval_nominal <= pval_nominal_threshold_perm,] # 以 permutation 的 pvalue 作为 pQTL 显著性的阈值
length(unique(hh1$gene_id)) 
hh2<-hh1[qval_perm<=0.05,]
nrow(hh2); length(unique(hh2$gene_id)); length(unique(hh2$variant_id)) # 1623 Sig pair; 34 Sig protein; 1623 Sig snp. 最后用这个数据做接下来的分析
hh2<-hh2[order(gene_id,pval_nominal,-abs(slope),slope_se)]
hh2[,is_LeadpQTL:=1:.N,by=.(gene_id)]; hh2[,is_LeadpQTL:=ifelse(is_LeadpQTL==1,1,0),by=.(gene_id)] # 这里是每个 protein 选一个最显著的，后续用每个 loci 选一个最显著的，见下
names(hh2)<-c("protein_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","num_var","pval_perm_LeadpQTL","qval_protein","pval_nominal_threshold","is_LeadpQTL")
hh2<-hh2[,c("protein_id","variant_id","is_LeadpQTL","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","pval_perm_LeadpQTL","qval_protein","pval_nominal_threshold")]
hh3 <- tidyr::extract(hh2, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F)
fwrite(hh3,"~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.SigPairs_pt.txt",sep="\t",col.names = T,quote = F) 

## manhattan plot ####
library(ggmanh); library(SeqArray) 
pp<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header = T)
pp2<- tidyr::extract(pp, col = 'variant_id', into = c('chr', 'pos','ref','alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); pp2<-as.data.table(pp2)
pp2$chr<-factor(pp2$chr,levels = 1:22); pp2$pos<-as.numeric(pp2$pos)
pt<-max(pp2[qval<=0.05]$pval_nominal); pt # 5.85902e-07
pp2<-merge(pp2, name_id_reviewd, by.x="gene_id",by.y="pro_id")
pp2[,lable_name:=ifelse(qval<=0.05, gene_name, "")] # for lable sig gene_name
pp2[,color_class:=ifelse(qval<=0.05,"sig","ns")]; highlight_colormap <- c("ns" = "grey", "sig" = "maroon")
manhattan_plot(pp2, pval.colname = "pval_nominal", chr.colname = "chr", pos.colname = "pos", y.label = "-log10(P)", thin=T, # thin 可以稀疏要画的点，加快画图速度
               signif = pt, point.size=2, label.font.size=4,
               label.colname="lable_name", 
               highlight.colname = "color_class", highlight.col = highlight_colormap, color.by.highlight = TRUE) + mytheme + theme(legend.position="none")

## proportion variance explained (PVE) by lead pQTL ####
SNP_PVE1 <- function(beta, MAF, se_beta, N) { # ref: PMID: 33067605; PMID: 36797296
  numerator <- 2 * (beta^2) * MAF * (1 - MAF) }

SNP_PVE2 <- function(beta, MAF, se_beta, N) { # ref https://www.researchgate.net/post/How_to_determine_the_percent_phenotypic_variation_explained_PVE_by_a_selected_SNP; https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0120758
  numerator <- 2 * (beta^2) * MAF * (1 - MAF)
  denominator <- numerator + ((se_beta^2) * 2 * N * MAF * (1 - MAF))
  SNP_PVE <- numerator / denominator ; return(SNP_PVE) }

sample_N<-186
hh1<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.SigPairs_pt.rsID.csv", header = T) # rsID
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh2[,pve_snp1:=SNP_PVE1(slope, maf, slope_se, sample_N)]
hh2[,pve_snp2:=SNP_PVE2(slope, maf, slope_se, sample_N)] 

hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[!duplicated(protein_id)] 
nrow(hh3[pve_snp2>=0.2]); length(unique(hh3[pve_snp2>=0.2]$protein_id)) # 8 pGenes, 8 loci
ggplot(hh3, aes(x=100*pve_snp2)) + 
  geom_histogram(bins=20, fill="steelblue",color='grey90', linewidth=0.5) +
  stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), vjust = -0.5, size = 6)+
  labs(x ="Proportion of variance explained (%)", y="Numbers of lead pQTL") + mytheme


## case plot
# SNP_PVE1 算的most PVE as example (GSTM1) to plot. 趋势画出来结果不好，暂时不放这个例子。
# 放SNP_PVE2 算的most PVE as example (ALDH2 P05091) to plot. 趋势画出来挺好，但是这个 SNP 不好, 但能解释。同时画出第二高的MRPL50 Q8N5N7备用
qqnorm<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt",header = T)
melt(qqnorm[,-1:-3],id.vars = "pro_id")->qq3; names(qq3)<-c("pro_id", "sample", "value")

sel_pro<-"P05091"; sel_snp<-hh3[protein_id == sel_pro][which.max(pve_snp2), variant_id]
# in linux: 将这个 SNP 写入 tmp.bed; bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o tmp.vcf
myvcf<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/tmp.vcf",skip = 490,header = T)
my_sel2<-grep("ID|NC",names(myvcf)); myvcf[,..my_sel2]->myvcf
melt.data.table(myvcf,id.vars = "ID")->myvcf2
names(myvcf2)<-c("id","sample","format"); myvcf2[,gt:=gsub(":.+", "",format)]
myvcf2[,ds:=as.numeric(sapply(strsplit(format,split = ":"), "[", 2))] # add dosage info

qq4<-merge(qq3[pro_id==sel_pro], myvcf2[id==sel_snp, .(id,sample, gt,ds)], by="sample")
qq4[,c("ref","alt"):=.(strsplit(id,":")[[1]][3], strsplit(id,":")[[1]][4])]
qq4[,gt2:=ifelse(gt=="0|0",paste(ref,ref,sep="/"),ifelse(gt=="1|1",paste(alt,alt,sep="/"),paste(ref,alt,sep="/")))]
qq4$gt2<-factor(qq4$gt2,levels = c(unique(paste(qq4$ref,qq4$ref,sep="/")),unique(paste(qq4$ref,qq4$alt,sep="/")),unique(paste(qq4$alt,qq4$alt,sep="/"))))
table(qq4$gt) # 0|0 0|1 1|0 1|1 = 161   8  15   2

max_inten<- 1.1*max(qq4$value,na.rm = T); min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), size=2.5, alpha = 0.8,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = gt2,group=factor(gt2)), outlier.shape = NA, fill = '#00000000', width = 0.7, lwd=1.1, position = position_dodge(width = 0.85))+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  annotate("text", x=2, y=max_inten, size=8, 
           label= paste0("Beta=", round(hh3[protein_id==sel_pro & variant_id==sel_snp, slope],2), 
                         ", P-value=", hh3[protein_id==sel_pro & variant_id==sel_snp, pval_nominal]
                        
                        , ", PVE=", round(100*hh3[protein_id==sel_pro & variant_id==sel_snp, pve_snp2], 2), "%"
                        ))+
  labs(x = 'Genotype', y = 'Normalized protein abundance', 
       title = paste(name_id_reviewd[pro_id==sel_pro, gene_name], hh1[variant_id==sel_snp]$rsID, sep=" ~ "),color = NULL)+
  mytheme + theme(legend.text = element_text(size = 16),legend.position = "none",legend.title=element_blank())

## effect siz vs MAF ####
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[!duplicated(protein_id)] # 按照 conditional, 每个 pGene 其实就一个独立信号，选最显著的 SNP 作为lead pQTL
cor.test(hh3$maf, abs(hh3$slope), method = "spearman") # P=8.061e-07, rho = -0.733 

ggplot(hh3, aes(x = maf, y = abs(slope))) + geom_point(size=2, color="steelblue") + 
  #geom_errorbar(aes(ymin = abs(slope)-slope_se, ymax = abs(slope)+slope_se),width=0.01)+
  #geom_smooth(method='loess', se=T, span = 1, formula= y~x, col="maroon", linetype="dashed") + 
  geom_smooth(method='lm', se=T, span = 1, formula= y~x, col="maroon", linetype="dashed") + 
  labs(x = 'Minor allele frequency', y = 'Effect size (absolute values)') +
  coord_cartesian(ylim = c(0, 2)) + mytheme

## effect siz vs distance to TSS ####
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[!duplicated(protein_id)] # 按照 conditional, 每个 pGene 其实就一个独立信号，选最显著的 SNP 作为lead pQTL
length(unique(hh3$genome_loci)); length(unique(hh3$protein_id)) # 34 loci; 34 protein
nrow(hh3[abs(tss_distance)<=1e5]); nrow(hh3) # 28/34=82.35%; 2e5: 31/34=91.18%
ggplot(hh3, aes(x = tss_distance/1e6, y = abs(slope))) + geom_point(size=3, color="steelblue") + 
  # geom_errorbar(aes(ymin = abs(slope)-slope_se, ymax = abs(slope)+slope_se),width=0.01)+
  # geom_smooth(method='lm', se=F,formula= y~x,col="red",linetype="dashed") + 
  coord_cartesian(xlim = c(-0.5, 0.5)) + labs(x = 'Distance to TSS (MB)', y = 'Effect size (absolute values)') +
  mytheme

nrow(hh2[abs(tss_distance)<=1e5]); nrow(hh2) # 1086/1623=66.91%; 2e5: 1329/1623=81.88%
ggplot(hh2, aes(x = tss_distance/1e6, y = abs(slope))) + geom_point(size=2, color="steelblue", alpha=0.4) + 
  # geom_errorbar(aes(ymin = abs(slope)-slope_se, ymax = abs(slope)+slope_se),width=0.01)+geom_smooth(method='lm', se=F,formula= y~x,col="red",linetype="dashed") + 
  coord_cartesian(xlim = c(-1, 1)) + labs(x = 'Distance to TSS (MB)', y = 'Effect size (absolute values)') + mytheme

## conditional QTL analysis by GCTA COJO ####
jj<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_merged_GCTA_COJO.cma.cojo", header = T)

hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh3<-merge(hh2, jj[,.(SNP, p, n, freq_geno, bC, bC_se, pC)], by.x="variant_id", by.y="SNP", all.x=T) 
nrow(hh3[pC<=pval_nominal]) # 0

## trans-pQTL identification ####
xx<-fread("~/Desktop/省皮/project/pQTL/trans_qtl/skin_protein_all_qqnorm.txt.QTLtools.gz.all.trans.adjust.hits.txt.gz", header = F)
xx1<-fread("~/Desktop/省皮/project/pQTL/trans_qtl/skin_protein_all_qqnorm.txt.QTLtools.gz.all.trans.adjust.hits.sig005.txt.gz", header = F) #这是在 protein 水平上矫正过的FDR<0.05的
nrow(xx1) # P23497(SP100)~20:5746925:C:G(rs6085261) ; P23497(SP100)~20:5746954:C:T(rs6085262)

## conditional QTL analysis by QTLtools ####
con<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.conditional.txt",header = F)
names(con)<-c("gene_id","chr_gene","start_gene","end_gene","strand","num_var","tss_distance","variant_id","chr_var","start_var","end_var","signal_rank","pval_forward","slope_forward","top_forward","pf_pt","pval_reverse","slope_reverse","top_reverse","pr_pt")
table(con$signal_rank) # 全是 0，说明没有额外的独立 pQTL

#### compared with other tissue pQTL by conditional analysis ####
## protein species and abundance compare ####
name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

tmp1<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/medRxiv_5Tissues_ProAbund.xlsx", sheet = 2); tmp1<-as.data.table(tmp1); tmp1[,tissue:="colon"]
tmp2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/medRxiv_5Tissues_ProAbund.xlsx", sheet = 3); tmp2<-as.data.table(tmp2); tmp2[,tissue:="heart"]
tmp3<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/medRxiv_5Tissues_ProAbund.xlsx", sheet = 4); tmp3<-as.data.table(tmp3); tmp3[,tissue:="liver"]
tmp4<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/medRxiv_5Tissues_ProAbund.xlsx", sheet = 5); tmp4<-as.data.table(tmp4); tmp4[,tissue:="lung"]
tmp5<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/medRxiv_5Tissues_ProAbund.xlsx", sheet = 6); tmp5<-as.data.table(tmp5); tmp5[,tissue:="thyroid"]
aa<-rbind(tmp1, tmp2, tmp3, tmp4, tmp5); aa[,Ensembl_ID:=gsub("[.]\\d+","", Ensembl_ID)]
aa1<-merge(aa, unique(name_id_reviewd[,1:2]), by.x="Ensembl_ID", by.y="gene_id")
aa1<-aa1[,.(pro_id, Ref_Intensity, tissue)]; names(aa1)[2]<-"intensity"
table(aa1$tissue) # colon   heart   liver    lung thyroid = 7735    7415    8615    9009    8357 

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
mm1 <- predict(fit, mm) # 已经是log2 转化了
mm2<-as.data.table(mm1); mm2$pro_id<-rownames(mm1)
mm2_long<-melt(mm2, id.vars = "pro_id"); names(mm2_long)<-c("pro_id","sample","intensity_log")
mm2_long[, intensity:=2^intensity_log] # 将 log2 值恢复原样，以进行样本间的 aggregate
mm3<-mm2_long[,.(intensity_median=median(intensity,na.rm = T), intensity_mean=mean(intensity,na.rm = T)), by=.(pro_id)] # 同一 tissue 的不同 rep 之间 aggregate
mm3[,c("intensity_median_log", "intensity_mean_log"):=.(ifelse(is.na(intensity_median), NA, log10(intensity_median)), ifelse(is.na(intensity_mean), NA, log10(intensity_mean)))]
mm4<-mm3[,.(pro_id, intensity_median_log)]; mm4$tissue<-"skin"; names(mm4)[2]<-"intensity"

aa2<-rbind(aa1, mm4); table(aa2$tissue) # colon   heart   liver    lung    skin thyroid = 7735    7415    8615    9009    5956    8357
# 由于 skin 和其他组织的 scale 并不一致，这里添加一个系数，使之 scale 一致便于画图
aa2[,intensity_mean:= mean(intensity[tissue!="skin"], na.rm=T), by=.(pro_id)]
aa2[, intensity_skin:=intensity[tissue=="skin"], by=.(pro_id)]
aa2[,fc:=intensity_mean/intensity_skin, by=.(pro_id)];
tmp<-unique(aa2[,.(pro_id, fc)]); summary(tmp$fc) # median=10.949, mean=11.009
aa2[,intensity2:=ifelse(tissue=="skin", intensity*10.949, intensity)]
# aa2[, intensity_zs := scale(intensity), by = .(pro_id)] # 每个蛋白内部z-score
aa4<-dcast(aa2, pro_id~tissue, value.var = "intensity2")
aa4<-aa4[,.(pro_id, skin, heart, colon, thyroid, liver, lung)]

cols <- setdiff(names(aa4), "pro_id"); n <- length(cols) # 除了第一列的所有列
overlap_mat <- matrix(NA_integer_, nrow = n, ncol = n, dimnames = list(cols, cols)) # 构造一个矩阵存放结果
for (i in seq_along(cols)) {
  for (j in seq_along(cols)) {
    overlap_mat[i, j] <- sum(!is.na(aa4[[cols[i]]]) & !is.na(aa4[[cols[j]]]))
  }
}

library("UpSetR") # protein species overlap
listInput<-list(aa4[!is.na(skin)]$pro_id, aa4[!is.na(heart)]$pro_id, aa4[!is.na(colon)]$pro_id, aa4[!is.na(thyroid)]$pro_id, aa4[!is.na(liver)]$pro_id, aa4[!is.na(lung)]$pro_id)
names(listInput)<-c("Skin", "Heart", "Colon", "Thyroid", "Liver", "Lung")
upset(fromList(listInput), order.by = "freq",sets = rev(names(listInput)),keep.order=T,mainbar.y.label = "Intersection of quantified proteins", nintersects = 10, 
      sets.x.label = "Quantified proteins", point.size = 3, line.size = 1.2, text.scale = c(2, 2, 2, 2, 1.5, 3))

library(corrplot) # protein abundance
mycor<-cor(aa4[,-1], method = "spearman", use = "pairwise.complete.obs") # spearman 0.73, pearson 0.78
# tmp<-cor(jj[,-1], method = "pearson", use = "pairwise.complete.obs")
# corrplot(tmp, method = "circle",type = "upper", addCoef.col="black", cl.cex=1.5, number.digits=2, number.cex =0.5,col.lim=c(min(tmp),1), is.cor=F, order = "original", tl.col = "black", tl.srt = 45, tl.cex=0.8, mar = c(0,0,0,0))
corrplot(mycor, method = "circle",type = "lower", addCoef.col="black", cl.cex=1.5, number.digits=2, number.cex =2, col.lim=c(min(mycor),1), is.cor=F, order = "original", tl.col = "black", tl.srt = 45, tl.cex=0.8, mar = c(0,0,0,0))
mycor2<-as.data.table(mycor); mycor2[, tissue:=rownames(mycor)]

library(viridis); library(ggpointdensity)
aa4_long <- melt(as.data.table(aa4), id.vars = c("pro_id", "skin"), variable.name = "tissue", value.name = "proteome")
ggplot(aa4_long, aes(x = skin, y = proteome)) +
  geom_pointdensity(size = 1, alpha = 0.5) + scale_color_viridis(option = "D") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Skin proteome", y = "Other tissue proteome", color = NULL) +
  facet_wrap(~ tissue, nrow=1, scales = "free_y") + 
  mytheme

## overall effect size correlation ####
## in linux
ee1<-fread("~/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq/blood_pQTL_combined_all_corrected.csv", header = T)
tmp<-fread("~/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq//blood_pQTL_combined_all_corrected.hg38.bed"); tmp[,id:=paste(V4, V5, sep="_")]
ee2<-merge(ee1, tmp[,c(3, 4, 7)], by="id"); names(ee2)[ncol(ee2)-1]<-"pos_hg38"; names(ee2)[ncol(ee2)]<-"pro_id"
ee2[, pos_id_hg38:=paste(chr, pos_hg38, ref, alt, sep = ":")]; ee2[, id_hg38:=paste(pro_id, pos_id_hg38, sep="_")]
names(ee2)[13]<-"cohort"; ee2[,study:="plasma"]

ff2<-fread("~/projects//pQTL/other_pQTL_datasets/getx_5tissue/medRxiv_5Tissues_pQTL.csv", header = T) 
ff2[,variant_id:=paste(chr, pos_hg38, ref_allele, effect_allele, sep=":")]
ff2[,id_hg38:=paste(pro_id, variant_id, sep = "_")]
tmp1<-ee2[,.(pro_id, pos_id_hg38, id_hg38, study, corrected_beta, beta_se, pvalue)]
names(tmp1)<-c("pro_id", "variant_id_hg38", "id_hg38", "study", "beta", "beta_se", "pvalue")
tmp2<-ff2[,.(pro_id, variant_id, id_hg38, study, beta, beta_se, pvalue)]; 
names(tmp2)<-c("pro_id", "variant_id_hg38", "id_hg38", "study", "beta", "beta_se", "pvalue")
ef<-rbind(tmp1, tmp2)

pp<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz", header=T)
pp[,id_hg38:=paste(gene_id, variant_id, sep="_")]
ef2<-merge(pp, ef[,3:7], by="id_hg38")
fwrite(ef2, "projects/pQTL/other_pQTL_datasets/getx_5tissue/medRxiv_5Tissues_pQTL_Overlap_skin.csv", col.names=T, row.names=F, sep=",", quote=F)

ef2<-fread("~/Desktop/省皮/project/pQTL/plasma_pQTL/medRxiv_5Tissues_pQTL_Overlap_skin.csv", header = T)
ef2[study=="plasma"]$study<-"Plasma"
nrow(ef2); length(unique(ef2$gene_id)); length(unique(ef2$variant_id)) # 1883; 1082; 1786
table(ef2$study) # Colon   Heart   Liver    Lung  plasma  Thyroid = 178     185     185     102    1074     159 

## 画图
pt_list <- c(1e0, 1e-1, 5e-2, 1e-2, 1e-3, 1e-4, 1e-5) 
result <- matrix(NA, nrow = length(pt_list), ncol = 5)
colnames(result) <- c("pt", "n_pair", "n_same_direction", "pearson_r", "p_value")
for (i in 1:length(pt_list)) { # the consistent number under different pvalue cutoff
  pt <- pt_list[i]
  tmp <- ef2[pval_nominal < pt & pvalue < pt, ]
  result[i, 1]<-pt
  result[i, 2] <- nrow(tmp)
  result[i, 3] <- sum(tmp$slope * tmp$beta >= 0)
  cor_res <- cor.test(tmp$slope, tmp$beta, method = "pearson")
  result[i, 4] <- as.numeric(cor_res$estimate)
  result[i, 5] <- cor_res$p.value
  }

result<-data.table(result); result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)], id.vars = "pt"); names(mydt2)<-c("pt","type","value")
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)] #为了在一个图中画两个y坐标系，需把要画的统一至同一尺度下
mydt2$type<-factor(mydt2$type,levels = c("percent","pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = pt_list)
ggplot(mydt2[pt%in%c("1","0.05", "0.001", 1e-5)], aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(45, 100)) + labs(x = 'P-value thresholds')+
  scale_fill_manual(values = c("percent" = "indianred", "pearson_r" = "steelblue"),
                    labels = c("percent" = "Percent of same direction (%)", "pearson_r" = "Pearson correlation"))+
  mytheme  + theme(legend.position="top",legend.direction = "horizontal", legend.text=element_text(size=16,face = "bold"), legend.title=element_blank())

library(viridis); library(ggpointdensity); library(ggpubr)
pt <- 0.05; tmp<-ef2[pval_nominal<pt & pvalue<pt]
tmp1<-tmp[,.(n = .N, R = cor.test(slope, beta, method = "pearson")$estimate, P = cor.test(slope, beta, method = "pearson")$`p.value`), by=.(study)]
ggplot(tmp, mapping = aes(x = slope, y = beta)) + # 分组织画图
  geom_point(color = "black", size = 1.5, alpha = 0.8) + 
  geom_hline(yintercept=0, color="black",linetype="dashed", linewidth=0.4) + geom_vline(xintercept=0, color="black",linetype="dashed", linewidth=0.4)+
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + labs(x = 'Effect size of skin pQTL', y = 'Effect size of other pQTL', color = NULL)+
  facet_wrap(~study, nrow=1, scales = "free")+
  mytheme + theme(axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1,color="black"))  


## plasma_pQTL hg19 to hg38 ####
# linux; pwd: /sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq
# less blood_pQTL_combined_all_corrected.csv |tail -n +2|sed 's/,/\t/g'|cut -f 15|sed -e 's/:/\t/g' -e 's/_/\t/g'|awk -v OFS="\t" '{print $2, $3-1, $3, $1, $2":"$3":"$4":"$5,"+"}' > blood_pQTL_combined_all_corrected.hg19.bed 
# ~/software/liftOver blood_pQTL_combined_all_corrected.hg19.bed ~/reference/hg19ToHg38.nochr.over.chain blood_pQTL_combined_all_corrected.hg38.bed tmp.hg38.unmapped

## plink ld ####
# pwd: /sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional
# 先前是只取 pGene的显著 pQTL。这里和先前不一样，是取出 pGene 的所有 pQTL snp。
# less TableS2_skin_pQTL_allSig_pairs.csv |grep -v '#'|tail -n +2 |cut -d ',' -f 1|sort -u > tmp
# zcat skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz | head -1 > skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt
# zcat skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz | tail -n +2 | grep -F -f tmp >> skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt
pp<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt", header=T)
names(pp)[1:2]<-c("pro_id","variant_id_hg38"); length(unique(pp$variant_id_hg38)) # 100471

ee1<-fread("~/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq/blood_pQTL_combined_all_corrected.csv", header = T) 
tmp<-fread("~/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq//blood_pQTL_combined_all_corrected.hg38.bed"); tmp[,id:=paste(V4, V5, sep="_")]
ee2<-merge(ee1, tmp[,c(3, 4, 7)], by="id"); names(ee2)[ncol(ee2)-1]<-"pos_hg38"; names(ee2)[ncol(ee2)]<-"pro_id"
ee2[, pos_id_hg38:=paste(chr, pos_hg38, ref, alt, sep = ":")]; ee2[, id_hg38:=paste(pro_id, pos_id_hg38, sep="_")]
names(ee2)[13]<-"cohort"; ee2[,study:="plasma"]

ff2<-fread("~/projects//pQTL/other_pQTL_datasets/getx_5tissue/medRxiv_5Tissues_pQTL.csv", header = T) 
ff2[,variant_id:=paste(chr, pos_hg38, ref_allele, effect_allele, sep=":")]
ff2[,id_hg38:=paste(pro_id, variant_id, sep = "_")]

tmp1<-ee2[,.(pro_id, pos_id_hg38, id_hg38, study)]; names(tmp1)<-c("pro_id", "variant_id_hg38", "id_hg38", "study")
tmp2<-ff2[,.(pro_id, variant_id, id_hg38, study)]; names(tmp2)<-c("pro_id", "variant_id_hg38", "id_hg38", "study")
ef<-rbind(tmp1, tmp2) # 28824; 28630; 5657; 14840; 6
nrow(ef); length(unique(ef$id_hg38)); length(unique(ef$pro_id)); length(unique(ef$variant_id_hg38)); length(unique(ef$study)) 

common_pro <- intersect(ef$pro_id, pp$pro_id); length(common_pro) # 28/34 pGenes共同被鉴定到
ef_sub <- ef[pro_id %in% common_pro, .(pro_id, variant_ef = variant_id_hg38)]
pp_sub <- pp[pro_id %in% common_pro, .(pro_id, variant_pp = variant_id_hg38)]
pairs <- merge(pp_sub, ef_sub, by = "pro_id", allow.cartesian = TRUE) # 每个蛋白内部，SNP 全排列
pairs<-unique(pairs); nrow(pairs) # 593767

all_snps <- unique(c(pairs$variant_pp, pairs$variant_ef))
length(all_snps); length(unique(pp_sub$variant_pp)); length(unique(ef_sub$variant_ef)) # 83365; 83253; 183
# plink2 --bfile /sh2/home/sunyuanqiang/reference/G1000_hg38/all_pop/ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged --extract all_snps_forLD.txt --make-bed --out /sh2/home/sunyuanqiang/reference/G1000_hg38/all_pop/ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged.snps_forLD
# https://zzz.bwh.harvard.edu/plink/ld.shtml: 
# plink --bfile ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged.snps_forLD --extract all_snps_forLD.txt --r2 gz --ld-window-kb 1000000 --ld-window 100000 --ld-window-r2 0.5 --out all_snps_forLD

ld<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/all_snps_forLD.ld.gz", header=T); ld<-ld[,c(3, 6, 7)]
pairs[, pair_id := paste(variant_pp, variant_ef, sep="_")]
ld[,pair_id:=paste(SNP_A, SNP_B, sep="_")]
ld[,pair_id2:=paste(SNP_B, SNP_A, sep="_")]
res<-merge(pairs, ld[,c(3,4)], by.x="pair_id", by.y="pair_id", all.x=T)
res2<-merge(res, ld[,c(3,5)], by.x="pair_id", by.y="pair_id2", all.x=T)

res2[,R2:=ifelse(!is.na(R2.x), R2.x, ifelse(!is.na(R2.y), R2.y, 0))]
res2[,R2:=ifelse(variant_pp==variant_ef, 1, R2)] # plink 计算 ld 时，相等的SNP 之间不会计算，赋值为 1
res2<-res2[, c(1:4, 7)]
length(unique(res2$pro_id)); nrow(res2); nrow(pairs) # 28; 593767=593767
fwrite(res2, "~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_6Tissue_pGenes_all_ld_r2_full.csv", col.names=T, row.names=F, sep="\t", quote=F)

## conditional by COJO ####
res2<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_6Tissue_pGenes_all_ld_r2_full.csv", header=T)
res2[,pro_variant_pp:=paste(pro_id, variant_pp,sep="_")]; res2[,pro_variant_ef:=paste(pro_id, variant_ef,sep="_")]
pp[,pro_variant_pp :=paste(pro_id, variant_id_hg38,sep="_")]
res3<-merge(res2, pp[,.(pro_variant_pp, maf, pval_nominal, slope)])
res4<-merge(res3, ef[,.(id_hg38, study)], by.x="pro_variant_ef", by.y = "id_hg38", all.x = TRUE, allow.cartesian = TRUE); res4<-res4[order(study, pro_id, -R2)]
res5<-res4[R2>=0.8, 4:11]; nrow(res5); length(unique(res5$pro_id)) # 4483; 27
# 切割 conditional SNP
outdir <- "/sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional"
res5[, {
  filename <- paste(pro_id, variant_pp, variant_ef, study, sep = "_")
  filepath <- file.path(outdir, filename)
  writeLines(variant_pp, filepath)
  NULL
}, by = seq_len(nrow(res5))]

# ll *:*|wc -l # 4483
# 切割 skin pQTL for common pGenes
outdir <- "/sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional"
N_val <- 186
pp_sub2[, {
  tmp <- tstrsplit(variant_id_hg38, ":", fixed = TRUE)
  chr_pos <- paste0(tmp[[1]], ":", tmp[[2]], ":", tmp[[3]], ":", tmp[[4]])
  alt <- tmp[[4]]; ref <- tmp[[3]]
  # 重新组织成COJO输入格式
  df <- data.table(SNP = chr_pos, A1 = alt, A2 = ref, fre = maf, b = slope, # A1 是 effect allele, A2是ref allele,  A1/A2不能颠倒。freq是 A1 的频率 
                   se = slope_se, p = format(pval_nominal, scientific = TRUE), N = N_val)
  filename <- paste0(unique(pro_id), "_skin"); filepath <- file.path(outdir, filename) # 文件名, 存贮路径
  fwrite(df, filepath, sep = "\t", quote = FALSE, col.names = TRUE) # 写文件
  NULL}, 
  by = pro_id
]

if (1) { #最初跑的时候设置--diff-freq 0.5，--cojo-collinear 0.8，会有一些 badsnps 和共线性过高的 snps没法 condition，后续取消限制
  # case test
  # /sh2/home/sunyuanqiang/software/gcta/gcta64 --bfile /sh2/home/sunyuanqiang/projects/pQTL/formal/conditional_COJO/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.plink.nochr --cojo-file Q9Y6N5_skin --cojo-cond Q9Y6N5_15:45676237:T:C_15:45676237:T:C_Thyroid --diff-freq 1 --cojo-wind 100000 --cojo-collinear 0.99 --out Q9Y6N5_15:45676237:T:C_15:45676237:T:C_Thyroid.COJO.conditiona
  # 批量运行
  # 生成批量运行的脚本：bash generate_COJO_sh
  # Warning: 找不到 P01860 对应的 cojo-cond 文件; 这是因为生成cojo-cond 文件，要求 other_SNP在 skin_SNP中必须有R2>0.8的。P01860虽然在 other tissue也被检测到，但是other_SNP与skin_SNP没有R2的配对
  # summary(res2[pro_id=="P01860" & variant_ef=="14:105707924:C:G"]$R2); summary(res2[pro_id=="P01860" & variant_ef=="14:105773531:C:G"]$R2)都为 0
  # bash generate_COJO_sh会产生run_cojo_tasks.sh脚本文件
  # 添加 & 和 wait: awk '{print $0 " &"; if (NR % 10 == 0) print "wait"} END {if (NR % 10 != 0) print "wait"}' run_cojo_tasks.sh > run_cojo_tasks_parallel.sh
  # nohup bash run_cojo_tasks_parallel.sh & 生成的结果 move 到~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res里，diff-freq 0.5，--cojo-collinear 0.8的结果在其下的DifFre_05_R2_08目录下
}

## process COJO-conditional res ####
## 整合多个蛋白的单个文件到一个文件
library(data.table); library(stringr)
myfiles <- list.files(path="~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res", pattern = "COJO\\.conditional\\.cma\\.cojo$", full.names = TRUE)
length(myfiles) # 4483; "/sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res/O95479_1:9231607:T:C_1:9246790:T:C_plasma.COJO.conditional.cma.cojo"
cc <- rbindlist(lapply(myfiles, function(f) { # 整合多个蛋白的单个文件到一个文件
  dt <- fread(f, header=T); fname <- basename(f) 
  parts <- str_split(fname, pattern = "_", simplify = TRUE) # 用下划线拆分
  dt[, pro_id := parts[1]]; dt[, conditional_SNP_skin := parts[2]]; dt[, conditional_SNP_other := parts[3]]; dt[, study := str_split(parts[4], pattern = "\\.", simplify = TRUE)[1]]   # 添加四列
  return(dt)
}))
nrow(cc) # 13156982
cc[, bC:=ifelse(is.na(bC), 0, bC)]; cc[, bC_se:=ifelse(is.na(bC_se), 1, bC_se)]; cc[, pC:=ifelse(is.na(pC), 1, pC)] 

myfiles2 <- list.files(path="~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res", pattern = "COJO\\.conditional\\.given\\.cojo$", full.names = TRUE)
length(myfiles2) # 4483; "/sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res/O95479_1:9231607:T:C_1:9246790:T:C_plasma.COJO.conditional.given.cojo"

res_list <- lapply(myfiles2, function(f) { # 批量处理
  fname <- basename(f); parts <- strsplit(fname, "_")[[1]]
  pro_id <- parts[1]; cond_skin <- parts[2]; cond_other <- parts[3]; study <- sub("\\.COJO.*", "", parts[4])
  dt <- fread(f)
  dt[, `:=`(n = 186, freq_geno = freq, bC = 0, bC_se = 1, pC = 1, pro_id = pro_id, # conditional snp本身pC设为 1，bC设为 0
            conditional_SNP_skin = cond_skin, conditional_SNP_other = cond_other, study = study)]
  return(dt)
})

cc2 <- rbindlist(res_list, fill = TRUE); nrow(cc2) # 4483
identical(cc2$SNP, cc2$conditional_SNP_skin) # TRUE
cc3<-rbind(cc, cc2)

cc3[, pair_id:=paste(conditional_SNP_skin, conditional_SNP_other, sep="_")]; nrow(cc3) # 13161465
res2<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_6Tissue_pGenes_all_ld_r2_full.csv", header=T); nrow(res2) # 593767
length(unique(cc3$pair_id)); length(unique(res2$pair_id)) # 4019, 593767
cc4<-merge(cc3, unique(res2[,.(pair_id, R2)]), by="pair_id", all.x=T) # 添加 R2 信息

pp<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt", header=T)
pp1<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pp2<-merge(pp, pp1[,.(gene_id, pval_nominal_threshold)], by="gene_id") # 对每个 pGene 添加 pvalue threshold

pp3<-merge(pp2[,-3:-5], cc4[,.(pro_id, SNP, conditional_SNP_skin, conditional_SNP_other, R2, bC, bC_se, pC, study)], by.x=c("gene_id","variant_id"), by.y=c("pro_id","SNP"), all.x=T)
pp3[,bC:=ifelse(is.na(bC), slope, bC)]; pp3[,bC_se:=ifelse(is.na(bC_se), slope_se, bC_se)]; pp3[,pC:=ifelse(is.na(pC), pval_nominal, pC)]
pp3[,conditional_SNP_skin:=ifelse(is.na(conditional_SNP_skin), "none", conditional_SNP_skin)]; pp3[,conditional_SNP_other:=ifelse(is.na(conditional_SNP_other), "none", conditional_SNP_other)]
pp3[,study:=ifelse(is.na(study), "none", study)]; pp3[,R2:=ifelse(is.na(R2), 2, R2)]
nrow(pp3[is.na(pC)]); table(pp3$study) # 0; Colon   Heart   Liver    Lung    none  plasma Thyroid

pp3 <- tidyr::extract( pp3, col = 'variant_id', into = c('chr', 'pos', 'ref', 'alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); pp3<-as.data.table(pp3)
pp3 <- tidyr::extract( pp3, col = 'conditional_SNP_skin', into = c('conditional_SNP_skin_chr', 'conditional_SNP_skin_pos', 'conditional_SNP_skin_ref', 'conditional_SNP_skin_alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); pp3<-as.data.table(pp3)
pp3[, pos:=as.numeric(pos)]; pp3[, conditional_SNP_skin_pos:=as.numeric(conditional_SNP_skin_pos)] # 便于后面画 manhantton 图

# 添加 condition SNP在other tissue的统计量
ee1<-fread("~/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq/blood_pQTL_combined_all_corrected.csv", header = T) 
tmp<-fread("~/projects/pQTL/other_pQTL_datasets/plasma_pQTL_syq/blood_pQTL_combined_all_corrected.hg38.bed"); tmp[,id:=paste(V4, V5, sep="_")]
ee2<-merge(ee1, tmp[,c(3, 4, 7)], by="id"); names(ee2)[ncol(ee2)-1]<-"pos_hg38"; names(ee2)[ncol(ee2)]<-"pro_id"
ee2[, pos_id_hg38:=paste(chr, pos_hg38, ref, alt, sep = ":")]; ee2[, id_hg38:=paste(pro_id, pos_id_hg38, sep="_")]
names(ee2)[13]<-"cohort"; ee2[,study:="plasma"]
ee2<-ee2[, .(study, pro_id, pos_id_hg38, corrected_beta, pvalue)]; names(ee2)[4]<-"beta"
ff2<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/medRxiv_5Tissues_pQTL.csv", header = T) 
ff2[,pos_id_hg38:=paste(chr, pos_hg38, ref_allele, effect_allele, sep=":")]
ff2<-ff2[, .(study, pro_id, pos_id_hg38, beta, pvalue)]
ef2<-rbind(ee2, ff2); ef2<-unique(ef2)
names(ef2)[4:5]<-c("conditional_SNP_other_beta", "conditional_SNP_other_pvalue")
ef2[,id_tmp:=paste(study, pro_id, pos_id_hg38, sep="_")]

pp3[,id_tmp:=paste(study, gene_id, conditional_SNP_other, sep="_")]
pp4<-merge(pp3, ef2[,4:6], by="id_tmp", all.x=T) # 此前conditional_SNP_skin/other/study 填充为 none 的 snp，conditional_SNP_other_beta/pvalue都为 NA
pp4[, conditional_SNP_other_beta:=ifelse(is.na(conditional_SNP_other_beta), 0, conditional_SNP_other_beta)]
pp4[, conditional_SNP_other_pvalue:=ifelse(is.na(conditional_SNP_other_pvalue), 1, conditional_SNP_other_pvalue)] # NA分别填补 beta 和 p 为 0 和 1
pp4<-pp4[, id_tmp:=NULL]
fwrite(pp4, "~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.condition.csv", col.names=T, row.names=F, sep=",", quote=F)

## condition 分析画图 ####
library(ggplot2); library(patchwork)
pp4<-fread("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.condition.csv", header=T)
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
pp4<-merge(pp4, unique(name_id_reviewd[,c(1,4)]), by.x="gene_id", by.y="pro_id", all.x=T)

length(unique(pp4$gene_id)); length(unique(pp4[pval_nominal<=pval_nominal_threshold]$gene_id)) # 34=34
length(unique(ef2$pro_id)) # 5657
length(intersect(unique(pp4$gene_id), unique(ef2$pro_id))) # 28个pGenes 在其他 tissue 中检测到了
length(setdiff(unique(pp4$gene_id), unique(ef2$pro_id))) # 6个pGenes 在其他 tissue 中没检测到, unique_pro

## 统计每个蛋白的每个 tissue condition 前后的最显著的P值
aa2 <- pp4[, { # 后半部分（每 gene_id × tissue × conditional_SNP_skin分组)
  # 当前  gene_id × tissue × conditional_SNP_skin 下最小 pval_nominal
  min_pN <- min(pval_nominal, na.rm = TRUE)
  row_min_pN <- .SD[pval_nominal == min_pN][1]
  
  # skin_least_snp_id 对应的最大 pC 
  df_least <- .SD[variant_id == row_min_pN$variant_id]
  max_pC <- max(df_least$pC, na.rm = TRUE)
  row_max_pC <- df_least[pC == max_pC]
  if (nrow(row_max_pC) > 1) row_max_pC <- row_max_pC[which.max(R2)]
  row_max_pC <- row_max_pC[1]
  
  # 当前 gene × tissue × conditional_SNP_skin 下最小 pC
  min_pC <- min(pC, na.rm = TRUE)
  row_min_pC <- .SD[pC == min_pC]
  if (nrow(row_min_pC) > 1) row_min_pC <- row_min_pC[which.max(R2)]
  row_min_pC <- row_min_pC[1]
  
  vals <- c( # 动态生成列名和值
    row_min_pN$variant_id, row_min_pN$pval_nominal, row_min_pN$slope, row_min_pN$pval_nominal_threshold,
    row_max_pC$pC, row_max_pC$bC, row_max_pC$conditional_SNP_other,
    row_min_pC$pC, row_min_pC$conditional_SNP_other,  row_min_pC$bC, row_min_pC$R2)
  
  names(vals) <- c(
    "skin_least_snp_id", "skin_least_snp_pN", "skin_least_snp_bN", "skin_least_snp_SigPt",
    "skin_least_snp_AftCond_most_pC", "skin_least_snp_AftCond_most_bC", "skin_least_snp_AftCond_most_ConSNP_other",
    "skin_AftCond_least_pC", "skin_AftCond_least_pC_ConSNP_other", "skin_AftCond_least_pC_bC", "skin_AftCond_least_pC_ConSNP_R2")
  
  as.list(vals) }, by = .(gene_id, study, conditional_SNP_skin)]

aa2$skin_AftCond_least_pC<-as.numeric(aa2$skin_AftCond_least_pC); 
aa2$skin_AftCond_least_pC_ConSNP_R2<-as.numeric(aa2$skin_AftCond_least_pC_ConSNP_R2)
aa2<-aa2[order(gene_id, study, -skin_AftCond_least_pC, -skin_AftCond_least_pC_ConSNP_R2)]
aa3 <- aa2[, .SD[1], by = .(gene_id, study)] # 每个gene × tissue中，取每个condition_SNP_skin的最小 pC中的最大 pC， 作为该gene × tissue的代表性pC

aa31<-unique(aa3[,.(gene_id, skin_least_snp_pN, skin_least_snp_bN, skin_least_snp_SigPt)]) # 取出 condition前 skin 中的 lead SNP的 P 和 beta
aa31_long<-melt(aa31, id.vars = "gene_id")
aa32<-aa3[,.(gene_id, study, skin_AftCond_least_pC, skin_AftCond_least_pC_bC)] # 取出各个组织condition后的ConSNP_skin、 P 和 beta
aa32_wide <- dcast(aa32, gene_id ~ study, value.var = setdiff(names(aa32), c("gene_id", "study")))
aa33<-melt.data.table(aa32_wide, id.vars = "gene_id") # 先 dcast 再 melt 是为了全排列所有的 gene × tissue情况，防止有些 gene 只有少数几种 tissue 的数据
aa4<-rbind(aa31_long, aa33); aa4<-aa4[order(gene_id, variable)]
 # 此时仍为 NA 的tissue证明 condition 的时候没有 它没有对应的 SNP hit 到，所以填补为 condition 之前的 p 和 beta 
aa4[, value_imp := ifelse(grepl("^skin_AftCond_least_pC_[^_]+$", variable) & is.na(value), value[variable == "skin_least_snp_pN"][1], 
                          ifelse(grepl("skin_AftCond_least_pC_bC_", variable) & is.na(value), value[variable == "skin_least_snp_bN"][1], value)), by = gene_id]
aa4[, study := fifelse(grepl("AftCond", as.character(variable)), sapply(strsplit(as.character(variable), "_"), tail, 1), "Skin")] # 补充 tissue 信息
aa4[, type := ifelse(grepl("skin_least_snp_pN$|skin_AftCond_least_pC_[^_]+$", variable), "pvalue",
                     ifelse(grepl("skin_least_snp_bN$|skin_AftCond_least_pC_bC_", variable), "beta",
                            ifelse(grepl("^skin_least_snp_SigPt$", variable), "sig_pt", NA_character_)))] # 命名统计量类型
aa4<-aa4[order(gene_id, study, type)]; aa4$value_imp<-as.numeric(aa4$value_imp)

aa6 <- dcast(aa4, gene_id + study ~ type, value.var = "value_imp")
aa6[, sig_pt:=ifelse(is.na(sig_pt), sig_pt[study=="Skin"], sig_pt), by=.(gene_id)]; 
aa6[study == "plasma", study := "Plasma"]
names(aa6)[2]<-"tissue"; aa6<-aa6[tissue!="none"]; aa6[, condition:=ifelse(tissue=="Skin", "BeforeCondition", "AfterCondition")]
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
aa6<-merge(aa6, unique(name_id_reviewd[,c(1,4)]), by.x="gene_id", by.y="pro_id")
fwrite(aa6, "~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/TableSX_skin_pGenes_condition_by6Tissue.csv", row.names = F, col.names = T, )

## 本地气泡图
aa6<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableSX_skin_pGenes_condition_by6Tissue.csv", header = T)
aa6[, sig_mark := ifelse(pvalue <= sig_pt, "*", "")] # 星号标注
aa6[, logp_bin := cut(-log10(pvalue), breaks = seq(0, ceiling(max(-log10(pvalue), na.rm = TRUE))+3, by = 4),
                      include.lowest = TRUE, right = FALSE )] # pvalue每两个数量级一个层级，便于画图
n <- length(levels(aa6$logp_bin)); values <- seq(2, by = 3, length.out = n) # 用于画点的大小
 # 把 tissue 和 gene 排序
aa6[, sig_num_tissue:=sum(pvalue <= sig_pt), by=.(gene_id)]; 
aa6[, sig_num_gene:=sum(pvalue <= sig_pt), by=.(tissue)]
gene_order <- unique(aa6[tissue == "Skin"][order(sig_num_tissue, -pvalue), gene_name])
tissue_order<-unique(aa6[order(-sig_num_gene), tissue])
aa6$tissue<-factor(aa6$tissue, levels = tissue_order); aa6$gene_name<-factor(aa6$gene_name, levels = gene_order)

if (1) { # 纵向图
  ggplot(aa6, aes(x = tissue, y = gene_name)) +
    geom_point(aes(size = logp_bin, color = beta)) +
    scale_size_manual(values = values, drop = FALSE) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B22222", midpoint = 0, na.value = "grey80") +
    geom_text(aes(label = sig_mark), vjust = 0.8, size = 10, color = "black") +
    labs(x = NULL, y = NULL, size = NULL, color = NULL) +
    scale_x_discrete(position = "top")+
    theme_minimal(base_size = 14) +
    theme(axis.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(face = "bold", size = 12, color = "black", angle = 0, vjust = 0),
          axis.text.y = element_text(face = "bold.italic", size = 12, color = "black")
          # , legend.position = "bottom", legend.box = "horizontal"
          # ,legend.position = "none"
    ) 
  # + guides(size = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1), color = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = 15))
}


## 本地画case manhatton
# 结论
length(unique(aa6[sig_num_tissue==7,]$gene_name)) # 12/34=35.3; 若在所有 tissue (other+skin=7) 中pC 的最小值都显著，则证明没有被 condition掉，是 skin-unique的
unique_pro<-setdiff(unique(aa6$gene_id), unique(ef2$pro_id)) # 6个是因为 other tissue没报告；Q9NWV4(CZIB)、Q15369(ELOC)、Q5T440(IBA57)、P07476(IVL)、Q8N5N7(MRPL50)、Q9HBL7(PLGRKT)
 # 剩余 6 个是other tissue报告了这个基因，但是没有 condition 成功: ACADSB/ALDH2/FCGR3A/GAA/GSTM1/IGHG3
tmp<-unique(aa6[,.(tissue, sig_num_gene)]); tmp[,conditioned_num:= 34 - sig_num_gene] # 每个 tissue 中有多少 pGene被condition掉了
tmp2<-unique(aa6[,.(gene_name, sig_num_tissue)]); tmp2[,conditioned_num:= 7 - sig_num_tissue] # 每个 pGene在多少tissue被condition掉了

# 画图： skin/other tissue都检测到，并且 condition 成功的，CPNE1, leadpQTL=20:35762424:G:A, skin_p=6.78446e-12, pC=2.28673e-03
sel_pro<-"CPNE1"
for (i in 1:length(sel_pro)) {
  # On linux
  tmp<-pp4[gene_name==sel_pro[1]]; 
  tmp<-tmp[order(gene_id, pval_nominal, -R2, -pC)] # tmp<-tmp[order(gene_id, pval_nominal, -pC, -R2)] #这样排序可能会选出距离峰值太远的conditional SNP
  conditional_SNP_skin_most <- tmp[1, conditional_SNP_skin]
  conditional_SNP_other_most <- tmp[1, conditional_SNP_other]
  conditional_SNP_skin_most_study <- tmp[1, study]
  tmp2 <- tmp[conditional_SNP_skin==conditional_SNP_skin_most & study == conditional_SNP_skin_most_study & conditional_SNP_other == conditional_SNP_other_most]; 
  tmp2<-tmp2[order(pos)]
  summary(tmp2$pval_nominal); summary(tmp2$pC) # 0.0000002 0.3063260 0.6069860 0.5531334 0.8049790 0.9998380; 0.005706 0.320231 0.562920 0.534607 0.727564 1.000000 
  fwrite(tmp2, paste0("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/Condition_manhaton_", sel_pro,".csv"), col.names = T, row.names = F, sep=",", quote = F)
  # On local
  sel_pro<-"CPNE1"
  tmp2<-fread(paste0("~/Desktop/省皮/project/pQTL/plasma_pQTL/conditional/Condition_manhaton_", sel_pro,".csv"), header = T)
  x_limits <- range(-log10(tmp2$pval_nominal), -log10(tmp2$pC), na.rm = TRUE) # -log10(P) 的最小值和最大值
  conditional_SNP_skin_most <- unique(tmp2$conditional_SNP_skin)
  conditional_SNP_other_most <-  unique(tmp2$conditional_SNP_other)
  conditional_SNP_skin_most_study <-  unique(tmp2$study)
  condition_snp_pos <- unique(tmp2[conditional_SNP_skin==conditional_SNP_skin_most]$conditional_SNP_skin_pos)
  condition_snp_p <- unique(tmp2[variant_id==conditional_SNP_skin_most]$pval_nominal)
  
  p2 <- ggplot(tmp2, aes(y = pos, x = -log10(pC))) +
    geom_point(aes(color = slope), size = 0.6, alpha = 0.8) + 
    geom_point(aes(x = -log10(condition_snp_p), y = condition_snp_pos), shape = 4, color = "firebrick", size = 4, alpha = 1) + # 标注 condition的 SNP
    geom_segment(aes(x = 0, xend = -log10(condition_snp_p), y = condition_snp_pos, yend = condition_snp_pos), linetype = "dashed", color = "black", linewidth = 0.1)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = range(tmp2$slope, na.rm = TRUE)) +
    scale_x_reverse(limits = rev(x_limits)) +
    scale_y_continuous(position = "right", labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = NULL, color = "Effect size") +
    mytheme + 
    theme(legend.position = c(0.07, 0.07), legend.justification = c("left", "bottom"), legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key.size = unit(0.8, "lines"), legend.title = element_text(size = 14), legend.text = element_text(size = 13),
          axis.title.y = element_blank(), axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank(), axis.line.y.left = element_blank(),axis.text.y.right = element_blank(), axis.ticks.y.right = element_line(),axis.line.y.right = element_line())
  
  p1 <- ggplot(tmp2, aes(y = pos, x = -log10(pval_nominal))) +
    geom_point(aes(color = slope), size = 0.6, alpha = 0.8) + 
    geom_point(aes(x = -log10(condition_snp_p), y = condition_snp_pos), shape = 4, color = "firebrick", size = 4, alpha = 1) + # 标注 condition的 SNP
    geom_text(aes(x = -log10(condition_snp_p), y = condition_snp_pos, label = paste(conditional_SNP_skin_most,  conditional_SNP_skin_most_study, sep="_")), hjust = 1.2, vjust = -12, color = "black", size = 3)+
    geom_segment(aes(x = 0, xend = -log10(condition_snp_p), y = condition_snp_pos, yend = condition_snp_pos), linetype = "dashed", color = "black", linewidth = 0.1)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = range(tmp2$slope, na.rm = TRUE), guide = "none") +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(position = "left",  labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = NULL) +
    mytheme + 
    theme(axis.title.y = element_blank(), axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(),
          axis.ticks.y.left = element_line(), axis.text.y.left = element_text(),axis.line.y.left = element_line())
  
  print(p2 | p1 + plot_layout(widths = c(1, 0)))  
}

for (i in 1:length(sel_pro)) { # 配合论文草稿用rs12481228 20:35630751:G:C作 condition SNP, in plasma
  # On linux
  tmp<-pp4[gene_name==sel_pro[1] & conditional_SNP_skin=="20:35630751:G:C" & conditional_SNP_other=="20:35630751:G:C" & study=="plasma"]
  tmp2<-tmp[order(pos)]
  summary(tmp2$pval_nominal); summary(tmp2$pC) # 0.0000002 0.3063260 0.6069860 0.5531334 0.8049790 0.9998380; 0.005706 0.320231 0.562920 0.534607 0.727564 1.000000 
  fwrite(tmp2, paste0("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/Condition_manhaton_", sel_pro,".csv"), col.names = T, row.names = F, sep=",", quote = F)
  # On local
  sel_pro<-"CPNE1"
  tmp2<-fread(paste0("~/Desktop/省皮/project/pQTL/plasma_pQTL/conditional/Condition_manhaton_", sel_pro,".csv"), header = T)
  x_limits <- range(-log10(tmp2$pval_nominal), -log10(tmp2$pC), na.rm = TRUE) # -log10(P) 的最小值和最大值
  conditional_SNP_skin_most <- unique(tmp2$conditional_SNP_skin)
  conditional_SNP_other_most <-  unique(tmp2$conditional_SNP_other)
  conditional_SNP_skin_most_study <-  unique(tmp2$study)
  condition_snp_pos <- unique(tmp2[conditional_SNP_skin==conditional_SNP_skin_most]$conditional_SNP_skin_pos)
  condition_snp_p <- unique(tmp2[variant_id==conditional_SNP_skin_most]$pval_nominal)
  
  p2 <- ggplot(tmp2, aes(y = pos, x = -log10(pC))) +
    geom_point(aes(color = slope), size = 0.6, alpha = 0.8) + 
    geom_point(aes(x = -log10(condition_snp_p), y = condition_snp_pos), shape = 4, color = "firebrick", size = 4, alpha = 1) + # 标注 condition的 SNP
    geom_segment(aes(x = 0, xend = -log10(condition_snp_p), y = condition_snp_pos, yend = condition_snp_pos), linetype = "dashed", color = "black", linewidth = 0.1)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = range(tmp2$slope, na.rm = TRUE)) +
    scale_x_reverse(limits = rev(x_limits)) +
    scale_y_continuous(position = "right", labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = NULL, color = "Effect size") +
    mytheme + 
    theme(legend.position = c(0.07, 0.07), legend.justification = c("left", "bottom"), legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key.size = unit(0.8, "lines"), legend.title = element_text(size = 14), legend.text = element_text(size = 13),
          axis.title.y = element_blank(), axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank(), axis.line.y.left = element_blank(),axis.text.y.right = element_blank(), axis.ticks.y.right = element_line(),axis.line.y.right = element_line())
  
  p1 <- ggplot(tmp2, aes(y = pos, x = -log10(pval_nominal))) +
    geom_point(aes(color = slope), size = 0.6, alpha = 0.8) + 
    geom_point(aes(x = -log10(condition_snp_p), y = condition_snp_pos), shape = 4, color = "firebrick", size = 4, alpha = 1) + # 标注 condition的 SNP
    geom_text(aes(x = -log10(condition_snp_p), y = condition_snp_pos, label = paste(conditional_SNP_skin_most,  conditional_SNP_skin_most_study, sep="_")), hjust = 1.2, vjust = -12, color = "black", size = 3)+
    geom_segment(aes(x = 0, xend = -log10(condition_snp_p), y = condition_snp_pos, yend = condition_snp_pos), linetype = "dashed", color = "black", linewidth = 0.1)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = range(tmp2$slope, na.rm = TRUE), guide = "none") +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(position = "left",  labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = NULL) +
    mytheme + 
    theme(axis.title.y = element_blank(), axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(),
          axis.ticks.y.left = element_line(), axis.text.y.left = element_text(),axis.line.y.left = element_line())
  
  print(p2 | p1 + plot_layout(widths = c(1, 0)))  
}


# 画图： skin/other tissue都检测到，但没有 condition 成功的。ALDH2, leadpQTL=20:35762424:G:A, skin_p=6.78446e-12, pC=2.28673e-03
sel_pro<-"ALDH2"
for (i in 1:length(sel_pro)) {
  # On linux
  tmp<-pp4[gene_name==sel_pro[1]]; 
  tmp<-tmp[order(gene_id, pval_nominal, -R2, -pC)] # tmp<-tmp[order(gene_id, pval_nominal, -pC, -R2)] #这样排序可能会选出距离峰值太远的conditional SNP
  conditional_SNP_skin_most <- tmp[1, conditional_SNP_skin]
  conditional_SNP_other_most <- tmp[1, conditional_SNP_other]
  conditional_SNP_skin_most_study <- tmp[1, study]
  tmp2 <- tmp[conditional_SNP_skin==conditional_SNP_skin_most & study == conditional_SNP_skin_most_study & conditional_SNP_other == conditional_SNP_other_most]; 
  tmp2<-tmp2[order(pos)]
  summary(tmp2$pval_nominal); summary(tmp2$pC) # 0.0000002 0.3063260 0.6069860 0.5531334 0.8049790 0.9998380; 0.005706 0.320231 0.562920 0.534607 0.727564 1.000000 
  fwrite(tmp2, paste0("~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/Condition_manhaton_", sel_pro,".csv"), col.names = T, row.names = F, sep=",", quote = F)
  # On local
  sel_pro<-"ALDH2"
  tmp2<-fread(paste0("~/Desktop/省皮/project/pQTL/plasma_pQTL/conditional/Condition_manhaton_", sel_pro,".csv"), header = T)
  x_limits <- range(-log10(tmp2$pval_nominal), -log10(tmp2$pC), na.rm = TRUE) # -log10(P) 的最小值和最大值
  conditional_SNP_skin_most <- unique(tmp2$conditional_SNP_skin)
  conditional_SNP_other_most <-  unique(tmp2$conditional_SNP_other)
  conditional_SNP_skin_most_study <-  unique(tmp2$study)
  condition_snp_pos <- unique(tmp2[conditional_SNP_skin==conditional_SNP_skin_most]$conditional_SNP_skin_pos)
  condition_snp_p <- unique(tmp2[variant_id==conditional_SNP_skin_most]$pval_nominal)
  
  p2 <- ggplot(tmp2, aes(y = pos, x = -log10(pC))) +
    geom_point(aes(color = slope), size = 0.8, alpha = 0.8) + 
    geom_point(aes(x = -log10(condition_snp_p), y = condition_snp_pos), shape = 4, color = "firebrick", size = 4, alpha = 1) + # 标注 condition的 SNP
    geom_segment(aes(x = 0, xend = -log10(condition_snp_p), y = condition_snp_pos, yend = condition_snp_pos), linetype = "dashed", color = "black", linewidth = 0.1)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = range(tmp2$slope, na.rm = TRUE)) +
    scale_x_reverse(limits = rev(x_limits)) +
    scale_y_continuous(position = "right", labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = NULL, color = "Effect size") +
    mytheme + 
    theme(legend.position = c(0.07, 0.07), legend.justification = c("left", "bottom"), legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key.size = unit(0.8, "lines"), legend.title = element_text(size = 14), legend.text = element_text(size = 13),
          axis.title.y = element_blank(), axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank(), axis.line.y.left = element_blank(),axis.text.y.right = element_blank(), axis.ticks.y.right = element_line(),axis.line.y.right = element_line())
  
  p1 <- ggplot(tmp2, aes(y = pos, x = -log10(pval_nominal))) +
    geom_point(aes(color = slope), size = 0.8, alpha = 0.8) + 
    geom_point(aes(x = -log10(condition_snp_p), y = condition_snp_pos), shape = 4, color = "firebrick", size = 4, alpha = 1) + # 标注 condition的 SNP
    geom_text(aes(x = -log10(condition_snp_p), y = condition_snp_pos, label = paste(conditional_SNP_skin_most,  conditional_SNP_skin_most_study, sep="_")), hjust = 1.2, vjust = -12, color = "black", size = 3)+
    geom_segment(aes(x = 0, xend = -log10(condition_snp_p), y = condition_snp_pos, yend = condition_snp_pos), linetype = "dashed", color = "black", linewidth = 0.1)+
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0, limits = range(tmp2$slope, na.rm = TRUE), guide = "none") +
    scale_x_continuous(limits = x_limits) + 
    scale_y_continuous(position = "left",  labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    labs(x = NULL, y = NULL) +
    mytheme + 
    theme(axis.title.y = element_blank(), axis.ticks.y.right = element_blank(), axis.text.y.right = element_blank(), axis.line.y.right = element_blank(),
          axis.ticks.y.left = element_line(), axis.text.y.left = element_text(),axis.line.y.left = element_line())
  
  print(p2 | p1 + plot_layout(widths = c(1, 0)))  
}

######################### skin RNA ######################### 
######## Preapre RNA matrix ########
## QC of skin RNA matrix ####
rr<-fread("/Volumes/disk2s1/pQTL/2023_10_wangzhenzhen_GNS2305142_gaoy_analysis_report_lncRNA_gene_count/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # a total table contain raw cnt, log2Normalized_cnt and FPKM value
sel_col<-grep("gene_ID|.+SRNA.+fpkm",names(rr),value=T); rr[,..sel_col]->rr_rpkm; names(rr_rpkm)<-gsub("SRNA_fpkm", "",names(rr_rpkm))
sel_col<-grep("gene_ID|.+SRNA$",names(rr),value=T); rr[,..sel_col]->rr_cnt; names(rr_cnt)<-gsub("SRNA","",names(rr_cnt))
names(rr_rpkm)[1]<-"gene_id";names(rr_cnt)[1]<-"gene_id"
rm(rr) # large memory load
identical(names(rr_rpkm),names(rr_cnt)) # TRUE

## select smps of pheno (RNA) overlapped with geno samoles (150 smp)
my_sel<-names(rr_rpkm)%in%c("gene_id",setdiff(sel_smp$V1,c("NC285","NC060"))) # skin这里要额外排除两个 smp
rr_rpkm[,..my_sel]->rr_rpkm1; rr_cnt[,..my_sel]->rr_cnt1; ncol(rr_rpkm1)-1 # 150 smps

## select genes of autosomes and non-predicting
tss_all<-fread("~/Desktop/reference/ensembl_allGene_pos_biomart.txt",header = T) # 最开始的 tss 信息只是针对 pro-coding
tss_all[,strand:=ifelse(strand==1,"+","-")]
tss_all2<-tss_all[,.(chr=chromosome,start=transcription_start_site-1,end=transcription_start_site,id=gene_id,strand=strand)]
tss_all2<-tss_all2[order(-id,chr,start),]
tss_all2[,sel_tss:=ifelse(length(unique(chr))>1,start[1],ifelse(strand=="+",min(start),max(start))),by=.(id)]
tss_all3<-tss_all2[start==sel_tss,.SD,by=.(id)];
tss_all3<-tss_all3[!duplicated(id),];
tss_all3<-tss_all3[id!="",.(chr,start,end,id,strand)]

my_sel2<- tss_all3[paste0("chr",chr)%in%sel_chr$V1,id]; length(my_sel2) # 59578 autosome genes
rr_rpkm1[gene_id%in%my_sel2,]->rr_rpkm1; rr_cnt1[gene_id%in%my_sel2,]->rr_cnt1; nrow(rr_rpkm1) # 57327 genes
rr_rpkm1[grep("MSTR",gene_id,invert = T),]->rr_rpkm1; rr_cnt1[grep("MSTR",gene_id,invert = T),]->rr_cnt1

## Adjustment:  RPKM >0.1 in at least 10 smps; summed cnt > 5 for each gene. Adapted from ref: GTEx PMID: 32913072: >10smps, >6cnt total
smp_num <- ncol(rr_rpkm1)-1
rr_cnt1[,ave_cnt:=rowMeans(as.matrix(rr_cnt1[,2:151]))]; rr_cnt1[,sum_cnt:=rowSums(as.matrix(rr_cnt1[,2:151]))]
melt(rr_rpkm1,id.vars = "gene_id")->rr_rpkm2
names(rr_rpkm2)<-c("gene_id","sample","value")
rr_rpkm2[,rpkm01_gene:=sum(value>0.1),by=.(gene_id)]
length(unique(rr_rpkm2[rpkm01_gene > (smp_num/2 -1), gene_id])) # 30429/28022/22600 for 10/20/half_sample
length(unique(rr_rpkm2[gene_id%in%rr_cnt1[ave_cnt>5,gene_id],gene_id])) # 27605/23772/20968 for ave_ant 2/5/10
rr_rpkm3<-rr_rpkm2[rpkm01_gene >= 10]; length(unique(rr_rpkm3$gene_id)) # 30775
rr_rpkm3<-rr_rpkm3[gene_id%in%rr_cnt1[sum_cnt >= 6,gene_id],]; length(unique(rr_rpkm3$gene_id)) # 30775

#RPKM >0.1 in > half smps; average cnt > 3 in each smp. 
#rr_rpkm3<-rr_rpkm2[rpkm01_gene > (smp_num/2 -1)]; length(unique(rr_rpkm3$gene_id)) # 19616
#rr_rpkm3<-rr_rpkm3[gene_id%in%rr_cnt1[ave_cnt>3,gene_id],]; length(unique(rr_rpkm3$gene_id)) # 19029

## Normalization of skin RNA matrix  ####
#和蛋白不一样的是，不存在缺失值，只有 0，本来不需要填补,但是后面 log0会有异常值，在那加 psedo cnt,不如这里直接填补
if (1) { rr_rpkm3[,min:=min(value[value!=0],na.rm = T),by=.(gene_id)]
  rr_rpkm3[,value:=ifelse(value==0,min/2,value)] }
rr_rpkm3[,value_log:=ifelse(is.na(value),NA,log2(value))] # log transformed
rr_rpkm3[,value_log_scale:=as.numeric(scale(value_log,center = T,scale = T)),by=.(gene_id)] # scale for each protein across sample
rr_rpkm3[,value_log_scale_rank:=rank(value_log_scale,na.last="keep",ties.method = "average"),by=.(sample)] # rank within each sample across  proteins
rr_rpkm3[,value_log_scale_rank_qnorm:= qnorm((value_log_scale_rank-0.5)/sum(!is.na(value_log_scale_rank))),by=.(sample)] # quantile normalized (rank inverse normalized) within each smp
rr_rpkm4<-dcast(rr_rpkm3,gene_id~sample,value.var = "value_log_scale_rank_qnorm") 

## add TSS info to RNA matrix ####
tmp<-merge(rr_rpkm4,tss_all3[,.(chr,start,end,id)],by.x = "gene_id",by.y="id"); nrow(tmp); ncol(tmp)-4 # 30775 genes and 150smp
if (1) { # for inspect the proteins without TSS info and manually add them to tss original file iteratively
  tmp2<-merge(rr_rpkm4,tss_all2[,.(chr,start,end,id)],by.x = "gene_id",by.y="id",all.x=T) 
  tmp2[is.na(tmp2$start),]->tmp3; nrow(tmp3) # 0个基因没匹配上坐标
}
my_sel2<-c(ncol(tmp)-2,ncol(tmp)-1,ncol(tmp),1,2:(ncol(tmp)-3))
rr_rpkm5<-tmp[,..my_sel2] # adjust column order to bed format
rr_rpkm5[,chr:=paste0("chr",chr)]
rr_rpkm5<-rr_rpkm5[order(chr,start)] # 后续index 需要排序
rr_rpkm6<-copy(rr_rpkm5); rr_rpkm6<-rr_rpkm6[chr%in%sel_chr$V1]; dim(rr_rpkm6) # 选择常染色体上的基因。这步可省略，QC那一步已经完成
names(rr_rpkm6)[1]<-"#Chr";
fwrite(rr_rpkm6, "~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt", sep = "\t",col.names = T)

######## skin RNA covariats processing ########
library(findPC)
pca<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca",header = T)
names(pca)[1]<-"PC_id"
pca[,PC_id:=paste0("geno_",gsub("hg38.+svd_", "",PC_id))]
pca1<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt.forPCA.gz.pca",header = T)
names(pca1)[1]<-"PC_id"
pca1[,PC_id:=paste0("pheno_",gsub("skin.+svd_", "",PC_id))]

# select K PCs for use
pve<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca_stats",header = F,nrows = 3)
pve<-transpose(pve,keep.names = "PC",make.names = "V1");pve[,PC:=1:.N];pve[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve$PC,pve$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve$prop_var,number = c(30,60,100),method = 'all',aggregate = 'median',figure = T)
sel_genoPCs<-5

pve1<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt.forPCA.gz.pca_stats",header = F,nrows = 3)
pve1<-transpose(pve1,keep.names = "PC",make.names = "V1");pve1[,PC:=1:.N];pve1[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve1$PC,pve1$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve1$prop_var,number = c(30,60,100,150),method = 'all',aggregate = 'median',figure = T)
sel_phenoPCs<-10

transpose(pca,keep.names = "sample", make.names = "PC_id")->tmp
transpose(pca1,keep.names = "sample", make.names = "PC_id")->tmp1
merge(tmp1[,1:(sel_phenoPCs+1)],tmp[,1:(sel_genoPCs+1)],by="sample")->tmp2 
merge(tmp2,sample_info[,.(No,gender,age,Other_disease, sample_pos)],by.x="sample",by.y="No",all.x = T)->tmp3 # 皮肤需考虑采样部位; all geno_pheno overlapped sample should be retained. The missing values should be encoded as NA;手动尽量填补NA 
tmp3[,gender:=ifelse(gender=="male",1,0)]
Other_disease<-tmp3$Other_disease; Other_disease <- model.matrix(~Other_disease); Other_disease<-dplyr::as_tibble(Other_disease)
sample_pos<-tmp3$sample_pos; sample_pos <- model.matrix(~sample_pos); sample_pos<-dplyr::as_tibble(sample_pos)

tmp3<-cbind(tmp3,Other_disease[,-1], sample_pos[,-1])
tmp3[,Other_disease:=NULL]; tmp3[,sample_pos:=NULL]; str(tmp3)

transpose(tmp3,keep.names = "PC_id", make.names = "sample")->tmp4
tmp4[is.na(tmp4)]<-"NA" #  The missing values should be encoded as NA 
fwrite(tmp4,"~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt.gz.12pheno_5geno.PCs",col.names = T,quote = F,sep = "\t")

#### eQTL charactering ####
## 计算显著的 SNP-RNA pair 数目 ####
# all pirs太大，在服务器跑
library(data.table)
name_id_reviewd<-fread("~/reference/uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header=T)
length(unique(nn$gene_id)); length(unique(nn$variant_id))
pp<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
length(unique(pp$gene_id)); length(unique(pp$variant_id))
names(pp)[-1]<-paste0(names(pp)[-1],"_perm")

hh<-merge(nn,pp[,.(gene_id,num_var_perm,pval_beta_perm,qval_perm,pval_nominal_threshold_perm)],by="gene_id")
hh1<-hh[pval_nominal <= pval_nominal_threshold_perm,]
hh2<-hh1[qval_perm<=0.05,]
nrow(hh2); length(unique(hh2$gene_id)); length(unique(hh2$variant_id)) # 24788 Sig pair; 282 Sig RNA,少于permu结果直接用qval<0.05算的284，因为ENSG00000272031、ENSG00000224106两个基因虽然满足 qval<0.05，但是不满足 nominal p < pt，不知道为啥，后续用 282 的结果; 20999 Sig snp. 最后用这个数据做接下来的分析
hh2<-hh2[order(gene_id,pval_nominal,-abs(slope),slope_se)]
hh2[,is_LeadpQTL:=1:.N,by=.(gene_id)]; hh2[,is_LeadpQTL:=ifelse(is_LeadpQTL==1,1,0),by=.(gene_id)];
names(hh2)<-c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","num_var","pval_beta","qval","pval_nominal_threshold","is_LeadpQTL")
hh2<-hh2[,c("gene_id","variant_id","is_LeadpQTL","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","pval_beta","qval","pval_nominal_threshold")]
hh3 <- tidyr::extract(hh2, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F)
fwrite(hh3,"~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.SigPairs_pt.txt",sep="\t",col.names = T,quote = F) 
fwrite(hh3,"~/projects/pQTL/formal/TableS3_skin_eQTL_allSig_pairs.csv",sep=",",col.names = T,quote = F) # then download on local and remove in server

## Counts of overlapped self eGenes vs GTEx eGenes ####
library(VennDiagram)
pp<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
names(pp)[1]<-"pro_id"
length(pp$pro_id); length(unique(pp$pro_id)); nrow(pp[qval<=0.05]) # 30470; 30470; 284

gtex<-fread("/Volumes/disk3s1/pQTL/GTEx_Analysis_v8_eQTL/Skin_Sun_Exposed_Lower_leg.v8.egenes.txt.gz") # # GTEx 605 samples; hg38 coordinates;"contain data for all genes tested (only top variant-gene pair). to obtain the list of eGenes, select the rows with 'qval' ≤ 0.05."
length(gtex$gene_id); length(unique(gtex$gene_id)); nrow(gtex[qval<=0.05]) # 25196; 25196; 16967
names(gtex)[1]<-"pro_id"; gtex[, pro_id:=gsub("[.]\\d+", "", pro_id)]; gtex<-gtex[order(pval_nominal)]
# All tested genes
venn.plot <-venn.diagram(list(unique(pp$pro_id),unique(gtex$pro_id)), inverted=T, resolution=600, filename=NULL,lwd=1.5,cex=1.5,cat.fontface=2,na="remove",cat.pos=180,cat.cex=1.2,fill=rainbow(2),category.names=c("This study","GTEx"))
grid.newpage(); grid.draw(venn.plot) # 3735; 21461; 9009
# All eGenes
venn.plot <-venn.diagram(list(unique(pp[qval<0.05]$pro_id),unique(gtex[qval<0.05]$pro_id)),resolution=600,filename=NULL,lwd=1.5,cex=1.5,cat.fontface=2,na="remove",cat.pos=c(315,45),cat.cex=1.2,fill=rainbow(2),category.names=c("This study","GTEx"))
grid.newpage(); grid.draw(venn.plot) # 16750, 217, 67
# eGenes among the overlapped tested genes
total<-merge(pp, gtex, by="pro_id")
venn.plot <-venn.diagram(list(unique(total[qval.x<0.05]$pro_id),unique(total[qval.y<0.05]$pro_id)),resolution=600,filename=NULL,lwd=1.5,cex=1.5,cat.fontface=2,na="remove",cat.pos=c(315,45),cat.cex=1.2,fill=rainbow(2),category.names=c("This study","GTEx"))
grid.newpage(); grid.draw(venn.plot) # 14706, 217, 9

plot(-log10(total$pval_nominal.x),-log10(total$pval_nominal.y),pch=20,col=gray(0.5,0.5),xlim=c(0,20),ylim=c(0,300),cex=1,xlab = "This study",ylab = "GTEx",las=1)
lines(0:300,0:300,col="red",lty=5); 
cor.test(-log10(total$pval_nominal.x),-log10(total$pval_nominal.y),method = "pearson") # 0.307

## Effect size of self eQTL vs gtex eQTL ####
# 太大了，在服务器上跑
nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T)
names(nn)[1]<-"pro_id"; nn<-nn[order(pro_id,pval_nominal,-abs(slope),slope_se)] #当Pvalue 也一致时，选 slope 绝对值较大的，如果还一样则选 slop_se较小的

gtex_nn<-fread("~/reference/GTEx_Analysis_v8_eQTL/Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs.txt.gz")  # 2414653 sig pairs
gtex_nn[,variant_id:=gsub("_b38","",gsub("chr","",variant_id))]; gtex_nn[,variant_id:=gsub("_",":",variant_id)]
gtex_nn[,gene_id:=gsub("[.]\\d+","",gene_id)]; names(gtex_nn)[2]<-"pro_id"
gtex_nn<-gtex_nn[order(gene_id,pval_nominal)]

total<-merge(nn, gtex_nn, by=c("pro_id","variant_id")); nrow(total) # 1478275
fwrite(total,"~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.OverlapGTEx.txt",sep="\t",quote=F) # gzip skin_RNA_all_qqnorm.txt.gz.allpairs.OverlapGTEx.txt 

pt_list <- c(1e-2, 1e-4, 1e-6, 1e-8)
result <- matrix(NA, nrow = 4, ncol = 5)
colnames(result) <- c("pt", "n_pair", "n_same_direction", "pearson_r", "p_value")
for (i in 1:length(pt_list)) {
  pt <- pt_list[i]
  tmp <- total[pval_nominal.x < pt & pval_nominal.y < pt, ]
  result[i, 1]<-pt
  result[i, 2] <- nrow(tmp)
  result[i, 3] <- sum(tmp$slope.x * tmp$slope.y >= 0)
  cor_res <- cor.test(tmp$slope.x, tmp$slope.y, method = "pearson")
  result[i, 4] <- as.numeric(cor_res$estimate)
  result[i, 5] <- cor_res$p.value
}

# 在本地画图
result<-data.table(pt=c(1e-2, 1e-4, 1e-6, 1e-8), n_pair=c(157655, 35725, 11741,3219),
                   n_same_direction=c(151162, 35292, 11718, 3219),pearson_r=c(0.8567206, 0.9016582, 0.9111776, 0.9428516),
                   p_value=c(0, 0, 0, 0)) # from linux R 
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)],id.vars = "pt"); names(mydt2)<-c("pt","type","value")
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)] #为了在一个图中画两个y坐标系，需把要画的统一至同一尺度下
mydt2$type<-factor(mydt2$type,levels = c("percent","pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = c(1e-02, 1e-04, 1e-06, 1e-08))
ggplot(mydt2, aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(85, 100)) + labs(x = 'P-value thresholds')+
  scale_fill_manual(values = c("percent" = "indianred", "pearson_r" = "steelblue"),
                    labels = c("percent" = "Percent of same direction (%)", "pearson_r" = "Pearson correlation"))+
  mytheme  + theme(legend.position="top",legend.direction = "horizontal", legend.title=element_blank())

oo<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt.gz.allpairs.OverlapGTEx.txt.gz", header = T); nrow(oo) # 1478275
pt <- 1e-04
ggplot(oo[pval_nominal.x<pt & pval_nominal.y<pt], mapping = aes(x = slope.x, y = slope.y)) +
  geom_pointdensity(size = 0.8) + scale_color_viridis() +
  labs(x = 'Effect size in our study', y = 'Effect size in GTEx', color = NULL)+
  # stat_cor(method = "pearson", digits=3,size=6)+
  coord_cartesian(xlim = c(-2, 2),ylim = c(-2, 2)) +
  geom_hline(yintercept=0, color="black",linetype="dashed", size=0.4)+geom_vline(xintercept=0, color="black",linetype="dashed", size=0.4)+
  geom_smooth(method='lm', se=T, formula= y~x, col="red", linetype="dashed") +
  mytheme + theme(legend.position = c(0.25, 0.9), legend.direction = "horizontal", legend.key.width = unit(1.2, "cm"), legend.text = element_text(size = 12)) # legend.position是相对位置, 范围必须在 [0,1] 之间

######## pQTL and eQTL replication ########
## Abundance correlation of RNA  vs protein ####
rr<-fread("/Volumes/disk3s1/pQTL/2023_10_wangzhenzhen_GNS2305142_gaoy_analysis_report_lncRNA_gene_count/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # a total table contain raw cnt, log2Normalized_cnt and FPKM value
sel_col<-grep("gene_ID|.+SRNA.+fpkm",names(rr),value=T); rr[,..sel_col]->rr_rpkm; names(rr_rpkm)<-gsub("SRNA_fpkm","",names(rr_rpkm))
sel_col<-grep("gene_ID|.+SRNA$",names(rr),value=T); rr[,..sel_col]->rr_cnt; names(rr_cnt)<-gsub("SRNA","",names(rr_cnt))
names(rr_rpkm)[1]<-"gene_id";names(rr_cnt)[1]<-"gene_id"
identical(names(rr_rpkm),names(rr_cnt)) # TRUE
rr_rpkm<-rr_rpkm[grep('ENSG',gene_id)]; rr_cnt<-rr_cnt[grep('ENSG',gene_id)] # 去掉公司预测的基因
identical(rr_rpkm$gene_id, rr_cnt$gene_id); nrow(rr_cnt) # TRUE; 60564

non_zero_sample_cols <- colSums(rr_rpkm[, -1, with = FALSE]) != 0
rr_rpkm <- rr_rpkm[, c(TRUE, non_zero_sample_cols), with = FALSE] # 去掉全为0的样本列
rr_rpkm <- rr_rpkm[rowSums(rr_rpkm[, -1, with = FALSE]) != 0]; dim(rr_rpkm) # 去掉全为0的基因行; 52080 gene; 208-1 smp

## defining expressed gene:  RPKM >0.1 in at least 10 smps; summed cnt > 5 for each gene. Adapted from ref: GTEx 29022597
rr_cnt[, `:=`(mean_expr = rowMeans(.SD), sum_expr = rowSums(.SD)), .SDcols = !'gene_id'] # # 计算每个基因的 cnt 的平均值和总和（不包括 gene_id 列）
melt(rr_rpkm, id.vars = "gene_id")->rr_rpkm2; names(rr_rpkm2)<-c("gene_id","sample","value")
rr_rpkm2[, rpkm01_gene:=sum(value>0.1), by=.(gene_id)]
rr_rpkm3 <- rr_rpkm2[rpkm01_gene >= 10]; length(unique(rr_rpkm3$gene_id)) # 33359；RPKM >0.1 in at least 10 smps
rr_rpkm3 <- rr_rpkm3[gene_id %in% rr_cnt[sum_expr >= 6, gene_id],]; length(unique(rr_rpkm3$gene_id)); length(unique(rr_rpkm3$sample)) # summed cnt > 5 for each gene; 33359,207

## rho for overlapped of samples and genes
name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
length(unique(name_id_reviewd$pro_id)); length(unique(name_id_reviewd$gene_id)) # 42433; 32814

qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; names(qq3)<-gsub("SP", ", ", names(qq3))
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
summary(jj1$rho) # Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  -0.332780 -0.007551  0.073984  0.086145  0.163403  0.806744 
summary(jj1[protein_log_mean>0 & rna_log_mean>0]$rho) # -0.331574 -0.009725  0.070044  0.081667  0.157217  0.806744 
summary(jj1[rho_p<0.05]$rho) # Min.   1st Qu.    Median      Mean   3rd Qu.      -0.3328  0.1607  0.2078  0.2044  0.2811  0.8067
nrow(jj1); nrow(jj1[rho>0]); nrow(jj1[rho<0]); nrow(jj1[rho_p<0.05 & rho>0]); nrow(jj1[rho_p<0.05 & rho<0]); nrow(jj1[rho_p>0.05]) # 5860; 4277; 1583; 1795; 201; 3864

ggplot(jj1, aes(x = "", y = rho)) + # violin
  geom_violin(width = 0.65, adjust = 1, fill = "#F35859", color=NA, linewidth=1.2) +
  geom_boxplot(width = 0.2, outlier.size = 1, outlier.alpha = 1, fatten = 2, lwd=1.2) +
  labs(x = NULL, y = "Spearman rho")+
  coord_cartesian(ylim = c(-0.4, 1))+
  mytheme  + theme(legend.position="none")

# overall
cor.test(jj1$protein_log_median, jj1$rna_log_median, method = "spearman") # 0.289; 0.210 for pearson
cor.test(jj1$protein_log_mean, jj1$rna_log_mean, method = "spearman") # 0.293; 0.23 for pearson
library(viridis); library(ggpointdensity); library(ggpubr)
nrow(jj1[protein_log_mean>0 & rna_log_mean>0]) # 5127
ggplot(data = jj1[protein_log_mean>0 & rna_log_mean>0], mapping = aes(x = protein_log_mean, y = rna_log_mean)) +
  geom_pointdensity() + scale_color_viridis() +
  labs(x = 'Protein abundance', y = 'RNA expression', color = NULL)+
  geom_smooth(method='lm', se=T,formula= y~x,col="red", linetype="dashed")+  
  # stat_cor(method = "spearman", label.x = 10, label.y = -4, digits=3,size=6, label.sep = "\n")+
  coord_cartesian(xlim = c(0, 11), ylim = c(0, 7))+
  mytheme 

## Count of QTL overlap ####
library(data.table); library(VennDiagram)
# eQTL 数据太大， 在服务器上跑
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
pro_nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); names(pro_nn)[1]<-"pro_id"
names(pro_nn)[3:9]<-paste0(names(pro_nn)[3:9], "_pqtl")

rna_nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); nrow(rna_nn) # 93,916,515 pairs
rna_nn1<-merge(rna_nn, name_id_reviewd, by="gene_id", all.x=T)
names(rna_nn1)[3:9]<-paste0(names(rna_nn1)[3:9], "_eqtl")

nn<-merge(pro_nn, rna_nn1, by=c("pro_id", "variant_id"))
fwrite(nn, "~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt", col.names=T, row.names=F, sep="\t", quote=F)
# gzip skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt 

# 置换检验 overlap 数目的显著性。在服务器上
pt<-1e-2
pqtl_pairs <- unique(paste(nn[pval_nominal_pqtl < pt]$pro_id, nn[pval_nominal_pqtl < pt]$variant_id, sep = "_"))
pqtl_n <- length(pqtl_pairs) # 175905
eqtl_pairs <- unique(paste(nn[pval_nominal_eqtl < pt]$pro_id, nn[pval_nominal_eqtl < pt]$variant_id, sep = "_"))
eqtl_n <- length(eqtl_pairs) # 210203 
true_overlap <- length(intersect(pqtl_pairs, eqtl_pairs)) # 7657
background_pairs <- unique(paste(nn$pro_id, nn$variant_id, sep = "_")); length(background_pairs) # 15384030

n_perm <- 1000; set.seed(123)
perm_overlap <- numeric(n_perm) # 初始化一个长度为 n_perm 的数值型向量 perm_overlap，并用 0 填充初始值。
for (i in 1:n_perm) {
  random_pqtl <- sample(background_pairs, pqtl_n)
  random_eqtl <- sample(background_pairs, eqtl_n)
  perm_overlap[i] <- length(intersect(random_pqtl, random_eqtl))
}
mean(perm_overlap) # 2404.629
empirical_p <- (sum(perm_overlap >= true_overlap) + 1) / (length(perm_overlap) + 1) # 9.9e-4
enrichment <- true_overlap / mean(perm_overlap) # 7657/2404.629 = 3.184275
fwrite(data.frame(overlap = perm_overlap), "~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.permutation.txt", col.names=T, row.names=F, sep="\t", quote=F)
tmp<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.permutation.txt")
ggplot(tmp, aes(x=overlap)) + 
  geom_histogram(bins=50, fill="steelblue",color='grey90', linewidth=0.5) +
  # geom_vline(xintercept = 7657, color = "red", linetype = "solid", linewidth = 1)+ #这个画在图上会导致左边压缩太厉害。后续手动标记
  labs(x = "Overlap count under permutations", y = "Frequency")+
  mytheme

## correaltion of beta ####
nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt.gz", header = T); nrow(nn); length(unique(nn$pro_id)); length(unique(nn$gene_id)) # 15443226; 5052; 5073
pt_list <- c(1e-2, 1e-3, 1e-4, 1e-6) 
result <- matrix(NA, nrow = 4, ncol = 5)
colnames(result) <- c("pt", "n_pair", "n_same_direction", "pearson_r", "p_value")
for (i in 1:length(pt_list)) { # the consistent number under different pvalue cutoff
  pt <- pt_list[i]
  tmp <- nn[pval_nominal_pqtl < pt & pval_nominal_eqtl < pt, ]
  result[i, 1]<-pt
  result[i, 2] <- nrow(tmp)
  result[i, 3] <- sum(tmp$slope_pqtl * tmp$slope_eqtl >= 0)
  cor_res <- cor.test(tmp$slope_pqtl, tmp$slope_eqtl, method = "pearson")
  result[i, 4] <- as.numeric(cor_res$estimate)
  result[i, 5] <- cor_res$p.value
}

# 在本地画图
result<-data.table(pt=c(0.05, 1e-2, 1e-4, 1e-6), n_pair=c(60765, 7660, 826, 132),
                   n_same_direction=c(43350, 6767, 826, 132),pearson_r=c(0.4405654, 0.7236292, 0.8143522, 0.9464664),
                   p_value=c(0.000000e+00, 0.000000e+00, 6.702596e-197, 1.070058e-65)) # from linux R 
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)], id.vars = "pt"); names(mydt2)<-c("pt","type","value")
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)] #为了在一个图中画两个y坐标系，需把要画的统一至同一尺度下
mydt2$type<-factor(mydt2$type,levels = c("percent","pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = c(0.05, 1e-2, 1e-4, 1e-6))
ggplot(mydt2, aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(45, 100)) + labs(x = 'P-value thresholds')+
  scale_fill_manual(values = c("percent" = "indianred", "pearson_r" = "steelblue"),
                    labels = c("percent" = "Percent of same direction (%)", "pearson_r" = "Pearson correlation"))+
  mytheme  + theme(legend.position="top",legend.direction = "horizontal", legend.text=element_text(size=16,face = "bold"), legend.title=element_blank())

library(viridis); library(ggpointdensity); library(ggpubr)
oo<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", header = T); nrow(oo) # 1478275
pt <- 1e-02
ggplot(oo[pval_nominal_pqtl<pt & pval_nominal_eqtl<pt], mapping = aes(x = slope_pqtl, y = slope_eqtl)) +
  geom_pointdensity(size = 1.5) + scale_color_viridis() +
  labs(x = 'Effect size of pQTL', y = 'Effect size of eQTL', color = NULL)+
  # stat_cor(method = "pearson", digits=3,size=6)+
  coord_cartesian(xlim = c(-1.2, 1.2),ylim = c(-1.2, 1.2)) +
  geom_hline(yintercept=0, color="black",linetype="dashed", linewidth=0.4) + geom_vline(xintercept=0, color="black",linetype="dashed", linewidth=0.4)+
  geom_smooth(method='lm', se=T, formula= y~x, col="red", linetype="dashed") +
  mytheme + theme(legend.position = c(0.25, 0.9), legend.direction = "horizontal", legend.key.width = unit(1.2, "cm"), legend.text = element_text(size = 14)) # legend.position是相对位置, 范围必须在 [0,1] 之间

## eQTL and pQTL opposite direction ####
oo<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", header = T); nrow(oo) # 1478275
pp1<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pp2<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
oo<-merge(oo, pp1[,.(gene_id, pval_nominal_threshold)], by="gene_id"); names(oo)[ncol(oo)]<-"pval_nominal_threshold_eqtl"
oo<-merge(oo, pp2[,.(gene_id, pval_nominal_threshold)], by.x = "pro_id", by.y ="gene_id"); names(oo)[ncol(oo)]<-"pval_nominal_threshold_pqtl"
oo[,p_sum:= -log10(pval_nominal_pqtl) -log10(pval_nominal_eqtl)]

oo1<-oo[slope_pqtl*slope_eqtl<0,]; oo1<-oo1[order(p_sum)]; nrow(oo1) # 893
nrow(oo1[pval_nominal_pqtl <= pval_nominal_threshold_pqtl & pval_nominal_eqtl <= pval_nominal_threshold_eqtl]) # 0 没有显著相反的点

library(ggplot2); library(scales)  # 用于色彩渐变和轴标签
ggplot(oo1, aes(x = -log10(pval_nominal_pqtl), y = -log10(pval_nominal_eqtl), color = slope_pqtl)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_hline(yintercept = max(-log10(oo1$pval_nominal_threshold_eqtl)), linetype = "dashed", color = "black") +
  geom_vline(xintercept = max(-log10(oo1$pval_nominal_threshold_pqtl)), linetype = "dashed", color = "black") + #标记显著性阈值
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "pQTL effect size", guide = guide_colorbar(title.position = "top", barwidth = 8, barheight = 1)) +
  labs(x = "-log10(P-value) of pQTL", y = "-log10(P-value) of eQTL", color = "pQTL effect size") +
  theme_minimal(base_size = 18) +  theme(legend.position = "top", legend.direction = "horizontal", legend.title.align = 0.5,  legend.title = element_text(size = 14), legend.text  = element_text(size = 14))

## RegulomeDB_rank compare ####
regulome<-fread("/Volumes/ElementsSE_syq/reference/RegulomeDB_rank.tsv", header = T) # https://regulomedb.org/regulome-search/
regulome<-regulome[,1:6]; regulome[, chrom:=gsub("chr", "", chrom)]; names(regulome)[c(1,3,6)]<-c("chr_var", "pos_var", "score")

pro<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T); pro[, chr_var:=as.character(chr_var)]
pro1<-merge(pro, regulome[, c(1, 3:6)], by=c("chr_var","pos_var")); nrow(pro); nrow(pro1) # 1623; 1622
pro1[,ranking2:= as.integer(sub("^(\\d).*", "\\1", ranking))]
summary(pro1$score) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 0.0000  0.3975  0.5544  0.4983  0.5544  1.0000 
nrow(pro1[ranking2<=2]); nrow(pro1[ranking2>2]) # 1259; 363; /1622=77.6%; 22.4%
pro1<-pro1[order(protein_id, pval_nominal)]; pro2<-pro1[!duplicated(protein_id)]; nrow(pro2) # 34,只考虑 lead pQTL
nrow(pro2[ranking2<=2]); nrow(pro2[ranking2>2]) # 23; 11; /34=67.6%; 32.4%。比例比 all_SigPair低

rna<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS3_skin_eQTL_allSig_pairs.csv", skip=3, header = T); rna[, chr_var:=as.character(chr_var)]
rna1<-rna[!duplicated(variant_id)] # 20999/24788
rna1<-merge(rna1, regulome[, c(1, 3:6)], by=c("chr_var","pos_var")); nrow(rna1) # 20938/20999
rna1[,ranking2:= as.integer(sub("^(\\d).*", "\\1", ranking))]
summary(rna1$score) # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 0.0000  0.2625  0.5139  0.4739  0.5544  1.0000
nrow(rna1[ranking2<=2]); nrow(rna1[ranking2>2]) # 13049; 7889; /20938=62.3%; 37.7%
nrow(rna1[is_LeadeQTL==1]); nrow(rna1[is_LeadeQTL==1 & ranking2<=2]); nrow(rna1[is_LeadeQTL==1 & ranking2>2]) # 884; 577; 307=65.3%; 34.7%。与 pqtl 保持一致，不用

# 背景variants 的富集比例。在服务器
regulome<-fread("~/reference/RegulomeDB_rank.tsv", header = T)
regulome<-regulome[,1:6]; regulome[, chrom:=gsub("chr", "", chrom)]; names(regulome)[c(1,3,6)]<-c("chr_var", "pos_var", "score")
regulome[,ranking2:= as.integer(sub("^(\\d).*", "\\1", ranking))]
class(regulome$chr_var); class(regulome$pos_var) # "character"; "integer"

pro_nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T)
pro_nn <- tidyr::extract(pro_nn, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); pro_nn<-as.data.table(pro_nn)
pro_nn<-pro_nn[order(variant_id, pval_nominal)]; pro_nn1<-pro_nn[!duplicated(variant_id)]
length(pro_nn$variant_id); length(pro_nn1$variant_id) # 15,669,687; 3,461,609
class(pro_nn1$chr_var); class(pro_nn1$pos_var)
pro_nn1[, pos_var := as.integer(pos_var)]
pro_nn2<-merge(pro_nn1, regulome[, -2], by=c("chr_var","pos_var")); nrow(pro_nn2) # 3,455,459
summary(pro_nn2$score) # 0.0000  0.1955  0.5139  0.4372  0.5896  1.0000   # 与 sig_pair没有差异，就不化 score 的 boxplot 了
nrow(pro_nn2[ranking2<=2]); nrow(pro_nn2[ranking2>2]) # 1,515,102; 1,940,357; /3455459=43.85%; 56.15%
chisq.test(matrix(c(1259, 363, 1515102, 1940357), byrow = T, nrow = 2))$`p.value` # 77.6%/43.85%=1.77, 4.767061e-165

pro_nn<-pro_nn[order(gene_id, pval_nominal)]; pro_nn2<-pro_nn[!duplicated(gene_id)]; nrow(pro_nn2) # 5152, 只考虑 lead variant
pro_nn2[, pos_var := as.integer(pos_var)]
pro_nn3<-merge(pro_nn2, regulome[, -2], by=c("chr_var","pos_var")); nrow(pro_nn3) # 5103
nrow(pro_nn3[ranking2<=2]); nrow(pro_nn3[ranking2>2]) # 2678; 2425; /5103=52.5%; 47.5%
chisq.test(matrix(c(23, 11, 2678, 2425), byrow = T, nrow = 2))$`p.value` # 67.6%/52.5%=1.77, 0.11; 0.086 for fisher test


rna_nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T)
rna_nn <- tidyr::extract(rna_nn, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); rna_nn<-as.data.table(rna_nn)
rna_nn<-rna_nn[order(variant_id, pval_nominal)]; rna_nn1<-rna_nn[!duplicated(variant_id)]
length(rna_nn$variant_id); length(rna_nn1$variant_id) # 93,916,515; 4,390,105
class(rna_nn1$chr_var); class(rna_nn1$pos_var)
rna_nn1$pos_var<-as.numeric(rna_nn1$pos_var)
rna_nn2<-merge(rna_nn1, regulome[, -2], by=c("chr_var","pos_var")); nrow(rna_nn2) # 4382618
summary(rna_nn2$score) #  0.0000  0.1841  0.5139  0.4185  0.5877  1.0000   # 与 sig_pair没有差异，就不化 score 的 boxplot 了
nrow(rna_nn2[ranking2<=2]); nrow(rna_nn2[ranking2>2]) # 1671518; 2711100; /4382618=38.14%; 61.86%
chisq.test(matrix(c(13049, 7889, 1671518, 2711100), byrow = T, nrow = 2))$`p.value` # 0

# stack plot
aa<-data.table(qtl=c("pQTL","pQTL","bg","bg"), replication=c("Regulome","No regulome","Regulome","No regulome"), num=c(1259, 363, 1515102, 1940357)) # consider all pQTLs
aa<-data.table(qtl=c("pQTL","pQTL","bg","bg"), replication=c("Regulome","No regulome","Regulome","No regulome"), num=c(23, 11, 2678, 2425)) # consider only lead pQTLs/variants
aa[, prop := 100*num / sum(num), by = qtl]
aa$qtl<-factor(aa$qtl,levels = c( "bg", "pQTL"))
aa$replication<-factor(aa$replication,levels = c("Regulome","No regulome"))
aa[, fill_group := paste(qtl, replication, sep = "_")]

ggplot(aa, aes(x = qtl, y = prop, fill = fill_group)) + # 展示比例
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#F4EDCA", "#C4961A", "#C3D7A4", "#52854C")) +
  # scale_fill_manual(values = c("mistyrose","indianred","lightsteelblue","steelblue")) +
  labs(x = '', y = 'Proportion of functional variants (%)') +
  scale_x_discrete(labels = c("bg"="Background", "pQTL"="pQTLs"))+
  coord_flip() +
  mytheme + theme(legend.position = "none")

## pGenes replicated in eGenes ####
# 服务器。‘replication’ defining as lead pQTL variants in has nominal PV <0.01 with consistent effect direction in eQTL. Vice versa. ref: PMID: 32773033 
# p<0.01: pGene: 13/33=39.4%; lead pQTL: 45/118=38.1%; 
# p<0.05: pGene: 16/33=48.5%; lead pQTL: 49/118=41.5%
pp<-fread("~/projects/pQTL/formal/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header=T)
names(pp)[1]<-"pro_id"
nrow(pp); length(unique(pp$pro_id)); length(unique(pp$variant_id)); length(unique(pp[is_LeadpQTL==1]$variant_id)) # 1623; 34; 1623; 120

rna_nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); nrow(rna_nn)
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt", header = T)
rna_nn1<-merge(rna_nn, name_id_reviewd[,1:2], by="gene_id")

pp1<-merge(pp, rna_nn1, by=c("pro_id", "variant_id"))
nrow(pp1); length(unique(pp1$pro_id)); length(unique(pp1$variant_id)); length(unique(pp1[is_LeadpQTL==1]$variant_id)) # 1613; 33; 1613; 118

pp2<-pp1[is_LeadpQTL==1] # 取出 lead pQTL
nrow(pp2); length(unique(pp2$pro_id)); length(unique(pp2$variant_id)); length(unique(pp2[is_LeadpQTL==1]$variant_id)) # 118; 33; 118; 118
nrow(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]); length(unique(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]$variant_id)) # 45; 13; 45 
# pGene: 13/33=39.4%; lead pQTL: 45/118=38.1%

## eGenes replicated in pGenes ####
# p<0.01: eGene: 11/42=26.2%; lead pQTL: 47/139=33.8%
# p<0.05: eGene: 16/42=38.1%; lead pQTL: 56/139=40.3%
pp<-fread("~/projects/pQTL/formal/TableS3_skin_eQTL_allSig_pairs.csv", skip=3, header=T)
nrow(pp); length(unique(pp$gene_id)); length(unique(pp$variant_id)); length(unique(pp[is_LeadeQTL==1]$variant_id)) # 24788; 282; 20999; 890
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt", header = T)
pp1<-merge(pp, name_id_reviewd, by="gene_id") # 把 gene_id添加到 pro_id
nrow(pp1); length(unique(pp1$gene_id)); length(unique(pp1$variant_id)); length(unique(pp1[is_LeadeQTL==1]$variant_id)) # 12154; 152; 11801; 507.有一些non-coding eGene没有对应的 pro_id

pro_nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); nrow(pro_nn) # 15669687
names(pro_nn)[1]<-"pro_id"

pp2<-merge(pp1, pro_nn, by=c("pro_id","variant_id"))
nrow(pp2); length(unique(pp2$pro_id)); length(unique(pp2$variant_id)); length(unique(pp2[is_LeadeQTL==1]$variant_id)) # 3549; 42; 3548; 139

pp3<-pp2[is_LeadeQTL==1] # 取出 lead pQTL
nrow(pp3); length(unique(pp3$pro_id)); length(unique(pp3$variant_id)); length(unique(pp3[is_LeadeQTL==1]$variant_id)) # 140; 42; 139; 139
nrow(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]); length(unique(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]$variant_id)) # 47; 11; 47
# eGene: 11/42=26.2%; lead pQTL: 47/139=33.8%

nrow(pp3[pval_nominal.y<=0.05 & slope.x*slope.y>0]); length(unique(pp3[pval_nominal.y<=0.05 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp3[pval_nominal.y<=0.05 & slope.x*slope.y>0]$variant_id)) # 57; 16; 56
# eGene: 16/42=38.1%; lead pQTL: 56/139=40.3%

## stacked plot 
aa<-data.table(qtl=c("pQTL","pQTL","eQTL","eQTL"),replication=c("no eQTL","eQTL","no pQTL","pQTL"),num=c(20,13,31,11))
aa[, prop := 100*num / sum(num), by = qtl]
aa$qtl<-factor(aa$qtl,levels = c("pQTL", "eQTL"))
aa$replication<-factor(aa$replication,levels = c("no eQTL","eQTL","no pQTL","pQTL"))

ggplot(aa, aes(x = qtl, y = prop, fill = replication)) + # 展示比例
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#F4EDCA", "#C4961A", "#C3D7A4", "#52854C")) +
  # scale_fill_manual(values = c("mistyrose","indianred","lightsteelblue","steelblue")) +
  labs(x = '', y = 'Proportion of replicated genes (%)') +
  scale_x_discrete(labels = c("eQTL"="eGenes", "pQTL"="pGenes"))+
  mytheme + theme(legend.position = "none")

ggplot(aa, aes(x=qtl, y=num, fill=replication)) + # 展示数目
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c("#F4EDCA","#C4961A", "#C3D7A4", "#52854C"))+
  labs(x = '', y = 'Number of genes')+
  scale_x_discrete(labels = c("eQTL"="eGene", "pQTL"="pGene"))+
  mytheme + theme(legend.position="none")

## colocolization of eQTL and pQTL ####
conda activate r-coloc
library(data.table); library(coloc)
if (1) {
  peqtl<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt.gz",header = T)
  nrow(peqtl); length(unique(peqtl$pro_id)); length(unique(peqtl$variant_id)) # 15,443,226 pairs; 5052 proteins; 3,444,802 variants
  peqtl[,id:=paste(pro_id, variant_id, sep="_")]; peqtl[,id2:=paste(gene_id, variant_id, sep="_")]
  
  p1 = 1e-4; p2 = 1e-4; p12 = 1e-5
  
  sel_gene<-unique(peqtl$pro_id); length(sel_gene) # 5052；以 proID 作为循环变量，因为有的geneID会对应多个proID，导致 snps 会出现重复（同一基因的不同 proID 具有的相同SNP）；虽然这样也会有一个 pro_id对应多个gene_id,但会比反过来的情况少
  pph4 <- NULL; pph4_snp <- list() # 存每个蛋白的整体共定位 summary（PPH0~PPH4）; 存每个蛋白的 SNP-level 共定位结果
  for (i in 1:length(sel_gene)) { # 5min
    tmp<-peqtl[pro_id==sel_gene[i],]
    tmp<-tmp[order(variant_id, pval_nominal_pqtl, pval_nominal_eqtl)]; tmp<-tmp[!duplicated(variant_id)] # 一个 pro_id对应多个gene_id,  merge 时variant会重复，去重防止报错
    my.coloc<-coloc.abf(dataset1=list(beta=tmp$slope_pqtl, varbeta=(tmp$slope_se_pqtl)^2, snp=tmp$variant_id, type='quant', N=186, MAF=tmp$maf_pqtl), 
                        dataset2=list(beta=tmp$slope_eqtl, varbeta=(tmp$slope_se_eqtl)^2, snp=tmp$variant_id, type='quant', N=150, MAF=tmp$maf_eqtl),
                        p1 = p1, p2 = p2, p12 = p12) # 先pQTL再eQTL
    # my.coloc是list，包含三个结果：summary 是每个基因的 PPH0~4；results是个 table，每行是一个测试 SNP 的情况，SNP.PP.H4该SNP是共享因果变异的概率；priors
    
    summary_row <- data.table(pro_id = sel_gene[i], t(my.coloc$summary)); pph4 <- rbind(pph4, summary_row)  # 存每个蛋白的整体共定位 summary
    if (!is.null(my.coloc$results)) { snp_res <- as.data.table(my.coloc$results); snp_res[, pro_id := sel_gene[i]]; pph4_snp[[i]] <- snp_res } # 存每个蛋白的 SNP-level 共定位结果
  }
  
  i; length(sel_gene) # 5052=5052，确保运行完毕没有中断
  pph4_snp_all <- rbindlist(pph4_snp, use.names = TRUE, fill = TRUE)
  length(unique(peqtl$pro_id)); nrow(pph4) # 5052=5052
  length(unique(peqtl$variant_id)); length(unique(pph4_snp_all$snp)) # 3444802=3444802
  
  pph4[,c("P1","P2","P3","P4"):=.(min(1,nsnps*p1), min(1,nsnps*p2), min(1,nsnps*(nsnps-1)*p1*p2), min(1,nsnps*p12)),by=.(pro_id)] # p1,p2,p12是每个SNP的先验概率，乘以此区域内的 SNP 数目，表示整个区域的先验概率。
  pph4[,c("relative_P4","relative_PP4"):=.(P4/(P4+P3),PP.H4.abf/(PP.H4.abf+PP.H3.abf)),by=.(pro_id)]
  nrow(pph4); head(pph4)
  fwrite(pph4,"~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt",col.names = T,sep="\t",quote = F)
  fwrite(pph4_snp_all,"~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_snp.txt",col.names = T,sep="\t",quote = F)
  # gzip skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_snp.txt 
  
  pph4<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt", header = T)
  nrow(pph4); nrow(pph4[PP.H0.abf>0.7]); nrow(pph4[PP.H1.abf>0.7]); nrow(pph4[PP.H2.abf>0.7]); nrow(pph4[PP.H3.abf>0.7]); nrow(pph4[PP.H4.abf>0.7])
  # 5052, 633, 7, 7, 0, 12
}

## plot PPH heatmap
library(ComplexHeatmap); library(circlize)
pph4<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt", header = T)
nrow(pph4); nrow(pph4[PP.H0.abf>0.7]); nrow(pph4[PP.H1.abf>0.7]); nrow(pph4[PP.H2.abf>0.7]); nrow(pph4[PP.H3.abf>0.7]); nrow(pph4[PP.H4.abf>0.7]) # 5052, 633, 7, 7, 0, 12

name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
pph4<-merge(pph4[,1:7], unique(name_id_reviewd[,c(1,4)]), by="pro_id"); nrow(pph4) # 5090.pro_id有重复的，对应多个 gene_name
names(pph4)<-gsub(".abf", "",names(pph4))
nrow(pph4[PP.H0>=0.7,]); nrow(pph4[PP.H1>=0.7,]); nrow(pph4[PP.H2>=0.7,]); nrow(pph4[PP.H3>=0.7,]); nrow(pph4[PP.H4>=0.7,]) 
# 650,7,7,0,12

hh2<-rbind(pph4[order(-PP.H1)][PP.H1>=0.7],pph4[order(-PP.H2)][PP.H2>=0.7],pph4[order(-PP.H3)][PP.H3>=0.7],pph4[order(-PP.H4)][PP.H4>=0.7])
fwrite(hh2[,c(8,1:7)], "~/Desktop/省皮/project/pQTL/manuscript/TableS4_skin_pQTL_eQTL.coloc.csv",col.names = T,sep=",",quote = F)
tmp<-as.data.frame(hh2); rownames(tmp)<-tmp$gene_name; mat <- as.matrix(tmp[, 3:7])

Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, col=colorRampPalette(c("white","firebrick"))(50),
        row_names_side = "left", column_names_side = "top", border = TRUE, column_names_rot = 45, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 12), legend_height = unit(4, "cm"), legend_width  = unit(0.6, "cm"), title = "PP", legend_direction = "vertical"),
        row_names_gp = gpar(fontsize = 14, fontface = "bold.italic"), column_names_gp = gpar(fontsize = 14, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) { grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey80", fill = NA, lwd = 0.5))}
        )

## Case plot: locus zoom/protein abundance~SNP boxplot. 这里只有 F12 的例子。
# F12 Protein (P00748)~5:177409531:A:G (rs1801020)
qqnorm<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt",header = T)
qq3<-melt(qqnorm[,-1:-3],id.vars = "pro_id"); names(qq3)<-c("pro_id", "sample", "value")
sel_pro<-"P00748"; sel_snp<- "5:177409531:A:G"
# in linux: 将这个 SNP 写入 tmp.bed; bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o tmp.vcf
myvcf<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/tmp.vcf",skip = 490,header = T)
my_sel2<-grep("ID|NC",names(myvcf)); myvcf[,..my_sel2]->myvcf
melt.data.table(myvcf,id.vars = "ID")->myvcf2
names(myvcf2)<-c("id","sample","format"); myvcf2[,gt:=gsub(":.+","",format)]
myvcf2[,ds:=as.numeric(sapply(strsplit(format,split = ":"), "[", 2))] # add dosage info

qq4<-merge(qq3[pro_id==sel_pro], myvcf2[id==sel_snp, .(id,sample, gt,ds)], by="sample")
qq4[,c("ref","alt"):=.(strsplit(id,":")[[1]][3], strsplit(id,":")[[1]][4])]
qq4[,gt2:=ifelse(gt=="0|0",paste(ref,ref,sep="/"),ifelse(gt=="1|1",paste(alt,alt,sep="/"),paste(ref,alt,sep="/")))]
qq4$gt2<-factor(qq4$gt2,levels = c(unique(paste(qq4$ref,qq4$ref,sep="/")),unique(paste(qq4$ref,qq4$alt,sep="/")),unique(paste(qq4$alt,qq4$alt,sep="/"))))
table(qq4$gt) # 0|0 0|1 1|0 1|1 = 85 47  33   21

max_inten<- 1.1*max(qq4$value,na.rm = T); min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), size=2.5, alpha = 0.8,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = gt2, group=factor(gt2)), outlier.shape = NA, fill = '#00000000', width = 0.7, lwd=1.1, position = position_dodge(width = 0.85))+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  labs(x = 'Genotype', y = 'Normalized protein abundance',color = NULL)+
  mytheme + theme(legend.text = element_text(size = 16),legend.position = "none",legend.title=element_blank())

# F12 mRNA (ENSG00000131187)~5:177409531:A:G (rs1801020)
qqnorm<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt",header = T)
qq3<-melt(qqnorm[,-1:-3],id.vars = "gene_id"); names(qq3)<-c("pro_id", "sample", "value")
sel_pro<-"ENSG00000131187"; sel_snp<- "5:177409531:A:G"
# in linux: 将这个 SNP 写入 tmp.bed; bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o tmp.vcf
myvcf<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/tmp.vcf",skip = 490,header = T)
my_sel2<-grep("ID|NC",names(myvcf)); myvcf[,..my_sel2]->myvcf
melt.data.table(myvcf,id.vars = "ID")->myvcf2
names(myvcf2)<-c("id","sample","format"); myvcf2[,gt:=gsub(":.+","",format)]
myvcf2[,ds:=as.numeric(sapply(strsplit(format,split = ":"), "[", 2))] # add dosage info

qq4<-merge(qq3[pro_id==sel_pro], myvcf2[id==sel_snp, .(id,sample, gt,ds)], by="sample")
qq4[,c("ref","alt"):=.(strsplit(id,":")[[1]][3], strsplit(id,":")[[1]][4])]
qq4[,gt2:=ifelse(gt=="0|0",paste(ref,ref,sep="/"),ifelse(gt=="1|1",paste(alt,alt,sep="/"),paste(ref,alt,sep="/")))]
qq4$gt2<-factor(qq4$gt2,levels = c(unique(paste(qq4$ref,qq4$ref,sep="/")),unique(paste(qq4$ref,qq4$alt,sep="/")),unique(paste(qq4$alt,qq4$alt,sep="/"))))
table(qq4$gt) # 0|0 0|1 1|0 1|1 = 71 38  24 17

max_inten<- 1.1*max(qq4$value,na.rm = T); min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), size=2.5, alpha = 0.8,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = gt2, group=factor(gt2)), outlier.shape = NA, fill = '#00000000', width = 0.7, lwd=1.1, position = position_dodge(width = 0.85))+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  labs(x = 'Genotype', y = 'Normalized RNA expression',color = NULL)+
  mytheme + theme(legend.text = element_text(size = 16),legend.position = "none",legend.title=element_blank())