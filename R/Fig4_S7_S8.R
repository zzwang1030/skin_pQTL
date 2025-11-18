######## common data ########
library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 24) + theme(
  axis.text = element_text(color = 'black',size=20), strip.background = element_blank(), 
  plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
name_id_reviewd2<-unique(name_id_reviewd[,c(1,4)]); name_id_reviewd2<-name_id_reviewd2[!duplicated(pro_id)]
length(unique(name_id_reviewd2$pro_id)); length(name_id_reviewd2$pro_id) # 42433 = 42433

## transcriptional start sites info
tss<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/ensembl_uniprot_id_pos_biomart.txt",header = T) 
tss[,strand:=ifelse(strand==1,"+","-")]
tss2<-tss[,.(chr=chromosome,start=transcription_start_site-1,end=transcription_start_site,id=uniprot_id,gene_id=gene_id,strand=strand)]
tss2<-tss2[order(-id,chr,start),]
tss2[,sel_tss:=ifelse(length(unique(chr))>1,start[1],ifelse(strand=="+",min(start),max(start))),by=.(id)]
tss3<-tss2[start==sel_tss,.SD,by=.(id)]; tss3<-tss3[!duplicated(id),]; tss3<-tss3[id!="",.(chr,start,end,id,gene_id,strand)]

# Extract geno samples from VCF. 198 smps: bcftools query -l hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf|sort > tmp_sample
sel_smp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/tmp_sample",header = F)
# Extracting autosomes from VCF. bcftools query -f '%CHROM\t%POS\n' hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz|awk '{n[$1]++} END {for(i in n) print i, n[i]}'|cut -d ' ' -f 1 | grep -v 'chrX'|grep -v 'chrY' | awk '{ if(length($1)<6) {print $1}}' > sel_chr
sel_chr<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/sel_chr",header = F)

#### skin eQTL mapping ######################### 
## QC of skin RNA matrix ####
rr<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # a total table contain raw cnt, log2Normalized_cnt and FPKM value; The complete file is too large, so only a demo was provided here.
# rr<-fread("~/Desktop/省皮/project/pQTL/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # complete data
sel_col<-grep("gene_ID|.+SRNA.+fpkm",names(rr),value=T); rr[,..sel_col]->rr_rpkm; names(rr_rpkm)<-gsub("SRNA_fpkm", "",names(rr_rpkm))
sel_col<-grep("gene_ID|.+SRNA$",names(rr),value=T); rr[,..sel_col]->rr_cnt; names(rr_cnt)<-gsub("SRNA","",names(rr_cnt))
names(rr_rpkm)[1]<-"gene_id";names(rr_cnt)[1]<-"gene_id"
rm(rr) # large memory load
identical(names(rr_rpkm),names(rr_cnt)) # TRUE
## select smps of pheno (RNA) overlapped with geno samoles (150 smp)
my_sel<-names(rr_rpkm)%in%c("gene_id",setdiff(sel_smp$V1,c("NC285","NC060"))) # Two more samples were excluded due to their low alignment quality.
rr_rpkm[,..my_sel]->rr_rpkm1; rr_cnt[,..my_sel]->rr_cnt1; ncol(rr_rpkm1)-1 # 150 smps
## select genes of autosomes and non-predicting
tss_all<-fread("~/Desktop/reference/ensembl_allGene_pos_biomart.txt",header = T) # download from ensembl biomart: all genes position info
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

##  RPKM >0.1 in at least 10 smps; summed cnt > 5 for each gene. Adapted from ref: GTEx PMID: 32913072: >10smps, >6cnt total
smp_num <- ncol(rr_rpkm1)-1
rr_cnt1[,ave_cnt:=rowMeans(as.matrix(rr_cnt1[,2:151]))]; rr_cnt1[,sum_cnt:=rowSums(as.matrix(rr_cnt1[,2:151]))]
melt(rr_rpkm1,id.vars = "gene_id")->rr_rpkm2
names(rr_rpkm2)<-c("gene_id","sample","value")
rr_rpkm2[,rpkm01_gene:=sum(value>0.1),by=.(gene_id)]
length(unique(rr_rpkm2[rpkm01_gene > (smp_num/2 -1), gene_id])) # 30429/28022/22600 for 10/20/half_sample
length(unique(rr_rpkm2[gene_id%in%rr_cnt1[ave_cnt>5,gene_id],gene_id])) # 27605/23772/20968 for ave_ant 2/5/10
rr_rpkm3<-rr_rpkm2[rpkm01_gene >= 10]; length(unique(rr_rpkm3$gene_id)) # 30775
rr_rpkm3<-rr_rpkm3[gene_id%in%rr_cnt1[sum_cnt >= 6,gene_id],]; length(unique(rr_rpkm3$gene_id)) # 30775

## Normalization of skin RNA matrix  ####
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
rr_rpkm6<-copy(rr_rpkm5); rr_rpkm6<-rr_rpkm6[chr%in%sel_chr$V1]; dim(rr_rpkm6)
names(rr_rpkm6)[1]<-"#Chr";
fwrite(rr_rpkm6, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt", sep = "\t",col.names = T)

## skin RNA covariats processing ########
# PCA method by QTLtools and select K PCs by findPC package
library(findPC)
pca<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca",header = T)
names(pca)[1]<-"PC_id"
pca[,PC_id:=paste0("geno_",gsub("hg38.+svd_", "",PC_id))]
pca1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.forPCA.gz.pca",header = T)
names(pca1)[1]<-"PC_id"
pca1[,PC_id:=paste0("pheno_",gsub("skin.+svd_", "",PC_id))]

# select K PCs for use
pve<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca_stats",header = F,nrows = 3)
pve<-transpose(pve,keep.names = "PC",make.names = "V1");pve[,PC:=1:.N];pve[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve$PC,pve$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve$prop_var,number = c(30,60,100),method = 'all',aggregate = 'median',figure = T)
sel_genoPCs<-5 

pve1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.forPCA.gz.pca_stats",header = F,nrows = 3)
pve1<-transpose(pve1,keep.names = "PC",make.names = "V1");pve1[,PC:=1:.N];pve1[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve1$PC,pve1$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve1$prop_var,number = c(30,60,100,150),method = 'all',aggregate = 'median',figure = T)
sel_phenoPCs<-10 

transpose(pca,keep.names = "sample", make.names = "PC_id")->tmp
transpose(pca1,keep.names = "sample", make.names = "PC_id")->tmp1
merge(tmp1[,1:(sel_phenoPCs+1)],tmp[,1:(sel_genoPCs+1)],by="sample")->tmp2 
sample_info<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Sample_info.csv", header = T)
merge(tmp2,sample_info[,.(No,gender,age,Sun_exposure)],by.x="sample",by.y="No",all.x = T)->tmp3 
# convert catagory variables to dummy variables
tmp3[,gender:=ifelse(gender=="male",1,0)]
tmp3[,Sun_exposure:=ifelse(Sun_exposure=="Sun_exposed",1,0)]

transpose(tmp3,keep.names = "PC_id", make.names = "sample")->tmp4
tmp4[is.na(tmp4)]<-"NA" #  The missing values should be encoded as NA 
fwrite(tmp4,"~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.12pheno_5geno.PCs",col.names = T,quote = F,sep = "\t")

# eQTL calling: ~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/pQTL_formal_code.docx

#### eQTL charactering ####
## 计算显著的 SNP-RNA pair 数目 ####
library(data.table)
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

# nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header=T); length(unique(nn$gene_id)); length(unique(nn$variant_id)) # 30470 RNAs; 4390105 snps
nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header=T) # The complete dataset is very large; a smaller demo dataset is provided here.
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
length(unique(pp$gene_id)); length(unique(pp$variant_id)) # 30470 RNAs; 29162 snps
names(pp)[-1]<-paste0(names(pp)[-1],"_perm")
nrow(pp[qval_perm<=0.05]) # 284 Sig eGenes 

hh<-merge(nn,pp[,.(gene_id,num_var_perm,pval_beta_perm,qval_perm,pval_nominal_threshold_perm)],by="gene_id")
hh1<-hh[pval_nominal <= pval_nominal_threshold_perm,]
hh2<-hh1[qval_perm<=0.05,]
hh2<-hh2[order(gene_id,pval_nominal,-abs(slope),slope_se)]
hh2[,is_LeadpQTL:=1:.N,by=.(gene_id)]; hh2[,is_LeadpQTL:=ifelse(is_LeadpQTL==1,1,0),by=.(gene_id)];
names(hh2)<-c("gene_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","num_var","pval_beta","qval","pval_nominal_threshold","is_LeadpQTL")
hh2<-hh2[,c("gene_id","variant_id","is_LeadpQTL","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","pval_beta","qval","pval_nominal_threshold")]
hh3 <- tidyr::extract(hh2, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F)
fwrite(hh3,"~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS3_skin_eQTL_allSig_pairs.csv",sep=",",col.names = T,quote = F) 

## Effect size of self eQTL vs gtex eQTL ####
# The all-pair data is too large, so here is a demo showing the sparsified data.
# nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T)
nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header = T)
nrow(nn); nrow(nn[pval_nominal<0.01,]); nrow(nn[pval_nominal<0.05,]) # 93916515 rows; p<0.01 (1262089); p<0.05(5307716)
names(nn)[1]<-"pro_id"; nn<-nn[order(pro_id,pval_nominal,-abs(slope),slope_se)] 

gtex_nn<-fread("~/reference/GTEx_Analysis_v8_eQTL/Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs.txt.gz")  # download from GTEx online; 2414653 sig pairs
gtex_nn[,variant_id:=gsub("_b38","",gsub("chr","",variant_id))]; gtex_nn[,variant_id:=gsub("_",":",variant_id)]
gtex_nn[,gene_id:=gsub("[.]\\d+","",gene_id)]; names(gtex_nn)[2]<-"pro_id"
gtex_nn<-gtex_nn[order(gene_id,pval_nominal)]

total<-merge(nn, gtex_nn, by=c("pro_id","variant_id")); nrow(total) # 1478275
fwrite(total,"~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.allpairs.OverlapGTEx.txt",sep="\t",quote=F) # gzip skin_RNA_all_qqnorm.txt.gz.allpairs.OverlapGTEx.txt 
total<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.allpairs.OverlapGTEx.txt.gz",header = T)
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

# plot
result<-data.table(pt=c(1e-2, 1e-4, 1e-6, 1e-8), n_pair=c(157655, 35725, 11741,3219),
                   n_same_direction=c(151162, 35292, 11718, 3219),pearson_r=c(0.8567206, 0.9016582, 0.9111776, 0.9428516),
                   p_value=c(0, 0, 0, 0)) # from linux R 
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)],id.vars = "pt"); names(mydt2)<-c("pt","type","value")
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)] 
mydt2$type<-factor(mydt2$type,levels = c("percent","pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = c(1e-02, 1e-04, 1e-06, 1e-08))
ggplot(mydt2, aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(85, 100)) + labs(x = 'P-value thresholds')+
  scale_fill_manual(values = c("percent" = "indianred", "pearson_r" = "steelblue"),
                    labels = c("percent" = "Percent of same direction (%)", "pearson_r" = "Pearson correlation"))+
  mytheme  + theme(legend.position="top",legend.direction = "horizontal", legend.title=element_blank())

pt <- 1e-04
ggplot(total[pval_nominal.x<pt & pval_nominal.y<pt], mapping = aes(x = slope.x, y = slope.y)) +
  geom_pointdensity(size = 0.8) + scale_color_viridis() +
  labs(x = 'Effect size in our study', y = 'Effect size in GTEx', color = NULL)+
  # stat_cor(method = "pearson", digits=3,size=6)+
  coord_cartesian(xlim = c(-2, 2),ylim = c(-2, 2)) +
  geom_hline(yintercept=0, color="black",linetype="dashed", size=0.4)+geom_vline(xintercept=0, color="black",linetype="dashed", size=0.4)+
  geom_smooth(method='lm', se=T, formula= y~x, col="red", linetype="dashed") +
  mytheme + theme(legend.position = c(0.25, 0.9), legend.direction = "horizontal", legend.key.width = unit(1.2, "cm"), legend.text = element_text(size = 12)) 

#### pQTL and eQTL replication ########
## Count of eQTL&pQTL overlap ####
library(data.table); library(VennDiagram)
# The all-pair data is too large, so here is a demo showing the sparsified data.
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

# pro_nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); names(pro_nn)[1]<-"pro_id"
pro_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header = T); names(pro_nn)[1]<-"pro_id"
nrow(pro_nn); length(unique(pro_nn$pro_id)); length(unique(pro_nn$variant_id)) # 15,669,687; 5152; 3,461,609
names(pro_nn)[3:9]<-paste0(names(pro_nn)[3:9], "_pqtl")

# rna_nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); nrow(rna_nn) # 93,916,515 pairs
rna_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header = T); names(pro_nn)[1]<-"pro_id"
rna_nn1<-merge(rna_nn, name_id_reviewd, by="gene_id", all.x=T)
names(rna_nn1)[3:9]<-paste0(names(rna_nn1)[3:9], "_eqtl")
nrow(rna_nn1); length(unique(rna_nn1$pro_id)); length(unique(rna_nn1$variant_id)) # 47,485,953; 15395; 4,197,272

nn<-merge(pro_nn, rna_nn1, by=c("pro_id", "variant_id"))
nrow(nn); length(unique(nn$pro_id)); length(unique(nn$variant_id)) # 15,443,226; 5052; 3,444,802
fwrite(nn, "~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt", col.names=T, row.names=F, sep="\t", quote=F)
# gzip skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt 

# Permutation test for the significance of the number of overlaps。
pt<-1e-2
pqtl_pairs <- unique(paste(nn[pval_nominal_pqtl < pt]$pro_id, nn[pval_nominal_pqtl < pt]$variant_id, sep = "_"))
pqtl_n <- length(pqtl_pairs) # 175905
eqtl_pairs <- unique(paste(nn[pval_nominal_eqtl < pt]$pro_id, nn[pval_nominal_eqtl < pt]$variant_id, sep = "_"))
eqtl_n <- length(eqtl_pairs) # 210203 
true_overlap <- length(intersect(pqtl_pairs, eqtl_pairs)) # 7657
background_pairs <- unique(paste(nn$pro_id, nn$variant_id, sep = "_")); length(background_pairs) # 15384030

n_perm <- 1000; set.seed(123)
perm_overlap <- numeric(n_perm) 
for (i in 1:n_perm) {
  random_pqtl <- sample(background_pairs, pqtl_n)
  random_eqtl <- sample(background_pairs, eqtl_n)
  perm_overlap[i] <- length(intersect(random_pqtl, random_eqtl))
}
mean(perm_overlap) # 2404.629
empirical_p <- (sum(perm_overlap >= true_overlap) + 1) / (length(perm_overlap) + 1) # 9.9e-4
enrichment <- true_overlap / mean(perm_overlap) # 7657/2404.629 = 3.184275
fwrite(data.frame(overlap = perm_overlap), "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.permutation.txt", col.names=T, row.names=F, sep="\t", quote=F)
tmp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.permutation.txt")
ggplot(tmp, aes(x=overlap)) + 
  geom_histogram(bins=50, fill="steelblue",color='grey90', linewidth=0.5) +
  labs(x = "Overlap count under permutations", y = "Frequency")+ mytheme

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

# plot
result<-data.table(pt=c(0.05, 1e-2, 1e-4, 1e-6), n_pair=c(60765, 7660, 826, 132),
                   n_same_direction=c(43350, 6767, 826, 132),pearson_r=c(0.4405654, 0.7236292, 0.8143522, 0.9464664),
                   p_value=c(0.000000e+00, 0.000000e+00, 6.702596e-197, 1.070058e-65)) # from linux R 
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)], id.vars = "pt"); names(mydt2)<-c("pt","type","value")
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)] 
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
oo<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", header = T); nrow(oo) # 1478275
pt <- 1e-02
ggplot(oo[pval_nominal_pqtl<pt & pval_nominal_eqtl<pt], mapping = aes(x = slope_pqtl, y = slope_eqtl)) +
  geom_pointdensity(size = 1.5) + scale_color_viridis() +
  labs(x = 'Effect size of pQTL', y = 'Effect size of eQTL', color = NULL)+
  coord_cartesian(xlim = c(-1.2, 1.2),ylim = c(-1.2, 1.2)) +
  geom_hline(yintercept=0, color="black",linetype="dashed", linewidth=0.4) + geom_vline(xintercept=0, color="black",linetype="dashed", linewidth=0.4)+
  geom_smooth(method='lm', se=T, formula= y~x, col="red", linetype="dashed") +
  mytheme + theme(legend.position = c(0.25, 0.9), legend.direction = "horizontal", legend.key.width = unit(1.2, "cm"), legend.text = element_text(size = 14))

## eQTL and pQTL opposite direction ####
oo<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", header = T); nrow(oo) # 1478275
pp1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pp2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
oo<-merge(oo, pp1[,.(gene_id, pval_nominal_threshold)], by="gene_id"); names(oo)[ncol(oo)]<-"pval_nominal_threshold_eqtl"
oo<-merge(oo, pp2[,.(gene_id, pval_nominal_threshold)], by.x = "pro_id", by.y ="gene_id"); names(oo)[ncol(oo)]<-"pval_nominal_threshold_pqtl"
oo[,p_sum:= -log10(pval_nominal_pqtl) -log10(pval_nominal_eqtl)]

oo1<-oo[slope_pqtl*slope_eqtl<0,]; oo1<-oo1[order(p_sum)]; nrow(oo1) # 893
nrow(oo1[pval_nominal_pqtl <= pval_nominal_threshold_pqtl & pval_nominal_eqtl <= pval_nominal_threshold_eqtl]) 

library(ggplot2); library(scales)  
ggplot(oo1, aes(x = -log10(pval_nominal_pqtl), y = -log10(pval_nominal_eqtl), color = slope_pqtl)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_hline(yintercept = max(-log10(oo1$pval_nominal_threshold_eqtl)), linetype = "dashed", color = "black") +
  geom_vline(xintercept = max(-log10(oo1$pval_nominal_threshold_pqtl)), linetype = "dashed", color = "black") + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "pQTL effect size", guide = guide_colorbar(title.position = "top", barwidth = 8, barheight = 1)) +
  labs(x = "-log10(P-value) of pQTL", y = "-log10(P-value) of eQTL", color = "pQTL effect size") +
  theme_minimal(base_size = 18) +  theme(legend.position = "top", legend.direction = "horizontal", legend.title.align = 0.5,  legend.title = element_text(size = 14), legend.text  = element_text(size = 14))

## pGenes replicated in eGenes ####
# ‘replication’ defining as lead pQTL variants in has nominal PV <0.01 with consistent effect direction in eQTL. Vice versa. ref: PMID: 32773033 
# p<0.01: pGene: 13/33=39.4%; lead pQTL: 45/118=38.1%; 
# p<0.05: pGene: 16/33=48.5%; lead pQTL: 49/118=41.5%
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header=T)
names(pp)[1]<-"pro_id"
nrow(pp); length(unique(pp$pro_id)); length(unique(pp$variant_id)); length(unique(pp[is_LeadpQTL==1]$variant_id)) # 1623; 34; 1623; 120

# The all-pair data is too large, so here is a demo showing the sparsified data.
# rna_nn<-fread("~/projects/pQTL/formal/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); nrow(rna_nn)
rna_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header = T); nrow(rna_nn)
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
rna_nn1<-merge(rna_nn, name_id_reviewd[,1:2], by="gene_id")

pp1<-merge(pp, rna_nn1, by=c("pro_id", "variant_id"))
nrow(pp1); length(unique(pp1$pro_id)); length(unique(pp1$variant_id)); length(unique(pp1[is_LeadpQTL==1]$variant_id)) # 1613; 33; 1613; 118
pp2<-pp1[is_LeadpQTL==1] # lead pQTL
nrow(pp2); length(unique(pp2$pro_id)); length(unique(pp2$variant_id)); length(unique(pp2[is_LeadpQTL==1]$variant_id)) # 118; 33; 118; 118
nrow(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]); length(unique(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]$variant_id)) # 45; 13; 45 
# pGene: 13/33=39.4%; lead pQTL: 45/118=38.1%

## eGenes replicated in pGenes ####
# p<0.01: eGene: 11/42=26.2%; lead pQTL: 47/139=33.8%
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS3_skin_eQTL_allSig_pairs.csv", skip=3, header=T)
nrow(pp); length(unique(pp$gene_id)); length(unique(pp$variant_id)); length(unique(pp[is_LeadeQTL==1]$variant_id)) # 24788; 282; 20999; 890
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
pp1<-merge(pp, name_id_reviewd, by="gene_id") # 把 gene_id添加到 pro_id
nrow(pp1); length(unique(pp1$gene_id)); length(unique(pp1$variant_id)); length(unique(pp1[is_LeadeQTL==1]$variant_id)) # 12154; 152; 11801; 507

# The all-pair data is too large, so here is a demo showing the sparsified data.
# pro_nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); nrow(pro_nn) # 15669687
pro_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header = T); nrow(pro_nn) # 15669687
names(pro_nn)[1]<-"pro_id"

pp2<-merge(pp1, pro_nn, by=c("pro_id","variant_id"))
nrow(pp2); length(unique(pp2$pro_id)); length(unique(pp2$variant_id)); length(unique(pp2[is_LeadeQTL==1]$variant_id)) # 3549; 42; 3548; 139

pp3<-pp2[is_LeadeQTL==1] # lead eQTL
nrow(pp3); length(unique(pp3$pro_id)); length(unique(pp3$variant_id)); length(unique(pp3[is_LeadeQTL==1]$variant_id)) # 140; 42; 139; 139
nrow(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]); length(unique(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]$variant_id)) # 47; 11; 47
# eGene: 11/42=26.2%; lead pQTL: 47/139=33.8%

## stacked plot 
aa<-data.table(qtl=c("pQTL","pQTL","eQTL","eQTL"),replication=c("no eQTL","eQTL","no pQTL","pQTL"),num=c(20,13,31,11))
aa[, prop := 100*num / sum(num), by = qtl]
aa$qtl<-factor(aa$qtl,levels = c("pQTL", "eQTL"))
aa$replication<-factor(aa$replication,levels = c("no eQTL","eQTL","no pQTL","pQTL"))

ggplot(aa, aes(x = qtl, y = prop, fill = replication)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#F4EDCA", "#C4961A", "#C3D7A4", "#52854C")) +
  labs(x = '', y = 'Proportion of replicated genes (%)') +
  scale_x_discrete(labels = c("eQTL"="eGenes", "pQTL"="pGenes"))+
  mytheme + theme(legend.position = "none")

## colocolization of eQTL and pQTL ####
library(data.table); library(coloc)
peqtl<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz",header = T) # The all-pair data is too large, so here is a demo showing the sparsified data.
nrow(peqtl); length(unique(peqtl$pro_id)); length(unique(peqtl$variant_id)) # 15,443,226 pairs; 5052 proteins; 3,444,802 variants
peqtl[,id:=paste(pro_id, variant_id, sep="_")]; peqtl[,id2:=paste(gene_id, variant_id, sep="_")]

p1 = 1e-4; p2 = 1e-4; p12 = 1e-5
sel_gene<-unique(peqtl$pro_id); length(sel_gene) 
pph4 <- NULL; pph4_snp <- list() 
for (i in 1:length(sel_gene)) { 
  tmp<-peqtl[pro_id==sel_gene[i],]
  tmp<-tmp[order(variant_id, pval_nominal_pqtl, pval_nominal_eqtl)]; tmp<-tmp[!duplicated(variant_id)] 
  my.coloc<-coloc.abf(dataset1=list(beta=tmp$slope_pqtl, varbeta=(tmp$slope_se_pqtl)^2, snp=tmp$variant_id, type='quant', N=186, MAF=tmp$maf_pqtl), 
                      dataset2=list(beta=tmp$slope_eqtl, varbeta=(tmp$slope_se_eqtl)^2, snp=tmp$variant_id, type='quant', N=150, MAF=tmp$maf_eqtl),
                      p1 = p1, p2 = p2, p12 = p12) 
  summary_row <- data.table(pro_id = sel_gene[i], t(my.coloc$summary)); pph4 <- rbind(pph4, summary_row)  
  if (!is.null(my.coloc$results)) { snp_res <- as.data.table(my.coloc$results); snp_res[, pro_id := sel_gene[i]]; pph4_snp[[i]] <- snp_res } 
}

pph4_snp_all <- rbindlist(pph4_snp, use.names = TRUE, fill = TRUE)
pph4[,c("P1","P2","P3","P4"):=.(min(1,nsnps*p1), min(1,nsnps*p2), min(1,nsnps*(nsnps-1)*p1*p2), min(1,nsnps*p12)),by=.(pro_id)] 
pph4[,c("relative_P4","relative_PP4"):=.(P4/(P4+P3),PP.H4.abf/(PP.H4.abf+PP.H3.abf)),by=.(pro_id)]
nrow(pph4); nrow(pph4[PP.H0.abf>0.7]); nrow(pph4[PP.H1.abf>0.7]); nrow(pph4[PP.H2.abf>0.7]); nrow(pph4[PP.H3.abf>0.7]); nrow(pph4[PP.H4.abf>0.7]) # 5052, 633, 7, 7, 0, 12
fwrite(pph4,"~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt",col.names = T,sep="\t",quote = F)

## plot PPH heatmap
library(ComplexHeatmap); library(circlize)
pph4<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt", header = T)
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
pph4<-merge(pph4[,1:7], unique(name_id_reviewd[,c(1,4)]), by="pro_id"); nrow(pph4) 
names(pph4)<-gsub(".abf", "",names(pph4))
nrow(pph4[PP.H0>=0.7,]); nrow(pph4[PP.H1>=0.7,]); nrow(pph4[PP.H2>=0.7,]); nrow(pph4[PP.H3>=0.7,]); nrow(pph4[PP.H4>=0.7,]) # 650,7,7,0,12
hh2<-rbind(pph4[order(-PP.H1)][PP.H1>=0.7],pph4[order(-PP.H2)][PP.H2>=0.7],pph4[order(-PP.H3)][PP.H3>=0.7],pph4[order(-PP.H4)][PP.H4>=0.7])
tmp<-as.data.frame(hh2); rownames(tmp)<-tmp$gene_name; mat <- as.matrix(tmp[, 3:7])

Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE, col=colorRampPalette(c("white","firebrick"))(50),
        row_names_side = "left", column_names_side = "top", border = TRUE, column_names_rot = 45, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"), labels_gp = gpar(fontsize = 12), legend_height = unit(4, "cm"), legend_width  = unit(0.6, "cm"), title = "PP", legend_direction = "vertical"),
        row_names_gp = gpar(fontsize = 14, fontface = "bold.italic"), column_names_gp = gpar(fontsize = 14, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) { grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey80", fill = NA, lwd = 0.5))}
)

## Case plot: locus zoom/protein abundance~SNP boxplot. Only show F12 as an example here。The approach is similar for other genes.
# F12 Protein (P00748)~5:177409531:A:G (rs1801020)
qqnorm<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt",header = T)
qq3<-melt(qqnorm[,-1:-3],id.vars = "pro_id"); names(qq3)<-c("pro_id", "sample", "value")
sel_pro<-"P00748"; sel_snp<- "5:177409531:A:G"
# write this SNP into tmp.bed, then bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o tmp.vcf
myvcf<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/F12_5_177409531_A_G.vcf",skip = 490,header = T)
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
qqnorm<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_RNA_all_qqnorm.txt",header = T)
qq3<-melt(qqnorm[,-1:-3],id.vars = "gene_id"); names(qq3)<-c("pro_id", "sample", "value")
sel_pro<-"ENSG00000131187"; sel_snp<- "5:177409531:A:G"
# write this SNP into tmp.bed, then bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o tmp.vcf
myvcf<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/F12_5_177409531_A_G.vcf",skip = 490,header = T)
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




