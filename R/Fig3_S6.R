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

#### skin pQTL mapping ########
## QC of skin protein matrix ####
sample_id<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/sample_id.csv",header = T)

qq2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro.xls",header  = T)
names(qq2)[6:ncol(qq2)]<-gsub("\\[\\d+\\] ","",names(qq2)[6:ncol(qq2)])
names(qq2)[6:ncol(qq2)]<-sample_id$names1[match(names(qq2)[6:ncol(qq2)],sample_id$names3)] # rename sample names
names(qq2)[1]<-"pro_id"
qq2 <- as.data.table(qq2 %>% tidyr::separate_rows(pro_id, sep = ";")); nrow(qq2) # 5979
qq2<-qq2[order(pro_id,PG.Qvalue),];qq2<-qq2[!duplicated(pro_id)] # remove duplicated proteins based on Qvalue, 5979
qq2[,-2:-5]->qq2;names(qq2)<-gsub("SP","",names(qq2)) # remove suffix of sample names to match genotype sample
## quality control
my_sel <- names(qq2)%in%c("pro_id",sel_smp$V1)
qq3 <- qq2[,..my_sel]; ncol(qq3)-1 # select pheno samples overlapped with geno samoles (186 = 249 overlap 198)
qq3<-melt(qq3, id.vars = "pro_id"); names(qq3)<-c("pro_id","sample","value")
qq3[, NA_num_pro:=sum(is.na(value)), by=.(pro_id)]; qq3[,NA_num_smp:=sum(is.na(value)),by=.(sample)] # NA counts for each protein and sample
qq3<-qq3[NA_num_pro < ceiling((ncol(qq2)-1)*0.5)] # remove proteins with NA in more than half of samples
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
qq5[,chr:=paste0("chr",chr)]; qq5<-qq5[order(chr,start)]
qq6<-copy(qq5); qq6<-qq6[chr%in%sel_chr$V1]; dim(qq6) # 5192 proteins and 186 smp
names(qq6)[1]<-"#Chr"; 
fwrite(qq6, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt", sep = "\t",col.names = T)
#  PCA of protein matrix  was performed in: ~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/pQTL_formal_code.docx

## skin protein covariats processing ########
# PCA method by QTLtools and select K PCs by findPC package 
library(findPC)
pca<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca",header = T)
names(pca)[1]<-"PC_id"
pca[,PC_id:=paste0("geno_",gsub("hg38.+svd_","",PC_id))]

pca1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.forPCA.gz.pca",header = T)
names(pca1)[1]<-"PC_id"
pca1[,PC_id:=paste0("pheno_",gsub("skin.+svd_","",PC_id))]

# select K PCs for use
pve<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.pca_stats",header = F,nrows = 3)
pve<-transpose(pve,keep.names = "PC",make.names = "V1");pve[,PC:=1:.N];pve[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve$PC,pve$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve$prop_var,number = c(30,60,100),method = 'all',aggregate = 'median',figure = T)
sel_genoPCs<-5 

pve1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.forPCA.gz.pca_stats",header = F,nrows = 3)
pve1<-transpose(pve1,keep.names = "PC",make.names = "V1"); pve1[,PC:=1:.N]; pve1[,c("sd","prop_var","cumm_prop"):=.(as.numeric(sd),as.numeric(prop_var),as.numeric(cumm_prop))]
plot(pve1$PC,pve1$prop_var,xlab="PC index",ylab="PVE",pch=20,col=gray(0.5,0.5),las=1)
findPC(sdev = pve1$prop_var,number = c(40,60,100,150),method = 'all',aggregate = 'median',figure = T)
sel_phenoPCs<-10 

transpose(pca,keep.names = "sample", make.names = "PC_id")->tmp
transpose(pca1,keep.names = "sample", make.names = "PC_id")->tmp1
merge(tmp1[,1:(sel_phenoPCs+1)],tmp[,1:(sel_genoPCs+1)],by="sample")->tmp2 
sample_info<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/Sample_info.csv", header = T)
merge(tmp2,sample_info[,.(No,gender,age,Sun_exposure)],by.x="sample",by.y="No",all.x = T)->tmp3 
# convert catagory variables to dummy variables
tmp3[,gender:=ifelse(gender=="male",1,0)]
tmp3[,Sun_exposure:=ifelse(Sun_exposure=="Sun_exposed",1,0)]

tmp3<-cbind(tmp3,Other_disease[,-1], sample_pos[,-1])
tmp3[,Other_disease:=NULL]; tmp3[,sample_pos:=NULL]; str(tmp3)

transpose(tmp3,keep.names = "PC_id", make.names = "sample")->tmp4
tmp4[is.na(tmp4)]<-"NA" #  The missing values should be encoded as NA 
fwrite(tmp4,"~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt.gz.12pheno_5geno.PCs",col.names = T,quote = F,sep = "\t")
#  pQTL mapping performed in: ~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/pQTL_formal_code.docx

#### skin pQTL charactering ####
## Calculate the number of significant SNP-protein pairs ####
library(data.table)
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

# nn<-fread("~/projects/pQTL/formal/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header=T); nrow(nn) # 15,669,687; # The complete dataset is too large to include here, so only a smaller demo dataset is provided.
nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header=T)
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
names(pp)[-1]<-paste0(names(pp)[-1],"_perm"); nrow(pp[qval_perm<=0.05]) # 34 Sig proteins 

hh<-merge(nn,pp[,.(gene_id, num_var_perm, pval_perm_perm, qval_perm, pval_nominal_threshold_perm)],by="gene_id")
hh1<-hh[pval_nominal <= pval_nominal_threshold_perm,] # The p-value of permutation is used as the threshold for the significance of pQTL.
hh2<-hh1[qval_perm<=0.05,]
nrow(hh2); length(unique(hh2$gene_id)); length(unique(hh2$variant_id)) # 1623 Sig pair; 34 Sig protein; 1623 Sig snp. We will use this data for the next analysis.
hh2<-hh2[order(gene_id,pval_nominal,-abs(slope),slope_se)]
hh2[,is_LeadpQTL:=1:.N,by=.(gene_id)]; hh2[,is_LeadpQTL:=ifelse(is_LeadpQTL==1,1,0),by=.(gene_id)] 
names(hh2)<-c("protein_id","variant_id","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","num_var","pval_perm_LeadpQTL","qval_protein","pval_nominal_threshold","is_LeadpQTL")
hh2<-hh2[,c("protein_id","variant_id","is_LeadpQTL","tss_distance","ma_samples","ma_count","maf","pval_nominal","slope","slope_se","pval_perm_LeadpQTL","qval_protein","pval_nominal_threshold")]
hh3 <- tidyr::extract(hh2, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F)
fwrite(hh3,"~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.SigPairs_pt.txt",sep="\t",col.names = T,quote = F) 

## manhattan plot ####
# The all-pair data is too large, so here is a demo showing the sparsified data.
library(ggmanh); library(SeqArray) 
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
nn3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.thin.txt.gz",header = T)
pp2<- tidyr::extract(nn3, col = 'variant_id', into = c('chr', 'pos','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F)
pp2<-as.data.table(pp2)
pp2$chr<-factor(pp2$chr,levels = 1:22); pp2$pos<-as.numeric(pp2$pos)
pt<-max(pp[qval<=0.05]$pval_nominal) # 5.85902e-07
pp2<-merge(pp2, name_id_reviewd, by.x="gene_id", by.y="pro_id")
pp2<-pp2[order(gene_id, pval_nominal)]
pp2[, p_order:=1:.N,by=.(gene_id)]
pp2[, lable_name:=ifelse(p_order==1, gene_name, ""), by=.(gene_id)] # for lable lead_snp for each gene
pp2[, lable_name2:=ifelse(pval_nominal<=pt, lable_name, "")] # for lable sig gene_name

pp2[,color_class:=ifelse(pval_nominal<=pt,"sig","ns")]; highlight_colormap <- c("ns" = "grey", "sig" = "maroon")
manhattan_plot(pp2, pval.colname = "pval_nominal", chr.colname = "chr", pos.colname = "pos",y.label = "-log10(P)", thin=T, # thin 可以稀疏要画的点，加快画图速度
               signif = pt, point.size=1.2, label.font.size=4,
               highlight.colname = "color_class", highlight.col = highlight_colormap, color.by.highlight = TRUE) + mytheme + theme(legend.position="none")

## proportion variance explained (PVE) by lead pQTL ####
SNP_PVE2 <- function(beta, MAF, se_beta, N) { 
  numerator <- 2 * (beta^2) * MAF * (1 - MAF)
  denominator <- numerator + ((se_beta^2) * 2 * N * MAF * (1 - MAF))
  SNP_PVE <- numerator / denominator ; return(SNP_PVE) }

sample_N<-186
hh1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.SigPairs_pt.rsID.csv", header = T) # rsID
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh2[,pve_snp2:=SNP_PVE2(slope, maf, slope_se, sample_N)]

hh2<-hh2[order(protein_id, pval_nominal, -slope)]; hh3<-hh2[!duplicated(protein_id)]
summary(hh3$pve_snp2)  # 0.1275  0.1408  0.1531  0.1723  0.1930  0.3091 
nrow(hh3[pve_snp2>=0.2]); length(unique(hh3[pve_snp2>=0.2]$protein_id)) # 8 pGenes, 8 loci
ggplot(hh3, aes(x=100*pve_snp2)) + 
  geom_histogram(bins=20, fill="steelblue",color='grey90', linewidth=0.5) +
  stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), vjust = -0.5, size = 6)+
  labs(x ="Proportion of variance explained (%)", y="Numbers of lead pQTL") + mytheme

## case plot: ALDH2 P05091
qqnorm<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt",header = T)
melt(qqnorm[,-1:-3],id.vars = "pro_id")->qq3; names(qq3)<-c("pro_id", "sample", "value")

sel_pro<-"P05091"; sel_snp<-hh3[protein_id == sel_pro][which.max(pve_snp2), variant_id]
# Write this SNP to tmp.bed; bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o ALDH2_12_112298314_A_G.vcf
myvcf<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/ALDH2_12_112298314_A_G.vcf",skip = 490,header = T)
my_sel2<-grep("ID|NC",names(myvcf)); myvcf[,..my_sel2]->myvcf
melt.data.table(myvcf,id.vars = "ID")->myvcf2
names(myvcf2)<-c("id","sample","format"); myvcf2[,gt:=gsub(":.+", "",format)]
myvcf2[,ds:=as.numeric(sapply(strsplit(format,split = ":"), "[", 2))] # add dosage info

qq4<-merge(qq3[pro_id==sel_pro], myvcf2[id==sel_snp, .(id,sample, gt,ds)], by="sample")
qq4[,c("ref","alt"):=.(strsplit(id,":")[[1]][3], strsplit(id,":")[[1]][4])]
qq4[,gt2:=ifelse(gt=="0|0",paste(ref,ref,sep="/"),ifelse(gt=="1|1",paste(alt,alt,sep="/"),paste(ref,alt,sep="/")))]
qq4$gt2<-factor(qq4$gt2,levels = c(unique(paste(qq4$ref,qq4$ref,sep="/")),unique(paste(qq4$ref,qq4$alt,sep="/")),unique(paste(qq4$alt,qq4$alt,sep="/"))))

max_inten<- 1.1*max(qq4$value,na.rm = T); min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), size=2.5, alpha = 0.8,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = gt2,group=factor(gt2)), outlier.shape = NA, fill = '#00000000', width = 0.7, lwd=1.1, position = position_dodge(width = 0.85))+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  annotate("text", x=2, y=max_inten, size=8, 
           label= paste0("Beta=", round(hh3[protein_id==sel_pro & variant_id==sel_snp, slope],2), 
                         ", P-value=", hh3[protein_id==sel_pro & variant_id==sel_snp, pval_nominal]
                         , ", PVE=", round(100*hh3[protein_id==sel_pro & variant_id==sel_snp, pve_snp2], 2), "%"))+
  labs(x = 'Genotype', y = 'Normalized protein abundance', title = paste(name_id_reviewd[pro_id==sel_pro, gene_name], hh1[variant_id==sel_snp]$rsID, sep=" ~ "),color = NULL)+
  mytheme + theme(legend.text = element_text(size = 16),legend.position = "none",legend.title=element_blank())

## effect siz vs MAF ####
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[!duplicated(protein_id)] 
cor.test(hh3$maf, abs(hh3$slope), method = "spearman") # P=8.061e-07, rho = -0.733 

ggplot(hh3, aes(x = maf, y = abs(slope))) + geom_point(size=2, color="steelblue") + 
  geom_smooth(method='lm', se=T, span = 1, formula= y~x, col="maroon", linetype="dashed") + 
  labs(x = 'Minor allele frequency', y = 'Effect size (absolute values)') +
  coord_cartesian(ylim = c(0, 2)) + mytheme

## effect siz vs distance to TSS ####
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[!duplicated(protein_id)] 
nrow(hh3[abs(tss_distance)<=1e5]); nrow(hh3) # 28/34=82.35%; 2e5: 31/34=91.18%
ggplot(hh3, aes(x = tss_distance/1e6, y = abs(slope))) + geom_point(size=3, color="steelblue") + 
  coord_cartesian(xlim = c(-0.5, 0.5)) + labs(x = 'Distance to TSS (MB)', y = 'Effect size (absolute values)') +
  mytheme

## feature annotation for pQTLs by VEP ####
# vep run in linux: vep -i skin_protein_all_qqnorm.txt.gz.SigPairs_pt.forVEP.vcf -o skin_protein_all_qqnorm.txt.gz.SigPairs_pt.forVEP.res.txt --gtf /sh2/home/sunyuanqiang/reference/Homo_sapiens.GRCh38.110.sorted.gtf.gz --fasta /sh2/home/sunyuanqiang/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tab --force_overwrite --show_ref_allele --sift b --canonical
ss7 <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/VEP_annotation_all_pQTLs.csv", header=T)

bb<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) 
bb<-bb[order(protein_id, pval_nominal, -slope)];   bb[,id:=paste(protein_id, variant_id, sep="_")]
bb2<-bb[!duplicated(protein_id)] 
ss8<-ss7[id %in% bb2$id]; nrow(ss8) # 34, extract lead pQTL
tmp1<- as.data.table(table(ss8$Consequence)); names(tmp1)<-c("effect","cnt"); tmp2<-tmp1
tmp2[,ratio:=round(cnt*100 / sum(cnt), 1)]; tmp2[, effect2:=gsub(" variant","",gsub("_", "", effect))]
tmp2<-tmp2[order(ratio,decreasing = T)]; tmp2$effect2 <- factor(tmp2$effect2,levels = tmp2$effect2)

# feature annotation for lead SNP of all genes: vep -i skin_protein_all_qqnorm.txt.gz.permu.genes.forVEP.vcf -o skin_protein_all_qqnorm.txt.gz.permu.genes.forVEP.res.txt --gtf /sh2/home/sunyuanqiang/reference/Homo_sapiens.GRCh38.110.sorted.gtf.gz --fasta /sh2/home/sunyuanqiang/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tab --force_overwrite --show_ref_allele --sift b --canonical
ss7<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/VEP_annotation_leadSNP_allGenes.csv", header = T)
tmp_bg1<- as.data.table(table(ss7$Consequence)); names(tmp_bg1)<-c("effect","cnt") 
tmp_bg2<-tmp_bg1[effect!="splice_acceptor_variant" &  effect!="splice_polypyrimidine_tract_variant"]; tmp_bg2[effect=="splice_region_variant"]$cnt<-3
tmp_bg2[,ratio:=round(cnt*100 / sum(cnt), 1)]; tmp_bg2[, effect2:=gsub(" variant","",gsub("_", "", effect))]

# lead pQTLs of pGenes vs lead pQTLs of all genes
tmp3<-rbind(tmp2[,-4], tmp_bg2[,-4]); tmp3$Type<-c(rep("pGenes", 5), rep("Other genes", 8))
tmp3$effect<-factor(tmp3$effect, levels = rev(c("upstream_gene_variant", "downstream_gene_variant", "intron_variant", "missense_variant", "3_prime_UTR_variant", "splice_region_variant", "synonymous_variant", "5_prime_UTR_variant")))

library(viridis)
ggplot(tmp3, aes(x = Type, y = ratio, fill = effect)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = NULL, y = "Proportion of functional region", fill = "Effect") +
  scale_fill_viridis_d(option = "F", name = "Effect")+
  mytheme + theme(axis.text.x = element_text(angle = 45, face="bold", hjust = 1))

## conditional QTL analysis by QTLtools ####
# conditional analysis refer to https://qtltools.github.io/qtltools/. Code see pQTL_formal_code.docx, output skin_protein_all_qqnorm.conditional.txt
con<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.conditional.txt",header = F)
names(con)<-c("gene_id","chr_gene","start_gene","end_gene","strand","num_var","tss_distance","variant_id","chr_var","start_var","end_var","signal_rank","pval_forward","slope_forward","top_forward","pf_pt","pval_reverse","slope_reverse","top_reverse","pr_pt")
table(con$signal_rank) # All zeros indicate that there is no additional independent pQTL.
bb<-con[top_reverse==1]; bb<-bb[order(gene_id,signal_rank)]
identical(length(bb$gene_id),length(unique(bb$gene_id))) # TRUE，Consistent with the above, this indicates that no protein has two independent pQTLs.

## trans-pQTL identification ####
# Using the approximate permutation method of QTLtools trans on Linux, code see in pQTL_formal_code.docx, output  skin_protein_all_qqnorm.txt.QTLtools.gz.all.trans.adjust.hits.sig005.txt.gz
xx1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.QTLtools.gz.all.trans.adjust.hits.sig005.txt.gz", header = F) # 这是在 protein 水平上矫正过的FDR<0.05的
nrow(xx1) # 3: P23497(SP100)~20:5746925:C:G(rs6085261) ; P23497(SP100)~20:5746954:C:T(rs6085262)

## RegulomeDB_rank annotation ####
regulome<-fread("/Volumes/ElementsSE_syq/reference/RegulomeDB_rank.tsv", header = T) # The file is too large, please download it here: https://regulomedb.org/regulome-search/
regulome<-regulome[,1:6]; regulome[, chrom:=gsub("chr", "", chrom)]; names(regulome)[c(1,3,6)]<-c("chr_var", "pos_var", "score")

pro<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableS2_skin_pQTL_allSig_pairs.csv", skip=3, header = T) ; pro[, chr_var:=as.character(chr_var)]
pro1<-merge(pro, regulome[, c(1, 3:6)], by=c("chr_var","pos_var"))
pro1[,ranking2:= as.integer(sub("^(\\d).*", "\\1", ranking))]
nrow(pro1[ranking2<=2]); nrow(pro1[ranking2>2]) # 1259; 363; /1622=77.6%; 22.4%

# enrichment proportion of all test variants (background variants)
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
nrow(pro_nn2[ranking2<=2]); nrow(pro_nn2[ranking2>2]) # 1,515,102; 1,940,357; /3455459=43.85%; 56.15%
chisq.test(matrix(c(1259, 363, 1515102, 1940357), byrow = T, nrow = 2))$`p.value` 

# stack plot
aa<-data.table(qtl=c("pQTL","pQTL","bg","bg"), replication=c("Regulome","No regulome","Regulome","No regulome"), num=c(1259, 363, 1515102, 1940357)) # consider all pQTLs
aa<-data.table(qtl=c("pQTL","pQTL","bg","bg"), replication=c("Regulome","No regulome","Regulome","No regulome"), num=c(23, 11, 2678, 2425)) # consider only lead pQTLs/variants
aa[, prop := 100*num / sum(num), by = qtl]
aa$qtl<-factor(aa$qtl,levels = c( "bg", "pQTL"))
aa$replication<-factor(aa$replication,levels = c("Regulome","No regulome"))
aa[, fill_group := paste(qtl, replication, sep = "_")]

ggplot(aa, aes(x = qtl, y = prop, fill = fill_group)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#F4EDCA", "#C4961A", "#C3D7A4", "#52854C")) +
  # scale_fill_manual(values = c("mistyrose","indianred","lightsteelblue","steelblue")) +
  labs(x = '', y = 'Proportion of functional variants (%)') +
  scale_x_discrete(labels = c("bg"="Background", "pQTL"="pQTLs"))+
  coord_flip() +
  mytheme + theme(legend.position = "none")

