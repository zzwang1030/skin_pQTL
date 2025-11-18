#### compared with other tissue pQTL ####
## common data ####
library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 24) + theme(
  axis.text = element_text(color = 'black',size=20), strip.background = element_blank(), 
  plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

## other pQTL datasets of multiple tissues from medRxiv https://doi.org/10.1101/2025.01.10.25320181
## protein species and abundance compare ####
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

tmp1<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_ProAbund.xlsx", sheet = 2); tmp1<-as.data.table(tmp1); tmp1[,tissue:="colon"]
tmp2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_ProAbund.xlsx", sheet = 3); tmp2<-as.data.table(tmp2); tmp2[,tissue:="heart"]
tmp3<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_ProAbund.xlsx", sheet = 4); tmp3<-as.data.table(tmp3); tmp3[,tissue:="liver"]
tmp4<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_ProAbund.xlsx", sheet = 5); tmp4<-as.data.table(tmp4); tmp4[,tissue:="lung"]
tmp5<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_ProAbund.xlsx", sheet = 6); tmp5<-as.data.table(tmp5); tmp5[,tissue:="thyroid"]
aa<-rbind(tmp1, tmp2, tmp3, tmp4, tmp5); aa[,Ensembl_ID:=gsub("[.]\\d+","", Ensembl_ID)]
aa1<-merge(aa, unique(name_id_reviewd[,1:2]), by.x="Ensembl_ID", by.y="gene_id")
aa1<-aa1[,.(pro_id, Ref_Intensity, tissue)]; names(aa1)[2]<-"intensity"
table(aa1$tissue) # colon   heart   liver    lung thyroid = 7735    7415    8615    9009    8357 

qq2 <-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename.csv",header  = T)
ff<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/20230514_153410_M-GSGC0290947_TISSUE_Report_pro_rename_filter.csv",header  = T)
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- qq3 %>% tidyr::separate_rows(PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] 
mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn); fit <- vsn2(mm); mm1 <- predict(fit, mm) 
mm2<-as.data.table(mm1); mm2$pro_id<-rownames(mm1)
mm2_long<-melt(mm2, id.vars = "pro_id"); names(mm2_long)<-c("pro_id","sample","intensity_log")
mm2_long[, intensity:=2^intensity_log] 
mm3<-mm2_long[,.(intensity_median=median(intensity,na.rm = T), intensity_mean=mean(intensity,na.rm = T)), by=.(pro_id)] 
mm3[,c("intensity_median_log", "intensity_mean_log"):=.(ifelse(is.na(intensity_median), NA, log10(intensity_median)), ifelse(is.na(intensity_mean), NA, log10(intensity_mean)))]
mm4<-mm3[,.(pro_id, intensity_median_log)]; mm4$tissue<-"skin"; names(mm4)[2]<-"intensity"

aa2<-rbind(aa1, mm4); table(aa2$tissue) # colon   heart   liver    lung    skin thyroid = 7735    7415    8615    9009    5956    8357
# Since the quantification scale of skin is not consistent with that of other tissues, a coefficient is added here to make the scale consistent for easier plotting.
aa2[,intensity_mean:= mean(intensity[tissue!="skin"], na.rm=T), by=.(pro_id)]
aa2[, intensity_skin:=intensity[tissue=="skin"], by=.(pro_id)]
aa2[,fc:=intensity_mean/intensity_skin, by=.(pro_id)];
tmp<-unique(aa2[,.(pro_id, fc)]); summary(tmp$fc) # median=10.949, mean=11.009
aa2[,intensity2:=ifelse(tissue=="skin", intensity*10.949, intensity)]
aa4<-dcast(aa2, pro_id~tissue, value.var = "intensity2")
aa4<-aa4[,.(pro_id, skin, heart, colon, thyroid, liver, lung)]

cols <- setdiff(names(aa4), "pro_id"); n <- length(cols) 
overlap_mat <- matrix(NA_integer_, nrow = n, ncol = n, dimnames = list(cols, cols)) 
for (i in seq_along(cols)) {
  for (j in seq_along(cols)) {
    overlap_mat[i, j] <- sum(!is.na(aa4[[cols[i]]]) & !is.na(aa4[[cols[j]]]))
  }
}

library(corrplot) # protein abundance
mycor<-cor(aa4[,-1], method = "spearman", use = "pairwise.complete.obs") # spearman 0.73, pearson 0.78
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
ee1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/blood_pQTL_combined_all_corrected.csv", header = T)
tmp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/blood_pQTL_combined_all_corrected.hg38.bed"); tmp[,id:=paste(V4, V5, sep="_")]
ee2<-merge(ee1, tmp[,c(3, 4, 7)], by="id"); names(ee2)[ncol(ee2)-1]<-"pos_hg38"; names(ee2)[ncol(ee2)]<-"pro_id"
ee2[, pos_id_hg38:=paste(chr, pos_hg38, ref, alt, sep = ":")]; ee2[, id_hg38:=paste(pro_id, pos_id_hg38, sep="_")]
names(ee2)[13]<-"cohort"; ee2[,study:="plasma"]

ff2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_pQTL.csv", header = T) 
ff2[,variant_id:=paste(chr, pos_hg38, ref_allele, effect_allele, sep=":")]
ff2[,id_hg38:=paste(pro_id, variant_id, sep = "_")]
tmp1<-ee2[,.(pro_id, pos_id_hg38, id_hg38, study, corrected_beta, beta_se, pvalue)]
names(tmp1)<-c("pro_id", "variant_id_hg38", "id_hg38", "study", "beta", "beta_se", "pvalue")
tmp2<-ff2[,.(pro_id, variant_id, id_hg38, study, beta, beta_se, pvalue)]; 
names(tmp2)<-c("pro_id", "variant_id_hg38", "id_hg38", "study", "beta", "beta_se", "pvalue")
ef<-rbind(tmp1, tmp2)

pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.allpairs.thin.txt.gz", header=T) # a demo for complete skin_protein_all_qqnorm.txt.gz.allpairs.txt
pp[,id_hg38:=paste(gene_id, variant_id, sep="_")]
ef2<-merge(pp, ef[,3:7], by="id_hg38")
fwrite(ef2, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_pQTL_Overlap_skin.csv", col.names=T, row.names=F, sep=",", quote=F)
ef2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_pQTL_Overlap_skin.csv", header = T)
ef2[study=="plasma"]$study<-"Plasma"
nrow(ef2); length(unique(ef2$gene_id)); length(unique(ef2$variant_id)) # 1883; 1082; 1786
table(ef2$study) # Colon   Heart   Liver    Lung  plasma  Thyroid = 178     185     185     102    1074     159 

## plot
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
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)] 
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

## plink ld ####
# pwd: /sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional
# extract all snps for all pGenes
# less TableS2_skin_pQTL_allSig_pairs.csv |grep -v '#'|tail -n +2 |cut -d ',' -f 1|sort -u > tmp
# zcat skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz | head -1 > skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt
# zcat skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz | tail -n +2 | grep -F -f tmp >> skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt

pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt", header=T)
names(pp)[1:2]<-c("pro_id","variant_id_hg38"); length(unique(pp$variant_id_hg38)) # 100471

ee1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code//blood_pQTL_combined_all_corrected.csv", header = T) 
tmp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code//blood_pQTL_combined_all_corrected.hg38.bed"); tmp[,id:=paste(V4, V5, sep="_")]
ee2<-merge(ee1, tmp[,c(3, 4, 7)], by="id"); names(ee2)[ncol(ee2)-1]<-"pos_hg38"; names(ee2)[ncol(ee2)]<-"pro_id"
ee2[, pos_id_hg38:=paste(chr, pos_hg38, ref, alt, sep = ":")]; ee2[, id_hg38:=paste(pro_id, pos_id_hg38, sep="_")]
names(ee2)[13]<-"cohort"; ee2[,study:="plasma"]

ff2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_pQTL.csv", header = T) 
ff2[,variant_id:=paste(chr, pos_hg38, ref_allele, effect_allele, sep=":")]
ff2[,id_hg38:=paste(pro_id, variant_id, sep = "_")]

tmp1<-ee2[,.(pro_id, pos_id_hg38, id_hg38, study)]; names(tmp1)<-c("pro_id", "variant_id_hg38", "id_hg38", "study")
tmp2<-ff2[,.(pro_id, variant_id, id_hg38, study)]; names(tmp2)<-c("pro_id", "variant_id_hg38", "id_hg38", "study")
ef<-rbind(tmp1, tmp2) # 28824; 28630; 5657; 14840; 6
nrow(ef); length(unique(ef$id_hg38)); length(unique(ef$pro_id)); length(unique(ef$variant_id_hg38)); length(unique(ef$study)) 

common_pro <- intersect(ef$pro_id, pp$pro_id); length(common_pro) # 28/34 pGenes were co-identified
ef_sub <- ef[pro_id %in% common_pro, .(pro_id, variant_ef = variant_id_hg38)]
pp_sub <- pp[pro_id %in% common_pro, .(pro_id, variant_pp = variant_id_hg38)]
pairs <- merge(pp_sub, ef_sub, by = "pro_id", allow.cartesian = TRUE) 
pairs<-unique(pairs); nrow(pairs) # 593767

all_snps <- unique(c(pairs$variant_pp, pairs$variant_ef))
length(all_snps); length(unique(pp_sub$variant_pp)); length(unique(ef_sub$variant_ef)) # 83365; 83253; 183
# plink2 --bfile /sh2/home/sunyuanqiang/reference/G1000_hg38/all_pop/ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged --extract all_snps_forLD.txt --make-bed --out /sh2/home/sunyuanqiang/reference/G1000_hg38/all_pop/ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged.snps_forLD
# plink --bfile ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged.snps_forLD --extract all_snps_forLD.txt --r2 gz --ld-window-kb 1000000 --ld-window 100000 --ld-window-r2 0.5 --out all_snps_forLD

ld<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/all_snps_forLD.ld.gz", header=T); ld<-ld[,c(3, 6, 7)]
pairs[, pair_id := paste(variant_pp, variant_ef, sep="_")]
ld[,pair_id:=paste(SNP_A, SNP_B, sep="_")]
ld[,pair_id2:=paste(SNP_B, SNP_A, sep="_")]
res<-merge(pairs, ld[,c(3,4)], by.x="pair_id", by.y="pair_id", all.x=T)
res2<-merge(res, ld[,c(3,5)], by.x="pair_id", by.y="pair_id2", all.x=T)

res2[,R2:=ifelse(!is.na(R2.x), R2.x, ifelse(!is.na(R2.y), R2.y, 0))]
res2[,R2:=ifelse(variant_pp==variant_ef, 1, R2)] # equal SNPs are not counted and are assigned a value of 1.
res2<-res2[, c(1:4, 7)]
length(unique(res2$pro_id)); nrow(res2); nrow(pairs) # 28; 593767=593767
fwrite(res2, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_6Tissue_pGenes_all_ld_r2_full.csv", col.names=T, row.names=F, sep="\t", quote=F)

## conditional by COJO ####
res2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_6Tissue_pGenes_all_ld_r2_full.csv", header=T)
res2[,pro_variant_pp:=paste(pro_id, variant_pp,sep="_")]; res2[,pro_variant_ef:=paste(pro_id, variant_ef,sep="_")]
pp[,pro_variant_pp :=paste(pro_id, variant_id_hg38,sep="_")]
res3<-merge(res2, pp[,.(pro_variant_pp, maf, pval_nominal, slope)])
res4<-merge(res3, ef[,.(id_hg38, study)], by.x="pro_variant_ef", by.y = "id_hg38", all.x = TRUE, allow.cartesian = TRUE); res4<-res4[order(study, pro_id, -R2)]
res5<-res4[R2>=0.8, 4:11]; nrow(res5); length(unique(res5$pro_id)) # 4483; 27
# split conditional SNP
outdir <- "/sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional"
res5[, {
  filename <- paste(pro_id, variant_pp, variant_ef, study, sep = "_")
  filepath <- file.path(outdir, filename)
  writeLines(variant_pp, filepath)
  NULL
}, by = seq_len(nrow(res5))]

# ll *:*|wc -l # 4483
# split skin pQTL for common pGenes
outdir <- "/sh2/home/sunyuanqiang/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional"
N_val <- 186
pp_sub2[, {
  tmp <- tstrsplit(variant_id_hg38, ":", fixed = TRUE)
  chr_pos <- paste0(tmp[[1]], ":", tmp[[2]], ":", tmp[[3]], ":", tmp[[4]])
  alt <- tmp[[4]]; ref <- tmp[[3]]
  # Reorganize into COJO input format
  df <- data.table(SNP = chr_pos, A1 = alt, A2 = ref, fre = maf, b = slope, # A1 is the effect allele, A2 is the reference allele. freq is the frequency of A1.
                   se = slope_se, p = format(pval_nominal, scientific = TRUE), N = N_val)
  filename <- paste0(unique(pro_id), "_skin"); filepath <- file.path(outdir, filename)
  fwrite(df, filepath, sep = "\t", quote = FALSE, col.names = TRUE) 
  NULL}, 
  by = pro_id
]

if (1) { 
  # case test
  # /sh2/home/sunyuanqiang/software/gcta/gcta64 --bfile /sh2/home/sunyuanqiang/projects/pQTL/formal/conditional_COJO/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.plink.nochr --cojo-file Q9Y6N5_skin --cojo-cond Q9Y6N5_15:45676237:T:C_15:45676237:T:C_Thyroid --diff-freq 1 --cojo-wind 100000 --cojo-collinear 0.99 --out Q9Y6N5_15:45676237:T:C_15:45676237:T:C_Thyroid.COJO.conditiona
  # Batch execution, generating a script for batch execution: `bash generate_COJO_sh`
  # bash generate_COJO_sh will generate the run_cojo_tasks.sh script file
  # awk '{print $0 " &"; if (NR % 10 == 0) print "wait"} END {if (NR % 10 != 0) print "wait"}' run_cojo_tasks.sh > run_cojo_tasks_parallel.sh
  # nohup bash run_cojo_tasks_parallel.sh & 
}

## process COJO-conditional res ####
##Integrate multiple protein files into one file
library(data.table); library(stringr)
myfiles <- list.files(path="~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res", pattern = "COJO\\.conditional\\.cma\\.cojo$", full.names = TRUE)
cc <- rbindlist(lapply(myfiles, function(f) { 
  dt <- fread(f, header=T); fname <- basename(f) 
  parts <- str_split(fname, pattern = "_", simplify = TRUE) 
  dt[, pro_id := parts[1]]; dt[, conditional_SNP_skin := parts[2]]; dt[, conditional_SNP_other := parts[3]]; dt[, study := str_split(parts[4], pattern = "\\.", simplify = TRUE)[1]]
  return(dt)
}))
# The NA values in these three columns are due to excessive collinearity with the conditional SNP (>0.8), so they are directly padded to 1 in this pC. conditional_SNP_skin/other/study is not NA and does not require padded.
cc[, bC:=ifelse(is.na(bC), 0, bC)]; cc[, bC_se:=ifelse(is.na(bC_se), 1, bC_se)]; cc[, pC:=ifelse(is.na(pC), 1, pC)] 
# The above file does not include information about the conditional SNP itself, so add it. Set pC to 1 and bC to 0; conditional_SNP_skin/other/study is not NA, so no padding is needed.
myfiles2 <- list.files(path="~/projects/pQTL/other_pQTL_datasets/getx_5tissue/conditional/res", pattern = "COJO\\.conditional\\.given\\.cojo$", full.names = TRUE)

res_list <- lapply(myfiles2, function(f) { 
  fname <- basename(f); parts <- strsplit(fname, "_")[[1]]
  pro_id <- parts[1]; cond_skin <- parts[2]; cond_other <- parts[3]; study <- sub("\\.COJO.*", "", parts[4])
  dt <- fread(f)
  dt[, `:=`(n = 186, freq_geno = freq, bC = 0, bC_se = 1, pC = 1, pro_id = pro_id, # Conditional SNP itself has pC set to 1 and bC set to 0.
            conditional_SNP_skin = cond_skin, conditional_SNP_other = cond_other, study = study)]
  return(dt)
})

cc2 <- rbindlist(res_list, fill = TRUE); nrow(cc2) # 4483
identical(cc2$SNP, cc2$conditional_SNP_skin) # TRUE
cc3<-rbind(cc, cc2)

# Add R2 information for conditional_SNP_skin/other
cc3[, pair_id:=paste(conditional_SNP_skin, conditional_SNP_other, sep="_")]; nrow(cc3) # 13161465
res2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_6Tissue_pGenes_all_ld_r2_full.csv", header=T); nrow(res2) # 593767
length(unique(cc3$pair_id)); length(unique(res2$pair_id)) # 4019, 593767
cc4<-merge(cc3, unique(res2[,.(pair_id, R2)]), by="pair_id", all.x=T) # add R2 
summary(cc4$R2) # 0.8005  0.9123  0.9443  0.9422  0.9726  2.0000 

# 将 pGene 的原始统计量与 conditional后的统计量，添加
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt", header=T)
pp1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pp2<-merge(pp, pp1[,.(gene_id, pval_nominal_threshold)], by="gene_id") 

pp3<-merge(pp2[,-3:-5], cc4[,.(pro_id, SNP, conditional_SNP_skin, conditional_SNP_other, R2, bC, bC_se, pC, study)], by.x=c("gene_id","variant_id"), by.y=c("pro_id","SNP"), all.x=T)
length(unique(pp$gene_id)); length(unique(pp3$gene_id)) # 34=34
length(unique(pp$variant_id)); length(unique(pp3$variant_id)) # 100471=100471, 
length(unique(pp3[is.na(pC)]$gene_id)); length(unique(pp3[is.na(pC)]$variant_id)) #  7 pGenes have NA. Of these 7, 6 are pGenes not detected in other tissues, and 1 is a pGene detected in other tissues but whose paired SNP R2 > 0.8.
# There are still 18262 SNPs without post-conditioning statistics, partly because 7 pGenes were completely undetectable; so why do some SNPs belonging to the same protein, under the same conditional SNP, show up as NA while others are conditional? Because these SNPs are bad SNPs with significantly different freq. However, after setting `--diff-freq 1`, this type of SNP no longer exists.
# At this point, bC, bC_se, and pC are all filled with the original slope, slop_se, and pval_nominal; conditional_SNP_skin/other and study are filled with none, and R2 is filled with 2.
pp3[,bC:=ifelse(is.na(bC), slope, bC)]; pp3[,bC_se:=ifelse(is.na(bC_se), slope_se, bC_se)]; pp3[,pC:=ifelse(is.na(pC), pval_nominal, pC)]
pp3[,conditional_SNP_skin:=ifelse(is.na(conditional_SNP_skin), "none", conditional_SNP_skin)]; pp3[,conditional_SNP_other:=ifelse(is.na(conditional_SNP_other), "none", conditional_SNP_other)]
pp3[,study:=ifelse(is.na(study), "none", study)]; pp3[,R2:=ifelse(is.na(R2), 2, R2)]
nrow(pp3[is.na(pC)]); table(pp3$study) # 0; Colon   Heart   Liver    Lung    none  plasma Thyroid

pp3 <- tidyr::extract( pp3, col = 'variant_id', into = c('chr', 'pos', 'ref', 'alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); pp3<-as.data.table(pp3)
pp3 <- tidyr::extract( pp3, col = 'conditional_SNP_skin', into = c('conditional_SNP_skin_chr', 'conditional_SNP_skin_pos', 'conditional_SNP_skin_ref', 'conditional_SNP_skin_alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); pp3<-as.data.table(pp3)
pp3[, pos:=as.numeric(pos)]; pp3[, conditional_SNP_skin_pos:=as.numeric(conditional_SNP_skin_pos)] # 便于后面画 manhantton 图

# Add the statistic of condition SNP in other tissues
ee1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/blood_pQTL_combined_all_corrected.csv", header = T) 
tmp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/blood_pQTL_combined_all_corrected.hg38.bed"); tmp[,id:=paste(V4, V5, sep="_")]
ee2<-merge(ee1, tmp[,c(3, 4, 7)], by="id"); names(ee2)[ncol(ee2)-1]<-"pos_hg38"; names(ee2)[ncol(ee2)]<-"pro_id"
ee2[, pos_id_hg38:=paste(chr, pos_hg38, ref, alt, sep = ":")]; ee2[, id_hg38:=paste(pro_id, pos_id_hg38, sep="_")]
names(ee2)[13]<-"cohort"; ee2[,study:="plasma"]
ee2<-ee2[, .(study, pro_id, pos_id_hg38, corrected_beta, pvalue)]; names(ee2)[4]<-"beta"
ff2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/medRxiv_5Tissues_pQTL.csv", header = T) 
ff2[,pos_id_hg38:=paste(chr, pos_hg38, ref_allele, effect_allele, sep=":")]
ff2<-ff2[, .(study, pro_id, pos_id_hg38, beta, pvalue)]
ef2<-rbind(ee2, ff2); ef2<-unique(ef2)
names(ef2)[4:5]<-c("conditional_SNP_other_beta", "conditional_SNP_other_pvalue")
ef2[,id_tmp:=paste(study, pro_id, pos_id_hg38, sep="_")]

pp3[,id_tmp:=paste(study, gene_id, conditional_SNP_other, sep="_")]
pp4<-merge(pp3, ef2[,4:6], by="id_tmp", all.x=T) 
pp4[, conditional_SNP_other_beta:=ifelse(is.na(conditional_SNP_other_beta), 0, conditional_SNP_other_beta)]
pp4[, conditional_SNP_other_pvalue:=ifelse(is.na(conditional_SNP_other_pvalue), 1, conditional_SNP_other_pvalue)] # NA is filled with beta and p as 0 and 1 respectively.
pp4<-pp4[, id_tmp:=NULL]
fwrite(pp4, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.condition.csv", col.names=T, row.names=F, sep=",", quote=F)

## condition analysis and plot ####
library(ggplot2); library(patchwork)
pp4<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.condition.csv", header=T) # here is the demo for complete skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.condition.csv, which is too large
name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
pp4<-merge(pp4, unique(name_id_reviewd[,c(1,4)]), by.x="gene_id", by.y="pro_id", all.x=T)

length(unique(pp4$gene_id)); length(unique(pp4[pval_nominal<=pval_nominal_threshold]$gene_id)) # 34=34
length(unique(ef2$pro_id)) # 5657
length(intersect(unique(pp4$gene_id), unique(ef2$pro_id))) # 28 pGenes were detected in other tissues.
length(setdiff(unique(pp4$gene_id), unique(ef2$pro_id))) # Six pGenes were not detected in other tissues, skin-unique protein

## Calculate the most significant p-value for each protein before and after each tissue condition.
# Grouped by gene × tissue × condition_SNP_skin, that is, for each gene × tissue, the minimum pC is taken for each condition_SNP_skin, and then the largest pC is selected as the representative pC for that gene × tissue, and compared with Pt. 
# This results in 2 fewer specific values than the initial method of directly taking the minimum pC. 
aa2 <- pp4[, { # The latter half (grouped per gene_id × tissue × conditional_SNP_skin)
  # The minimum pval_nominal value under the current gene_id × tissue × conditional_SNP_skin
  min_pN <- min(pval_nominal, na.rm = TRUE)
  row_min_pN <- .SD[pval_nominal == min_pN][1]
  
  # Maximum pC corresponding to skin_least_snp_id
  df_least <- .SD[variant_id == row_min_pN$variant_id]
  max_pC <- max(df_least$pC, na.rm = TRUE)
  row_max_pC <- df_least[pC == max_pC]
  if (nrow(row_max_pC) > 1) row_max_pC <- row_max_pC[which.max(R2)]
  row_max_pC <- row_max_pC[1]
  
  # Minimum pC under current gene × tissue × conditional_SNP_skin
  min_pC <- min(pC, na.rm = TRUE)
  row_min_pC <- .SD[pC == min_pC]
  if (nrow(row_min_pC) > 1) row_min_pC <- row_min_pC[which.max(R2)]
  row_min_pC <- row_min_pC[1]
  
  vals <- c( 
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
aa3 <- aa2[, .SD[1], by = .(gene_id, study)] # Within each gene × tissue, take the largest pC among the smallest pCs of each condition_SNP_skin, and use it as the representative pC for that gene × tissue.

aa31<-unique(aa3[,.(gene_id, skin_least_snp_pN, skin_least_snp_bN, skin_least_snp_SigPt)]) # Remove the P and beta of the lead SNP in the skin before condition removal
aa31_long<-melt(aa31, id.vars = "gene_id")
aa32<-aa3[,.(gene_id, study, skin_AftCond_least_pC, skin_AftCond_least_pC_bC)] # Extract ConSNP_skin, P, and beta after each tissue condition.
aa32_wide <- dcast(aa32, gene_id ~ study, value.var = setdiff(names(aa32), c("gene_id", "study")))
aa33<-melt.data.table(aa32_wide, id.vars = "gene_id") 
aa4<-rbind(aa31_long, aa33); aa4<-aa4[order(gene_id, variable)]
# At this point, the tissue, which is still NA, does not have a corresponding SNP hit when proving the condition. Therefore, it is filled with p and beta before the condition.
aa4[, value_imp := ifelse(grepl("^skin_AftCond_least_pC_[^_]+$", variable) & is.na(value), value[variable == "skin_least_snp_pN"][1], 
                          ifelse(grepl("skin_AftCond_least_pC_bC_", variable) & is.na(value), value[variable == "skin_least_snp_bN"][1], value)), by = gene_id]
aa4[, study := fifelse(grepl("AftCond", as.character(variable)), sapply(strsplit(as.character(variable), "_"), tail, 1), "Skin")] 
aa4[, type := ifelse(grepl("skin_least_snp_pN$|skin_AftCond_least_pC_[^_]+$", variable), "pvalue",
                     ifelse(grepl("skin_least_snp_bN$|skin_AftCond_least_pC_bC_", variable), "beta",
                            ifelse(grepl("^skin_least_snp_SigPt$", variable), "sig_pt", NA_character_)))]
aa4<-aa4[order(gene_id, study, type)]; aa4$value_imp<-as.numeric(aa4$value_imp)

aa6 <- dcast(aa4, gene_id + study ~ type, value.var = "value_imp")
aa6[, sig_pt:=ifelse(is.na(sig_pt), sig_pt[study=="Skin"], sig_pt), by=.(gene_id)]; 
aa6[study == "plasma", study := "Plasma"]
names(aa6)[2]<-"tissue"; aa6<-aa6[tissue!="none"]; aa6[, condition:=ifelse(tissue=="Skin", "BeforeCondition", "AfterCondition")]
name_id_reviewd<-fread("~/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
aa6<-merge(aa6, unique(name_id_reviewd[,c(1,4)]), by.x="gene_id", by.y="pro_id")
fwrite(aa6, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableSX_skin_pGenes_condition_by6Tissue.csv", row.names = F, col.names = T, )

## plot
aa6<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/code/TableSX_skin_pGenes_condition_by6Tissue.csv", header = T)
aa6[, sig_mark := ifelse(pvalue <= sig_pt, "*", "")]
aa6[, logp_bin := cut(-log10(pvalue), breaks = seq(0, ceiling(max(-log10(pvalue), na.rm = TRUE))+3, by = 4),
                      include.lowest = TRUE, right = FALSE )] 
n <- length(levels(aa6$logp_bin)); values <- seq(2, by = 3, length.out = n) 
aa6[, sig_num_tissue:=sum(pvalue <= sig_pt), by=.(gene_id)]; 
aa6[, sig_num_gene:=sum(pvalue <= sig_pt), by=.(tissue)]
gene_order <- unique(aa6[tissue == "Skin"][order(sig_num_tissue, -pvalue), gene_name])
tissue_order<-unique(aa6[order(-sig_num_gene), tissue])
aa6$tissue<-factor(aa6$tissue, levels = tissue_order); aa6$gene_name<-factor(aa6$gene_name, levels = gene_order)

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


# results
length(unique(aa6[sig_num_tissue==7,]$gene_name)) # 12/34=35.3; If the minimum value of pC is significant across all tissues (other + skin = 7), it proves that the condition has not been eliminated and is skin-unique.
unique_pro<-setdiff(unique(aa6$gene_id), unique(ef2$pro_id)) # Six were due to other tissues not being reported; Q9NWV4 (CZIB), Q15369 (ELOC), Q5T440 (IBA57), P07476 (IVL), Q8N5N7 (MRPL50), Q9HBL7 (PLGRKT).
# The remaining 6 are other tissues that reported this gene, but no condition was successfully established: ACADSB/ALDH2/FCGR3A/GAA/GSTM1/IGHG3
tmp<-unique(aa6[,.(tissue, sig_num_gene)]); tmp[,conditioned_num:= 34 - sig_num_gene] # How many pGenes were conditional in each tissue?
tmp2<-unique(aa6[,.(gene_name, sig_num_tissue)]); tmp2[,conditioned_num:= 7 - sig_num_tissue] # How many tissues were conditioned for each pGene?
