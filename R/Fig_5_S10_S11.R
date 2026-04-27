library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 24) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
)

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

#### protein abundance cross-tissue compare: Fig. S10A ####
aa1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/5Tissues_ProAbund.csv", header = T) # data from doi: https://doi.org/10.1101/2025.01.10.25320181
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
mm3<-mm2_long[,.(intensity_median=median(intensity,na.rm = T), 
                 intensity_mean=mean(intensity,na.rm = T)), 
              by=.(pro_id)] 

mm3[,c("intensity_median_log", "intensity_mean_log"):=
      .(ifelse(is.na(intensity_median), NA, log10(intensity_median)), 
        ifelse(is.na(intensity_mean), NA, log10(intensity_mean)))
    ]

mm4<-mm3[,.(pro_id, intensity_median_log)]; 
mm4$tissue<-"skin"; 
names(mm4)[2]<-"intensity"

aa2<-rbind(aa1, mm4)
# Since the scale of skin is not consistent with that of other tissues, a coefficient is added here to make the scale consistent for easier plotting.
aa2[,intensity_mean:= mean(intensity[tissue!="skin"], na.rm=T), by=.(pro_id)]
aa2[, intensity_skin:=intensity[tissue=="skin"], by=.(pro_id)]
aa2[,fc:=intensity_mean/intensity_skin, by=.(pro_id)];
tmp<-unique(aa2[,.(pro_id, fc)]); 
summary(tmp$fc) # median=10.949, mean=11.009
aa2[,intensity2:=ifelse(tissue=="skin", intensity*10.949, intensity)]
aa4<-dcast(aa2, pro_id~tissue, value.var = "intensity2")
aa4<-aa4[,.(pro_id, skin, heart, colon, thyroid, liver, lung)]

library(viridis); 
library(ggpointdensity)
aa4_long <- melt(as.data.table(aa4), 
                 id.vars = c("pro_id", "skin"), 
                 variable.name = "tissue", 
                 value.name = "proteome")
ggplot(aa4_long, aes(x = skin, y = proteome)) +
  geom_pointdensity(size = 1, alpha = 0.5) + 
  scale_color_viridis(option = "D") +
  geom_smooth(method = "lm", 
              se = TRUE, 
              color = "red", 
              linetype = "dashed", 
              linewidth = 1) +
  labs(x = "Skin proteome", y = "Other tissue proteome", color = NULL) +
  facet_wrap(~ tissue, nrow=1, scales = "free_y") + 
  mytheme

#### cross-tissue beta correlation: Fig. S10B & S10C ####
ef2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/5Tissues_pQTL_Overlap_skin.csv", header = T)
ef2[study=="plasma"]$study<-"Plasma"

## plot Fig. S10C
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
result<-data.table(result); 
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)], id.vars = "pt"); 
names(mydt2)<-c("pt","type","value")
mydt2[,scale_value:=ifelse(type=="pearson_r", value*100, value)]
mydt2$type<-factor(mydt2$type,levels = c("percent","pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = pt_list)
ggplot(mydt2[pt%in%c("1","0.05", "0.001", 1e-5)], aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + 
  labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(45, 100)) + 
  labs(x = 'P-value thresholds')+
  scale_fill_manual(values = c("percent" = "indianred", 
                               "pearson_r" = "steelblue"
                               ),
                    labels = c("percent" = "Percent of same direction (%)", 
                               "pearson_r" = "Pearson correlation"
                               )
                    )+
  mytheme + 
  theme(legend.position="top", 
        legend.direction = "horizontal", 
        legend.text=element_text(size=16,face = "bold"), 
        legend.title=element_blank())

# plot Fig. S10B
library(viridis); 
library(ggpointdensity); 
library(ggpubr)
pt <- 1; 
tmp<-ef2[pval_nominal<pt | pvalue<pt] 
tmp1<-tmp[,.(n = .N, 
             R = cor.test(slope, beta, method = "pearson")$estimate, 
             P = cor.test(slope, beta, method = "pearson")$`p.value`), 
          by=.(study)
          ]
ggplot(tmp, mapping = aes(x = slope, y = beta)) + 
  geom_point(color = "black", size = 0.6, alpha = 0.3) + 
  geom_hline(yintercept=0, color="black",linetype="dashed", linewidth=0.4) + 
  geom_vline(xintercept=0, color="black",linetype="dashed", linewidth=0.4)+
  geom_smooth(method = "lm", se = TRUE, color = "indianred",  linewidth= 1.2) + 
  labs(x = 'Effect size of skin pQTL', y = 'Effect size of other pQTL', color = NULL)+
  facet_wrap(~study, nrow=1, scales = "free")+
  mytheme + 
  theme(axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(angle = 45, 
                                   vjust = 0.95, 
                                   hjust=1, 
                                   color="black"
                                   )
        )  

#### heatmap: skin-specific expression: Fig. S10D ####
aa1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/5Tissues_ProAbund.csv", header = T) # data from doi: https://doi.org/10.1101/2025.01.10.25320181
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv",header = T)

qq4<-melt.data.table(qq3, id.vars = "pro_id"); 
names(qq4)<-c("pro_id","sample","intensity")
qq4[,intensity_log := ifelse(is.na(intensity), NA, log2(intensity))]; 
qq4[,intensity_log_median := median(intensity_log, na.rm = T), by=.(pro_id)]
qq5<-unique(qq4[,.(pro_id, intensity_log_median)]); 
qq5$tissue<-"skin"; names(qq5)[2]<-"intensity"

aa2<-rbind(aa1, qq5); 
aa2<-merge(aa2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
skin_uniq<-c("IVL", "CZIB", "RANGAP1", "IBA57","MRPL50", "PLGRKT") # 6 pGene reported only in skin
aa3<-aa2[gene_name%in%skin_uniq];
aa3[, intensity_mean:= mean(intensity[tissue!="skin"], na.rm=T), by=.(gene_name)]
aa3[, intensity_skin:=intensity[tissue=="skin"], by=.(gene_name) ]
aa3[,fc:=intensity_mean/intensity_skin]; 
summary(unique(aa3$fc)) # median=3.1, mean=3.2
aa3[,intensity2:=ifelse(tissue=="skin", intensity*3.1, intensity)]

aa4<-dcast(aa3, gene_name~tissue, value.var = "intensity2")
library(pheatmap) 
mat <- as.data.frame(aa4); 
rownames(mat) <- mat$gene_name; mat$gene_name <- NULL
mat <- mat[, c("skin", "colon", "liver", "lung", "thyroid", "heart")]
row_medians <- apply(mat, 1, median, na.rm = TRUE); 
mat <- mat[order(-row_medians), ] 
tissue_labels <- c("Skin", "Colon", "Liver", "Lung", "Thyroid", "Heart")
tissue_labels_bold <- parse(text = paste0("bold('", tissue_labels, "')"))
gene_labels_italic <- parse(text = paste0("italic('", rownames(mat), "')"))

pheatmap(t(as.matrix(mat)), 
         na_col = "grey", 
         cluster_rows = F,
         cluster_cols = F, 
         labels_col = gene_labels_italic, 
         angle_col = 45, 
         labels_row =  tissue_labels_bold,
         fontsize = 16, 
         fontsize_row = 16, 
         fontsize_col = 17, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)


#### Replication: LD aligment and conditional analysis ####
## plink ld ####
# less TableS4_skin_pQTL_allSig_pairs.csv |grep -v '#'|tail -n +2 |cut -d ',' -f 1|sort -u > tmp # 36 pGenes
# zcat skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz | head -1 > skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt
# zcat skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz | tail -n +2 | grep -F -f tmp >> skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt
# plink2 --bfile ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink --extract all_snps_forLD.txt --make-bed --out ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged.snps_forLD
# pairwise LD: plink --bfile ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.EUR.pos_id.plink.merged.snps_forLD --extract all_snps_forLD.txt --r2 gz --ld-window-kb 1000000 --ld-window 100000 --ld-window-r2 0.5 --out all_snps_forLD

ld<-fread("all_snps_forLD.ld.gz", header=T); 
ld<-ld[,c(3, 6, 7)]
pairs[, pair_id := paste(variant_pp, variant_ef, sep="_")]
ld[,pair_id:=paste(SNP_A, SNP_B, sep="_")]
ld[,pair_id2:=paste(SNP_B, SNP_A, sep="_")]
res<-merge(pairs, ld[,c(3,4)], by.x="pair_id", by.y="pair_id", all.x=T)
res2<-merge(res, ld[,c(3,5)], by.x="pair_id", by.y="pair_id2", all.x=T)

res2[,R2:=ifelse(!is.na(R2.x), R2.x, ifelse(!is.na(R2.y), R2.y, 0))]
res2[,R2:=ifelse(variant_pp==variant_ef, 1, R2)] # equal SNPs are assigned 1.
res2<-res2[, c(1:4, 7)]

## conditional by COJO ####
res2[,pro_variant_pp:=paste(pro_id, variant_pp,sep="_")]; 
res2[,pro_variant_ef:=paste(pro_id, variant_ef,sep="_")]
pp[,pro_variant_pp :=paste(pro_id, variant_id_hg38,sep="_")]
res3<-merge(res2, pp[,.(pro_variant_pp, maf, pval_nominal, slope)])
res4<-merge(res3, ef[,.(id_hg38, study)], by.x="pro_variant_ef", by.y = "id_hg38", all.x = TRUE, allow.cartesian = TRUE); 
res4<-res4[order(study, pro_id, -R2)]
res5<-res4[R2>=0.8, 4:11]

# Cutting conditional SNP
outdir <- "/5Tissues/"
res5[, {
  filename <- paste(pro_id, variant_pp, variant_ef, study, sep = "_")
  filepath <- file.path(outdir, filename)
  writeLines(variant_pp, filepath)
  NULL
  }, 
by = seq_len(nrow(res5))
]

# Cutting skin pQTL for common pGenes
outdir <- "/5Tissues/"
N_val <- 186
pp_sub2 <- pp[pro_id %in% common_pro,]
pp_sub2[, {
  tmp <- tstrsplit(variant_id_hg38, ":", fixed = TRUE)
  chr_pos <- paste0(tmp[[1]], ":", tmp[[2]], ":", tmp[[3]], ":", tmp[[4]])
  alt <- tmp[[4]]; ref <- tmp[[3]]
  # Reorganize into COJO input format
  df <- data.table(SNP = chr_pos, 
                   A1 = alt, 
                   A2 = ref, 
                   fre = maf, 
                   b = slope, 
                   se = slope_se, 
                   p = format(pval_nominal, scientific = TRUE), 
                   N = N_val)
  
  filename <- paste0(unique(pro_id), "_skin"); 
  filepath <- file.path(outdir, filename) 
  fwrite(df, filepath, sep = "\t", quote = FALSE, col.names = TRUE)
  NULL
  }, 
  by = pro_id 
  ]

if (1) { # in shell
  # see generate_COJO_sh # This will generate the run_cojo_tasks.sh script file.
  # awk '{print $0 " &"; if (NR % 10 == 0) print "wait"} END {if (NR % 10 != 0) print "wait"}' run_cojo_tasks.sh > run_cojo_tasks_parallel.sh
  # nohup bash run_cojo_tasks_parallel.sh 
}

## process COJO-conditional res ####
library(stringr)
myfiles <- list.files(path="/5Tissues/conditional_res", pattern = "COJO\\.conditional\\.cma\\.cojo$", full.names = TRUE)
cc <- rbindlist(lapply(myfiles, function(f) { # Integrate multiple protein files into one file
  dt <- fread(f, header=T); 
  fname <- basename(f) 
  parts <- str_split(fname, pattern = "_", simplify = TRUE) 
  dt[, pro_id := parts[1]]; 
  dt[, conditional_SNP_skin := parts[2]]; 
  dt[, conditional_SNP_other := parts[3]]; 
  dt[, study := str_split(parts[4], pattern = "\\.", simplify = TRUE)[1]]
  return(dt)
  })
)

cc[, bC:=ifelse(is.na(bC), 0, bC)]; cc[, bC_se:=ifelse(is.na(bC_se), 1, bC_se)]; cc[, pC:=ifelse(is.na(pC), 1, pC)] 

myfiles2 <- list.files(path="/5Tissues/conditional_res", pattern = "COJO\\.conditional\\.given\\.cojo$", full.names = TRUE)
res_list <- lapply(myfiles2, function(f) { # Batch processing
  fname <- basename(f); 
  parts <- strsplit(fname, "_")[[1]]
  pro_id <- parts[1]; 
  cond_skin <- parts[2]; 
  cond_other <- parts[3]; 
  study <- sub("\\.COJO.*", "", parts[4])
  dt <- fread(f)
  dt[, `:=`(n = 186, 
            freq_geno = freq, 
            bC = 0, 
            bC_se = 1, 
            pC = 1, 
            pro_id = pro_id, 
            conditional_SNP_skin = cond_skin, conditional_SNP_other = cond_other, study = study)
     ]
  return(dt)
  }
)

cc2 <- rbindlist(res_list, fill = TRUE)
cc3<-rbind(cc, cc2)

# Add R2 information for conditional_SNP_skin/other
cc3[, pair_id := paste(conditional_SNP_skin, conditional_SNP_other, sep="_")];
res2<-fread("/5Tissues/skin_6Tissue_pGenes_all_ld_r2_full.csv", header=T)
cc4<-merge(cc3, unique(res2[,.(pair_id, R2)]), by="pair_id", all.x=T) 

# Add the original statistics of pGene and the conditional statistics.
pp<-fread("/5Tissues/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.txt", header=T)
pp1<-fread("/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pp2<-merge(pp, pp1[,.(gene_id, pval_nominal_threshold)], by="gene_id") 

pp3<-merge(pp2[,-3:-5], cc4[,.(pro_id, SNP, conditional_SNP_skin, conditional_SNP_other, R2, bC, bC_se, pC, study)], by.x=c("gene_id","variant_id"), by.y=c("pro_id","SNP"), all.x=T)
pp3[,bC:=ifelse(is.na(bC), slope, bC)]; 
pp3[,bC_se:=ifelse(is.na(bC_se), slope_se, bC_se)]; 
pp3[,pC:=ifelse(is.na(pC), pval_nominal, pC)]
pp3[,conditional_SNP_skin:=ifelse(is.na(conditional_SNP_skin), "none", conditional_SNP_skin)]; 
pp3[,conditional_SNP_other:=ifelse(is.na(conditional_SNP_other), "none", conditional_SNP_other)]
pp3[,study:=ifelse(is.na(study), "none", study)]; 
pp3[,R2:=ifelse(is.na(R2), 2, R2)]

pp3 <- tidyr::extract( pp3, col = 'variant_id', into = c('chr', 'pos', 'ref', 'alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); 
pp3<-as.data.table(pp3)
pp3 <- tidyr::extract( pp3, col = 'conditional_SNP_skin', into = c('conditional_SNP_skin_chr', 'conditional_SNP_skin_pos', 'conditional_SNP_skin_ref', 'conditional_SNP_skin_alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); 
pp3<-as.data.table(pp3)
pp3[, pos:=as.numeric(pos)]; 
pp3[, conditional_SNP_skin_pos:=as.numeric(conditional_SNP_skin_pos)] 

## Condition Analysis and Graphing: Fig. 5B ####
pp4<-fread("/5Tissues/conditional_res/skin_protein_all_qqnorm.txt.gz.pGenes.allpairs.condition.csv", header=T)
pp4<-merge(pp4, unique(name_id_reviewd[,c(1,4)]), by.x="gene_id", by.y="pro_id", all.x=T)

# Grouped by gene × tissue × condition_SNP_skin, that is, for each gene × tissue, the minimum pC is taken for each condition_SNP_skin, 
# and then the largest pC is selected as the representative pC for that gene × tissue, and compared with Pt. This results in 2 fewer specific values than the initial method of directly taking the minimum pC.
aa2 <- pp4[, { 
  # Minimum pval_nominal under current gene_id × tissue × conditional_SNP_skin
  min_pN <- min(pval_nominal, na.rm = TRUE)
  row_min_pN <- .SD[pval_nominal == min_pN][1]
  
  # Maximum PC corresponding to skin_least_snp_id
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
    row_min_pN$variant_id, 
    row_min_pN$pval_nominal, 
    row_min_pN$slope, 
    row_min_pN$pval_nominal_threshold,
    row_max_pC$pC, 
    row_max_pC$bC, 
    row_max_pC$conditional_SNP_other,
    row_min_pC$pC, 
    row_min_pC$conditional_SNP_other,  
    row_min_pC$bC, row_min_pC$R2)
  
  names(vals) <- c("skin_least_snp_id", 
                   "skin_least_snp_pN", 
                   "skin_least_snp_bN", 
                   "skin_least_snp_SigPt", 
                   "skin_least_snp_AftCond_most_pC", 
                   "skin_least_snp_AftCond_most_bC", 
                   "skin_least_snp_AftCond_most_ConSNP_other", 
                   "skin_AftCond_least_pC", 
                   "skin_AftCond_least_pC_ConSNP_other", 
                   "skin_AftCond_least_pC_bC", 
                   "skin_AftCond_least_pC_ConSNP_R2"
                   )
  as.list(vals) 
  }, 
  by = .(gene_id, study, conditional_SNP_skin)
  ]

aa2$skin_AftCond_least_pC <- as.numeric(aa2$skin_AftCond_least_pC); 
aa2$skin_AftCond_least_pC_ConSNP_R2 <- as.numeric(aa2$skin_AftCond_least_pC_ConSNP_R2)
aa2<-aa2[order(gene_id, study, -skin_AftCond_least_pC, -skin_AftCond_least_pC_ConSNP_R2)]
aa3 <- aa2[, .SD[1], by = .(gene_id, study)] # Within each gene × tissue, take the largest pC among the smallest pCs of each condition_SNP_skin, and use it as the representative pC for that gene × tissue.

aa31<-unique(aa3[,.(gene_id, skin_least_snp_pN, skin_least_snp_bN, skin_least_snp_SigPt)]) # Remove the P and beta of the lead SNP in the skin before condition removal
aa31_long<-melt(aa31, id.vars = "gene_id")
aa32<-aa3[,.(gene_id, study, skin_AftCond_least_pC, skin_AftCond_least_pC_bC)] # Extract ConSNP_skin, P, and beta after retrieving each tissue condition
aa32_wide <- dcast(aa32, gene_id ~ study, value.var = setdiff(names(aa32), c("gene_id", "study")))
aa33<-melt.data.table(aa32_wide, id.vars = "gene_id") 
aa4<-rbind(aa31_long, aa33); 
aa4<-aa4[order(gene_id, variable)]
aa4[, value_imp := ifelse(grepl("^skin_AftCond_least_pC_[^_]+$", variable) & is.na(value), 
                          value[variable == "skin_least_snp_pN"][1], 
                          ifelse(grepl("skin_AftCond_least_pC_bC_", variable) & is.na(value), 
                                 value[variable == "skin_least_snp_bN"][1], value)), 
    by = gene_id]
aa4[, study := fifelse(grepl("AftCond", as.character(variable)), 
                       sapply(strsplit(as.character(variable), "_"), tail, 1), "Skin")]
aa4[, type := ifelse(grepl("skin_least_snp_pN$|skin_AftCond_least_pC_[^_]+$", variable), "pvalue",
                     ifelse(grepl("skin_least_snp_bN$|skin_AftCond_least_pC_bC_", variable), "beta",
                            ifelse(grepl("^skin_least_snp_SigPt$", variable), "sig_pt", NA_character_)))]
aa4<-aa4[order(gene_id, study, type)]; aa4$value_imp<-as.numeric(aa4$value_imp)

aa6 <- dcast(aa4, gene_id + study ~ type, value.var = "value_imp")
aa6[, sig_pt:=ifelse(is.na(sig_pt), sig_pt[study=="Skin"], sig_pt), by=.(gene_id)]; 
aa6[study == "plasma", study := "Plasma"]
names(aa6)[2]<-"tissue"; 
aa6<-aa6[tissue!="none"]; 
aa6[, condition:=ifelse(tissue=="Skin", "BeforeCondition", "AfterCondition")]
aa6<-merge(aa6, unique(name_id_reviewd[,c(1,4)]), by.x="gene_id", by.y="pro_id")
fwrite(aa6, "/TableSX_skin_pGenes_condition_by6Tissue.csv", row.names = F, col.names = T, )

## Bubble chart
library(patchwork)
aa6<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS8_skin_pGenes_replicated_by6Tissue.csv", header = T)
aa6[, sig_mark := ifelse(pvalue <= sig_pt, "*", "")] 
aa6[, logp_bin := cut(-log10(pvalue), 
                      breaks = seq(0, ceiling(max(-log10(pvalue), na.rm = TRUE))+3, by = 4),
                      include.lowest = TRUE, 
                      right = FALSE
                      )
    ] 
n <- length(levels(aa6$logp_bin)); 
values <- seq(2, by = 3, length.out = n) 

aa6[, sig_num_tissue:=sum(pvalue <= sig_pt), by=.(gene_id)]; 
aa6[, sig_num_gene:=sum(pvalue <= sig_pt), by=.(tissue)]
gene_order <- unique(aa6[tissue == "Skin"][order(sig_num_tissue, -pvalue), gene_name])
tissue_order<-unique(aa6[order(-sig_num_gene), tissue])
aa6$tissue<-factor(aa6$tissue, levels = tissue_order); 
aa6$gene_name<-factor(aa6$gene_name, levels = gene_order)

ggplot(aa6, aes(x = tissue, y = gene_name)) +
  geom_point(aes(size = logp_bin, color = beta)) +
  scale_size_manual(values = values, drop = FALSE) +
  scale_color_gradient2(low = "#2166AC", 
                        mid = "white", 
                        high = "#B22222",
                        midpoint = 0, 
                        na.value = "grey80"
                        ) +
  geom_text(aes(label = sig_mark), vjust = 0.8, size = 10, color = "black") +
  labs(x = NULL, y = NULL, size = NULL, color = NULL) +
  scale_x_discrete(position = "top")+
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 12, color = "black", angle = 0, vjust = 0),
        axis.text.y = element_text(face = "bold.italic", size = 12, color = "black")
  ) 


