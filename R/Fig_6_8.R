library(data.table); library(ggplot2)
library(ggmanh); library(SeqArray)
library(qvalue); library(VennDiagram)

mytheme <- theme_classic(base_size = 26) + 
  theme(
    axis.text = element_text(color = 'black',size=20), 
    strip.background = element_blank(), 
    plot.title = element_text(hjust = 0.5), 
    plot.subtitle = element_text(hjust = 0.5), 
    legend.text=element_text(size=9,face = "italic")
)

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)


#### PWAS analysis ####
## The pipeline for the proteome-wide association analysis (PWAS) is provided in *PWAS_fusion_code.docx*.
#### colocalization pipeline ####
library(data.table); library(coloc)
## leprosy
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz") # demo
gwas<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/leprosy_GWAS_metadata_lwc_20240204.pos_id.demo.txt", header=T)
pp_gwas<-merge(pp, gwas, by.x="variant_id", by.y="pos_id")

sel_gene<-unique(pp_gwas$gene_id); 
pph4 <- NULL;
for (i in 1:length(sel_gene)) { 
  tmp<-pp_gwas[gene_id==sel_gene[i],]
  tmp<-tmp[order(variant_id, pval_nominal)]; 
  tmp<-tmp[!duplicated(variant_id)] 
  my.coloc<-coloc.abf(dataset1=list(beta=tmp$slope, 
                                    varbeta=(tmp$slope_se)^2, 
                                    snp=tmp$variant_id, 
                                    type='quant', 
                                    N=186, 
                                    MAF=tmp$maf), 
                      dataset2=list(beta=tmp$beta, 
                                    varbeta=(tmp$se)^2, 
                                    snp=tmp$variant_id, 
                                    type='cc', 
                                    N=15022, 
                                    s=0.52,  
                                    MAF=tmp$coded_af)
                      )
  summary_row <- data.table(pro_id = sel_gene[i], t(my.coloc$summary)); 
  pph4 <- rbind(pph4, summary_row)
}

## Atopic dermatitis (AD) and psoriasis were analyzed using a coloc pipeline following the same workflow as described above for leprosy.
## The full GWAS summary statistics for AD and psoriasis were obtained from the FinnGen project: https://finngen.gitbook.io/documentation/data-download.

#### MR pipeline ####
library(data.table); library(dplyr)
library(TwoSampleMR);library(MendelianRandomization)
library(MRPRESSO); library(mr.raps); library(purrr)
## leprosy eas
## input data 
skin_protein_permu <- fread('~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz')
skin_protein_allpairs <- fread('~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz')
leprosy_GWAS <- fread('~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/leprosy_GWAS_metadata_lwc_20240204.pos_id.demo.txt') 

skin_protein_allpairs_less005 <- merge(skin_protein_allpairs,
                                       skin_protein_permu[,c('gene_id','qval','pval_nominal_threshold')],
                                       by = 'gene_id',
                                       all.y = T
)%>%
  filter(
    pval_nominal <= pval_nominal_threshold
  )  %>% 
  mutate(
    chr = sapply(strsplit(variant_id,':'),'[',1),
    pos = sapply(strsplit(variant_id,':'),'[',2)
  ) %>%
  mutate_at(
    .vars = c('chr','pos'),
    .funs = as.integer
  )

skin_pro_less005_leprosy_data <- merge(skin_protein_allpairs_less005,
                                       leprosy_GWAS,
                                       by.x=c('chr','pos'),
                                       by.y=c('chr','pos_hg38'),
                                       suffixes=c('.skin_pro','.leprosy_gwas'))%>%
  mutate(
    skin_pro_ref_A1 = sapply(strsplit(variant_id,':'),'[',3),
    skin_pro_oth_A2 = sapply(strsplit(variant_id,':'),'[',4),
    Ref = ifelse(slope > 0,skin_pro_ref_A1,skin_pro_oth_A2),
    Oth = ifelse(slope > 0,skin_pro_oth_A2,skin_pro_ref_A1),
    x_beta = ifelse(slope > 0,slope,-1*slope),
    x_se = slope_se,
    x_p = pval_nominal,
    chr_pos = paste0(chr,'_',pos),
    OR = exp(slope) 
  ) 

skin_pro_less005_leprosy_to_clump <- skin_pro_less005_leprosy_data %>%
  select('chr',
         'chr_pos',
         'pos',
         'x_p', 
         'skin_pro_ref_A1', 
         'OR', 
         'x_beta', 
         'x_se', 
         'gene_id'
         )

temp_dat <- skin_pro_less005_leprosy_to_clump
temp_dat <- as.data.frame(temp_dat)
for (i in unique(temp_dat$gene_id)) {
  
  file_name <- i
  
  single_skin_pro_data <- subset(temp_dat,gene_id == i)
  single_skin_pro_data <- single_skin_pro_data %>%
    select('chr', 
           'chr_pos', 
           'pos', 
           'x_p', 
           'skin_pro_ref_A1', 
           'OR', 
           'x_beta', 
           'x_se'
           )
  
  colnames(single_skin_pro_data) <- c("CHR","SNP","BP","P", "A1","OR","BETA","SE")

  fwrite(single_skin_pro_data, 
         file = paste0("./MR_QTLtool_nominal_less005/", file_name, "_leprosy.txt"), 
         sep = '\t',
         quote = F,
         row.names = FALSE)
}


## plink  to clump
cd ./MR_QTLtool_nominal_less005/
  for protein_trait in $(ls *_leprosy.txt|cut -d '.' -f 1)
do
# Linkage disequilibrium (LD) relationships were inferred using 1000 Genomes Project data from East Asian populations. The complete BIM/BED/FAM files were not provided. You should download the full 1000 Genomes VCF data and generated plink bim/fam/bed files.
plink --bfile g1000_eas_hg38.pos_id \
--clump ./MR_QTLtool_nominal_less005/${protein_trait}.txt \
--clump-p1 1 \
--clump-r2 0.001 \
--clump-kb 1000 \
--out ./MR_QTLtool_nominal_less005/${protein_trait}_clump

done

## extract clumped data
skin_protein_less005_clump <- NULL
skin_protein_less005_clump <- as.data.frame(skin_protein_less005_clump)

for (i in unique(temp_dat$gene_id)) { 
  eval(parse(text = paste0(i,"_clump <- fread('./MR_QTLtool_nominal_less005/",i,"_leprosy_clump.clumped')")))
  eval(parse(text = paste0(i,"_clump$gene_id <- i")))
  eval(parse(text = paste0("skin_protein_less005_clump <- rbind(",
                           i,"_clump,skin_protein_less005_clump
  )")))
} 

skin_protein_less005_iv_leprosy <- merge(skin_pro_less005_leprosy_data,
                                         skin_protein_less005_clump,
                                         by.x = c('gene_id','chr_pos'),
                                         by.y = c('gene_id','SNP'),
                                         all.y = T) %>%
  mutate(
    y_beta =  ifelse(allele_major == Ref,beta_minor,-1*beta_minor),
    y_se = se,
    y_p = P_value
  ) %>%
  select(
    "chr",'pos','gene_id','variant_id','tss_distance','rsID',
    'Ref','Oth','maf','x_beta','x_se','x_p',
    'y_beta','y_se','y_p','ma_samples'
  )

skin_protein_less005_iv_leprosy  <- skin_protein_less005_iv_leprosy %>%
  group_by(gene_id) %>%
  mutate(
    k = n(),
    max_N = max(ma_samples),
    F_statistic = (x_beta / x_se)^2
  ) %>%
  ungroup() %>%
  as.data.frame() %>%
  rename(
    rsids = 'rsID'
  )%>%
  filter(F_statistic >= 5)

run_mr_per_gene <- function(dat_gene) {
  dat_gene <- as.data.frame(dat_gene)
  gene <- unique(dat_gene$gene_id)
  snp_n <- nrow(dat_gene)
  
  res_list <- list(gene_id = gene, n_SNP = snp_n)
  
  # if only one iv → Wald ratio
  if (snp_n == 1) {
    beta_X <- dat_gene$x_beta
    beta_Y <- dat_gene$y_beta
    se_X <- dat_gene$x_se
    se_Y <- dat_gene$y_se
    
    wald_ratio <- beta_Y / beta_X
    se_wald <- sqrt((se_Y^2 / beta_X^2) + (beta_Y^2 * se_X^2) / beta_X^4)
    z_stat <- wald_ratio / se_wald
    p_value <- 2 * pnorm(-abs(z_stat))
    
    res_list$wald_beta <- wald_ratio
    res_list$wald_se <- se_wald
    res_list$wald_pval <- p_value
    
  } else if (snp_n >= 2) {
    # Create an MRInput object (used for most methods)
    mr_input_obj <- tryCatch({
      mr_input(
        bx = dat_gene$x_beta,
        bxse = dat_gene$x_se,
        by = dat_gene$y_beta,
        byse = dat_gene$y_se,
        snps = dat_gene$variant_id
      )
    }, error = function(e) NULL)
    
    # if two or more ivs → IVW
    ivw <- tryCatch(mr_ivw(mr_input_obj, model = "random", robust = FALSE), 
                    error = function(e) NULL)
    if (!is.null(ivw)) {
      res_list$ivw_beta <- ivw@Estimate
      res_list$ivw_se <- ivw@StdError
      res_list$ivw_pval <- ivw@Pvalue
    }
  }
  
  return(as.data.frame(res_list))
}

res_df <- skin_protein_less005_iv_leprosy %>%
  group_split(gene_id) %>%
  map_df(run_mr_per_gene)

## protein-gene transform
protein_gene <- fread(".\\skin_protein_data_gene_pos.tsv")
result_all_mr <- merge(
  res_df,
  protein_gene[,c("PG.Genes",
                  "PG.ProteinAccessions",
                  "PG.ProteinDescriptions")],
  by.x = 'gene_id',
  by.y = 'PG.ProteinAccessions'
  )

write.csv(result_all_mr, 
          file = ".//MR_results_summary_leprosy.csv", 
          row.names = FALSE,
          quote = F)

## AD and psoriasis were analyzed using a MR pipeline following the same workflow as described above for leprosy.
## The full GWAS summary statistics for AD and psoriasis were obtained from the FinnGen project: https://finngen.gitbook.io/documentation/data-download.

#### integration analysis ####
pp<-fread("all_disease_PWAS.csv", header = T); 
names(pp)[1]<-"pro_id"
mm<-fread("all_disease_MR.csv", header = T); 
names(mm)[1]<-"pro_id"
cc<-fread("all_disease_coloc.csv", header = T); 
names(cc)[1]<-"pro_id"

qq<-merge(merge(pp, mm, by=c("pro_id","disease"), all = T), cc, by=c("pro_id","disease"), all = T)
qq1<-merge(qq, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
qq1<- unique(qq1, by = c("gene_name", "disease")); 
qq1<-qq1[order(disease, TWAS.Q)]

## leprosy ####
## PWAS manhattan: Fig. 6A
library(ggmanh); 
library(SeqArray)
hh1<-hh[!is.na(TWAS.Q)]
hh1[,c("CHR","P0","TWAS.P") := .(as.numeric(CHR), as.numeric(P0), as.numeric(TWAS.P))]
pt<- max(hh1[TWAS.Q<=0.05]$TWAS.P) 

hh1[, lable_name:=ifelse(TWAS.Q<=0.05, gene_name, "")] 
hh1[,color_class:=ifelse(TWAS.Q<=0.05, "sig", "ns")]; 
highlight_colormap <- c("ns" = "grey", "sig" = "indianred")

manhattan_plot(hh1, 
               pval.colname = "TWAS.P", 
               chr.colname = "CHR", 
               pos.colname = "P0", 
               y.label = "-log10(P)", 
               thin=T, 
               signif = pt, 
               point.size=1, 
               label.font.size=5,
               label.colname="lable_name", 
               highlight.colname = "color_class", 
               highlight.col = highlight_colormap, 
               color.by.highlight = TRUE) + 
  coord_cartesian(ylim = c(0, 10)) + 
  mytheme + 
  theme(legend.position="none")

## MR: Fig. 6B 
hh<-qq1[disease=="lep_eas"]
hh[, OR := ifelse(!is.na(wald_beta), exp(wald_beta), exp(ivw_beta))]
hh[, OR_lower := ifelse(!is.na(wald_beta), exp(wald_beta - 1.96 * wald_se), exp(ivw_beta - 1.96 * ivw_se))]
hh[, OR_upper := ifelse(!is.na(wald_beta), exp(wald_beta + 1.96 * wald_se), exp(ivw_beta + 1.96 * ivw_se))]
hh[, OR_CI := sprintf("%.2f (%.2f–%.2f)", OR, OR_lower, OR_upper)]

library(ggpubr); 
library(patchwork)

jj <- hh[TWAS.Q<=0.05]
jj[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]; 
jj <- jj[order(MR_pval)]
jj$gene_name <- factor(jj$gene_name, levels = rev(jj$gene_name))

ggplot(jj, aes(x = OR, y = gene_name)) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), 
                 height = 0.2, 
                 linewidth = 1, 
                 color = "darkblue"
                 ) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.2, 1.8)) +
  xlab("OR (95% CI)") + 
  ylab(NULL) + 
  mytheme + 
  theme(axis.text.y = element_text(face = "bold.italic"))

## integration heatmap: Fig. 6C
library(patchwork)
paper_sup<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS9_PWAS_gene_other_evidence.csv", header = T)
tmp<-paper_sup[Disease=="Leprosy",.(Gene_symbol, `Publication (PMID)`)]; 
names(tmp)<-c("gene_name", "publication")
tmp1<-tmp[!duplicated(gene_name)]
tmp1[, publication:=ifelse(publication=="Not reported", 0, 1)]

hh<-qq1[disease=="lep_eas"]
hh[, MR_beta := ifelse(!is.na(wald_beta), wald_beta, ivw_beta)]; 
hh[,MR_beta_se:=ifelse(!is.na(wald_se), wald_se, ivw_se)]; 
hh[, MR_Z := ifelse(!is.na(wald_beta), wald_beta/wald_se, ivw_beta/ivw_se)]
hh[, MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]
hh1<-hh[,.(gene_name, TWAS.Z, TWAS.Q, MR_beta, MR_beta_se, MR_Z, MR_pval, PP.H4.abf_r1)]; 
hh2<-hh1[TWAS.Q<=0.05]; 
hh3<-merge(hh2, tmp1, by="gene_name", all.x = T); 
hh3<-hh3[order(TWAS.Z, decreasing = T)] 

plot_dt<-melt.data.table(hh3[, .(gene_name, TWAS.Z, MR_Z, PP.H4.abf_r1, publication)], 
                         id.vars = "gene_name", 
                         variable.name = "metric", 
                         value.name = "value")
plot_dt[, star := ""]
plot_dt[metric == "TWAS.Z" & hh3$TWAS.Q <= 0.05, star := "*"]
plot_dt[metric == "MR_Z" & hh3$MR_pval <= 0.05, star := "*"]
plot_dt[metric == "PP.H4.abf_r1" & hh3$PP.H4.abf_r1 >= 0.5, star := "*"]
plot_dt[metric == "publication" & hh3$publication == 1, star := "*"]
plot_dt[, gene_name := factor(gene_name, levels = hh2[order(TWAS.Z, decreasing = T)]$gene_name)]

make_tile_plot <- function(dt, title, fill_scale, show_y = FALSE) {
  ggplot(dt, aes(x = metric, y = gene_name, fill = value)) +
    geom_tile(color = "grey", linewidth = 0.6) +
    geom_text(aes(label = star), size = 10, fontface = "bold", vjust = 0.7) +
    scale_x_discrete(expand = c(0,0), 
                     labels = c("TWAS.Z" = "PWAS Z", 
                                "MR_Z" = "MR Z", 
                                "PP.H4.abf_r1" = "Coloc PPH4", 
                                "publication" = "Prev. Reported")) + 
    scale_y_discrete(expand = c(0,0)) +
    fill_scale + 
    theme( axis.title = element_blank(), 
           axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1, vjust = 1), 
           axis.ticks.x = element_blank(),
           axis.text.y = if (show_y) element_text(size = 20, face = "bold.italic") else element_blank(), axis.ticks.y = element_blank(),
           plot.title = element_text(hjust = 0.5, size = 10), 
           panel.grid = element_blank(),
           legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12, face = "bold"), 
           panel.spacing = unit(0, "pt"), 
           plot.margin = unit(c(0, 0, 0, 0), "pt"), 
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
  }

dt_twas<-plot_dt[metric == "TWAS.Z"]; 
dt_mr<-plot_dt[metric == "MR_Z"]; 
dt_coloc<-plot_dt[metric == "PP.H4.abf_r1"]; 
dt_paper<-plot_dt[metric == "publication"];
dt_paper$value <- factor(dt_paper$value, levels = c(0, 1))

p1 <- make_tile_plot(dt_twas, "TWAS Z", 
                     scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "PWAS Z"), 
                     show_y = TRUE)

p2 <- make_tile_plot(dt_mr, "MR Z", 
                     scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "MR Z")
                     )

p3 <- make_tile_plot(dt_coloc, "Coloc PP", 
                     scale_fill_gradient(low = "white", high = "#7b3294", name = "PP.H4")
                     )

p4 <- make_tile_plot(dt_paper, "Prev. Reported", 
                     scale_fill_manual(values = c("0" = "grey80", "1" = "firebrick"), labels = c("No reported", "Reported"), name = "Prev. Reported")
                     )

(p1 + p2 + p3 + p4) + plot_layout(nrow = 1, widths = c(1, 1, 1, 1), guides = "collect") & theme(legend.position = "right",plot.margin = unit(c(0, 0, 0, 0), "pt"))


## ALDH2 pQTL~rs671: Fig. 6E
qqnorm<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt",header = T)
qq3<-melt(qqnorm[,-1:-3],id.vars = "pro_id"); 
names(qq3)<-c("pro_id", "sample", "value")

sel_pro<-"P05091"; 
sel_snp<-"12:111803962:G:A" 
# In Linux: Write this SNP to tmp.bed; bcftools view -R tmp.bed hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.vcf.gz -O v -o tmp.vcf
myvcf<-fread("tmp.vcf", skip = 490, header = T)
my_sel2<-grep("ID|NC",names(myvcf)); 
myvcf<-myvcf[,..my_sel2]
myvcf2<-melt.data.table(myvcf,id.vars = "ID")
names(myvcf2)<-c("id","sample","format"); 
myvcf2[,gt:=gsub(":.+","",format)]
myvcf2[,ds:=as.numeric(sapply(strsplit(format,split = ":"), "[", 2))] 

qq4<-merge(qq3[pro_id==sel_pro], myvcf2[id==sel_snp, .(id,sample, gt,ds)], by="sample")
qq4[,c("ref","alt"):=.(strsplit(id,":")[[1]][3], strsplit(id,":")[[1]][4])]
qq4[,gt2:=ifelse(gt=="0|0", paste(ref,ref,sep="/"), 
                 ifelse(gt=="1|1",paste(alt,alt,sep="/"), paste(ref,alt,sep="/")))
    ]
qq4$gt2<-factor(qq4$gt2,levels = c(unique(paste(qq4$ref,qq4$ref,sep="/")),unique(paste(qq4$ref,qq4$alt,sep="/")),unique(paste(qq4$alt,qq4$alt,sep="/"))))

max_inten<- 1.1*max(qq4$value,na.rm = T); 
min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), 
              size=2.5, 
              alpha = 0.8, 
              position = position_jitter(width = 0.2)
              )+
  geom_boxplot(aes(color = gt2,group=factor(gt2)), 
               outlier.shape = NA, 
               fill = '#00000000', 
               width = 0.7, 
               lwd=1.1, 
               position = position_dodge(width = 0.85)
               )+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  scale_color_manual(values = c("#4E79A7", "#59A14F", "#F28E2B"))+
  labs(x = 'Genotype', y = 'Normalized protein abundance',color = NULL)+
  mytheme + 
  theme(legend.text = element_text(size = 16), 
        legend.position = "none",legend.title=element_blank()
        )

## AD ####
hh<-qq1[disease=="ad_eur"]

## PWAS manhattan: Fig. 8A for AD
library(ggmanh); 
library(SeqArray)
hh1<-hh[!is.na(TWAS.Q)]
hh1[,c("CHR","P0","TWAS.P"):=.(as.numeric(CHR),as.numeric(P0),as.numeric(TWAS.P))]
pt<- max(hh1[TWAS.Q<=0.05,]$TWAS.P)

hh1[, lable_name:=ifelse(TWAS.Q<=0.05, gene_name, "")]
hh1[, color_class := ifelse(TWAS.Q<=0.05 & TWAS.Z >= 0, "sig_up", 
                            ifelse(TWAS.Q<=0.05 & TWAS.Z <= 0, "sig_down", "normal")
                            )
    ]
highlight_colormap <- c(sig_up = "red", sig_down = "blue", normal = "gray60")

manhattan_plot(hh1, 
               pval.colname = "TWAS.P", 
               chr.colname = "CHR", 
               pos.colname = "P0", 
               y.label = "-log10(P)", 
               thin=T,
               signif = pt, 
               point.size=2, 
               label.font.size=5,
               label.colname="lable_name", 
               highlight.colname = "color_class", 
               highlight.col = highlight_colormap, 
               color.by.highlight = TRUE
               ) + 
  mytheme + 
  theme(legend.position="none")

## integration heatmap: Fig. 8B for AD
library(patchwork)
paper_sup<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS9_PWAS_gene_other_evidence.csv", header = T)
tmp<-paper_sup[Disease=="AD",.(Gene_symbol, `Publication (PMID)`)]; 
names(tmp)<-c("gene_name", "publication")
tmp1<-tmp[!duplicated(gene_name)]
tmp1[, publication:=ifelse(publication=="Not reported", 0, 1)]

hh<-qq1[disease=="ad_eur"]
hh[, MR_beta := ifelse(!is.na(wald_beta), wald_beta, ivw_beta)]; 
hh[,MR_beta_se:=ifelse(!is.na(wald_se), wald_se, ivw_se)]; 
hh[, MR_Z := ifelse(!is.na(wald_beta), wald_beta/wald_se, ivw_beta/ivw_se)]
hh[, MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]
hh1<-hh[,.(gene_name, TWAS.Z, TWAS.Q, MR_beta, MR_beta_se, MR_Z, MR_pval, PP.H4.abf_r1)]; 
hh2<-hh1[TWAS.Q<=0.05]; 
hh3<-merge(hh2, tmp1, by="gene_name", all.x = T); 
hh3<-hh3[order(TWAS.Z, decreasing = T)] 

plot_dt<-melt.data.table(hh3[, .(gene_name, TWAS.Z, MR_Z, PP.H4.abf_r1, publication)], 
                         id.vars = "gene_name", 
                         variable.name = "metric", 
                         value.name = "value")
plot_dt[, star := ""]
plot_dt[metric == "TWAS.Z" & hh3$TWAS.Q <= 0.05, star := "*"]
plot_dt[metric == "MR_Z" & hh3$MR_pval <= 0.05, star := "*"]
plot_dt[metric == "PP.H4.abf_r1" & hh3$PP.H4.abf_r1 >= 0.5, star := "*"]
plot_dt[metric == "publication" & hh3$publication == 1, star := "*"]
plot_dt[, gene_name := factor(gene_name, levels = hh2[order(TWAS.Z, decreasing = T)]$gene_name)]

make_tile_plot <- function(dt, title, fill_scale, show_y = FALSE) {
  ggplot(dt, aes(x = metric, y = gene_name, fill = value)) +
    geom_tile(color = "grey", linewidth = 0.6) +
    geom_text(aes(label = star), size = 10, fontface = "bold", vjust = 0.7) +
    scale_x_discrete(expand = c(0,0), 
                     labels = c("TWAS.Z" = "PWAS Z", 
                                "MR_Z" = "MR Z", 
                                "PP.H4.abf_r1" = "Coloc PPH4", 
                                "publication" = "Prev. Reported")) + 
    scale_y_discrete(expand = c(0,0)) +
    fill_scale + 
    theme( axis.title = element_blank(), 
           axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1, vjust = 1), 
           axis.ticks.x = element_blank(),
           axis.text.y = if (show_y) element_text(size = 20, face = "bold.italic") else element_blank(), axis.ticks.y = element_blank(),
           plot.title = element_text(hjust = 0.5, size = 10), 
           panel.grid = element_blank(),
           legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12, face = "bold"), 
           panel.spacing = unit(0, "pt"), 
           plot.margin = unit(c(0, 0, 0, 0), "pt"), 
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
}

dt_twas<-plot_dt[metric == "TWAS.Z"]; 
dt_mr<-plot_dt[metric == "MR_Z"]; 
dt_coloc<-plot_dt[metric == "PP.H4.abf_r1"]; 
dt_paper<-plot_dt[metric == "publication"];
dt_paper$value <- factor(dt_paper$value, levels = c(0, 1))

p1 <- make_tile_plot(dt_twas, "TWAS Z", 
                     scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "PWAS Z"), 
                     show_y = TRUE)

p2 <- make_tile_plot(dt_mr, "MR Z", 
                     scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "MR Z")
)

p3 <- make_tile_plot(dt_coloc, "Coloc PP", 
                     scale_fill_gradient(low = "white", high = "#7b3294", name = "PP.H4")
)

p4 <- make_tile_plot(dt_paper, "Prev. Reported", 
                     scale_fill_manual(values = c("0" = "grey80", "1" = "firebrick"), labels = c("No reported", "Reported"), name = "Prev. Reported")
)

(p1 + p2 + p3 + p4) + plot_layout(nrow = 1, widths = c(1, 1, 1, 1), guides = "collect") & theme(legend.position = "right",plot.margin = unit(c(0, 0, 0, 0), "pt"))


## psoriasis ####
hh<-qq1[disease=="pso_eur"]

## PWAS manhattan: Fig. 8A for psoriasis
library(ggmanh); 
library(SeqArray)
hh1<-hh[!is.na(TWAS.Q)]
hh1[,c("CHR","P0","TWAS.P"):=.(as.numeric(CHR),as.numeric(P0),as.numeric(TWAS.P))]
pt <- max(hh1[TWAS.Q<=0.05,]$TWAS.P)

hh1[, lable_name:=ifelse(TWAS.Q<=0.05, gene_name, "")]
hh1[, color_class := ifelse(TWAS.Q<=0.05 & TWAS.Z >= 0, "sig_up", 
                            ifelse(TWAS.Q<=0.05 & TWAS.Z <= 0, "sig_down", "normal")
)
]
highlight_colormap <- c(sig_up = "red", sig_down = "blue", normal = "gray60")

manhattan_plot(hh1, 
               pval.colname = "TWAS.P", 
               chr.colname = "CHR", 
               pos.colname = "P0", 
               y.label = "-log10(P)", 
               thin=T,
               signif = pt, 
               point.size=2, 
               label.font.size=5,
               label.colname="lable_name", 
               highlight.colname = "color_class", 
               highlight.col = highlight_colormap, 
               color.by.highlight = TRUE
) + 
  mytheme + 
  theme(legend.position="none")

## integration heatmap: Fig. 8B for psoriasis
library(patchwork)
paper_sup<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS9_PWAS_gene_other_evidence.csv", header = T)
tmp<-paper_sup[Disease=="Psoriasis",.(Gene_symbol, `Publication (PMID)`)]; 
names(tmp)<-c("gene_name", "publication")
tmp1<-tmp[!duplicated(gene_name)]
tmp1[, publication:=ifelse(publication=="Not reported", 0, 1)]

hh<-qq1[disease=="pso_eur"]
hh[, MR_beta := ifelse(!is.na(wald_beta), wald_beta, ivw_beta)]; 
hh[,MR_beta_se:=ifelse(!is.na(wald_se), wald_se, ivw_se)]; 
hh[, MR_Z := ifelse(!is.na(wald_beta), wald_beta/wald_se, ivw_beta/ivw_se)]
hh[, MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]
hh1<-hh[,.(gene_name, TWAS.Z, TWAS.Q, MR_beta, MR_beta_se, MR_Z, MR_pval, PP.H4.abf_r1)]; 
hh2<-hh1[TWAS.Q<=0.05]; 
hh3<-merge(hh2, tmp1, by="gene_name", all.x = T); 
hh3<-hh3[order(TWAS.Z, decreasing = T)] 

plot_dt<-melt.data.table(hh3[, .(gene_name, TWAS.Z, MR_Z, PP.H4.abf_r1, publication)], 
                         id.vars = "gene_name", 
                         variable.name = "metric", 
                         value.name = "value")
plot_dt[, star := ""]
plot_dt[metric == "TWAS.Z" & hh3$TWAS.Q <= 0.05, star := "*"]
plot_dt[metric == "MR_Z" & hh3$MR_pval <= 0.05, star := "*"]
plot_dt[metric == "PP.H4.abf_r1" & hh3$PP.H4.abf_r1 >= 0.5, star := "*"]
plot_dt[metric == "publication" & hh3$publication == 1, star := "*"]
plot_dt[, gene_name := factor(gene_name, levels = hh2[order(TWAS.Z, decreasing = T)]$gene_name)]

make_tile_plot <- function(dt, title, fill_scale, show_y = FALSE) {
  ggplot(dt, aes(x = metric, y = gene_name, fill = value)) +
    geom_tile(color = "grey", linewidth = 0.6) +
    geom_text(aes(label = star), size = 10, fontface = "bold", vjust = 0.7) +
    scale_x_discrete(expand = c(0,0), 
                     labels = c("TWAS.Z" = "PWAS Z", 
                                "MR_Z" = "MR Z", 
                                "PP.H4.abf_r1" = "Coloc PPH4", 
                                "publication" = "Prev. Reported")) + 
    scale_y_discrete(expand = c(0,0)) +
    fill_scale + 
    theme( axis.title = element_blank(), 
           axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1, vjust = 1), 
           axis.ticks.x = element_blank(),
           axis.text.y = if (show_y) element_text(size = 20, face = "bold.italic") else element_blank(), axis.ticks.y = element_blank(),
           plot.title = element_text(hjust = 0.5, size = 10), 
           panel.grid = element_blank(),
           legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12, face = "bold"), 
           panel.spacing = unit(0, "pt"), 
           plot.margin = unit(c(0, 0, 0, 0), "pt"), 
           panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
    )
}

dt_twas<-plot_dt[metric == "TWAS.Z"]; 
dt_mr<-plot_dt[metric == "MR_Z"]; 
dt_coloc<-plot_dt[metric == "PP.H4.abf_r1"]; 
dt_paper<-plot_dt[metric == "publication"];
dt_paper$value <- factor(dt_paper$value, levels = c(0, 1))

p1 <- make_tile_plot(dt_twas, "TWAS Z", 
                     scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "PWAS Z"), 
                     show_y = TRUE)

p2 <- make_tile_plot(dt_mr, "MR Z", 
                     scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "MR Z")
)

p3 <- make_tile_plot(dt_coloc, "Coloc PP", 
                     scale_fill_gradient(low = "white", high = "#7b3294", name = "PP.H4")
)

p4 <- make_tile_plot(dt_paper, "Prev. Reported", 
                     scale_fill_manual(values = c("0" = "grey80", "1" = "firebrick"), labels = c("No reported", "Reported"), name = "Prev. Reported")
)

(p1 + p2 + p3 + p4) + plot_layout(nrow = 1, widths = c(1, 1, 1, 1), guides = "collect") & theme(legend.position = "right",plot.margin = unit(c(0, 0, 0, 0), "pt"))


