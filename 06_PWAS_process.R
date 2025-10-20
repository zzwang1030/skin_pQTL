library(data.table); library(ggplot2)
library(ggmanh); library(SeqArray)
library(qvalue); library(VennDiagram)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=9,face = "italic") #设置图例字体的大小及字体https://ggplot2.tidyverse.org/articles/ggplot2-specs.html
  #,axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1,color="black"),
  #,axis.line = element_line(size = 1.2),axis.ticks=element_line(size=1.2),#坐标轴及 ticks 的粗度
)

#### data prepare ####
name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

gene_name<-fread("~/Desktop/reference/Homo_sapiens.GRCh38.110.geneID_trID_name.txt",header = T)
gene_name<-gene_name[order(gene_id, chr,start)]
gene_tss<-gene_name[!duplicated(gene_id)]
gene_name<-unique(gene_name[,.(gene_id,gene_name,gene_type)])

tss<-fread("~/Desktop/reference/ensembl_uniprot_id_pos_biomart.txt",header = T) # download from ensembl biomart and manually curated/add ~50 proteins by syq
pro_gene<-unique(tss[,.(gene_id,uniprot_id)])
pro_gene<-pro_gene[uniprot_id!=""]

###################### PWAS, MR and coloc ###################### 
#### integration analysis ####
pp<-fread("~/Desktop/省皮/project/pQTL/PWAS/all_disease_PWAS.csv", header = T); names(pp)[1]<-"pro_id"
mm<-fread("~/Desktop/省皮/project/pQTL/MR/all_disease_MR.csv", header = T); names(mm)[1]<-"pro_id"
cc<-fread("~/Desktop/省皮/project/pQTL/colocalization_disease/syq_coloc/all_disease_coloc.csv", header = T); names(cc)[1]<-"pro_id"
table(pp$disease); table(mm$disease); table(cc$disease)

qq<-merge(merge(pp, mm, by=c("pro_id","disease"), all = T), cc, by=c("pro_id","disease"), all = T)
qq1<-merge(qq, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
qq1<- unique(qq1, by = c("gene_name", "disease")); nrow(qq1); length(unique(qq1$pro_id)); length(unique(qq1$gene_name)) # 48634; 5149; 5186
qq1<-qq1[order(disease, TWAS.Q)]

#### lep_eas ####
hh<-qq1[disease=="lep_eas" & hsq_pv<0.05]

## PWAS manhattan ####
library(ggmanh); library(SeqArray)
hh1<-hh[!is.na(TWAS.Q)]
hh1[,c("CHR","P0","TWAS.P"):=.(as.numeric(CHR),as.numeric(P0),as.numeric(TWAS.P))]
pt<- max(hh1[TWAS.Q<=0.05,]$TWAS.P)*1.01; pt # ; 2.97e-4

hh1[, lable_name:=ifelse(TWAS.Q<=0.05, gene_name, "")] # for lable sig gene_name
hh1[,color_class:=ifelse(TWAS.Q<=0.05, "sig", "ns")]; highlight_colormap <- c("ns" = "grey", "sig" = "maroon")

manhattan_plot(hh1, pval.colname = "TWAS.P", chr.colname = "CHR", pos.colname = "P0", y.label = "-log10(P)", thin=T, # thin 可以稀疏要画的点，加快画图速度
               signif = pt, point.size=2, label.font.size=5,
               label.colname="lable_name", 
               highlight.colname = "color_class", highlight.col = highlight_colormap, color.by.highlight = TRUE) + 
  coord_cartesian(ylim = c(0, 10))+ mytheme + theme(legend.position="none")

## MR ####
hh<-qq1[disease=="lep_eas" & hsq_pv<0.05]
nrow(hh[TWAS.Q<=0.05]); hh[TWAS.Q<=0.05]$gene_name  # 11; "ALDH2","TSFM","DECR1","ANKRD13A,NAA25","APOBR","GBA1,HGFAC,DDB2,MPST,OARD1"   
nrow(hh[wald_pval<=0.05|ivw_pval<=0.05]) # 207
nrow(hh[TWAS.Q<=0.05 & (wald_pval<=0.05|ivw_pval<=0.05)]); # 4
hh[TWAS.Q<=0.05 & (wald_pval<=0.05|ivw_pval<=0.05)]$gene_name # "ALDH2" "NAA25" "APOBR" "DDB2"  


## MR forest plot 
library(ggplot2);library(ggpubr); library(patchwork)

hh[, OR := ifelse(!is.na(wald_beta), exp(wald_beta), exp(ivw_beta))]
hh[, OR_lower := ifelse(!is.na(wald_beta), exp(wald_beta - 1.96 * wald_se), exp(ivw_beta - 1.96 * ivw_se))]
hh[, OR_upper := ifelse(!is.na(wald_beta), exp(wald_beta + 1.96 * wald_se), exp(ivw_beta + 1.96 * ivw_se))]
hh[, OR_CI := sprintf("%.2f (%.2f–%.2f)", OR, OR_lower, OR_upper)]

jj<-hh[TWAS.Q<=0.05 ]
jj[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]; jj<-jj[order(MR_pval)] # 按照 MR 的 pvalue 进行排序
jj$gene_name<-factor(jj$gene_name, levels = rev(jj$gene_name))

ggplot(jj, aes(x = OR, y = gene_name)) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2, linewidth = 1, color = "darkblue") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.2, 1.8)) +
  xlab("OR (95% CI)") + ylab(NULL)+ mytheme + theme(axis.text.y = element_text(face = "bold.italic"))

tmp <- jj[, .(`Protein`=gene_name, `OR (95% CI)`=OR_CI, `MR Pvalue` = formatC(MR_pval, format = "e", digits = 2))]
ggtexttable(tmp, rows = NULL, theme = ttheme("light"))

## colocalization  ####
hh<-qq1[disease=="lep_eas"& hsq_pv<0.05]

## integration heatmap ####
library(patchwork); library(qvalue)
hh<-qq1[disease=="lep_eas" & hsq_pv<0.05]
hh[, MR_beta := ifelse(!is.na(wald_beta), wald_beta, ivw_beta)]; hh[,MR_beta_se:=ifelse(!is.na(wald_se), wald_se, ivw_se)]; hh[, MR_Z := ifelse(!is.na(wald_beta), wald_beta/wald_se, ivw_beta/ivw_se)]
hh[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]
hh1<-hh[,.(gene_name, TWAS.Z, TWAS.Q, MR_beta, MR_beta_se, MR_Z, MR_pval, PP.H4.abf_r1)]; 
hh2<-hh1[TWAS.Q<=0.05]; hh2<-hh2[order(TWAS.Z, decreasing = T)] #只在 PWAS-sig基因基础上画
fwrite(hh2, "~/Desktop/省皮/project/pQTL/manuscript/TableSX_integration_lep_eas.csv", col.names = T, row.names = F, sep=",", quote = F)

plot_dt<-melt.data.table(hh2[, .(gene_name, TWAS.Z, MR_Z, PP.H4.abf_r1)], id.vars = "gene_name", variable.name = "metric", value.name = "value")
plot_dt[, star := ""]
plot_dt[metric == "TWAS.Z" & hh2$TWAS.Q <= 0.05, star := "*"]
plot_dt[metric == "MR_Z" & hh2$MR_pval <= 0.05, star := "*"]
plot_dt[metric == "PP.H4.abf_r1" & hh2$PP.H4.abf_r1 >= 0.5, star := "*"]
plot_dt[, gene_name := factor(gene_name, levels = hh2[order(TWAS.Z, decreasing = T)]$gene_name)]

make_tile_plot <- function(dt, title, fill_scale, show_y = FALSE) {
  ggplot(dt, aes(x = metric, y = gene_name, fill = value)) +
    geom_tile(color = "grey", linewidth = 0.6) +
    geom_text(aes(label = star), size = 10, fontface = "bold") +
    scale_x_discrete(expand = c(0,0), labels = c("TWAS.Z" = "PWAS Z", "MR_Z" = "MR Z","PP.H4.abf_r1" = "Coloc PPH4")) + 
    scale_y_discrete(expand = c(0,0)) +
    fill_scale + 
    theme( axis.title = element_blank(), axis.text.x = element_text(size = 16, face = "bold"), axis.ticks.x = element_blank(),
           axis.text.y = if (show_y) element_text(size = 20, face = "bold.italic") else element_blank(), axis.ticks.y = element_blank(),
           plot.title = element_text(hjust = 0.5, size = 10), panel.grid = element_blank(),legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12, face = "bold"), panel.spacing = unit(0, "pt"), plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))
}

dt_twas<-plot_dt[metric == "TWAS.Z"]; dt_mr<-plot_dt[metric == "MR_Z"]; dt_coloc<-plot_dt[metric == "PP.H4.abf_r1"]
p1 <- make_tile_plot(dt_twas, "TWAS Z", scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "PWAS Z"), show_y = TRUE)
p2 <- make_tile_plot(dt_mr, "MR Z",   scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "MR Z"))
p3 <- make_tile_plot(dt_coloc,"Coloc PP", scale_fill_gradient(low = "white", high = "firebrick", name = "PP.H4"))

(p1 + p2 + p3) + plot_layout(nrow = 1, widths = c(1, 1, 1), guides = "collect") & theme(legend.position = "right",plot.margin = unit(c(0, 0, 0, 0), "pt"))

## ALDH2 case plot ####
## ALDH2 protein~rs671
qqnorm<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_protein_all_qqnorm.txt",header = T)
melt(qqnorm[,-1:-3],id.vars = "pro_id")->qq3; names(qq3)<-c("pro_id", "sample", "value")

sel_pro<-"P05091"; sel_snp<-"12:111803962:G:A" # ALDH2~rs671
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
table(qq4$gt) # 0|0 0|1 1|0 1|1 = 136   25  22   3

max_inten<- 1.1*max(qq4$value,na.rm = T); min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), size=2.5, alpha = 0.8,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = gt2,group=factor(gt2)), outlier.shape = NA, fill = '#00000000', width = 0.7, lwd=1.1, position = position_dodge(width = 0.85))+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  scale_color_manual(values = c("#4E79A7", "#59A14F", "#F28E2B"))+
  labs(x = 'Genotype', y = 'Normalized protein abundance',color = NULL)+
  mytheme + theme(legend.text = element_text(size = 16),legend.position = "none",legend.title=element_blank())

## ALDH2 RNA~rs671
qqnorm<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt",header = T)
names(qqnorm)[4]<-"pro_id"
melt(qqnorm[,-1:-3],id.vars = "pro_id")->qq3; names(qq3)<-c("pro_id", "sample", "value")

sel_pro<-"ENSG00000111275"; sel_snp<-"12:111803962:G:A" # ALDH2~rs671
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
table(qq4$gt) # 0|0 0|1 1|0 1|1 = 103   24  20   3

max_inten<- 1.1*max(qq4$value,na.rm = T); min_inten<- 1.1*min(qq4$value,na.rm = T)
ggplot(qq4, aes(x = gt2, y = value))+
  geom_jitter(aes(color = gt2), size=2.5, alpha = 0.8,position = position_jitter(width = 0.2))+
  geom_boxplot(aes(color = gt2,group=factor(gt2)), outlier.shape = NA, fill = '#00000000', width = 0.7, lwd=1.1, position = position_dodge(width = 0.85))+
  coord_cartesian(ylim = c(min_inten, max_inten)) + 
  scale_color_manual(values = c("#4E79A7", "#59A14F", "#F28E2B"))+
  labs(x = 'Genotype', y = 'Normalized mRNA abundance',color = NULL)+
  mytheme + theme(legend.text = element_text(size = 16),legend.position = "none",legend.title=element_blank())

## ALDH2 rs671: GWAS, maf
gg<-data.table(class=c("Control","Control","Control","Leprosy","Leprosy","Leprosy"),
               GT=c("AA","GA","GG","AA","GA","GG"), cnt=c(245, 2420, 5646, 285, 2251, 4748)) # from LWC
gg[,all_cnt_byClass:=sum(cnt), by=.(class)]; gg[,freq_byClass:= 100*cnt/all_cnt_byClass]
gg[,all_cnt_byGT:=sum(cnt), by=.(GT)]; gg[,freq_byGT:= 100*cnt/all_cnt_byGT]
gg$GT<-factor(gg$GT, levels = c("GG","GA","AA"))


ggplot(gg, aes(x = GT, y = freq_byGT, fill = class)) + # 展示比例
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("Control" = "#4575b4", "Leprosy" = "firebrick"))+
  coord_cartesian(ylim = c(40, 60)) + labs(x = 'Genotype', y = 'Percentage(%)') +
  scale_x_discrete(labels = c("GG"="G/G", "GA"="G/A","AA"="A/A"))+
  mytheme + theme(legend.title = element_blank(), legend.text = element_text(color = "black", face = "bold", size = 18), legend.key.size = unit(0.8, "cm"))


#### ad_eur ####
hh<-qq1[disease=="ad_eur" & hsq_pv<0.05]
## PWAS manhattan ####
library(ggmanh); library(SeqArray)
hh1<-hh[!is.na(TWAS.Q)]
hh1[,c("CHR","P0","TWAS.P"):=.(as.numeric(CHR),as.numeric(P0),as.numeric(TWAS.P))]
pt<- max(hh1[TWAS.Q<=0.05,]$TWAS.P)*1.01; pt # 7.56e-4

hh1[, lable_name:=ifelse(TWAS.Q<=0.05, gene_name, "")] # for lable sig gene_name
hh1[, color_class := ifelse(TWAS.Q<=0.05 & TWAS.Z >= 0, "sig_up", ifelse(TWAS.Q<=0.05 & TWAS.Z <= 0, "sig_down", "normal"))]
highlight_colormap <- c(sig_up = "red", sig_down = "blue", normal = "gray60") # 这个可以标注上下调，成簇不明显
# hh1[,color_class:=ifelse(TWAS.Q<=0.05, "sig", "ns")]; highlight_colormap <- c("ns" = "grey", "sig" = "maroon")

manhattan_plot(hh1, pval.colname = "TWAS.P", chr.colname = "CHR", pos.colname = "P0", y.label = "-log10(P)", thin=T, # thin 可以稀疏要画的点，加快画图速度
               signif = pt, point.size=2, label.font.size=5,
               label.colname="lable_name", 
               highlight.colname = "color_class", highlight.col = highlight_colormap, color.by.highlight = TRUE) + 
  # coord_cartesian(ylim = c(0, 10))+ 
  mytheme + theme(legend.position="none")

## MR ####
hh<-qq1[disease=="ad_eur"& hsq_pv<0.05]

## MR forest plot
library(ggplot2);library(ggpubr); library(patchwork)

hh[, OR := ifelse(!is.na(wald_beta), exp(wald_beta), exp(ivw_beta))]
hh[, OR_lower := ifelse(!is.na(wald_beta), exp(wald_beta - 1.96 * wald_se), exp(ivw_beta - 1.96 * ivw_se))]
hh[, OR_upper := ifelse(!is.na(wald_beta), exp(wald_beta + 1.96 * wald_se), exp(ivw_beta + 1.96 * ivw_se))]
hh[, OR_CI := sprintf("%.2f (%.2f–%.2f)", OR, OR_lower, OR_upper)]

jj<-hh[TWAS.Q<=0.05]
jj[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]; jj<-jj[order(MR_pval)] # 按照 MR 的 pvalue 进行排序
jj$gene_name<-factor(jj$gene_name, levels = rev(jj$gene_name))

ggplot(jj, aes(x = OR, y = gene_name)) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2, linewidth = 1, color = "darkblue") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.5, 1.6)) +
  xlab("OR (95% CI)") + ylab(NULL)+ mytheme + theme(axis.text.y = element_text(face = "bold.italic"))

tmp <- jj[, .(`Protein`=gene_name, `OR (95% CI)`=OR_CI, `MR Pvalue` = formatC(MR_pval, format = "e", digits = 2))]
ggtexttable(tmp, rows = NULL, theme = ttheme("light"))

## colocalization ####
hh<-qq1[disease=="ad_eur" & hsq_pv<0.05]

## integration heatmap ####
library(patchwork)
hh<-qq1[disease=="ad_eur" & hsq_pv<0.05]
hh[, MR_beta := ifelse(!is.na(wald_beta), wald_beta, ivw_beta)]; hh[,MR_beta_se:=ifelse(!is.na(wald_se), wald_se, ivw_se)]; hh[, MR_Z := ifelse(!is.na(wald_beta), wald_beta/wald_se, ivw_beta/ivw_se)]
hh[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]
hh1<-hh[,.(gene_name, TWAS.Z, TWAS.Q, MR_beta, MR_beta_se, MR_Z, MR_pval, PP.H4.abf_r1)]; hh2<-hh1[TWAS.Q<=0.05]; 
paper_sup<-readxl::read_excel("~/Desktop/省皮/project/pQTL/PWAS/PWAS_gene_paper_support.xlsx", sheet = 2,col_names = T); paper_sup<-as.data.table(paper_sup); names(paper_sup)[1]<-"gene_name"
paper_sup[,gene_name:=gsub("^\\d+:\\s*", "", gene_name)]; paper_sup[,paper2:=ifelse(is.na(paper), 0, 1)]
hh3<-merge(hh2, paper_sup[,.(gene_name, paper2)]); hh3<-hh3[order(TWAS.Z, decreasing = T)] 
hh3<-hh3[order(TWAS.Z, decreasing = T)] #只在 PWAS-sig基因基础上画
fwrite(hh3, "~/Desktop/省皮/project/pQTL/manuscript/TableSX_integration_ad_eur.csv", col.names = T, row.names = F, sep=",", quote = F)

plot_dt<-melt.data.table(hh3[, .(gene_name, TWAS.Z, MR_Z, PP.H4.abf_r1, paper2)], id.vars = "gene_name", variable.name = "metric", value.name = "value")
plot_dt[, star := ""]
plot_dt[metric == "TWAS.Z" & hh3$TWAS.Q <= 0.05, star := "*"]
plot_dt[metric == "MR_Z" & hh3$MR_pval <= 0.05, star := "*"]
plot_dt[metric == "PP.H4.abf_r1" & hh3$PP.H4.abf_r1 >= 0.5, star := "*"]
plot_dt[metric == "paper2" & hh3$paper2 == 1, star := "*"]
plot_dt[, gene_name := factor(gene_name, levels = hh2[order(TWAS.Z, decreasing = T)]$gene_name)]

make_tile_plot <- function(dt, title, fill_scale, show_y = FALSE) {
  ggplot(dt, aes(x = metric, y = gene_name, fill = value)) +
    geom_tile(color = "grey", linewidth = 0.6) +
    geom_text(aes(label = star), size = 10, fontface = "bold", vjust = 0.7) +
    scale_x_discrete(expand = c(0,0), labels = c("TWAS.Z" = "PWAS Z", "MR_Z" = "MR Z","PP.H4.abf_r1" = "Coloc PPH4", "paper2" = "Prev. Reported")) + 
    scale_y_discrete(expand = c(0,0)) +
    fill_scale + 
    theme( axis.title = element_blank(), axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1, vjust = 1), axis.ticks.x = element_blank(),
           axis.text.y = if (show_y) element_text(size = 20, face = "bold.italic") else element_blank(), axis.ticks.y = element_blank(),
           plot.title = element_text(hjust = 0.5, size = 10), panel.grid = element_blank(),legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12, face = "bold"), panel.spacing = unit(0, "pt"), plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))
}

dt_twas<-plot_dt[metric == "TWAS.Z"]; dt_mr<-plot_dt[metric == "MR_Z"]; dt_coloc<-plot_dt[metric == "PP.H4.abf_r1"]; dt_paper<-plot_dt[metric == "paper2"]; dt_paper$value <- factor(dt_paper$value, levels = c(0, 1))
p1 <- make_tile_plot(dt_twas, "TWAS Z", scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "PWAS Z"), show_y = TRUE)
p2 <- make_tile_plot(dt_mr, "MR Z",   scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "MR Z"))
p3 <- make_tile_plot(dt_coloc, "Coloc PP", scale_fill_gradient(low = "white", high = "firebrick", name = "PP.H4"))
p4 <- make_tile_plot(dt_paper, "Prev. Reported", scale_fill_manual(values = c("0" = "grey80", "1" = "firebrick"), labels = c("No reported", "Reported"), name = "Prev. Reported"))

(p1 + p2 + p3 + p4) + plot_layout(nrow = 1, widths = c(1, 1, 1, 1), guides = "collect") & theme(legend.position = "right",plot.margin = unit(c(0, 0, 0, 0), "pt"))
# 12 out 29 at least 2 methods excluding GBA1 have opposite effect direction between PWAS and MR; only 1 (AAMDC) supported by three methods: [ PMID ]36167494; 33772001

#### pso_eur ####
hh<-qq1[disease=="pso_eur" & hsq_pv<0.05]
## PWAS manhattan ####
library(ggmanh); library(SeqArray)
hh1<-hh[!is.na(TWAS.Q)]
hh1[,c("CHR","P0","TWAS.P"):=.(as.numeric(CHR),as.numeric(P0),as.numeric(TWAS.P))]
pt<- max(hh1[TWAS.Q<=0.05,]$TWAS.P)*1.01; pt # 4.38e-4

hh1[, lable_name:=ifelse(TWAS.Q<=0.05, gene_name, "")] # for lable sig gene_name
hh1[, color_class := ifelse(TWAS.Q<=0.05 & TWAS.Z >= 0, "sig_up", ifelse(TWAS.Q<=0.05 & TWAS.Z <= 0, "sig_down", "normal"))]
highlight_colormap <- c(sig_up = "red", sig_down = "blue", normal = "gray60") # 这个可以标注上下调，成簇不明显

manhattan_plot(hh1, pval.colname = "TWAS.P", chr.colname = "CHR", pos.colname = "P0", y.label = "-log10(P)", thin=T, # thin 可以稀疏要画的点，加快画图速度
               signif = pt, point.size=2, label.font.size=5,
               label.colname="lable_name", 
               highlight.colname = "color_class", highlight.col = highlight_colormap, color.by.highlight = TRUE) + 
  # coord_cartesian(ylim = c(0, 10))+ 
  mytheme + theme(legend.position="none")

## MR ####
hh<-qq1[disease=="pso_eur" & hsq_pv<0.05]

## MR forest plot
library(ggplot2);library(ggpubr); library(patchwork)

hh[, OR := ifelse(!is.na(wald_beta), exp(wald_beta), exp(ivw_beta))]
hh[, OR_lower := ifelse(!is.na(wald_beta), exp(wald_beta - 1.96 * wald_se), exp(ivw_beta - 1.96 * ivw_se))]
hh[, OR_upper := ifelse(!is.na(wald_beta), exp(wald_beta + 1.96 * wald_se), exp(ivw_beta + 1.96 * ivw_se))]
hh[, OR_CI := sprintf("%.2f (%.2f–%.2f)", OR, OR_lower, OR_upper)]

jj<-hh[TWAS.Q<=0.05]
jj[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]; jj<-jj[order(MR_pval)] # 按照 MR 的 pvalue 进行排序
jj$gene_name<-factor(jj$gene_name, levels = rev(jj$gene_name))

ggplot(jj, aes(x = OR, y = gene_name)) +
  geom_point(size = 3, color = "black") +
  geom_errorbarh(aes(xmin = OR_lower, xmax = OR_upper), height = 0.2, linewidth = 1, color = "darkblue") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.5, 1.6)) +
  xlab("OR (95% CI)") + ylab(NULL)+ mytheme + theme(axis.text.y = element_text(face = "bold.italic"))

tmp <- jj[, .(`Protein`=gene_name, `OR (95% CI)`=OR_CI, `MR Pvalue` = formatC(MR_pval, format = "e", digits = 2))]
ggtexttable(tmp, rows = NULL, theme = ttheme("light"))

## colocalization ####
hh<-qq1[disease=="pso_eur"& hsq_pv<0.05]
nrow(hh[TWAS.Q<=0.05]); hh[TWAS.Q<=0.05]$gene_name  # 16
nrow(hh[PP.H4.abf_r1>=0.5]); hh[PP.H4.abf_r1>=0.5]$gene_name # 3: "HNRNPAB" "KRT14"   "AP1M2"  
nrow(hh[TWAS.Q<=0.05 & PP.H4.abf_r1>=0.5]); # 1
hh[TWAS.Q<=0.05 & PP.H4.abf_r1>=0.5]$gene_name; # HNRNPAB

## integration heatmap ####
library(patchwork)
hh<-qq1[disease=="pso_eur" & hsq_pv<0.05]
hh[, MR_beta := ifelse(!is.na(wald_beta), wald_beta, ivw_beta)]; hh[,MR_beta_se:=ifelse(!is.na(wald_se), wald_se, ivw_se)]; hh[, MR_Z := ifelse(!is.na(wald_beta), wald_beta/wald_se, ivw_beta/ivw_se)]
hh[,MR_pval:=ifelse(!is.na(wald_pval), wald_pval, ivw_pval)]
hh1<-hh[,.(gene_name, TWAS.Z, TWAS.Q, MR_beta, MR_beta_se, MR_Z, MR_pval, PP.H4.abf_r1)]; hh2<-hh1[TWAS.Q<=0.05]; 
paper_sup<-readxl::read_excel("~/Desktop/省皮/project/pQTL/PWAS/PWAS_gene_paper_support.xlsx", sheet = 3, col_names = T); paper_sup<-as.data.table(paper_sup); names(paper_sup)[1]<-"gene_name"
paper_sup[,gene_name:=gsub("^\\d+:\\s*", "", gene_name)]; paper_sup[,paper2:=ifelse(is.na(paper), 0, 1)]
hh3<-merge(hh2, paper_sup[,.(gene_name, paper2)]); hh3<-hh3[order(TWAS.Z, decreasing = T)] 
fwrite(hh3, "~/Desktop/省皮/project/pQTL/manuscript/TableSX_integration_pso_eur.csv", col.names = T, row.names = F, sep=",", quote = F)

plot_dt<-melt.data.table(hh3[, .(gene_name, TWAS.Z, MR_Z, PP.H4.abf_r1, paper2)], id.vars = "gene_name", variable.name = "metric", value.name = "value")
plot_dt[, star := ""]
plot_dt[metric == "TWAS.Z" & hh3$TWAS.Q <= 0.05, star := "*"]
plot_dt[metric == "MR_Z" & hh3$MR_pval <= 0.05, star := "*"]
plot_dt[metric == "PP.H4.abf_r1" & hh3$PP.H4.abf_r1 >= 0.5, star := "*"]
plot_dt[metric == "paper2" & hh3$paper2 == 1, star := "*"]
plot_dt[, gene_name := factor(gene_name, levels = hh2[order(TWAS.Z, decreasing = T)]$gene_name)]

make_tile_plot <- function(dt, title, fill_scale, show_y = FALSE) {
  ggplot(dt, aes(x = metric, y = gene_name, fill = value)) +
    geom_tile(color = "grey", linewidth = 0.6) +
    geom_text(aes(label = star), size = 10, fontface = "bold", vjust = 0.7) +
    scale_x_discrete(expand = c(0,0), labels = c("TWAS.Z" = "PWAS Z", "MR_Z" = "MR Z","PP.H4.abf_r1" = "Coloc PPH4", "paper2" = "Prev. Reported")) + 
    scale_y_discrete(expand = c(0,0)) +
    fill_scale + 
    theme( axis.title = element_blank(), axis.text.x = element_text(size = 16, face = "bold", angle = 45, hjust = 1, vjust = 1), axis.ticks.x = element_blank(),
           axis.text.y = if (show_y) element_text(size = 20, face = "bold.italic") else element_blank(), axis.ticks.y = element_blank(),
           plot.title = element_text(hjust = 0.5, size = 10), panel.grid = element_blank(),legend.title = element_text(size = 14, face = "bold"),
           legend.text = element_text(size = 12, face = "bold"), panel.spacing = unit(0, "pt"), plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))
}

dt_twas<-plot_dt[metric == "TWAS.Z"]; dt_mr<-plot_dt[metric == "MR_Z"]; dt_coloc<-plot_dt[metric == "PP.H4.abf_r1"]; dt_paper<-plot_dt[metric == "paper2"]; dt_paper$value <- factor(dt_paper$value, levels = c(0, 1))
p1 <- make_tile_plot(dt_twas, "TWAS Z", scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "PWAS Z"), show_y = TRUE)
p2 <- make_tile_plot(dt_mr, "MR Z",   scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "MR Z"))
p3 <- make_tile_plot(dt_coloc, "Coloc PP", scale_fill_gradient(low = "white", high = "firebrick", name = "PP.H4"))
p4 <- make_tile_plot(dt_paper, "Prev. Reported", scale_fill_manual(values = c("0" = "grey80", "1" = "firebrick"), labels = c("No reported", "Reported"), name = "Prev. Reported"))

(p1 + p2 + p3 + p4) + plot_layout(nrow = 1, widths = c(1, 1, 1, 1), guides = "collect") & theme(legend.position = "right",plot.margin = unit(c(0, 0, 0, 0), "pt"))
# 7 out 16 at least 2 methods ; only 1 (HNRNPAB) supported by three methods: [ PMID ]15804702; 32697368

