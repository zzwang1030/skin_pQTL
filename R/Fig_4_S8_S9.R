library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 24) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
)

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

tss3 <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/ensembl_uniprot_id_pos_biomart.csv", header = T)

#### eQTL calling ####
# For the eQTL calling pipelines, see the shell script in pQTL_eQTL_calling.docx.
tmp<-fread("~/Desktop/省皮/project/pQTL/QTL_calling_formal/skin_RNA_all_qqnorm.txt")
fwrite(tmp[1:100, 1:100], "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.demo.txt", col.names = T, sep=",", quote = F, row.names = F)

#### eQTL charactering ####
## manhattan plot: Fig. S8A ####
library(ggmanh); 
library(SeqArray)
nn3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T) # demo
pp2<- tidyr::extract(nn3, col = 'variant_id', into = c('chr', 'pos','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F); 
pp2<-as.data.table(pp2)
pp2$chr<-factor(pp2$chr,levels = 1:22); 
pp2$pos<-as.numeric(pp2$pos)

gene_name<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/Homo_sapiens.GRCh38.110.geneID_name.txt",header = F) 
names(gene_name)<-c("gene_id","gene_name","gene_type"); 
pp2<-merge(pp2,gene_name[,-3],by="gene_id",all.x = T); 
pp2[,gene_name:=ifelse(is.na(gene_name), gene_id, gene_name)]
pp2<-pp2[order(gene_id, pval_nominal)]
pp2[,p_order:=1:.N, by=.(gene_id)]
pp2[,lable_name:=ifelse(p_order==1, gene_name, ""), by=.(gene_id)]

pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pt<-max(pp[qval<=0.05]$pval_nominal)
pp2[, lable_name2:=ifelse(pval_nominal<=pt, lable_name, "")] 
pp2[,color_class:=ifelse(pval_nominal<=pt,"sig","ns")]; 
highlight_colormap <- c("ns" = "grey", "sig" = "maroon")

manhattan_plot(pp2, 
               pval.colname = "pval_nominal", 
               chr.colname = "chr",
               pos.colname = "pos", 
               y.label = "-log10(P)", 
               thin=T,
               signif = 5e-8, 
               point.size=0.6, 
               label.font.size=3
               ) + 
  mytheme + 
  theme(legend.position="none")

## our eQTL vs gtex eQTL: Fig. S8B & S8C ####
library(viridis); 
library(ggpointdensity)

nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T) # demo
names(nn)[1]<-"pro_id"; 
nn<-nn[order(pro_id,pval_nominal,-abs(slope),slope_se)] 

gtex_nn<-fread("~/reference/GTEx_Analysis_v8_eQTL/Skin_Sun_Exposed_Lower_leg.v8.signif_variant_gene_pairs.txt.gz")  # data download from GTEx portal v8
nrow(gtex_nn); # 3244342
nrow(gtex_nn[pval_nominal<=pval_nominal_threshold]) # 3244342. The fact that they are the same indicates that they are SIG pairs filtered based on pval_nominal_threshold.
gtex_nn[,variant_id:=gsub("_b38","",gsub("chr","",variant_id))]; 
gtex_nn[,variant_id:=gsub("_",":",variant_id)]
gtex_nn[,gene_id:=gsub("[.]\\d+","",gene_id)]; 
names(gtex_nn)[2]<-"pro_id"
gtex_nn<-gtex_nn[order(gene_id,pval_nominal)]
total<-merge(nn, gtex_nn, by=c("pro_id","variant_id")); 

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

# plot Fig. S8C
result<-data.table(pt=result$pt, 
                   n_pair=result$n_pair,
                   n_same_direction=result$n_same_direction,
                   pearson_r=result$pearson_r,
                   p_value=result$p_value)
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)],id.vars = "pt"); 
names(mydt2)<-c("pt","type","value")
mydt2[, scale_value := ifelse(type=="pearson_r", value*100, value)] 
mydt2$type <- factor(mydt2$type, levels = c("percent","pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = c(1e-02, 1e-04, 1e-06, 1e-08))
ggplot(mydt2, aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + 
  labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(85, 100)) + 
  labs(x = 'P-value thresholds')+
  scale_fill_manual(values = c("percent" = "indianred", "pearson_r" = "steelblue"),
                    labels = c("percent" = "Percent of same direction (%)", "pearson_r" = "Pearson correlation"))+
  mytheme + 
  theme(legend.position="top", 
        legend.direction = "horizontal", 
        legend.title=element_blank())

# plot Fig. S8B
oo<-copy(total)
pt <- 1e-04
ggplot(oo[pval_nominal.x<pt & pval_nominal.x<pt], mapping = aes(x = slope.x, y = slope.y)) + 
  geom_pointdensity(size = 0.8) + 
  scale_color_viridis() +
  labs(x = 'Effect size in our study', y = 'Effect size in GTEx', color = NULL)+
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  geom_hline(yintercept=0, color="black",linetype="dashed") + 
  geom_vline(xintercept=0, color="black",linetype="dashed") +
  geom_smooth(method='lm', se=T, formula= y~x, col="red", linetype="dashed") +
  mytheme + 
  theme(legend.position = c(0.25, 0.9), 
        legend.direction = "horizontal", 
        legend.key.width = unit(1.2, "cm"), 
        legend.text = element_text(size = 12))

#### pQTL and eQTL replication #### 
## Count of QTL overlap: Fig. S8D ####
library(data.table); 
library(VennDiagram)

pro_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); # demo
names(pro_nn)[1]<-"pro_id"
names(pro_nn)[3:9] <- paste0(names(pro_nn)[3:9], "_pqtl")

rna_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); # demo
rna_nn1<-merge(rna_nn, unique(name_id_reviewd[,.(pro_id, gene_id)]), by="gene_id", all.x=T, allow.cartesian = TRUE)
names(rna_nn1)[3:9] <- paste0(names(rna_nn1)[3:9], "_eqtl")

nn<-merge(pro_nn, rna_nn1, by=c("pro_id", "variant_id"))
tmp<-nn[pval_nominal_pqtl < 0.01 & pval_nominal_eqtl < 0.01]
fwrite(tmp, "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", col.names=T, row.names=F, sep="\t", quote=F)

# Permutation of the number of overlaps
pt<-1e-2
pqtl_pairs <- unique(paste(nn[pval_nominal_pqtl < pt]$pro_id, nn[pval_nominal_pqtl < pt]$variant_id, sep = "_"))
pqtl_n <- length(pqtl_pairs)
eqtl_pairs <- unique(paste(nn[pval_nominal_eqtl < pt]$pro_id, nn[pval_nominal_eqtl < pt]$variant_id, sep = "_"))
eqtl_n <- length(eqtl_pairs)
true_overlap <- length(intersect(pqtl_pairs, eqtl_pairs)) 
background_pairs <- unique(paste(nn$pro_id, nn$variant_id, sep = "_")); 

n_perm <- 1000; 
set.seed(123)
perm_overlap <- numeric(n_perm)
for (i in 1:n_perm) {
  random_pqtl <- sample(background_pairs, pqtl_n)
  random_eqtl <- sample(background_pairs, eqtl_n)
  perm_overlap[i] <- length(intersect(random_pqtl, random_eqtl))
}
mean(perm_overlap) 
empirical_p <- (sum(perm_overlap >= true_overlap) + 1) / (length(perm_overlap) + 1) 
enrichment <- true_overlap / mean(perm_overlap)
fwrite(data.frame(overlap = perm_overlap), "skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.permutation.txt", col.names=T, row.names=F, sep="\t", quote=F)

tmp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.permutation.txt")
ggplot(tmp, aes(x=overlap)) +  
  geom_histogram(bins=50, 
                 fill="steelblue",
                 color='grey90', 
                 linewidth=0.5) +
  labs(x = "Overlap count under permutations", 
       y = "Frequency") + 
  mytheme




## correaltion of beta: Fig. 4A & 4B ####
nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt.gz", header = T); # demo
pt_list <- c(0.05, 1e-2, 1e-4, 1e-6) 
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

# plot Fig. 4B
result<-data.table(pt=result$pt, 
                   n_pair=result$n_pair,
                   n_same_direction=result$n_same_direction,
                   pearson_r=result$pearson_r,
                   p_value=result$p_value)
result[,percent := 100*n_same_direction/n_pair]

mydt2<-melt(result[,c(1, 4, 6)], id.vars = "pt"); 
names(mydt2)<-c("pt","type","value")
mydt2[, scale_value:=ifelse(type=="pearson_r", value*100, value)] 
mydt2$type <- factor(mydt2$type,levels = c("percent", "pearson_r")); 
mydt2$pt <- factor(mydt2$pt, levels = c(0.05, 1e-2, 1e-4, 1e-6))
ggplot(mydt2, aes(x=pt, y = scale_value, fill= type)) + 
  geom_col(position="dodge") + 
  scale_y_continuous(sec.axis = sec_axis(~ . / 100, name = "Pearson correlation")) + 
  labs(y="Percent of same direction (%)") +
  coord_cartesian(ylim = c(45, 100)) + 
  labs(x = 'P-value thresholds') +
  scale_fill_manual(values = c("percent" = "indianred", "pearson_r" = "steelblue"),
                    labels = c("percent" = "Percent of same direction (%)", "pearson_r" = "Pearson correlation")) + 
  mytheme + 
  theme(legend.position="top", 
        legend.direction = "horizontal", 
        legend.text=element_text(size=16,face = "bold"), 
        legend.title=element_blank())

# plot Fig. 4A
library(viridis); 
library(ggpointdensity); 
library(ggpubr)
oo<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", header = T); 
pt <- 1e-02
ggplot(oo[pval_nominal_pqtl<pt & pval_nominal_eqtl<pt], 
       mapping = aes(x = slope_pqtl, y = slope_eqtl)) +
  geom_pointdensity(size = 1.5) + 
  scale_color_viridis() +
  labs(x = 'Effect size of pQTL', 
       y = 'Effect size of eQTL', 
       color = NULL)+
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  geom_hline(yintercept=0, 
             color="black",
             linetype="dashed", 
             linewidth=0.4) + 
  geom_vline(xintercept=0, 
             color="black", 
             linetype="dashed", 
             linewidth=0.4) +
  geom_smooth(method='lm', 
              se=T, 
              formula= y~x, 
              col="red", 
              linetype="dashed") +
  mytheme + 
  theme(legend.position = c(0.25, 0.9), 
        legend.direction = "horizontal", 
        legend.key.width = unit(1.2, "cm"), 
        legend.text = element_text(size = 14))

## eQTL and pQTL opposite direction: Fig. S8E ####
oo<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.p001.txt.gz", header = T); 
pp1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
pp2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
oo<-merge(oo, pp1[,.(gene_id, pval_nominal_threshold)], by="gene_id"); 
names(oo)[ncol(oo)]<-"pval_nominal_threshold_eqtl"
oo<-merge(oo, pp2[,.(gene_id, pval_nominal_threshold)], by.x = "pro_id", by.y ="gene_id"); 
names(oo)[ncol(oo)]<-"pval_nominal_threshold_pqtl"
oo[,p_sum:= -log10(pval_nominal_pqtl) -log10(pval_nominal_eqtl)]

oo1<-oo[slope_pqtl*slope_eqtl<0,]; 
oo1<-oo1[order(p_sum)]; nrow(oo1) 
nrow(oo1[pval_nominal_pqtl <= pval_nominal_threshold_pqtl & pval_nominal_eqtl <= pval_nominal_threshold_eqtl]) 

library(scales)  
ggplot(oo1, aes(x = -log10(pval_nominal_pqtl), y = -log10(pval_nominal_eqtl), color = slope_pqtl)) +
  geom_point(size = 2, alpha = 0.8) + 
  geom_hline(yintercept = max(-log10(oo1$pval_nominal_threshold_eqtl)), 
             linetype = "dashed", 
             color = "black") +
  geom_vline(xintercept = max(-log10(oo1$pval_nominal_threshold_pqtl)), 
             linetype = "dashed", 
             color = "black") + 
  scale_color_gradient2(low = "blue", 
                        mid = "white", 
                        high = "red", 
                        midpoint = 0, 
                        name = "pQTL effect size", 
                        guide = guide_colorbar(title.position = "top", 
                                               barwidth = 8, 
                                               barheight = 1)) +
  labs(x = "-log10(P-value) of pQTL", 
       y = "-log10(P-value) of eQTL", 
       color = "pQTL effect size") + 
  theme_minimal(base_size = 18) +  
  theme(legend.position = "top", 
        legend.direction = "horizontal", 
        legend.title.align = 0.5,  
        legend.title = element_text(size = 14), 
        legend.text  = element_text(size = 14))


## pGenes replicated in eGenes ####
# 'replication' defining as lead pQTL variants in has nominal PV <0.01 with consistent effect direction in eQTL. Vice versa. ref: PMID: 32773033 
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS4_skin_pQTL_allSig_pairs.csv",  header=T);
names(pp)[1]<-"pro_id"
rna_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_RNA_all_qqnorm.txt.gz.allpairs.txt.gz",header = T) # demo
rna_nn1<-merge(rna_nn, unique(name_id_reviewd[,.(pro_id, gene_id)]), by="gene_id", all.x=T, allow.cartesian = TRUE)
pp1<-merge(pp, rna_nn1, by=c("pro_id", "variant_id"))
nrow(pp1); length(unique(pp1$pro_id)); length(unique(pp1$variant_id)); length(unique(pp1[is_LeadpQTL==1]$variant_id)) # 1752; 36; 1752; 35

pp2<-pp1[is_LeadpQTL==1] 
nrow(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]); 
length(unique(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp2[pval_nominal.y<=0.01 & slope.x*slope.y>0]$variant_id)) # 13; 13; 13 
# pGene: 13/35=37.14%;

## eGenes replicated in pGenes ####
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS5_skin_eQTL_allSig_pairs.csv", skip=3, header=T)
pp1<-merge(pp, unique(name_id_reviewd[,.(pro_id, gene_id)]), by="gene_id")
pro_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T); 
names(pro_nn)[1]<-"pro_id"
pp2<-merge(pp1, pro_nn, by=c("pro_id","variant_id"))
nrow(pp2); length(unique(pp2$pro_id)); length(unique(pp2$variant_id)); length(unique(pp2[is_LeadeQTL==1]$variant_id)) # 3549; 42; 3548; 139

pp3<-pp2[is_LeadeQTL==1] 
nrow(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]); 
length(unique(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]$pro_id)); 
length(unique(pp3[pval_nominal.y<=0.01 & slope.x*slope.y>0]$variant_id)) # 47; 11; 47
# p<0.01 eGene: 11/42=26.2%; 

## stacked plot: Fig. S4C ####
aa<-data.table(qtl=c("pQTL","pQTL","eQTL","eQTL"), 
               replication=c("no eQTL","eQTL","no pQTL","pQTL"), 
               num=c(22,13,31,11))
aa[, prop := 100*num / sum(num), by = qtl]
aa$qtl<-factor(aa$qtl,levels = c("pQTL", "eQTL"))
aa$replication<-factor(aa$replication,levels = c("no eQTL","eQTL","no pQTL","pQTL"))

ggplot(aa, aes(x = qtl, y = prop, fill = replication)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#F4EDCA", "#C4961A", "#C3D7A4", "#52854C")) +
  labs(x = '', y = 'Proportion of replicated genes (%)') +
  scale_x_discrete(labels = c("eQTL"="eGenes", "pQTL"="pGenes"))+
  mytheme + 
  theme(legend.position = "none")


## colocolization of eQTL and pQTL: Fig. 4D ####
library(coloc)
peqtl<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.txt.gz", header = T) # demo
peqtl[,id:=paste(pro_id, variant_id, sep="_")]; 
peqtl[,id2:=paste(gene_id, variant_id, sep="_")]

sel_gene<-unique(peqtl$pro_id); 
pph4 <- NULL; 
pph4_snp <- list() 
for (i in 1:length(sel_gene)) { 
  tmp<-peqtl[pro_id==sel_gene[i],]
  tmp<-tmp[order(variant_id, pval_nominal_pqtl, pval_nominal_eqtl)]; 
  tmp<-tmp[!duplicated(variant_id)] 
  my.coloc<-coloc.abf(dataset1=list(beta=tmp$slope_pqtl, 
                                    varbeta=(tmp$slope_se_pqtl)^2, 
                                    snp=tmp$variant_id, 
                                    type='quant', 
                                    N=186, 
                                    MAF=tmp$maf_pqtl
                                    ), 
                      dataset2=list(beta=tmp$slope_eqtl, 
                                    varbeta=(tmp$slope_se_eqtl)^2, 
                                    snp=tmp$variant_id, 
                                    type='quant', 
                                    N=150, 
                                    MAF=tmp$maf_eqtl)
                      ) 
  summary_row <- data.table(pro_id = sel_gene[i], t(my.coloc$summary)); 
  pph4 <- rbind(pph4, summary_row) 
}

fwrite(pph4,"skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt",col.names = T,sep="\t",quote = F)

## plot PPH heatmap
library(ComplexHeatmap); 
library(circlize)
pph4<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.Overlap_eQTL.colo2.pph4_gene.txt", header = T)
pph4<-merge(pph4[,1:7], unique(name_id_reviewd[,c(1,4)]), by="pro_id"); 
names(pph4)<-gsub(".abf", "",names(pph4))
hh2<-rbind(pph4[order(-PP.H1)][PP.H1>=0.7], pph4[order(-PP.H2)][PP.H2>=0.7], pph4[order(-PP.H3)][PP.H3>=0.7], pph4[order(-PP.H4)][PP.H4>=0.7])
tmp<-as.data.frame(hh2); 
rownames(tmp)<-tmp$gene_name; 
mat <- as.matrix(tmp[, 3:7])

Heatmap(mat, 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col=colorRampPalette(c("white","firebrick"))(50),
        row_names_side = "left", 
        column_names_side = "top", 
        border = TRUE, 
        column_names_rot = 45, 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 14, fontface = "bold"), 
                                    labels_gp = gpar(fontsize = 12), 
                                    legend_height = unit(4, "cm"), 
                                    legend_width  = unit(0.6, "cm"), 
                                    title = "PP", 
                                    legend_direction = "vertical"),
        row_names_gp = gpar(fontsize = 14, fontface = "bold.italic"), 
        column_names_gp = gpar(fontsize = 14, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) { 
          grid.rect(x = x, 
                    y = y, 
                    width = width, 
                    height = height,
                    gp = gpar(col = "grey80", fill = NA, lwd = 0.5))}
)


