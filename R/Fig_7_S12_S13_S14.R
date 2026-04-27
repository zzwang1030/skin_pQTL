library(data.table); library(ggplot2); library(scales)
std <- function(x) sd(x,na.rm = T)/sqrt(length(x[!is.na(x)]))
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=26),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=9,face = "italic")
  )

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

## macrophage marker expression: Fig. S12C ####
sel_marker<-data.table(gene_name = c("DCN","COL1A1","COL1A2","VIM",  
                                     "KRT5","KRT10","KRT1","KRT14",  
                                     "CD68","CD163"),
                       celltype = c(rep("Fibroblasts", 4), 
                                    rep("Keratinocyte", 4), 
                                    rep("Macrophage", 2))
                       );

## protein level
qq3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/proteome_matrix_imputed.csv"); 
mm<-as.matrix(qq3); 
rownames(mm)<-mm[,1]; 
mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) 
mm1 <- predict(fit, mm) 
mm2 <- as.data.table(mm1); 
mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); 
names(mm3)<-c("pro_id","sample","intensity")
mm4<-merge(mm3, unique(name_id_reviewd[,c(1,4)]), by = "pro_id")

all_comb <- CJ(sample = unique(mm4$sample), 
               gene_name = sel_marker$gene_name) 
res <- merge(all_comb, mm4, by = c("sample", "gene_name"), all.x = TRUE) 
res <- merge(res, sel_marker, by="gene_name", all.x=T)
res[is.na(intensity), intensity := 0] 

res1<-copy(res) 
res1$gene_name<-factor(res1$gene_name, 
                       levels = c("COL1A1", "VIM", "DCN", "COL1A2", 
                                  "KRT1", "KRT14", "KRT10", "KRT5", 
                                  "CD163", "CD68")
                       )
ggplot(res1, aes(x = gene_name, y = intensity)) +
  geom_boxplot(aes(color = celltype), 
               linewidth = 1, 
               median.linewidth = 2, 
               outlier.shape = NA) +
  geom_jitter( color="grey1", 
               width = 0.25, 
               alpha = 0.3, 
               size = 0.6) +
  scale_color_manual(values = c("Fibroblasts" = "#d95f02", 
                                "Keratinocyte" = "#E09411", 
                                "Macrophage" = "#1b9e77"))+
  labs(x = NULL, y = "Protein level") +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, face = "bold.italic"))



## RNA level
rr_rpkm1 <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/gene_sample_FPKM.csv",header = T)
rr_rpkm2<-melt(rr_rpkm1, id.vars = "gene_id"); 
names(rr_rpkm2)<-c("gene_id","sample","intensity")
rr_rpkm3<-merge(rr_rpkm2, unique(name_id_reviewd[,c(2,4)]), by = "gene_id")

all_comb <- CJ(sample = unique(rr_rpkm3$sample), gene_name = sel_marker$gene_name) 
res <- merge(all_comb, rr_rpkm3, by = c("sample", "gene_name"), all.x = TRUE) 
res <- merge(res, sel_marker, by="gene_name", all.x=T) 
res[is.na(intensity), intensity := 0] 

res1<-copy(res)
res1$gene_name<-factor(res1$gene_name, levels = c("COL1A1", "VIM", "DCN", "COL1A2", 
                                                  "KRT1", "KRT14", "KRT10", "KRT5", 
                                                  "CD163", "CD68")) 
ggplot(res1, aes(x = gene_name, y = log10(intensity+0.01))) +
  geom_boxplot(aes(color = celltype), 
               linewidth = 1, 
               fatten = 2, 
               outlier.shape = NA
               ) +
  geom_jitter( width = 0.25, 
               alpha = 0.3, 
               size = 0.6
               ) +
  scale_color_manual(values = c("Fibroblasts" = "#d95f02", 
                                "Keratinocyte" = "#E09411", 
                                "Macrophage" = "#1b9e77")
                     )+
  labs(x = NULL, y = "log10(FPKM) of RNA level") +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 16, face = "bold"), 
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, face = "bold.italic"))

## ALDH2-KD expression: Fig. 7B ####
ggplot(tmp2, aes(x = genes, y = exp_mean)) +
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           width = 0.6, 
           aes(fill = genes)
  ) +
  scale_fill_manual(values = c("NC" = "gray70", 
                               "ALDH2_siR1" = "#4682B4")
  ) + 
  geom_errorbar(aes(ymin = exp_mean - exp_mean_std, ymax = exp_mean + exp_mean_std),
                width = 0.15, 
                position = position_dodge(0.6), 
                color = "black", 
                linewidth = 0.6) +
  geom_jitter(aes(y = rel_exp), 
              width = 0.2, 
              size = 2, 
              color = "black", 
              alpha = 1)+
  labs(x = NULL, y = "Relative RNA expression") +
  scale_x_discrete(labels = c("NC" = "NC", "ALDH2_siR1" = "ALDH2 KD")) +
  mytheme + 
  theme(legend.position="none",legend.title=element_blank())

## ALDH2-KD cell CFU: Fig. 7C ####
ggplot(tmp2, aes(x = genes, y = CFU_mean/1e5)) + 
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           width = 0.6, 
           aes(fill = genes)
           ) +
  scale_fill_manual(values = c("NC" = "gray70", 
                               "ALDH2" = "#4682B4")
                    ) + 
  geom_errorbar(aes(ymin = (CFU_mean - CFU_mean_std)/1e5, 
                    ymax = (CFU_mean + CFU_mean_std)/1e5), 
                width = 0.15, 
                position = position_dodge(0.6), 
                color = "black", 
                linewidth = 0.6
                ) +
  geom_jitter(aes(y = CFU2/1e5), 
              width = 0.2, 
              size = 2, 
              color = "black", 
              alpha = 1
              )+
  labs(x = NULL, y = expression("CFU counts (×10"^5*")")) +
  scale_x_discrete(labels = c("NC" = "NC", "ALDH2" = "ALDH2 KD")) +
  mytheme + 
  theme(legend.position="none", 
        legend.title=element_blank())

## ALDH2-KO tail lesions: Fig. 7F ####
ggplot(hh, aes(x = Mice, y = 100*lesion_mean)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, aes(fill = Mice)) +
  scale_fill_manual(values = c("WT" = "gray70", "KO" = "#4682B4")) + 
  geom_errorbar(aes(ymin = 100*(lesion_mean - lesion_std), 
                    ymax = 100*(lesion_mean + lesion_std)), 
                width = 0.15, position = position_dodge(0.6), color = "black", linewidth = 0.8) +
  geom_jitter(aes(y = 100*tail_lesion), width = 0.2, size = 3, color = "black", alpha = 1)+
  coord_cartesian(ylim = c(20, 90))+
  labs(x = NULL, y = "Tail lesion (%)") +
  scale_x_discrete(labels = c("WT" = "WT", "KO" = "ALDH2 KO")) +
  mytheme + 
  theme(legend.position="none",legend.title=element_blank())

wilcox.test(hh[Mice=="WT"]$tail_lesion, hh[Mice=="KO"]$tail_lesion)$`p.value` # 0.0079

## ALDH2-KO tail CFU: Fig. 7H ####
ggplot(hh, aes(x = Mice, y = mean100/1e6)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, aes(fill = Mice)) +
  scale_fill_manual(values = c("WT" = "gray70", "KO" = "#4682B4")) + 
  geom_errorbar(aes(ymin = (mean100 - std100)/1e6, ymax = (mean100 + std100)/1e6), 
                width = 0.15, position = position_dodge(0.6), color = "black", linewidth = 0.8) +
  geom_jitter(aes(y = dilution100/1e6), width = 0.2, size = 3, color = "black", alpha = 1)+
  labs(x = NULL, y = expression("CFU counts (×10"^6*")")) +
  scale_x_discrete(labels = c("WT" = "WT", "KO" = "ALDH2 KO")) +
  mytheme+ theme(legend.position="none",legend.title=element_blank())

wilcox.test(hh[Mice=="WT"]$dilution100, hh[Mice=="KO"]$dilution1000)$`p.value` 

## ALDH2-KD engulfed cells: Fig. S13B ####
ggplot(tmp2, aes(x = genes, y = percent_mean)) +
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           width = 0.6, 
           aes(fill = genes)
           ) +
  scale_fill_manual(values = c("NC" = "gray70", 
                               "ALDH2" = "#4682B4")
                    ) + 
  geom_errorbar(aes(ymin = percent_mean - percent_mean_std, 
                    ymax = percent_mean + percent_mean_std),
                width = 0.15, 
                position = position_dodge(0.6), 
                color = "black", 
                linewidth = 0.6) +
  geom_jitter(aes(y = Statistic), 
              width = 0.2, 
              size = 2, 
              color = "black", 
              alpha = 1) +
  labs(x = NULL, y = "Percent of fluorescent cells") +
  scale_x_discrete(labels = c("NC" = "NC", "ALDH2" = "ALDH2 KD")) +
  mytheme + 
  theme(legend.position="none", 
        legend.title=element_blank())

## ALDH2-KO spleen mass: Fig. S13E ####
ggplot(hh, aes(x = Mice, y = spleen_mean)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, aes(fill = Mice)) +
  scale_fill_manual(values = c("WT" = "gray70", "KO" = "#4682B4")) + 
  geom_errorbar(aes(ymin = spleen_mean - spleen_std, 
                    ymax = spleen_mean + spleen_std), 
                width = 0.15, position = position_dodge(0.6), color = "black", linewidth = 0.8) +
  geom_jitter(aes(y = spleen_mass2), width = 0.2, size = 3, color = "black", alpha = 1)+
  coord_cartesian(ylim = c(0.3, 0.6)) +
  labs(x = NULL, y = "Spleen mass (g)") +
  scale_x_discrete(labels = c("WT" = "WT", "KO" = "ALDH2 KO")) +
  mytheme + 
  theme(legend.position="none",legend.title=element_blank())

## NAA25-KD tail lesion: Fig. S14C ####
ggplot(mydf, aes(x = group, y = value, fill = group)) +
  stat_summary(fun = mean, 
               geom = "bar", 
               width = 0.6, 
               color = NA) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               width = 0.15, 
               linewidth = 0.6) +
  geom_jitter(width = 0.08, 
              linewidth = 2) +
  scale_fill_manual(values = c("shCtrl" = "#BDBDBD", 
                               "shNAA25" = "#4C78A8")) +  
  labs(x = NULL, y = "Tail lesion (%)") + 
  coord_cartesian(ylim = c(0, 95)) +
  mytheme + 
  theme(legend.position = "none")

## NAA25-KD tail CFU: Fig. S14E ####
ggplot(mydf, aes(x = group, y = value, fill = group)) +
  stat_summary(fun = mean, 
               geom = "bar", 
               width = 0.6, 
               color = NA) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar", 
               width = 0.15, 
               size = 0.6) +
  geom_jitter(width = 0.08, 
              size = 2) +
  scale_fill_manual(values = c("shCtrl" = "#BDBDBD", 
                               "shNAA25" = "#4C78A8")) +  
  labs(x = NULL, y = expression(CFU~(x10^6))) + 
  mytheme + 
  theme(legend.position = "none")

