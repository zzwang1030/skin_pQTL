library(data.table); library(ggplot2)
mytheme <- theme_classic(base_size = 24) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5)
  )

name_id_reviewd<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)

tss3 <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/ensembl_uniprot_id_pos_biomart.csv", header = T)

#### pQTL calling ####
# For the pQTL calling pipelines, see the shell script in pQTL_eQTL_calling.docx.
#### pQTL charactering ####
## manhattan plot: Fig. 3A ####
library(ggmanh); 
library(SeqArray)
nn3<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz", header = T) # demo 
pp2<- tidyr::extract(nn3, col = 'variant_id', into = c('chr', 'pos','ref','alt'),regex = '(.+):(.+):(.+):(.+)', remove=F)
pp2<-as.data.table(pp2)
pp2$chr<-factor(pp2$chr,levels = 1:22); 
pp2$pos<-as.numeric(pp2$pos)
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header = T)
pt<-max(pp[qval<=0.05]$pval_nominal)

pp2<-merge(pp2, unique(name_id_reviewd[,.(pro_id, gene_name)]), by.x="gene_id", by.y="pro_id")
pp2<-pp2[order(gene_id, pval_nominal)]
pp2[, p_order:=1:.N,by=.(gene_id)]
pp2[, lable_name:=ifelse(p_order==1, gene_name, ""), by=.(gene_id)] 
pp2[, lable_name2:=ifelse(pval_nominal<=pt, lable_name, "")] 

pp2[,color_class:=ifelse(pval_nominal<=pt, "sig", "ns")];
highlight_colormap <- c("ns" = "grey90", "sig" = "indianred")
manhattan_plot(pp2, 
               pval.colname = "pval_nominal", 
               chr.colname = "chr", 
               pos.colname = "pos",
               y.label = "-log10(P)", 
               thin=T, 
               signif = pt, 
               point.size = 0.4, 
               label.font.size=4,
               label.colname="lable_name2", 
               highlight.colname = "color_class", 
               highlight.col = highlight_colormap, 
               color.by.highlight = TRUE) + 
  mytheme + 
  theme(legend.position="none") # Some labels couldn't be drawn in time, so I added some manually, but some were still omitted.

## proportion variance explained (PVE) by lead pQTL: Fig. 3C ####
SNP_PVE <- function(beta, MAF, se_beta, N) { 
  numerator <- 2 * (beta^2) * MAF * (1 - MAF)
  denominator <- numerator + ((se_beta^2) * 2 * N * MAF * (1 - MAF))
  SNP_PVE <- numerator / denominator ; return(SNP_PVE) }

sample_N<-186
hh1<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.SigPairs_pt.txt", header = T)
hh1[,pve_snp:=SNP_PVE(slope, maf, slope_se, sample_N)]

hh2<-hh1[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[!duplicated(protein_id)]
summary(hh3$pve_snp)  # 0.1205  0.1433  0.1582  0.1748  0.1936  0.3257 
ggplot(hh3, aes(x=100*pve_snp)) + 
  geom_histogram(bins=20, fill="steelblue", color='grey90', linewidth=0.5) +
  stat_bin(bins = 20, geom = "text", aes(label = after_stat(count)), vjust = -0.5, size = 6)+
  labs(x ="Proportion of variance explained (%) by lead pQTL", y="Numbers of lead pQTL") + 
  mytheme

## effect siz vs MAF: Fig. 3B ####
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS4_skin_pQTL_allSig_pairs.csv", header = T) 
hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[is_LeadpQTL==1]
cor.test(hh3$maf, abs(hh3$slope), method = "spearman") # P=1.205e-05, rho = -0.6597142

ggplot(hh3, aes(x = maf, y = abs(slope))) + 
  geom_point(size=2, color="steelblue") + 
  geom_smooth(method='lm', se=T, span = 1, formula= y~x, col="indianred", linetype="dashed") + 
  labs(x = 'Minor allele frequency', y = 'Effect size (absolute values)') +
  coord_cartesian(ylim = c(0, 2)) + 
  mytheme

## effect siz vs distance to TSS: Fig. 3D ####
hh2<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS4_skin_pQTL_allSig_pairs.csv", header = T) 
hh2<-hh2[order(protein_id, pval_nominal, -slope)]
hh3<-hh2[is_LeadpQTL==1] 
ggplot(hh3, aes(x = tss_distance/1e6, y = abs(slope))) + 
  geom_point(size=3, color="steelblue") + 
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  labs(x = 'Distance to TSS (MB)', y = 'Effect size (absolute values)') +
  mytheme

## feature annotation by VEP: Fig. 3E ####
# lead pQTL for pGene
ss7 <- readxl::read_excel("/Users/sunyuanqian/Library/CloudStorage/OneDrive-个人/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/annotation_VEP.xls", sheet = 1)
ss7<-as.data.table(ss7)
bb<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS4_skin_pQTL_allSig_pairs.csv", header = T)
bb<-bb[order(protein_id, pval_nominal, -slope)]; 
bb[,id:=paste(protein_id, variant_id, sep="_")]
bb2<-bb[!duplicated(protein_id)] 
ss8<-ss7[id %in% bb2$id];
tmp1<- as.data.table(table(ss8$Consequence)); 
names(tmp1)<-c("effect","cnt"); 
tmp2<-tmp1
tmp2[,ratio:=round(cnt*100 / sum(cnt), 1)]; 
tmp2[, effect2:=gsub(" variant","",gsub("_", "", effect))]

# lead SNP of all proteins
ss7 <- readxl::read_excel("/Users/sunyuanqian/Library/CloudStorage/OneDrive-个人/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/annotation_VEP.xls", sheet = 2)
ss7<-as.data.table(ss7)
tmp_bg1<- as.data.table(table(ss7$Consequence)); 
names(tmp_bg1)<-c("effect","cnt") 
tmp_bg2<-tmp_bg1[effect!="splice_acceptor_variant" &  effect!="splice_polypyrimidine_tract_variant"]; 
tmp_bg2[effect=="splice_region_variant"]$cnt<-3 # For simplicity, combine several splice-related terms.
tmp_bg2[,ratio:=round(cnt*100 / sum(cnt), 1)]; 
tmp_bg2[, effect2:=gsub(" variant","",gsub("_", "", effect))]

tmp3<-rbind(tmp2[,-4], tmp_bg2[,-4]); 
tmp3$Type<-c(rep("pGenes", 5), rep("Other genes", 8))
tmp3$effect<-factor(tmp3$effect, levels = rev(c("upstream_gene_variant", "downstream_gene_variant", "intron_variant", "missense_variant", "3_prime_UTR_variant", "splice_region_variant", "synonymous_variant", "5_prime_UTR_variant")))
library(viridis)
ggplot(tmp3, aes(x = Type, y = ratio, fill = effect)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = NULL, y = "Proportion of functional region", fill = "Effect") +
  scale_fill_viridis_d(option = "F", name = "Effect")+
  mytheme + 
  theme(axis.text.x = element_text(angle = 45, face="bold", hjust = 1))


## RegulomeDB_rank of pQTL: Fig. 3F ####
regulome<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/RegulomeDB_rank.tsv", header = T) # demo; full from https://regulomedb.org/regulome-search/
# pGene
pro<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/TableS4_skin_pQTL_allSig_pairs.csv", header = T) ; 
pro[, chr_var:=as.character(chr_var)]
pro1<-merge(pro, regulome[, c(1, 3:6)], by=c("chr_var","pos_var")); 
nrow(pro1) # 1764
pro1[,ranking2:= as.integer(sub("^(\\d).*", "\\1", ranking))]
nrow(pro1[ranking2<=2]);  # 1386/1764 = 78.6% 
nrow(pro1[ranking2>2]) # 378/1764 = 21.4% 

# Backgroud gene
pro_nn<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz",header = T)
pro_nn <- tidyr::extract(pro_nn, col = 'variant_id', into = c('chr_var', 'pos_var','ref','alt'), regex = '(.+):(.+):(.+):(.+)', remove=F); 
pro_nn<-as.data.table(pro_nn)
pro_nn<-pro_nn[order(variant_id, pval_nominal)]; 
pro_nn1<-pro_nn[!duplicated(variant_id)]
pro_nn1[, pos_var := as.integer(pos_var)]
pro_nn2<-merge(pro_nn1, regulome[, -2], by=c("chr_var","pos_var")); 
nrow(pro_nn2) # 3,455,459
pro_nn2[,ranking2:= as.integer(sub("^(\\d).*", "\\1", ranking))]
nrow(pro_nn2[ranking2<=2]);   # 1,515,102/3455459 = 43.85%
nrow(pro_nn2[ranking2>2]) # 1,940,357/3455459 = 56.15%

chisq.test(matrix(c(1386, 378, 1515102, 1940357), byrow = T, nrow = 2))$`p.value` #  p=1.808773e-189

# stack plot
aa<-data.table(qtl=c("pQTL","pQTL","bg","bg"), 
               replication=c("Regulome","No regulome","Regulome","No regulome"), 
               num=c(1386, 378, 1515102, 1940357)) 
aa[, prop := 100*num / sum(num), by = qtl]
aa$qtl <- factor(aa$qtl,levels = c( "bg", "pQTL"))
aa$replication <- factor(aa$replication,levels = c("Regulome","No regulome"))
aa[, fill_group := paste(qtl, replication, sep = "_")]

ggplot(aa, aes(x = qtl, y = prop, fill = fill_group)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#F4EDCA", "#C4961A", "#C3D7A4", "#52854C")) +
  labs(x = '', y = 'Proportion of functional variants (%)') +
  scale_x_discrete(labels = c("bg"="Background", "pQTL"="pQTLs"))+
  coord_flip() +
  mytheme + 
  theme(legend.position = "none")

#### pQTL fine-mapping by susieR: Fig. S6A ####
## 计算 ld R ####
# # Generate a list of snp for each gene for susieR to calculate ld R.
nn3 <- fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz")
nn3 <- nn3[ !is.na(gene_id) & !is.na(variant_id) & !is.na(slope) & !is.na(slope_se) & slope_se > 0] # 基础过滤

pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header=T)
length(unique(pp$gene_id)); length(unique(pp[qval<0.05,]$gene_id)) # 5152; 36
pgenes<-unique(pp[qval<0.05,]$gene_id)

dir.create("/sh2/home/sunyuanqiang/projects/pQTL/formal/susie_variant_lists", showWarnings = FALSE)
for (g in unique(nn3$gene_id)) {
  snps <- unique(nn3[gene_id == g, .(variant_id)])
  fwrite(snps, file = file.path("/sh2/home/sunyuanqiang/projects/pQTL/formal/UsingDataDuration/susie_variant_lists", paste0(g, ".snplist")), sep = "\t",  col.names = FALSE)
} 

# Calculate the LD correlation matrix for all pGenes in batches in shell
# mkdir -p susie_ld; for f in susie_variant_lists/*.snplist; do g=$(basename "$f" .snplist); echo "Processing $g"; plink --bfile /sh2/home/sunyuanqiang/projects/pQTL/formal/hg38_ASA_NH_196samples_474857SNPs_FWD_chr1_22_exclude_imputation_r8m5.newHeader.plink --extract "$f" --write-snplist --r square --out "cd/${g}"; done

# Add SNP names to the LDR matrix, perform quality checks
ld_dir <- "susie_ld"; 
snp_dir <- "susie_variant_lists"; 
out_dir <- "susie_R_tables"
dir.create(out_dir, showWarnings = FALSE)
ld_files <- list.files(ld_dir, pattern = "\\.ld$", full.names = TRUE)
for (f in ld_files) {
  g <- sub("\\.ld$", "", basename(f))
  snp_file <- file.path(snp_dir, paste0(g, ".snplist"))
  
  if (!file.exists(snp_file)) {
    message("Missing snplist for ", g)
    next
  }
  
  snps <- fread(snp_file, header = FALSE)$V1   # Read SNP order
  R <- as.matrix(read.table(f))   # read LD matrix
  
  if (nrow(R) != length(snps) || ncol(R) != length(snps)) {
    message("Dimension mismatch for ", g)
    next
  }   # Dimensional Inspection
  
  rownames(R) <- snps; colnames(R) <- snps   # Assign row and column names
  
  if (!isTRUE(all.equal(R, t(R)))) {
    stop("R matrix not symmetric for ", g)
  }   # Symmetry check
  
  if (!all(abs(diag(R) - 1) < 1e-6)) {
    stop("Diagonal not equal to 1 for ", g)
  }  # Diagonal check
  
  R_dt <- as.data.table(R, keep.rownames = "variant_id")
  fwrite(R_dt, file = file.path(out_dir, paste0(g, ".R.txt")), sep = "\t", quote = FALSE, na = "NA")
  message("Finished ", g)
}

## 运行susieR ####
sumstats_file <- "~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.allpairs.txt.gz"
R_dir <- "susie_R_tables"; 
out_dir <- "susie_results"
n_sample <- 186; 
cs_coverage <- 0.95; 
cs_min_abs_corr <- 0.5 
L_max <- 1 

ss_all <- fread(sumstats_file)
ss_all <- ss_all[!is.na(gene_id) & !is.na(variant_id) & !is.na(slope) & !is.na(slope_se) & slope_se > 0 ]
ss_all <- unique(ss_all, by = c("gene_id", "variant_id"))
ss_all[, z := slope / slope_se]

## Define a function: Read an R.txt file and restore it as a matrix
read_R_txt <- function(R_file) {
  R_dt <- fread(R_file)
  if (ncol(R_dt) < 3) { stop("R table has too few columns: ", R_file) }
  
  snps <- R_dt[[1]] 
  R_mat <- as.matrix(R_dt[, -1, with = FALSE])
  mode(R_mat) <- "numeric"
  rownames(R_mat) <- snps; colnames(R_mat) <- colnames(R_dt)[-1]
  
  if (!identical(rownames(R_mat), colnames(R_mat))) { 
    stop("Row names and column names of R do not match: ", R_file)
  }
  
  if (!isTRUE(all.equal(R_mat, t(R_mat)))) { 
    stop("R matrix is not symmetric: ", R_file)
  }
  
  if (!all(abs(diag(R_mat) - 1) < 1e-6)) { 
    stop("Diagonal of R is not 1: ", R_file)
    }
  
  return(R_mat)
}

## Define a function: run susie_rss on a single gene
run_susie_one_gene <- function(gene_name, 
                               ss_all, 
                               R_file, 
                               n_sample = 186, 
                               L_max = 1, 
                               cs_coverage = 0.95, 
                               cs_min_abs_corr = 0.5) {
  R <- read_R_txt(R_file)
  snps_R <- colnames(R)
  
  # Extract summary stats for this gene
  ss <- ss_all[gene_id == gene_name]
  if (nrow(ss) == 0) {
    stop("No summary stats found for gene: ", gene_name)
    }
  
  # Align according to the SNP order of R
  ss <- ss[match(snps_R, variant_id)]
  
  if (any(is.na(ss$variant_id))) {
    missing_snps <- snps_R[is.na(ss$variant_id)]
    stop("Some SNPs in R are missing from summary stats for ", gene_name, ": ", paste(head(missing_snps, 5), collapse = ", "))
  }
  
  if (!identical(ss$variant_id, snps_R)) { 
    stop("SNP order mismatch after matching for gene: ", gene_name) 
    }
  
  # perform susie_rss
  fit <- susie_rss(
    z = ss$z,
    R = R,
    n = n_sample,
    L = L_max,
    coverage = cs_coverage,
    min_abs_corr = cs_min_abs_corr,
    estimate_prior_variance = FALSE, 
    estimate_residual_variance = FALSE, 
    check_R = TRUE, check_z = FALSE, max_iter = 1000
  )
  
  # SNP-level result
  snp_res <- copy(ss)
  snp_res[, pip := fit$pip]
  snp_res[, cs_id := NA_integer_]
  
  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    for (k in seq_along(fit$sets$cs)) {
      idx <- fit$sets$cs[[k]]
      snp_res[idx, cs_id := k]
    }
  }
  
  # top PIP SNP
  top_idx <- which.max(snp_res$pip); 
  lead_idx <- which.min(snp_res$pval_nominal)
  
  gene_summary <- data.table(
    gene_id = gene_name, 
    n_snps = nrow(snp_res),
    n_cs = ifelse(is.null(fit$sets$cs), 0, length(fit$sets$cs)),
    top_pip_variant = snp_res$variant_id[top_idx],
    top_pip = snp_res$pip[top_idx],
    top_pip_pval = snp_res$pval_nominal[top_idx],
    lead_variant = snp_res$variant_id[lead_idx],
    lead_pval = snp_res$pval_nominal[lead_idx],
    top_equals_lead = snp_res$variant_id[top_idx] == snp_res$variant_id[lead_idx],
    max_pip = max(snp_res$pip, na.rm = TRUE),
    any_pip_gt_095 = any(snp_res$pip > 0.95, na.rm = TRUE)
  )
  
  list( fit = fit, snp_res = snp_res, gene_summary = gene_summary)
}

## Batch Execution
R_files <- list.files(R_dir, pattern = "\\.R\\.txt$", full.names = TRUE)
all_gene_summaries <- list(); 
all_snp_results <- list()

for (R_file in R_files) {
  gene_name <- sub("\\.R\\.txt$", "", basename(R_file))
  message("Running SuSiE for ", gene_name)
  
  res <- tryCatch(
    run_susie_one_gene(
      gene_name = gene_name,
      ss_all = ss_all,
      R_file = R_file,
      n_sample = n_sample,
      L_max = L_max, 
      cs_coverage = cs_coverage,
      cs_min_abs_corr = cs_min_abs_corr
    ),
    error = function(e) {
      message("Failed for ", gene_name, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(res)) next
  
  all_gene_summaries[[gene_name]] <- res$gene_summary
  all_snp_results[[gene_name]] <- res$snp_res
  fwrite(res$snp_res, file = file.path(out_dir, paste0(gene_name, ".snp_results.txt")), sep = "\t") 
} 

## Summary Table
gene_summary_table <- rbindlist(all_gene_summaries, use.names = TRUE, fill = TRUE)
snp_result_table <- rbindlist(all_snp_results, use.names = TRUE, fill = TRUE)
fwrite(snp_result_table, file = file.path(out_dir, "SuSiE_all_snp_results.txt"), sep = "\t")

## plot
ss<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/SuSiE_all_snp_results.txt.gz")
pp<-fread("~/Desktop/省皮/project/pQTL/manuscript/pQTL_MS_20251106_NatComm/revision1_202512/code/skin_protein_all_qqnorm.txt.gz.permu.genes.txt.gz",header = T) 
ss1<-merge(ss, unique(pp[,.(gene_id, pval_nominal_threshold)]), by="gene_id", all.x=T)
ss1<-merge(ss1, unique(name_id_reviewd[,.(pro_id, gene_name)]), by.x="gene_id", by.y="pro_id", all.x=T)
ss1$gene_id <- factor(ss1$gene_id, levels = unique(ss1$gene_id))
lab_map <- unique(ss1[, c("gene_id", "gene_name")]); 
lab_map <- setNames(lab_map$gene_name, lab_map$gene_id)
thres_df <- unique(ss1[, c("gene_id", "pval_nominal_threshold")]) 

ggplot(ss1, aes(x = -log10(pval_nominal), y = pip)) +
  geom_point(color = "grey8", alpha = 0.2, size = 0.5) + 
  geom_point(data = subset(ss1, !is.na(cs_id)), color = "indianred", alpha = 0.6, size = 1.2) +  
  geom_point( data = subset(ss1, pip > 0.9), shape = 8, size = 3, color = "indianred") +   
  geom_vline(data = thres_df, aes(xintercept = -log10(pval_nominal_threshold)), linetype = "dashed", linewidth = 0.3) +  
  geom_hline(yintercept = 0.9, linetype = "dashed", linewidth = 0.3) +   # PIP threshold
  facet_wrap( ~gene_id, nrow = 6, scales = "free_x", labeller = labeller(gene_id = lab_map) ) +
  labs( x = "-log10(P-value)", y = "PIP") + coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  theme( strip.text = element_text(face = "bold.italic", size = 12),
         axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
         axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12) )

#### conditional QTL analysis by GCTA COJO ####
# see the shell script in pQTL_eQTL_calling.docx and run_cojo_by_protein.sh.

#### trans-pQTL identification ####
# see the shell script in pQTL_eQTL_calling.docx.




