library(data.table); library(ggplot2)
name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
age_paper<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/Cell_proteome_aging/mmc2.xlsx", sheet = 2, skip = 2); age_paper<-as.data.table(age_paper) 
table(age_paper$Tissue)
mytheme <- theme_classic(base_size = 26) + theme(
  axis.text = element_text(color = 'black',size=20),
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  plot.subtitle = element_text(hjust = 0.5),
  legend.text=element_text(size=16, face = "bold") #设置图例字体的大小及字体https://ggplot2.tidyverse.org/articles/ggplot2-specs.html
  #,axis.text.x = element_text(angle = 45,vjust = 0.95, hjust=1,color="black"),
  #,axis.line = element_line(size = 1.2),axis.ticks=element_line(size=1.2),#坐标轴及 ticks 的粗度
)

#### protein and age ####
## proteome process ####
# 1/2 minimum impution, VSN normalizaiton, not lmFit, 同一年龄不同样本取平均
qq2<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename.csv"); qq2<-qq2[,-2]
ff<-fread("~/Desktop/省皮/project/pQTL/M-GSGC0290947数据-剔除离群值-20230822/组织剔除离群值/20230514_153410_M-GSGC0290947_TISSUE_Report(不填补空值)_去除离群值_pro_rename_filter.csv")
cols_to_keep <- colnames(qq2)[2:ncol(qq2)] %in% colnames(ff); selected_cols <- c(TRUE, cols_to_keep)
rows_to_keep <- qq2$PG.ProteinAccessions %in% ff$protein
qq3 <- qq2[rows_to_keep, ..selected_cols]; colnames(qq3)<-gsub("SP", "", colnames(qq3))
qq4<- tidyr::separate_rows(qq3, PG.ProteinAccessions, sep = ";"); qq4<-as.data.table(qq4)
qq5<-melt.data.table(qq4, id.vars = "PG.ProteinAccessions"); names(qq5)<-c("pro_id","sample","intensity")
qq5[, min:=min(intensity,na.rm = T), by=.(pro_id)]; qq5[,intensity:=ifelse(is.na(intensity), min/2, intensity)] # 填补 NA

mm<-dcast(qq5,formula = pro_id~sample, value.var = "intensity" )
mm<-as.matrix(mm); rownames(mm)<-mm[,1]; mm<-mm[,-1]
library(vsn)
fit <- vsn2(mm) # 建立 vsn 模型并归一化，不能提前 log2 转化
meanSdPlot(fit) # 标准差(sd)在不同均值区间(rank mean)近乎恒定，说明拟合很好
mm1 <- predict(fit, mm) # 已经是log2 转化了
mm2 <- as.data.table(mm1); mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); names(mm3)<-c("pro_id","sample","intensity")

## protein level vs age: rho ####
hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_省立正常对照_外送20200623.xls", sheet = 1) #里面含有采样时间信息
hh2<-as.data.table(hh2); names(hh2)[1:2]<-c("sample","date"); hh2<-hh2[,1:2]
hh2[,date2:=as.Date(date, format = "%Y-%m-%d")]; hh2[,date_duration:=as.numeric(date2 - min(date2, na.rm = TRUE))]
hh3<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS1_261Samples_info.csv", header = T)
hh3<-merge(hh3, hh2[,c(1,4)], by.x="No", by.y="sample", all.x=T); names(hh3)[3]<-"age"; hh3$age<-as.numeric(hh3$age)

jj<-merge(mm3, hh3, by.x="sample", by.y="No"); jj<-jj[order(pro_id, age)]
jj[, intensity2:=mean(intensity), by=.(pro_id, age)] 
jj1<-unique(jj[,.(pro_id, intensity2, age)])
jj1[, c("rho_age", "rho_p_age") := .(cor.test(intensity2, age, method = "spearman")$estimate, cor.test(intensity2, age, method = "spearman")$p.value), by = .(pro_id)]
fwrite(jj1, "~/Desktop/省皮/project/pQTL/manuscript/Protein_ageMean_cor.csv", col.names = T, row.names = F, sep=",", quote = F)
jj2<-unique(jj1[,.(pro_id, rho_age, rho_p_age)])
jj2<-merge(jj2, unique(name_id_reviewd[,c(1,4)]), by="pro_id")
fwrite(jj2[rho_p_age<0.05, c(4,2,3)], "~/Desktop/省皮/project/pQTL/manuscript/Table SX. Age-correlated_proteins.csv", col.names = T, row.names = F, sep=",", quote = F)
nrow(jj2); nrow(jj2[rho_p_age<0.05]); nrow(jj2[rho_p_age<0.05 & rho_age>0]); nrow(jj2[rho_p_age<0.05 & rho_age<0]) # 5956; 891; 461; 430

## protein level vs age: heatmap ####
library(ComplexHeatmap); library(circlize)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/Protein_ageMean_cor.csv", header = T)
jj_mat <- dcast(jj1, pro_id ~ age, value.var = "intensity2")
mat <- as.matrix(jj_mat[, -1]); rownames(mat) <- jj_mat$pro_id
mat_scaled <- t(scale(t(mat))) # z-score 标准化每一行（每个蛋白）
age_order <- sort(as.numeric(colnames(mat_scaled))); mat_scaled <- mat_scaled[, as.character(age_order)]; colnames(mat_scaled) # 列按照 age 大小排列

jj_pos<-jj2[rho_p_age<0.05 & rho_age>0][order(rho_p_age, -rho_age)] # 461
jj_neg<-jj2[rho_p_age<0.05 & rho_age<0][order(rho_p_age, rho_age)] # 430
jj_sig<-rbind(jj_pos, jj_neg)
tmp <- mat_scaled[match(jj_pos$pro_id, rownames(mat_scaled)), ] # 取出包含的行，并按顺序排列
tmp <- mat_scaled[match(jj_sig$pro_id, rownames(mat_scaled)), ] # 取出包含的行，并按顺序排列

age_labels <- rep("", ncol(tmp)); idx_to_label <- seq(1, ncol(tmp), by = 5)
age_labels[idx_to_label] <-  colnames(tmp)[idx_to_label]; age_labels # 每隔 10 个列标注一次
ht<-Heatmap(tmp, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, 
        column_title = "Age", column_title_side = "bottom", row_title = "Protein level",
        column_names_rot = 45, column_names_gp = gpar(fontface = "bold"), column_labels = age_labels, heatmap_legend_param = list(direction = "horizontal"))
draw(ht, heatmap_legend_side = "top")

## protein variance vs age ####
## 先试试同一年龄组不同样本取平均值的数据。结果挺好，就用这个
# CV 和 SD 结果类似，因为取了 log，所以展示SD
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/Protein_ageMean_cor.csv", header = T)
summary(jj1[pro_id=="A0A075B6H7"]$age) # 13.00   34.50   52.00   52.08   69.50   95.00 
jj1[, age_group := cut(age, breaks = c(13, 34.5, 52, 69.5, 95), labels = c("age1", "age2", "age3", "age4"), include.lowest = TRUE, right = TRUE)]
table(jj1[pro_id=="A0A075B6H7"]$age_group) # age1 age2 age3 age4  = 18   18   17   18 
jj1_cv <- jj1[, .(SD = sd(intensity2, na.rm = TRUE), CV = sd(intensity2, na.rm = TRUE) / mean(intensity2, na.rm = TRUE)), by = .(pro_id, age_group)]

tmp<-jj1_cv[age_group=="age1"|age_group=="age4"]
tmp$age_group<-factor(tmp$age_group,levels = c("age1","age4"))

boxplot(jj1_cv[age_group=="age1"]$SD, jj1_cv[age_group=="age2"]$SD, jj1_cv[age_group=="age3"]$SD, jj1_cv[age_group=="age4"]$SD, outline=F)
summary(jj1_cv[age_group=="age1"]$SD) # 0.06636 0.21710 0.33004 0.38470 0.49590 3.03636
summary(jj1_cv[age_group=="age4"]$SD) # 0.08077 0.25300 0.38217 0.44697 0.56839 2.78095
wilcox.test(jj1_cv[age_group=="age1"]$SD, jj1_cv[age_group=="age4"]$SD, paired = T)$`p.value` # 1.221947e-257
ggplot(tmp, aes(x=age_group, y=SD))+
  geom_line(aes(group = pro_id), color = "grey60", linewidth = 0.2, alpha = 0.2) +
  geom_jitter(aes(color = age_group, group = age_group), width = 0.25, size = 0.8, alpha = 0.2) +
  geom_boxplot(aes(group = age_group), width = 0.4,  outlier.shape = NA, fill = NA, color = "black", lwd = 1.1) +
  scale_color_manual(values = c("indianred", "steelblue"))+ # scale_fill_manual(values = alpha(c("#E04643", "#686F76"), 0.8))+
  labs(x = ' ', y = 'SD of protein level') + coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(labels = c("age1" = "Young", "age4" = "Old"))+
  mytheme + theme(legend.position="none")

## protein predict age model ####
## 测试了如下条件：同一年龄的样本组合并与否；intensity zscore与否; 蛋白纳入标准age_rho_p 1~1e-13; mlr3verse超参搜索与否
## 结论：同一年龄的样本组合并>不合并; zs>log2intensity; 蛋白纳入标准age_rho_p<1e-9, 纳入 8 个蛋白最好；mlr3verse超参搜索>固定alpha = 0.5，但在限定p<1e-9时没啥提升，并且耗时
library(glmnet); library(iml)
paper_model<-readxl::read_excel("~/Desktop/省皮/project/pQTL/plasma_pQTL/Cell_proteome_aging/mmc3.xlsx", sheet = 1, skip = 2); paper_model<-as.data.table(paper_model) 
length(unique(paper_model[tissue=="Skin"]$feature)) # 纳入了 56 个蛋白

library(mlr3verse); library(mlr3tuning); library(paradox); library(mlr3learners)
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/Protein_ageMean_cor.csv", header = T)
expr_mat <- dcast(jj1[rho_p_age<1e-9], age ~ pro_id, value.var = "intensity2_zs")  # 用 intensity2/ intensity2_zs
expr_mat <- na.omit(expr_mat)  # 去除含NA样本

# 0. LOOCV 和超参数调优的建模部分
y <- expr_mat$age; X <- as.data.frame(expr_mat[, !c("age")]) # 拆分 X 和 y
n <- nrow(X); pred_age <- rep(NA, n)

for (i in 1:n) { # ~15 min
  # 1. LOOCV划分
  train_idx <- setdiff(1:n, i); X_train <- X[train_idx, ]; y_train <- y[train_idx]; X_test <- X[i, , drop = FALSE]
  # 2. 构造 Task
  train_data <- data.table::data.table(X_train); train_data$age <- y_train
  task <- TaskRegr$new("age_prediction", backend = train_data, target = "age")
  # 3. 定义 learner (elastic net via glmnet)
  learner <- lrn("regr.glmnet", predict_type = "response")
  # 4. 搜索空间
  param_space <- ps(
    alpha = p_dbl(lower = 0, upper = 1),
    lambda = p_dbl(lower = -4, upper = 1, trafo = function(x) 10^x)) # 实际是1e-4~10,log-scale 搜索 lambda	采样更均匀、搜索更有效
  # 5. Tuning
  instance <- TuningInstanceBatchSingleCrit$new(task = task, learner = learner, resampling = rsmp("cv", folds = 5), 
                                                measure = msr("regr.mae"), search_space = param_space, terminator = trm("evals", n_evals = 50),  # 尝试 30 个组合
                                                )
  tuner <- tnr("random_search"); tuner$optimize(instance)
  # 6. 最佳参数建模
  learner$param_set$values <- instance$result_learner_param_vals
  learner$train(task)
  # 7. 预测
  pred_age[i] <- learner$predict_newdata(X_test)$response
}

# 这个不作为模型性能的最终数据。因为下述全体数据训练最终模型的性能数据更好
mae <- mean(abs(pred_age - y)); mae # 7.63
cor.test(y, pred_age, method = "pearson") # R=0.882, P=3.49e-24

# 1. 用全体数据训练最终模型
# 全数据构建 task
final_data <- data.table::data.table(X); final_data$age <- y
task_final <- TaskRegr$new("age_prediction", backend = final_data, target = "age")
# learner 和搜索空间
learner_final <- lrn("regr.glmnet", predict_type = "response")
param_space <- ps(
  alpha = p_dbl(lower = 0, upper = 1),
  lambda = p_dbl(lower = -4, upper = 1, trafo = function(x) 10^x)  # log-scale lambda
)
# 调参实例
instance_final <- TuningInstanceBatchSingleCrit$new(task = task_final, learner = learner_final, resampling = rsmp("cv", folds = 5), 
                                                    measure = msr("regr.mae"), search_space = param_space, terminator = trm("evals", n_evals = 30))
# 调参并设置最优参数
tuner <- tnr("random_search")
tuner$optimize(instance_final)
learner_final$param_set$values <- instance_final$result_learner_param_vals
# 训练最终模型
learner_final$train(task_final)
# 模型性能
pred_all <- learner_final$predict(task_final) 
pred_age_all <- pred_all$response
mae_all <- mean(abs(pred_age_all - y)) # MAE= 6.86
cor.test(pred_age_all, y, method = "pearson") # R=0.91; P=1.50e-27

tmp<-as.data.table(cbind(y, pred_age_all)); names(tmp)<-c("actual", "predict")
ggplot(tmp, aes(x = actual, y = predict)) + geom_point(color = "grey30", size = 1.5, alpha = 1) + geom_smooth(method = "lm", se = TRUE, color = "indianred",  size = 1.2) + 
  labs(x = "Actual age", y = "Predicted age") + coord_cartesian(xlim=c(13,100), ylim = c(13, 100))+ mytheme

# 2. 用 iml 包识别重要蛋白
library(iml)
# 构建 iml 的 Predictor 对象
predictor <- Predictor$new(model = learner_final, data = X, y = y) 
# 使用 permutation importance
imp <- FeatureImp$new(predictor, loss = "mae")  # or "rmse"
# 查看前几个重要蛋白
imp_dt_top <- head(imp$results[order(-imp$results$importance), ], 20); imp_dt_top<-as.data.table(imp_dt_top)
imp_dt_top<-merge(imp_dt_top, unique(name_id_reviewd[,c(1,4)]), by.x="feature", by.y="pro_id"); imp_dt_top<-imp_dt_top[order(-importance)]
fwrite(imp_dt_top, "~/Desktop/省皮/project/pQTL/manuscript/Table.XS Protein_age_model_pro_importance.csv", col.names = T, row.names = F, sep=",", quote = F)
# 3. 可视化 top features
imp_dt_top<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS3_Protein_age_model_importance.csv", header = T)
ggplot(imp_dt_top, aes(x = reorder(gene_name, importance), y = importance)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = importance.05, ymax = importance.95), linewidth = 1.1, width = 0.2, color = "steelblue") +
  coord_flip() + labs(x = "Gene", y = "Permutation Importance") +
  mytheme + theme(axis.text.y = element_text(face = "italic"),legend.position = "none")

## time-resolved DEPs ####
jj1<-fread("~/Desktop/省皮/project/pQTL/manuscript/Protein_ageMean_cor.csv", header = T)
jj2<-jj1[rho_p_age<0.05] # rho_p_age<1 for all protiens

# Step 0: 预设参数
window_size <- 20; step_size <- 5
age_min <- min(jj2$age); age_max <- max(jj2$age); age_windows <- seq(age_min, age_max - window_size, by = step_size)
# Step 1: 得到 protein list（所有蛋白 + 各组织蛋白总数）
all_proteins <- unique(jj2$pro_id); n_detected_proteins <- length(all_proteins)
time_resolved_DEPs2 <- list() # 初始化 list 储存每个窗口的 time-resolved DEPs
# Step 2: 滑窗循环
for (start_age in age_windows) {
  window_label <- paste0("Age_", start_age + window_size/2)
  end_age <- start_age + window_size
  dat_win <- jj2[age >= start_age & age <= end_age] # 提取当前窗口的数据
  if (length(unique(dat_win$age)) < 5) next   # 若窗口样本太少跳过

  # 对每个蛋白进行 spearman 相关分析
  dep_result <- dat_win[, {
    cor_test <- cor.test(intensity2, age, method = "spearman")
    list(rho = cor_test$estimate, p = cor_test$p.value)}, by = pro_id]
  
  dep_sig <- dep_result[p < 0.05, pro_id] # 筛选显著的 DEPs (p < 0.05)
  time_resolved_DEPs2[[window_label]] <- unique(dep_sig)   # 保存每个 window 的DEP结果
  }
# Step 3: 构建 cumulative DEP 曲线
cumulative_set <- c(); cumulative_counts2 <- data.frame()

for (i in seq_along(time_resolved_DEPs2)) {
  window_name <-  names(time_resolved_DEPs2)[i]
  current_prots <- time_resolved_DEPs2[[i]]
  
  cumulative_set <- union(cumulative_set, current_prots)
  count <- length(cumulative_set)
  
  cumulative_counts2 <- rbind(cumulative_counts2, data.frame(
    window = window_name,
    midpoint_age = as.numeric(gsub("Age_", "", window_name)),
    cumulative_unique_DEPs = count,
    cumulative_percent = count / n_detected_proteins * 100))
}

# Step 4: 画图
library(ggplot2); library(patchwork)

setDT(cumulative_counts2)
cumulative_counts2$count_DEPs<-as.numeric(sapply(time_resolved_DEPs2, length))
cumulative_counts2[, delta_DEPs := c(0, diff(cumulative_unique_DEPs))] # cumulative_unique_DEPs[1]
cumulative_counts2[, delta_cumperc := c(0, diff(cumulative_percent))] # cumulative_percent[1]

scale_factor <- max(cumulative_counts2$delta_DEPs) / max(cumulative_counts2$count_DEPs) # 计算缩放因子，让柱状图能显示在同一张图上
cumulative_counts2[midpoint_age==min(midpoint_age)]$delta_DEPs<-NA # 第一个时间点的 newly merged 不限制
ggplot(cumulative_counts2, aes(x = midpoint_age)) +
  geom_line(aes(y = delta_DEPs), color = "firebrick", size = 1, linetype = "solid") +
  geom_point(aes(y = delta_DEPs), color = "firebrick", size = 2) + # 左 y 轴: 新出现 ACPs (红色虚线+点)
  geom_col(aes(y = count_DEPs * scale_factor), fill = "gray50", alpha = 0.6) +  # 右 y 轴: 总 ACPs 柱状图, 缩放显示
  scale_x_continuous(breaks = cumulative_counts2$midpoint_age) +
  scale_y_continuous(name = "Newly emerged ACPs per window", sec.axis = sec_axis(~ . / scale_factor, name = "Total ACPs per window")) +
  labs(x = "Age") + theme_minimal(base_size = 18) +
  theme(text = element_text(face = "bold"),axis.title.y.left = element_text(color = "firebrick"), axis.text.y.left = element_text(color = "firebrick"), 
        axis.title.y.right = element_text(color = "gray50"), axis.text.y.right = element_text(color = "gray50"))


#### RNA与protein 的相关性 ####
## protein process ####
# 1/2 minum impution, VSN protein normalization, no lmFit, 同一age组取平均
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
meanSdPlot(fit) # 标准差(sd)在不同均值区间(rank mean)近乎恒定，说明拟合很好
mm1 <- predict(fit, mm) # 已经是log2 转化了
mm2 <- as.data.table(mm1); mm2$pro_id<-rownames(mm)
mm3<-melt.data.table(mm2, id.vars = "pro_id"); names(mm3)<-c("pro_id","sample","intensity")

## RNA process by FPKM ####
name_id_reviewd<-fread("~/Desktop/reference/human_uniprotkb_proteome_UP000005640_reviewed_2023_09_07_ID.txt",header = T)
length(unique(name_id_reviewd$pro_id)); length(unique(name_id_reviewd$gene_id)) # 42433; 32814 , pro_id:gene_id=n:1更多，把 gene_id转为 pro_id

rr<-fread("~/Desktop/省皮/project/pQTL/gene_sample_count_with_symbol_GO_KEGG_FPKM_islnc.xls") # a total table contain raw cnt, log2Normalized_cnt and FPKM value
sel_col<-grep("gene_ID|.+SRNA.+fpkm",names(rr),value=T); rr_rpkm<-rr[,..sel_col]; names(rr_rpkm)<-gsub("SRNA_fpkm","",names(rr_rpkm))
names(rr_rpkm)[1]<-"gene_id"; rr_rpkm<-rr_rpkm[grep('ENSG',gene_id)]; # 去掉公司预测的基因
rm(rr) # large memory load
rr_rpkm2<-melt(rr_rpkm, id.vars = "gene_id"); names(rr_rpkm2)<-c("gene_id","sample","value") # 这里不进行过滤，可在与蛋白 merge 后直接subset
rr_rpkm4<-merge(rr_rpkm2, unique(name_id_reviewd[,1:2]), by="gene_id")

## RNA vs protein: overall mean across sample ####
jj<-merge(mm3, rr_rpkm4[, .(pro_id, sample, value)], by=c("pro_id", "sample")); names(jj)[3:4]<-c("pro_raw","rna_raw")
jj<-jj[order(pro_id, sample, -pro_raw, -rna_raw)]; jj[,id:=paste(pro_id, sample, sep="_")]; jj<-jj[!duplicated(id)]
length(unique(mm3$pro_id)); length(unique(rr_rpkm4$pro_id)); length(unique(jj$pro_id)) # 5956; 19275; 5929
length(unique(mm3$sample)); length(unique(rr_rpkm4$sample)); length(unique(jj$sample)) # 248; 207; 199
jj[,c("pro_log","rna_log"):=.(pro_raw, log2(rna_raw+0.001))] # protein 已经是VSN log2 转化了；先 log，再 mean更稳健，减少极端值影响
jj[,c("pro_log_mean", "pro_log_median", "rna_log_mean", "rna_log_median"):=.(mean(pro_log), median(pro_log), mean(rna_log), median(rna_log)), by=.(pro_id)]

jj1<-unique(jj[, .(pro_id, pro_log_mean, pro_log_median, rna_log_mean, rna_log_median)])
cor.test(jj1$pro_log_mean, jj1$rna_log_mean, method = "spearman") # 0.269; 0.116 for pearson
cor.test(jj1$pro_log_median, jj1$rna_log_median, method = "spearman") # 0.272; 0.105 for pearson
tmp<-jj1[pro_log_mean>0 & rna_log_mean>0]; length(unique(tmp$pro_id)) # 5115
cor.test(tmp$pro_log_mean, tmp$rna_log_mean, method = "spearman") # 0.327; 0.398 for pearson
cor.test(tmp$pro_log_median, tmp$rna_log_median, method = "spearman") # 0.325; 0.395 for pearson

## RNA vs protein: each gene across sample by age group ####
hh2<-readxl::read_excel("~/Desktop/省皮/project/pQTL/DIA样本信息_省立正常对照_外送20200623.xls", sheet = 1) #里面含有采样时间信息
hh2<-as.data.table(hh2); names(hh2)[1:2]<-c("sample","date"); hh2<-hh2[,1:2]
hh2[,date2:=as.Date(date, format = "%Y-%m-%d")]; hh2[,date_duration:=as.numeric(date2 - min(date2, na.rm = TRUE))]
hh3<-fread("~/Desktop/省皮/project/pQTL/manuscript/TableS1_249smp_info.csv")
hh3<-merge(hh3, hh2[,c(1,4)], by.x="No", by.y="sample", all.x=T)
summary(hh3$age) # 13.00   38.00   52.00   50.14   62.00   95.00 
summary(hh3$date_duration) # 1.0    62.0   121.0   135.9   199.0   304.0 

jj<-merge(mm3, rr_rpkm4[, .(pro_id, sample, value)], by=c("pro_id", "sample")); names(jj)[3:4]<-c("pro_raw","rna_raw")
jj<-jj[order(pro_id, sample, -pro_raw, -rna_raw)]; jj[,id:=paste(pro_id, sample, sep="_")]; jj<-jj[!duplicated(id)]
length(unique(mm3$pro_id)); length(unique(rr_rpkm4$pro_id)); length(unique(jj$pro_id)) # 5956; 19275; 5929
length(unique(mm3$sample)); length(unique(rr_rpkm4$sample)); length(unique(jj$sample)) # 248; 207; 199
jj[,c("pro_log","rna_log"):=.(pro_raw, log2(rna_raw+0.001))] # protein 已经是VSN log2 转化了；先 log，再 mean更稳健，减少极端值影响
jj[,c("pro_log_mean", "pro_log_median", "rna_log_mean", "rna_log_median"):=.(mean(pro_log), median(pro_log), mean(rna_log), median(rna_log)), by=.(pro_id)]

jj<-merge(jj, hh3[, c(1:3,8)], by.x="sample", by.y="No"); jj<-jj[order(pro_id, age, -date_duration)]
jj[, age_group := cut(age, breaks = c(13, 38, 52, 62, 95), labels = c("age1", "age2", "age3", "age4"), include.lowest = TRUE, right = TRUE)]
jj[, date_group := cut(date_duration, breaks = c(1, 62, 121, 199, 304), labels = c("date1", "date2", "date3", "date4"), include.lowest = TRUE, right = TRUE)]
table(jj[pro_id=="A0FGR8"]$age_group) # age1 age2 age3 age4 = 59   51   44   45  
table(jj[pro_id=="A0FGR8"]$date_group) # date1 date2 date3 date4  = 52    49    40    58 
table(jj[pro_id=="A0FGR8"]$gender) # female   male = 72    127 
jj[, c("rho_age", "rho_p_age", "R_age", "R_p_age") := .(cor.test(pro_log, rna_log, method = "spearman")$estimate, cor.test(pro_log, rna_log, method = "spearman")$p.value,
                                                        cor.test(pro_log, rna_log, method = "pearson")$estimate,  cor.test(pro_log, rna_log, method = "pearson")$p.value), by = .(pro_id, age_group)]
jj[, c("rho_date", "rho_p_date", "R_date", "R_p_date") := .(cor.test(pro_log, rna_log, method = "spearman")$estimate, cor.test(pro_log, rna_log, method = "spearman")$p.value,
                                                            cor.test(pro_log, rna_log, method = "pearson")$estimate,  cor.test(pro_log, rna_log, method = "pearson")$p.value), by = .(pro_id, date_group)]
jj[, c("rho_sex", "rho_p_sex", "R_sex", "R_p_sex") := .(cor.test(pro_log, rna_log, method = "spearman")$estimate, cor.test(pro_log, rna_log, method = "spearman")$p.value,
                                                        cor.test(pro_log, rna_log, method = "pearson")$estimate,  cor.test(pro_log, rna_log, method = "pearson")$p.value), by = .(pro_id, gender)]

tmp<-jj[pro_log_mean>0 & rna_log_mean>0]; length(unique(tmp$pro_id)) # 5115
jj_age<-unique(jj[pro_id%in%tmp$pro_id, .(pro_id, age_group, rho_age, rho_p_age, R_age, R_p_age)]); jj_age<-jj_age[!is.na(rho_age)]

# count barplot
nrow(jj_age[age_group=="age1"]); nrow(jj_age[age_group=="age1" & rho_p_age<0.05]); nrow(jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age>0]); nrow(jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age<0]) # 5115; 1026; 940; 86 /1026=91.62 
nrow(jj_age[age_group=="age2"]); nrow(jj_age[age_group=="age2" & rho_p_age<0.05]); nrow(jj_age[age_group=="age2" & rho_p_age<0.05 & rho_age>0]); nrow(jj_age[age_group=="age2" & rho_p_age<0.05 & rho_age<0]) # 5115; 753; 662; 91 /753=87.92
nrow(jj_age[age_group=="age3"]); nrow(jj_age[age_group=="age3" & rho_p_age<0.05]); nrow(jj_age[age_group=="age3" & rho_p_age<0.05 & rho_age>0]); nrow(jj_age[age_group=="age3" & rho_p_age<0.05 & rho_age<0]) # 5115; 681; 568; 113 /681=83.41
nrow(jj_age[age_group=="age4"]); nrow(jj_age[age_group=="age4" & rho_p_age<0.05]); nrow(jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age>0]); nrow(jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age<0]) # 5115; 598; 526; 72 /598=87.96
chisq.test(matrix(c(5115, 940, 5115, 526), byrow = T, ncol = 2))$`p.value` # 6.069992e-24
mymat<-data.table(group = rep(c("Young", "Old"), each = 2), status = rep(c("Positive", "Negative"), times = 2),count = c(940, 86, 526, 72))
mymat$group<-factor(mymat$group, levels = c("Young", "Old")); mymat$status<-factor(mymat$status, levels = c("Negative", "Positive"))
ggplot(mymat, aes(x = group, y = count, fill = status)) + 
  geom_bar(stat = "identity", width = 0.7) + labs(x = NULL, y = "Number of genes") + 
  scale_fill_manual(values = c("Positive" = "indianred", "Negative" = "steelblue"), labels = c("Negative cor", "Positive cor")) +
  mytheme+ theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank())

# decoupling GO
tmp1<-jj_age[age_group=="age1" & rho_p_age<0.05 & rho_age>0] # positive in young
tmp2<-jj_age[age_group=="age4" & rho_p_age<0.05 & rho_age>0] # positive in old
old_loss<-setdiff(tmp1$pro_id, tmp2$pro_id); old_gain<-setdiff(tmp2$pro_id, tmp1$pro_id)
nrow(jj_age[age_group=="age1"]); nrow(tmp1); nrow(tmp2); length(age_loss); length(age_gain) # 5115; 940; 526; 687; 273

df <- data.frame(age_group = rep(c("Young", "Old"), each = 2), status = rep(c("Positive cor", "Positive loss in old"), 2), count = c(940, 0, 253, 687))
df$age_group <- factor(df$age_group, levels = c("Young", "Old"))
ggplot(df, aes(x = age_group, y = count, fill = status)) + geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual( values = c("Positive cor" = "indianred", "Positive loss in old" = "grey70"), name = NULL, labels = c("Couping in young", "Decoupling in Old")) +
  labs(x = NULL, y = "Number of Genes") + mytheme +
  theme(legend.position = "top",legend.direction = "horizontal", legend.text = element_text(size = 16))

library("clusterProfiler"); library("org.Hs.eg.db"); library(enrichplot)
enrich_go3 <- enrichGO( old_loss, # ECM/血管/没指望
                        OrgDb = org.Hs.eg.db, keyType = "UNIPROT", universe= unique(jj_age$pro_id), # background genes
                        ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
dotplot(enrich_go3, x='GeneRatio',showCategory = 20, font.size = 18); barplot(enrich_go3, showCategory = 20, x = "GeneRatio", color = "p.adjust")

tmp3<-as.data.table(enrich_go3@result); tmp3<-tmp3[order(p.adjust)]; tmp3<-tmp3[p.adjust<0.05] # 手动简练
openxlsx::write.xlsx(tmp3, file = "~/Desktop/省皮/project/pQTL/manuscript/GO_RNA_pro_decoupling_old.xlsx", sheetName = "Sheet1", rowNames = FALSE) 
tmp4<-readxl::read_excel("~/Desktop/省皮/project/pQTL/manuscript/GO_RNA_pro_decoupling_old.xlsx", sheet = 2); tmp4<-as.data.table(tmp4)
df_plot <- tmp4[order(p.adjust)]; df_plot$negLogP <- -log10(df_plot$p.adjust)
df_plot$Description <- factor(df_plot$Description, levels = rev(df_plot$Description[order(df_plot$negLogP, decreasing = TRUE)]))
ggplot(df_plot, aes(x = negLogP, y = Description, fill = Count)) + geom_bar(stat = "identity") +
  scale_fill_gradient(name = "Gene Count", low = "mistyrose", high = "brown4") +
  labs(x = "-log10(p.adj)", y = NULL) + mytheme + 
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.title = element_text(size = 14), legend.text = element_text(size = 14))

# rho compare ECDF
summary(jj_age[age_group=="age1"]$rho_age) # -0.48966 -0.02580  0.08925  0.09902  0.21189  0.76639 
summary(jj_age[age_group=="age4"]$rho_age) # -0.53715 -0.04427  0.06640  0.07167  0.18103  0.72714
wilcox.test(jj_age[age_group=="age1"]$rho_age, jj_age[age_group=="age4"]$rho_age, paired = TRUE) # 3.964501e-20,用这个
wilcox.test(jj_age[age_group=="age1"]$rho_age, jj_age[age_group=="age4"]$rho_age) # 2.081e-12
ks.test(jj_age[age_group=="age1"]$rho_age, jj_age[age_group=="age4"]$rho_age) # 1.086e-08
ggplot(jj_age[age_group %in% c("age1", "age4")], aes(x = rho_age, color = age_group)) + 
  stat_ecdf(geom = "step", linewidth = 1) +
  scale_color_manual(values = c("age1" = "indianred", "age4" = "steelblue"), labels = c("age1" = "Young", "age4" = "Old")) +
  labs(x = "Spearman rho", y = "Cumulative Distribution Fraction", color = NULL) + coord_cartesian(xlim = c(-0.3, 0.6)) +
  mytheme+theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.key.width = unit(2.5, "lines"), legend.text=element_text(size=18, face = "bold"))
