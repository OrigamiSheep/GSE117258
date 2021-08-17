library(tidyverse)

library(GEOquery)
GSEs <- getGEO('GSE117285',destdir = ".",
               AnnotGPL = F,
               getGPL = F) 
GSEs <- GSEs[[1]]
#如果是本地文件读取 GSEs <- getGEO(filename = "GSE117285_series_matrix.txt.gz",AnnotGPL = F, getGPL = F)
GSE_expr <- exprs(GSEs)
GSE_pd <- pData(GSEs)

#提取表达矩阵
View(GSE_pd)
GSE_targets <- GSE_pd %>%
  dplyr::select(title, geo_accession, characteristics_ch1) %>%
  mutate(group = if_else(str_sub(characteristics_ch1, 9, 9) %in% c("N", "o"),
                         "Normal", "EMPD"))
view(GSE_targets)

#normalization
boxplot(GSE_expr, las = 2, outline = FALSE, main = "GSE117285")

#数据转换
#log转换
ex <- GSE_expr
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
GSE_expr_log2 <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}
GSE_expr_log2 <- na.omit(GSE_expr_log2)

#注释探针
#注释基因
GPL17077_anno <- data.table::fread("GPL17077-17467.txt", skip = "ID")


library(dplyr)
library(tidyr)
GSE117285_probe2symbol <- GPL17077_anno %>%
  dplyr::select(ID, GENE_SYMBOL)

#探针转换
GSE117285_anno <- GSE_expr_log2 %>%
  as.data.frame() %>%
  dplyr::select(GSE_targets$geo_accession) %>%
  rownames_to_column(var = "ID") %>%
  inner_join(GSE117258_probe2symbol, ., by = "ID") %>%
  dplyr::select(-1) %>%
  as.data.frame() %>%
  aggregate(.~ GENE_SYMBOL, data = ., mean) %>%
  column_to_rownames(var = "GENE_SYMBOL")

save(GSE117285_anno,file = 'GSE117285_anno.Rdata')
##差异分析
{
  library(limma)
  #加载芯片数据&定义分组
  group <- factor(GSE_targets$group, levels = c("Normal", "EMPD"))
  expr <- GSE117285_anno
  contrast <- paste0(rev(levels(group)), collapse = "-")
  contrast
  
  #分组矩阵
  design <- model.matrix( ~ 0 + group)
  colnames(design) <- levels(group)
  #差异表达矩阵
  contrast.matrix <- makeContrasts(contrast, levels = design)
  contrast.matrix
  
  #limma差异表达分析
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  #查看结果
  GSE117285_DEG <- topTable(fit, coef = 1, n = Inf) %>%
    rownames_to_column(var = "symbol")
  #筛选|logFC|>0.5且p<0.05的结果
  GSE117285_DEG_p <- GSE117285_DEG %>%
    dplyr::filter(abs(logFC) > 2, P.Value < 0.05)
  head(GSE117285_DEG_p)
}
save(GSE117285_DEG, GSE117285_DEG_p, file = "GSE117285_DEG.Rdata")
write.csv(GSE117285_DEG_p, file = "data/GSE117285_DEG_p.csv")

##火山图与热图
{
  #火山图
  plot_df <- GSE117285_DEG
  #定义阈值
  plot_df$threshold = factor(ifelse(plot_df$P.Value < 0.05 & abs(plot_df$logFC) > 2,
                                    ifelse(plot_df$logFC > 2, 'Up', 'Down'), 'NS'),
                             levels = c('Up', 'Down', 'NS'))
  #绘图
  library(ggplot2)
  library(ggrepel)
  pdf("figures/Volcano.pdf", 12,8)
  ggplot(plot_df, aes(y = logFC, x= -log10(P.Value), color = threshold)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("#DC143C", "#143cdc", "#808080")) +
    xlab('-log10 (P-value)') +
    ylab('log2 (FoldChange)') +
    ggtitle("Volcano Plot") +
    geom_hline(yintercept = c(-2, 2), lty = 3, col = "black", lwd = 0.5) +
    geom_vline(xintercept = -log10(0.05), lty = 3, col = "black", lwd = 0.5) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  dev.off()
  
  #热图
  library(pheatmap)
  heat <- GSE117285_DEG
  heat <- arrange(heat, desc(abs(logFC)))
  rownames(heat) <- heat[, 1]
  heat <- heat[, -1]
  choose_gene = head(rownames(heat), 25)
  choose_matrix = GSE117285_anno[choose_gene,]
  choose_matrix = t(scale(t(choose_matrix)))   
  p = pheatmap(choose_matrix,filename = 'figures/DEG_top25_heatmap.pdf') 
}

##富集分析
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  
  input <- GSE117285_DEG_p
  #转换ID
  gene.df <- bitr(input$symbol,
                  fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"), 
                  OrgDb = org.Hs.eg.db)
  gene <- gene.df$ENTREZID
  
  
  #GO
  ego <- enrichGO(gene = gene,
                  OrgDb = org.Hs.eg.db,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  head(ego)
  #KEGG
  kk <- enrichKEGG(gene = gene,
                   organism = "hsa",
                   pvalueCutoff = 0.05)
  k <- as.data.frame(head(kk))
  #可视化
  #柱状图
  pdf("figures/GO_barplot.pdf", 20, 8)
  barplot(ego, showCategory = 10, color = "pvalue")
  dev.off()
  
  pdf("figures/KEGG_barplot.pdf", 20, 8)
  barplot(kk, showCategory = 10, color = "pvalue")
  dev.off()
  #气泡图
  pdf("figures/GO_dotplot.pdf", 20, 8)
  dotplot(ego, showCategory = 10, color = "pvalue")
  dev.off()
  
  pdf("figures/KEGG_dotplot.pdf", 20, 8)
  dotplot(kk, showCategory = 10, color = "pvalue")
  dev.off()
  #网络图
  pdf("GO_netplot.pdf", 20, 8)
  cnetplot(ego, categorySize = "pvalue",
           circular = TRUE, colorEdge = TRUE,
           node_label = "category")
  dev.off()
  
  pdf("KEGG_netplot.pdf", 10, 8)
  cnetplot(kk, categorySize = "pvalue",
           circular = TRUE, colorEdge = TRUE,
           node_label = "category")
  dev.off()
  
  
  
##免疫浸润分析
  #install_packages
  {### 安装estimate包
  library(utils)
  rforge <- "http://r-forge.r-project.org"#开发中的包
  if(!"estimate" %in% installed.packages()) 
    install.packages("estimate", repos=rforge, dependencies=TRUE)#装不上的话应该是网络问题
  
  ### 安装MCPcounter
  if(!"e1071" %in% installed.packages())
    install.packages("e1071")
  if(!"preprocessCore" %in% installed.packages()) 
    BiocManager::install("preprocessCore")
  
  if(!"curl" %in% installed.packages())
    install.packages("curl")
  if(!"devtools" %in% installed.packages())
    install.packages("devtools")
  
  # 从github上安装MCPcounter
  if(!"MCPcounter" %in% installed.packages()){
    library(devtools)
    install_github("ebecht/MCPcounter",ref="master", subdir="Source")
  }
  
  #如果使用上面的代码从GitHub安装失败，可以使用以下代码本地安装
  if(!"MCPcounter" %in% installed.packages())
    install.packages("MCPcounter_1.1.0.tar.gz", repos = NULL, type = "source")
  
  ### 安装GSVA
  if(!"GSVA"%in%installed.packages()) 
    BiocManager::install("GSVA")
  
  ### 安装pheatmap
  if(!"pheatmap"%in%installed.packages()) 
    BiocManager::install("pheatmap")
  
  if(!dir.exists("data")) dir.create("data")
  if(!dir.exists("figures")) dir.create("figures")}
 
  #estimate
  library(estimate)
  
  input.f <- GSE117285_anno
    input.f$ID <- rownames(input.f) 
    input.f <- select(input.f, ID, GSM3290431, GSM3290432, GSM3290433, GSM3290434)
    input.f <- input.f[-1, ]
    write.table(input.f, file = "input.txt", sep = "\t", row.names = F, quote = F)
    input.c <- select(input.f, ID, everything())
    write.table(input.c, file = "input_c.txt", sep = "\t", row.names = F, quote = F)

  
  filterCommonGenes(input.f="input.txt", output.f="data/estimate_input.gct", id="GeneSymbol")#过滤掉细胞中共有的基因
  
  ### 计算免疫浸润打分 ###
  estimateScore("data/estimate_input.gct", "data/estimate_score.gct", platform="affymetrix")#计算肿瘤微环境的评分 ##这里paltform可以选c("affymetrix", "agilent", "illumina")，但只有affymetrix计算结果中会有肿瘤纯度，才能画图
  scores <- read.table("data/estimate_score.gct", skip = 2, header = T)
  
  ### 可视化 ###
  # 绘制单个样本
  plotPurity(scores = "data/estimate_score.gct", samples="GSM3290431", platform="affymetrix", output.dir = "figures/")
  # 绘制多个样本
  plotPurity(scores="data/estimate_score.gct", platform="affymetrix", output.dir = "figures/")
  
  ###### CIBERSORT免疫浸润注释
  # 加载CIBERSORT函数
  source("CIBERSORT.R")
  
  # 运行CIBERSORT分析
  CIBERSORT_Result <- CIBERSORT("data/CIBERSORT_Ref.txt", #模型的参考矩阵文件路径
                                "input_c.txt", #要预测的混合样本表达矩阵
                                perm=1000, QN=TRUE)#perm，置换检验的次数，默认为0，大于100次时，可以计算显著性p值；QN，是否对输入的表达矩阵进行标准化
  #相当于不断计算1000次，看结果在1000次是否是随机的
  # 将结果写出到文件
  write.table(CIBERSORT_Result, file = "data/CIBERSORT_Result.txt",sep="\t", quote=F, col.names=T)
  
  # 提取CIBERSORT的打分结果
  CIBERSORT_Score <- CIBERSORT_Result[,1:(ncol(CIBERSORT_Result)-3)]#因为绘图不需要最后三列，所以去掉最后三列，取22列
  
  # 生成标准化函数normalize，构建一个function
  normalize <- function(x){# （x-x的最小值）/(X的最大值-X的最小值)
    return((x-min(x))/(max(x)-min(x)))}
  
  # 使用normalize标准化CIBERSORT打分
  CIBERSORT_Score <- normalize(CIBERSORT_Score)
  
  # 将处理后CIBERSORT的打分写出到文件
  write.table(CIBERSORT_Score, file = "data/CIBERSORT_Score.txt",sep="\t", quote=F, col.names=T)
  
  # 加载pheatmap包
  library(pheatmap)
  
  # 构建热图注释数据框
  annotation <- data.frame(rownames(CIBERSORT_Score))
  colnames(annotation) <- "sample"#设定列名
  rownames(annotation) <- rownames(CIBERSORT_Score)#设定行名
  annotation_row <- cbind(annotation, CIBERSORT_Result[,(ncol(CIBERSORT_Result)-2):ncol(CIBERSORT_Result)])#合并列，CIBERSORT_Result的最后三列
  
  # 绘制热图
  pheatmap(CIBERSORT_Score,
           show_colnames = T,
           cluster_rows = F,
           cluster_cols = F,#cluster_cols = F 对列不进行聚类
           annotation_row = annotation_row,#注释数据框放在行
           cellwidth=10,
           cellheight=40,
           fontsize=5,
           fontsize_row=5,
           scale="row",
           fontsize_col=3,
           filename = 'figures/CIBERSORT-heatmap.tiff')
  
  #可视化
  results <- read.table("data/CIBERSORT_Score.txt", sep = '\t', header = T) 
  results_p <- t(results)
  
  bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
  # breaks用来定义数值和颜色的对应关系。
  
  # Step4:将CIBERSORT_Result进行可视化
  #1）热图
  library(pheatmap)
  library(RColorBrewer)
  pdf("figures/pheatmap.pdf",12,8)
  pheatmap(
    results_p,
    breaks = bk,
    cluster_cols = F,
    scale = "row",
    cluster_row = F,
    border_color = NA,
    show_colnames = T,
    show_rownames = T,
    color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
              colorRampPalette(colors = c("white","red"))(length(bk)/2)
    ))
  dev.off()
  #调整参数让热图更加美观。
  
  #柱状图可视化细胞占比预测
  library(RColorBrewer)
  mypalette <- colorRampPalette(brewer.pal(8,"Set1"))###8种颜色
  cibersort_barplot <- results %>%
    gather(key = Cell_type,value = Proportion,1:22)
  #使用RColorBrewer包配置需要的色彩方案，使用gather函数中的key-value对应关系重建细胞名称和比例的对应关系并赋值给cibersort_barplot
  
  #cibersort_barplot$Patient1 <- factor(cibersort_barplot$Patient,
  #                                   levels = str_sort(unique(cibersort_barplot$Patient),
  #                                                      numeric = T))
  pdf("figures/barplot.pdf",12,8)
  ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
    geom_bar(position = "stack",stat = "identity") +
    labs(fill = "Cell Type",x = '',y = "Estiamted Proportion") + theme_bw() +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_y_continuous(expand = c(0.01,0)) +
    scale_fill_manual(values = mypalette(23))
  dev.off()
  #调整参数让柱状图更加美观。
  
  #直观箱线图
  pdf("figures/box.pdf", 12, 8)
  ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
    geom_boxplot(outlier.shape = 21,coulour = "black") + theme_bw() + 
    labs(x = "", y = "Estimated Proportion") +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = mypalette(23))
  dev.off()
    ##区分肿瘤及正常组
  cibersort_barplot_normal <- results[1:4, ] %>% 
    gather(key = Cell_type,value = Proportion,1:22)
  cibersort_barplot_tumor <- results[5:8, ] %>%
    gather(key = Cell_type,value = Proportion,1:22)
    
  ##肿瘤
  pdf("figures/box_tumor.pdf", 12, 8)
  ggplot(cibersort_barplot_tumor,aes(Cell_type,Proportion,fill = Cell_type)) + 
    geom_boxplot(outlier.shape = 21,coulour = "black") + theme_bw() + 
    labs(x = "", y = "Estimated Proportion") +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = mypalette(23))
  dev.off()
  
  ##正常
  pdf("figures/box_normal.pdf", 12, 8)
  ggplot(cibersort_barplot_normal,aes(Cell_type,Proportion,fill = Cell_type)) + 
    geom_boxplot(outlier.shape = 21,coulour = "black") + theme_bw() + 
    labs(x = "", y = "Estimated Proportion") +
    theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
    scale_fill_manual(values = mypalette(23))
  dev.off()
  
  
  #调整参数让柱状图更加美观。
  
  #corHeatmap
  rt=results
  flag <- apply(rt,2,function(x) sum(x == 0) < 
                  dim(rt)[1]/2)
  # 筛选出0值超过样本的一半的一些细胞
  rt <- rt[, which(flag)] %>%
    as.matrix() %>%
    t()
  # 留下在大部分样本中有所表达的细胞。
  
  library(corrplot)
  pdf("figures/corHeatmap.pdf",height=13,width=13)
  corrplot(corr=cor(rt),
           method = "color",
           order = "hclust",
           tl.col="black",
           addCoef.col = "black",
           number.cex = 0.8,
           col=colorRampPalette(c("blue", "white", "red"))(50),
  )
  dev.off()
  #筛选肿瘤做相关分析
  rt_tumor <- rt[-c(1:4), ]
   rt_tumor <- rt_tumor[, which(flag)] %>%
    as.matrix() 
  rt_tumor_p <- cor.mtest(rt_tumor)
  pdf("figures/corHeatmap_T.pdf",height=13,width=13)
  corrplot(corr=cor(rt_tumor),
           method = "color",
           order = "hclust",
           p.mat = rt_tumor_p$p,
           sig.level = 0.05,
           insig = "pch",
           tl.col="black",
           addCoef.col = "black",
           number.cex = 0.8,
           col=colorRampPalette(c("blue", "white", "red"))(50),
  )
  dev.off()

  
  
 ##两组对比
  library(reshape2)
  library(ggpubr)
  
  rt <- CIBERSORT_Score %>%
    as.data.frame()
  rt$type <- c("normal", "normal", "normal", "normal", "tumor", "tumor", "tumor", "tumor")  
  rt <- select(rt, type, everything())
  x <- colnames(rt)[1]
  colnames(rt)[1] <- "Type"  
  data <- melt(rt, id.vars = c("Type"))
  colnames(data) <- c("Type", "Lymphocytes", "Value")
  head(data)  
    #绘图boxplot
  p <- ggboxplot(data, x = "Lymphocytes", y = "Value", color = "Type",
                 ylab = "",
                 xlab = "",
                 legend.title = x,
                 palette = c("blue", "red"),
                 width = 0.6, add = "none")
  
  p = p + rotate_x_text(60)
  p1 = p + stat_compare_means(aes(group = Type),
                              method = "wilcox.test",
                              symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                              label = "p.signif")
    #输出
  pdf("boxplot.pdf", width = 12, height = 8)
  print(p1)
  dev.off()
  

##提取出IL-17信号通路富集基因，之后与免疫细胞进行相关分析
  immune <- subset(results, select = T.cells.gamma.delta) %>%
    t() %>%
    as.data.frame() %>%
    subset(., select = c(5:8)) %>%
    t() %>%
    as.data.frame()
  gene <- subset(GSE117285_anno, select = c(5:8)) %>%
    t() %>%
    as.data.frame() %>%
    subset(., select = c("DEFB4A", "S100A7A", "USP25", "CXCL1", "S100A7", "CXCL10", "MMP13", "S100A8", "MAPK6", "CCL11", "CCL20", "IFNG", "LCN2"))
  
  gene$T.cells.gamma.delta <- immune$T.cells.gamma.delta
  	
 data <- gene %>% 
   select(T.cells.gamma.delta, everything()) %>%
   gather(key = gene_cell,value = expression,1:14)  #宽数据
   
  #Visualize
 compare_means(expression ~ gene_cell, data = data, ref.group = "T.cells.gamma.delta",
               method = "wilcox.test")
  pdf("figures/compare_gene_cell.pdf", 20, 8)
  ggboxplot(data, x = "gene_cell", y = "expression",
            color = "gene_cell", add = "jitter") +
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = "T.cells.gamma.delta")
  dev.off()
  
  compare_means(expression ~ gene_cell, data = data,
                method = "anova")
  
  ##cor
  library(ggExtra)
  library(ggplot2)
  library(ggpubr)
  corT <- cor.test(data$T.cells.gamma.delta, data$LCN2, method = "spearman")  
  cor <- corT$estimate  
  pValue <- corT$p.value  

  #ssGSEA
  # 清空当前变量
  rm(list = ls())
  
  ##### 关于读入gmt格式文件的代码示例，
  ##### 请参考模块02的GSVA分析部分，
  ##### 这里演示读入自定义数据集。
  ## 读入细胞注释文件
  immune <- read.table("pan_cancer_immune.txt", sep="\t", 
                       header=T, check.names=F)
  
  # 提取细胞类型并去重
  CellType <- as.vector(unique(immune$`Cell type`))
  
  # 生成免疫细胞对应的基因标签的列表，用于ssGSEA分析
  geneset <- lapply(CellType, function(x){
    x <- as.vector(immune[immune$`Cell type`==x,1])
    return(x)
  })
  names(geneset)<- CellType
  
  ### 计算ssGSEA score
  # 加载GSVA包
  library(GSVA)
  # 读取输入表达矩阵
  input <- GSE117285_anno
  
  
  # 表达矩阵数据框转换为矩阵
  exp <- as.matrix(input)
  
  # 进行ssGSEA分析
  gsva_matrix<- gsva(exp, 
                     geneset,
                     method='ssgsea',#设置成ssgsea
                     kcdf='Gaussian',#高斯分布
                     abs.ranking=TRUE)
  
  # 设置一个标准化函数
  normalize <- function(x){
    return((x-min(x))/(max(x)-min(x)))}
  
  # 对ssGSEA结果进行标准化
  ssgseaScore <- normalize(gsva_matrix)
  # 将ssGSEA结果写出到文件
  write.table(ssgseaScore,file="data/ssgseaScore.txt",sep="\t", quote=F, col.names=T)
  
  # 构建热图注释信息
  annotation_col <- data.frame(colnames(ssgseaScore))
  colnames(annotation_col) <- "sample"
  rownames(annotation_col) <- colnames(ssgseaScore)
  
  # 绘制热图
  library(pheatmap)
  pheatmap(ssgseaScore,
           show_colnames = F,
           cluster_rows = F,
           cluster_cols = F,
           annotation_col = annotation_col,#注释信息在列
           cellwidth=5,
           cellheight=5,
           fontsize=5,
           gaps_row = c(12,20),
           filename = 'figures/ssgsea-heatmap.tiff')
  
  ## 将estimate和ssGSEA的结果合并展示在热图上
  # 读入estimate结果
  estimateScore <- read.table("data/estimate_score.gct", sep = "\t", skip = 2,header = T,row.names = 1)
  
  # 将estimate结果和其他注释信息合并
  estimateScore <- t(estimateScore)
  estimateScore <- estimateScore[-1,]#删掉第一行
  annotation_tumor <- annotation_col[-c(1:4),] %>%
    as.data.frame()
  colnames(annotation_tumor) <- "sample"
  rownames(annotation_tumor) <- annotation_tumor$sample
  annotation <- cbind(annotation_tumor, estimateScore)
  ssgseaScore_T <- as.data.frame(ssgseaScore) %>%
    select(. ,GSM3290431, GSM3290432, GSM3290433, GSM3290434)
  
  # 加载pheatmap软件包
  library(pheatmap)
  # 绘制热图
  pheatmap(ssgseaScore_T,annotation = annotation,
           cluster_cols = F,fontsize=5,fontsize_row=5,
           scale="row",show_colnames=F,
           fontsize_col=3,filename = 'figures/ssgsea-estimate-heatmap_T.tiff')
  
  
  
  ##两组对比
  library(reshape2)
  library(ggpubr)
  
  rt <- ssgseaScore %>%
    t() %>%
    as.data.frame()
  rt$type <- c("normal", "normal", "normal", "normal", "tumor", "tumor", "tumor", "tumor")  
  rt <- select(rt, type, everything())
  x <- colnames(rt)[1]
  colnames(rt)[1] <- "Type"  
  data <- melt(rt, id.vars = c("Type"))
  colnames(data) <- c("Type", "Lymphocytes", "Value")
  head(data)  
  #绘图boxplot
  p <- ggboxplot(data, x = "Lymphocytes", y = "Value", color = "Type",
                 ylab = "",
                 xlab = "",
                 legend.title = x,
                 palette = c("blue", "red"),
                 width = 0.6, add = "none")
  
  p = p + rotate_x_text(60)
  p1 = p + stat_compare_means(aes(group = Type),
                              method = "wilcox.test",
                              symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                              label = "p.signif")
  #输出
  pdf("figures/boxplot_ssgsea.pdf", width = 12, height = 8)
  print(p1)
  dev.off()
  
  ##CCR3&CCL11
  
  gene <- t(GSE117285_anno) %>%
    as.data.frame() %>%
    subset(., select = c("CCL11", "CCR3"))
  
  gene$type <- c("normal", "normal", "normal", "normal", "tumor", "tumor", "tumor", "tumor")
  gene <- select(gene, type, everything())
  x <- colnames(gene)[1]
  colnames(gene)[1] <- "Type"  
  data <- melt(gene, id.vars = c("Type"))
  colnames(data) <- c("Type", "Gene", "Value")
  
  #Visualize
  p <- ggboxplot(data, x = "Gene", y = "Value", color = "Type",
                 ylab = "",
                 xlab = "",
                 legend.title = x,
                 palette = c("blue", "red"),
                 width = 0.6, add = "none")
  
  p = p + rotate_x_text(60)
  p1 = p + stat_compare_means(aes(group = Type),
                              method = "wilcox.test",
                              symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                              label = "p.signif")
  #输出
  pdf("boxplot_gene.pdf", width = 12, height = 12)
  print(p1)
  dev.off()
