# interaction across gene


gene_ia = read.table("../Juicer/merge/ST.gene_IA", stringsAsFactors = F)
gene_ia = gene_ia[!duplicated(gene_ia$V1),]
rownames(gene_ia) = gene_ia$V1
for (i in Sys.glob("../Juicer/merge/*gene_IA")){
  ia = read.table(i, stringsAsFactors = F)
  ia = ia[!duplicated(ia$V1),]
  rownames(ia) = ia$V1
  sample = strsplit(basename(i), "[.]")[[1]][1]
  gene_ia[[sample]] = ia[rownames(gene_ia), 2] 
}
gene_ia = gene_ia[,c("LT","ST","CMP","CLP","MEP","GMP","MK","G")]
gene_ia = na.omit(gene_ia)
boxplot(log2(gene_ia+1))


gene_ia_norm = normalize.quantiles(as.matrix(gene_ia))
colnames(gene_ia_norm) = colnames(gene_ia)
rownames(gene_ia_norm) = rownames(gene_ia)
gene_ia_norm = as.data.frame(gene_ia_norm)



gene_ia = read.table("../Juicer/merge/ST.gene_IA", stringsAsFactors = F)
gene_ia = gene_ia[!duplicated(gene_ia$V1),]

expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = as.matrix(expr[,c(12,7,20,2,17,18,16,25)])
colnames(expr) = c("LT","ST","CMP","CLP","MEP","GMP","MK","G")


gene_ia_norm_expr = cbind(gene_ia_norm, expr[rownames(gene_ia_norm),])
gene_ia_norm_expr1 = gene_ia_norm_expr[,c(2,10)]
plot(gene_ia_norm_expr1[,1], log2(gene_ia_norm_expr1[,2]+1), xlim=c(0,10))

gene_ia_norm_expr1$expr = "0-1"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.1 <= 1, 3] = "0-1"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.1 > 1 & gene_ia_norm_expr1$ST.1 <= 10 , 3] = "1-10"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.1 > 10 & gene_ia_norm_expr1$ST.1 <= 30 , 3] = "10-30"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.1 > 30 & gene_ia_norm_expr1$ST.1 <= 80 , 3] = "30-80"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.1 > 80 , 3] = ">80"
gene_ia_norm_expr1$expr = factor(gene_ia_norm_expr1$expr, levels=c("0-1","1-10","10-30","30-80",">80"))
write.csv(gene_ia_norm_expr1, "../data/ST_gene_expr_interation.tab", row.names = T, quote = F)

ggplot(gene_ia_norm_expr1, aes(expr, ST, fill=expr))+
  geom_boxplot(size=1) + ylim(0,12)+
  theme_classic(base_size = 20) +
  xlab("gene expression (TPM)") + ylab("Normalized Intra-gene interaction") +
  theme(legend.position = "none")
ggsave("../figure/ST_Gene_interaction_expr.pdf", device = "pdf", width = 6, height = 5)


smoothScatter(log(gene_ia_norm_expr1[,1:2]+1))

pdf("../figure/ST_Gene_interaction_expr_scatter.pdf",5,5)
plot(gene_ia_norm_expr1[,1],log(gene_ia_norm_expr1[,2]+1),xlim=c(1,10),pch=20,xlab=("intra-gene interaction"),ylab=c("log(TPM+1)"),cex.lab=1.5, cex.axis=1.5)
reg1 = lm(log(gene_ia_norm_expr1[,2]+1)~gene_ia_norm_expr1[,1])
abline(reg1, col="darkred", lwd=2)
dev.off()


# G
gene_ia_norm_G = data.frame(gene_ia_norm[,c(2,8)])
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(6,9,25,27)]
#expr = expr[rowMeans(expr)>=1,]
gene_ia_norm_G$expr_ST = rowMeans(expr[rownames(gene_ia_norm_G),1:2])
gene_ia_norm_G$expr_G =  rowMeans(expr[rownames(gene_ia_norm_G),3:4])
gene_ia_norm_G$interaction = (gene_ia_norm_G$ST) - (gene_ia_norm_G$G)
gene_ia_norm_G[gene_ia_norm_G$interaction > 1, "interaction"] <- 1
gene_ia_norm_G[gene_ia_norm_G$interaction < -1, "interaction"] <- -1
gene_ia_norm_G$Gene = rownames(gene_ia_norm_G)
g=ggplot(gene_ia_norm_G, aes(log(expr_ST+1), log(expr_G+1), name=Gene))+
  geom_point(aes(col=interaction), shape=3, stroke=1.5) + xlim(0,6)+ylim(0,6)+
  xlab("log(ST expr.)")+ylab("log(G expr.)")+
  geom_abline(slope = 1, color="red", linetype="dashed", size=1)+
  theme_classic(base_size=15)+
  scale_colour_gradient2()
                                     

ggplotly(g) 


library(ggrepel)
g + geom_text_repel(
    data = gene_ia_norm_G[c("Syne1","App", "Lrrk2", "Lyst","Prkcb","Lyn","Abca13",
                            "Angpt1","Fut8","Eya2","Meis1", "Msi2","Erg"),],
    aes(label = Gene),
    size = 5
  )
ggsave("../figure/Gene_interaction_expr.pdf", device = "pdf", width = 6, height = 5)


# 
gene_ia_norm_G$interaction_change = "no-change"
gene_ia_norm_G[gene_ia_norm_G$interaction >=1, "interaction_change"] <- "decrease"
gene_ia_norm_G[gene_ia_norm_G$interaction <= -1, "interaction_change"] <- "increase"
gene_ia_norm_G$interaction_change = factor(gene_ia_norm_G$interaction_change, levels=c("no-change","increase", "decrease"))

ggplot(gene_ia_norm_G, aes(interaction_change, log2(expr_G/expr_ST), fill=interaction_change))+
  geom_boxplot(size=1) + 
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 0, color="red", linetype="dashed", size=1)+
  xlab("intra-gene interaction (ST vs. G)") + ylab("gene expr: log2(G/ST)") +
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))+
  scale_fill_manual( values = c("gray","red","blue"))

ggsave("../figure/ST_G_Gene_interaction_expr.pdf", device = "pdf", width = 5, height = 4)

gene_ia_norm_G$log2_G_ST = log2(gene_ia_norm_G$expr_G/gene_ia_norm_G$expr_ST)
gene_ia_norm_G1 = na.omit(gene_ia_norm_G)
gene_ia_norm_G1 = gene_ia_norm_G1[!is.infinite(gene_ia_norm_G1$log2_G_ST),]
t.test(gene_ia_norm_G1[gene_ia_norm_G1$interaction_change == "no-change", "log2_G_ST"],
       gene_ia_norm_G1[gene_ia_norm_G1$interaction_change == "decrease", "log2_G_ST"])


gene_ST = gene_ia_norm_G[gene_ia_norm_G$interaction == 1, "Gene"]
gene_G = gene_ia_norm_G[gene_ia_norm_G$interaction == -1, "Gene"]

library(goseq)
library(org.Mm.eg.db)

genes = rep(0, nrow(gene_ia_norm_G))
names(genes) = rownames(gene_ia_norm_G)
genes[names(genes) %in% gene_G] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_G=goseq(pwf,"mm10","geneSymbol")
GO_G = GO_G[GO_G$ontology == "BP",][6:1,]
GO_G$term = factor(GO_G$term, levels = GO_G$term)
ggplot(GO_G, aes(term,-log10(over_represented_pvalue) )) +
  geom_bar(stat = "identity", fill="black") + coord_flip()+
  theme_minimal() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/Gene_interaction_up.GO.pdf", device = "pdf", width = 4, height = 2)


genes = rep(0, nrow(gene_ia_norm_G))
names(genes) = rownames(gene_ia_norm_G)
genes[names(genes) %in% gene_ST] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_LT=goseq(pwf,"mm10","geneSymbol")
GO_LT = GO_LT[GO_LT$ontology == "BP",][6:1,]
GO_LT$term = factor(GO_LT$term, levels = GO_LT$term)
ggplot(GO_LT, aes(term,-log10(over_represented_pvalue) )) +
  geom_bar(stat = "identity", fill="black") + coord_flip()+
  theme_minimal() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/Gene_interaction_down.GO.pdf", device = "pdf", width = 5.5, height = 2)



# CLP CMP



# cell type specifc

gene_ia_norm_sig_gene = c()
for (i in 1:nrow(gene_ia_norm)){
  t.test(gene_ia_norm[1,], c(1:3,3,4,4,5,5))$p.value
  gene_ia_norm_sig_gene = c(names(a[a >= 1]))
}
gene_ia_norm_sig = na.omit(gene_ia_norm[apply(gene_ia_norm,1,sd)>1,])
gene_ia_norm_sig = na.omit(gene_ia_norm[apply(gene_ia_norm,1,max)-apply(gene_ia_norm,1,mean) > 1,])
gene_ia_norm_sig=gene_ia_norm_sig[rowMeans(gene_ia_norm_sig)>1,]
gene_ia_norm_sig[gene_ia_norm_sig>8] <- 8
ph = pheatmap(gene_ia_norm_sig, scale="none", cluster_cols = F)

expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = as.matrix(expr[,c(12,7,3,2,17,18,32,25)])
colnames(expr) = c("LT","ST","CMP","CLP","MEP","GMP","MK","G")

expr.ia = merge(expr, gene_ia_norm, by="row.names")
rownames(expr.ia) = expr.ia$Row.names
expr.ia = expr.ia[,-1]
expr.ia_sig = na.omit(expr.ia[apply(expr.ia[,1:8],1,max)-apply(expr.ia[,1:8],1,mean) > 30,])

breaksList = seq(-1.2, 1.2, by = 0.1)

png("../figure/Gene_interaction_expr_cellType.png",1000,1000, res=300)
ph = pheatmap((expr.ia_sig[,1:8]), scale="row", cluster_cols = F,clustering_method = "ward.D",show_rownames = F,border_color = NA,
              col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),breaks = breaksList)
dev.off()
png("../figure/Gene_interaction_IA_cellType.png",800,1000, res=300)
pheatmap((expr.ia_sig[ph$tree_row$labels[ph$tree_row$order],9:16]), scale="row", show_rownames = F,border_color = NA,
         col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),breaks = breaksList,
         cluster_cols = F,cluster_rows =F)
dev.off()

pheatmap((expr.ia_sig[sample(1:nrow(expr.ia_sig)),9:16]), scale="row", show_rownames = F,border_color = NA,
         col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),breaks = breaksList,
         cluster_cols = F,cluster_rows =F)

cor(as.numeric(as.matrix(expr.ia_sig[,1:8])), as.numeric(as.matrix(expr.ia_sig[ph$tree_row$labels[ph$tree_row$order],9:16])),method="spearman")

sample_c = c()
for(i in 1:1000){
  c = cor(as.numeric(as.matrix(expr.ia_sig[sample(1:nrow(expr.ia_sig)),1:8])), as.numeric(as.matrix(expr.ia_sig[sample(1:nrow(expr.ia_sig)),9:16])),method="spearman")
  sample_c = c(sample_c,c)
}

plot(density(sample_c))
abline(v=0.0245571, col = "red")



### ### ### ### ### ### ### ### ### ### ### ### ### 
### Gene associated domain
### ### ### ### ### ### ### ### ### ### ### ### ### 


gene_ia = read.table("../Juicer/merge/LT.gene_IA", stringsAsFactors = F)
gene_ia = gene_ia[!duplicated(gene_ia$V1),]
rownames(gene_ia) = gene_ia$V1
for (i in Sys.glob("../Juicer/merge/*gene_IA")){
  ia = read.table(i, stringsAsFactors = F)
  ia = ia[!duplicated(ia$V1),]
  rownames(ia) = ia$V1
  sample = strsplit(basename(i), "[.]")[[1]][1]
  gene_ia[[sample]] = 2*ia[rownames(gene_ia), 2] / (ia[rownames(gene_ia), 3]+ia[rownames(gene_ia), 4])
  #gene_ia[[sample]] = ia[rownames(gene_ia), 2] 
}
gene_ia = gene_ia[,c("LT","ST","MPP","CMP","CLP","MEP","GMP","MKP","MK","G")]
gene_ia = na.omit(gene_ia)
boxplot(log2(gene_ia+1))

# remove gene in bad regions
badGenes = read.table("/lustre/user/liclab/zhangc/proj/blood3D_3/homer_merge/badRegions.genes.bed", stringsAsFactors = F)$V1
gene_ia = gene_ia[!rownames(gene_ia) %in% badGenes, ]
gene_ia = gene_ia[rowMeans(gene_ia)!=0,]

# normalize by colsum
#col_Sum = apply(gene_ia, 2, sum)/nrow(gene_ia)
#col_Sum = col_Sum/col_Sum[1]
#gene_ia_norm = as.data.frame(t(apply(gene_ia, 1, function(x) x/col_Sum)))
#boxplot(log2(gene_ia_norm+1))
gene_ia_norm = gene_ia
genes = rownames(gene_ia_norm)

breaksList = seq(-2, 2, by = 0.1)
gene_ia_norm_sig = na.omit(gene_ia_norm[apply(gene_ia_norm,1,max)-apply(gene_ia_norm,1,mean) > 0.4,])

pheatmap((gene_ia_norm_sig), scale="row", cluster_cols = F,clustering_method = "ward.D",show_rownames = T,border_color = NA,
              col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),breaks = breaksList)

head(gene_ia_norm[order(gene_ia_norm$LT, decreasing = T),])

gene_ia_norm_LT = sort(gene_ia_norm$LT)
plot(1:nrow(gene_ia_norm),gene_ia_norm_LT)

dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  d
} 

plot_GAD_define <- function(sample){
  gene_ia_norm_sample = gene_ia_norm[order(gene_ia_norm[[sample]]), sample]
  names(gene_ia_norm_sample) = rownames(gene_ia_norm)[order(gene_ia_norm[[sample]])]
  d = c()
  for(i in 1:nrow(gene_ia_norm)){
    d = c(d, dist2d(c(i,gene_ia_norm_sample[i]), c(0,min(gene_ia_norm_sample)), c(nrow(gene_ia_norm),max(gene_ia_norm_sample))))
  }
  dd = which(d==max(d))
  x = 1:length(gene_ia_norm_sample)
  gene_ia_norm_sample1 = as.data.frame(cbind(x,gene_ia_norm_sample))
  colnames(gene_ia_norm_sample1) = c("sample","GAD")
  g = ggplot(gene_ia_norm_sample1, aes(sample, GAD))+
    geom_point(size=2,aes(color=ifelse( gene_ia_norm_sample1$sample>=dd, "darkred", "gray60"))) +
    scale_color_identity()+
    theme_classic(base_size = 20) +
    xlab("genes ranked by GAD score") + ylab("GAD score") + ggtitle(sample)+
    theme(legend.position = "none",
          axis.line = element_line(size = 1, colour = "black"),
          axis.text = element_text(colour = "black")) 
  gene_ia_norm_sample1$Gene = rownames(gene_ia_norm_sample1)
  g + geom_text_repel(
    data = tail(gene_ia_norm_sample1),
    aes(label = Gene),
    size = 6
  )
  #ggsave(paste0("../figure/GAD/GAD_",sample,"_define.boxplot.pdf"),device = "pdf",width = 5, height = 4.5)
  #return(g)
  return(names(gene_ia_norm_sample[dd:length(gene_ia_norm_sample)]))
}

ST_GAD = plot_GAD_define("ST")


library("gridExtra")
pl <- lapply(c("LT","ST","MPP","CMP","CLP","MEP","GMP","MKP","MK","G"), plot_GAD_define)
marrangeGrob(pl, nrow=2, ncol=5)


# GAD GWAS
GAD_GWAS = data.frame(row.names = rownames(gene_ia_norm),
                      GAD = rep("F",2202), stringsAsFactors = F)

GAD_GWAS[rownames(GAD_GWAS) %in% ST_GAD, 1] <- "T"
gene_info = read.table("/lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.protein_coding.bed", stringsAsFactors = F)
GAD_GWAS = cbind( gene_info[match(rownames(GAD_GWAS),gene_info$V4),],GAD_GWAS)
#write.table(GAD_GWAS, "../data/ST_GAD_GWAS.txt", row.names = F, quote = F, sep="\t", col.names = F)

GAD_GWAS = read.table("../data/ST_GAD_GWAS_mm10.bed", stringsAsFactors = F)
GAD_GWAS = unique(GAD_GWAS[,c(1:12,15)])
library(dplyr)
GAD_GWAS_count = GAD_GWAS %>%
  group_by(V4,V7) %>%
  summarize(count = n())
GAD_GWAS_count = as.data.frame(GAD_GWAS_count)

GAD_GWAS_count1 = data.frame(row.names = rownames(gene_ia_norm),
                             GAD = rep("F",2202), stringsAsFactors = F)
GAD_GWAS_count1[rownames(GAD_GWAS_count1) %in% ST_GAD, 1] <- "T"
GAD_GWAS_count1$count = 0

GAD_GWAS_count1[GAD_GWAS_count$V4, 2] <- GAD_GWAS_count$count



ggplot(GAD_GWAS_count1, aes(GAD, log2(count+1), fill=GAD)) +
  geom_boxplot() +
  #stat_summary(fun.y=mean, geom="point", shape=20, size=14, color="red", fill="red")+
  theme_classic(base_size = 20) +
  xlab("") + ylab("log2(#GWAS_SNP+1)") + 
  scale_fill_aaas() +
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) 
ggsave("../figure/GAD/GAD_GWAS_SNP_number.pdf",device = "pdf",width = 4, height = 4)

wilcox.test(GAD_GWAS_count1[GAD_GWAS_count1$GAD == "T", 2], GAD_GWAS_count1[GAD_GWAS_count1$GAD == "F", 2])
t.test(GAD_GWAS_count1[GAD_GWAS_count1$GAD == "T", 2], GAD_GWAS_count1[GAD_GWAS_count1$GAD == "F", 2])

GAD_GWAS_count2 = GAD_GWAS %>%
  group_by(V7,V15) %>%
  summarize(count = n())
GAD_GWAS_count2 = as.data.frame(GAD_GWAS_count2)
GAD_GWAS_count2$T1 = ifelse(GAD_GWAS_count2$V7, sum(GAD_GWAS_count2[GAD_GWAS_count2$V7==T, "count"]),  sum(GAD_GWAS_count2[GAD_GWAS_count2$V7==F, "count"]))
colnames(GAD_GWAS_count2)  = c("GAD", "Var1","Freq","T1")
GAD_GWAS_count2$ratio = 100 * GAD_GWAS_count2$Freq / GAD_GWAS_count2$T1

select_term = c("Body_mass_index","Obesity-related_traits","Neuroticism","White_blood_cell_count","Blood_protein_levels",
                "Mean_platelet_volume","Red_cell_distribution_width"," Red_blood_cell_count","Eosinophil_counts")
GAD_GWAS_count3 = GAD_GWAS_count2[GAD_GWAS_count2$Var1 %in% select_term,]
GAD_GWAS_count3$Var1 = factor(GAD_GWAS_count3$Var1, levels = rev(select_term))
GAD_GWAS_count3$GAD = ifelse(GAD_GWAS_count3$GAD, "GAD", "non GAD")
ggplot(GAD_GWAS_count3, aes(Var1, ratio, fill=(GAD)))+
  geom_bar(position = "dodge", stat = 'identity') +
  coord_flip() +
  xlab("GWAS SNPs ratio")+ylab("")+
  theme_classic(base_size=15)+
  scale_fill_aaas()
ggsave("../figure/GAD/GAD_GWAS_SNP_ratio.pdf",device = "pdf",width = 7, height = 3.5)






# RNA-seq
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = as.matrix(expr[,c(12,7,3,2,17,18,32,25)])
colnames(expr) = c("LT","ST","CMP","CLP","MEP","GMP","MK","G")

ST_GAD = plot_GAD_define("ST")
ST_GAD_expr = as.data.frame( expr[genes,"ST"])
ST_GAD_expr$gene = rownames(ST_GAD_expr)
ST_GAD_expr$GAD = ifelse(ST_GAD_expr$gene %in% ST_GAD,"within GAD","without GAD")
colnames(ST_GAD_expr) = c("expr","gene","GAD")
ggplot(ST_GAD_expr, aes(GAD, log(expr+1), fill=ifelse( GAD=="within GAD", "red3", "gray") ))+
  geom_boxplot(size=1) +
  scale_fill_identity()+
  theme_classic(base_size = 20) +
  xlab("") + ylab("log(TMP+1)") + 
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) 
ggsave("../figure/GAD/GAD_ST_expr.boxplot.pdf",device = "pdf",width = 4, height = 3.5)
t.test(ST_GAD_expr[ST_GAD_expr$GAD == "within GAD", 1], ST_GAD_expr[ST_GAD_expr$GAD != "within GAD", 1])

write.csv(ST_GAD_expr, "../data/ST_GAD_expr.tab", row.names = T, quote = F)


gene_length = read.table("/lustre/user/liclab/zhangc/Taolab/guan/rna-seq/reference/gene_position.protein_coding.bed")

ST_GAD_expr = ST_GAD_expr[order(ST_GAD_expr$GAD),]
write.csv(gene_length[match(rownames(ST_GAD_expr), gene_length$V4),], "../data/ST_GAD.tab", row.names = T, quote = F)


expr.ia = merge(expr, gene_ia_norm, by="row.names")
rownames(expr.ia) = expr.ia$Row.names
expr.ia = expr.ia[,-1]
expr.ia_sig = na.omit(expr.ia[apply(expr.ia[,1:8],1,max)-apply(expr.ia[,1:8],1,mean) > 30 & apply(expr.ia[,1:8],1,max)/apply(expr.ia[,1:8],1,mean) > 1.5,])

breaksList = seq(-1.5, 1.5, by = 0.1)

png("../figure/Gene_interaction_expr_cellType.png",1000,1000, res=300)
ph = pheatmap((expr.ia_sig[,1:8]), scale="row", cluster_cols = F,clustering_method = "ward.D",show_rownames = F,border_color = NA,
              col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),breaks = breaksList)
dev.off()
png("../figure/Gene_interaction_IA_cellType.png",800,1000, res=300)
pheatmap((expr.ia_sig[ph$tree_row$labels[ph$tree_row$order],9:16]), scale="row", show_rownames = F,border_color = NA,
         col=colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(length(breaksList)),breaks = breaksList,
         cluster_cols = F,cluster_rows =F)
dev.off()


# gene expression Vs. GAD

gene_ia_norm_expr1 = expr.ia[,c(2,10)]

gene_ia_norm_expr1$expr = "0-1"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.x <= 1, 3] = "0-1"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.x > 1 & gene_ia_norm_expr1$ST.x <= 10 , 3] = "1-10"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.x > 10 & gene_ia_norm_expr1$ST.x <= 30 , 3] = "10-30"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.x > 30 & gene_ia_norm_expr1$ST.x <= 80 , 3] = "30-80"
gene_ia_norm_expr1[gene_ia_norm_expr1$ST.x > 80 , 3] = ">80"
gene_ia_norm_expr1$expr = factor(gene_ia_norm_expr1$expr, levels=c("0-1","1-10","10-30","30-80",">80"))
write.csv(gene_ia_norm_expr1, "../data/ST_gene_expr_interation.tab", row.names = T, quote = F)

ggplot(gene_ia_norm_expr1, aes(expr, ST.y, fill=expr))+
  geom_boxplot(size=1) + 
  theme_classic(base_size = 20) +
  xlab("gene expression (TPM)") + ylab("GAD score") +
  theme(legend.position = "none")
ggsave("../figure/ST_Gene_interaction_expr.pdf", device = "pdf", width = 6, height = 5)


# G
gene_ia_norm_G = data.frame(gene_ia_norm[,c(2,8)])
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(6,9,25,27)]
#expr = expr[rowMeans(expr)>=1,]
gene_ia_norm_G$expr_ST = rowMeans(expr[rownames(gene_ia_norm_G),1:2])
gene_ia_norm_G$expr_G =  rowMeans(expr[rownames(gene_ia_norm_G),3:4])
gene_ia_norm_G$interaction = (gene_ia_norm_G$ST) - (gene_ia_norm_G$G)
gene_ia_norm_G[gene_ia_norm_G$interaction > 0.5, "interaction"] <- 0.5
gene_ia_norm_G[gene_ia_norm_G$interaction < -0.5, "interaction"] <- -0.5
gene_ia_norm_G$Gene = rownames(gene_ia_norm_G)
g=ggplot(gene_ia_norm_G, aes(log(expr_ST+1), log(expr_G+1), name=Gene))+
  geom_point(aes(col=interaction), shape=3, stroke=1.5) + xlim(0,6)+ylim(0,6)+
  xlab("log(ST expr.)")+ylab("log(G expr.)")+
  geom_abline(slope = 1, color="red", linetype="dashed", size=1)+
  theme_classic(base_size=15)+
  scale_colour_gradient2()


ggplotly(g) 


library(ggrepel)
g + geom_text_repel(
  data = gene_ia_norm_G[c("Syne1","App", "Lrrk2", "Lyst","Prkcb","Lyn","Abca13",
                          "Angpt1","Fut8","Eya2","Meis1", "Msi2","Erg", "Rftn1"),],
  aes(label = Gene),
  size = 5
)
ggsave("../figure/Gene_interaction_expr.pdf", device = "pdf", width = 6, height = 5)


# 
gene_ia_norm_G$interaction_change = "no-change"
gene_ia_norm_G[gene_ia_norm_G$interaction >=0.2, "interaction_change"] <- "decrease"
gene_ia_norm_G[gene_ia_norm_G$interaction <= -0.2, "interaction_change"] <- "increase"
gene_ia_norm_G$interaction_change = factor(gene_ia_norm_G$interaction_change, levels=c("no-change","increase", "decrease"))

ggplot(gene_ia_norm_G, aes(interaction_change, log2(expr_G/expr_ST), fill=interaction_change))+
  geom_boxplot(size=1) + 
  theme_classic(base_size = 20) +
  geom_hline(yintercept = 0, color="red", linetype="dashed", size=1)+
  xlab("intra-gene interaction (ST vs. G)") + ylab("gene expr: log2(G/ST)") +
  theme(legend.position = "none",
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))+
  scale_fill_manual( values = c("gray","red","blue"))

ggsave("../figure/ST_G_Gene_interaction_expr.pdf", device = "pdf", width = 5, height = 4)

gene_ia_norm_G$log2_G_ST = log2(gene_ia_norm_G$expr_G/gene_ia_norm_G$expr_ST)
gene_ia_norm_G1 = na.omit(gene_ia_norm_G)
gene_ia_norm_G1 = gene_ia_norm_G1[!is.infinite(gene_ia_norm_G1$log2_G_ST),]
t.test(gene_ia_norm_G1[gene_ia_norm_G1$interaction_change == "no-change", "log2_G_ST"],
       gene_ia_norm_G1[gene_ia_norm_G1$interaction_change == "decrease", "log2_G_ST"])


gene_ST = gene_ia_norm_G[gene_ia_norm_G$interaction == 0.5, "Gene"]
gene_G = gene_ia_norm_G[gene_ia_norm_G$interaction == -0.5, "Gene"]

library(goseq)
library(org.Mm.eg.db)

genes = rep(0, nrow(gene_ia_norm_G))
names(genes) = rownames(gene_ia_norm_G)
genes[names(genes) %in% gene_G] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_G=goseq(pwf,"mm10","geneSymbol")
GO_G = GO_G[GO_G$ontology == "BP",][6:1,]
GO_G$term = factor(GO_G$term, levels = GO_G$term)
ggplot(GO_G, aes(term,-log10(over_represented_pvalue) )) +
  geom_bar(stat = "identity", fill="black") + coord_flip()+
  theme_minimal() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/Gene_interaction_up.GO.pdf", device = "pdf", width = 4, height = 2)


genes = rep(0, nrow(gene_ia_norm_G))
names(genes) = rownames(gene_ia_norm_G)
genes[names(genes) %in% gene_ST] = 1
pwf=nullp(genes,"mm10","geneSymbol", plot.fit = FALSE)
GO_LT=goseq(pwf,"mm10","geneSymbol")
GO_LT = GO_LT[GO_LT$ontology == "BP",][6:1,]
GO_LT$term = factor(GO_LT$term, levels = GO_LT$term)
ggplot(GO_LT, aes(term,-log10(over_represented_pvalue) )) +
  geom_bar(stat = "identity", fill="black") + coord_flip()+
  theme_minimal() + xlab("") + ylab("-log10(p-value)") +
  theme(axis.text = element_text(size=10, color="black"))+
  scale_colour_brewer(palette="Dark2")
ggsave("../figure/Gene_interaction_down.GO.pdf", device = "pdf", width = 5.5, height = 2)




# CMP vs CLP

gene_ia_norm_CMP = data.frame(gene_ia_norm[,c(3,4)])
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(20,3,2,14)]
#expr = expr[rowMeans(expr)>=1,]
gene_ia_norm_CMP$expr_CMP = rowMeans(expr[rownames(gene_ia_norm_CMP),1:2])
gene_ia_norm_CMP$expr_CLP =  rowMeans(expr[rownames(gene_ia_norm_CMP),3:4])
gene_ia_norm_CMP$interaction = (gene_ia_norm_CMP$CMP) - (gene_ia_norm_CMP$CLP)
gene_ia_norm_CMP[gene_ia_norm_CMP$interaction > 0.5, "interaction"] <- 0.5
gene_ia_norm_CMP[gene_ia_norm_CMP$interaction < -0.5, "interaction"] <- -0.5
gene_ia_norm_CMP$Gene = rownames(gene_ia_norm_CMP)
g=ggplot(gene_ia_norm_CMP, aes(log(expr_CMP+1), log(expr_CLP+1), name=Gene))+
  geom_point(aes(col=interaction), shape=3, stroke=1.5) + xlim(0,6)+ylim(0,6)+
  xlab("log(ST expr.)")+ylab("log(G expr.)")+
  geom_abline(slope = 1, color="red", linetype="dashed", size=1)+
  theme_classic(base_size=15)+
  scale_colour_gradient2()
ggplotly(g) 


# MEP ST
gene_ia_norm_CMP = data.frame(gene_ia_norm[,c(2,5)])
expr = read.table("../../blood3D/RNA_seq/stringtie/gene.TPM.tab", row.names = 1, header = T)
expr.mat = normalize.quantiles(as.matrix(expr))
colnames(expr.mat) = colnames(expr)
rownames(expr.mat) = rownames(expr)
expr = as.data.frame(expr.mat)
expr = expr[,c(6,9,17,5)]
#expr = expr[rowMeans(expr)>=1,]
gene_ia_norm_CMP$expr_CMP = rowMeans(expr[rownames(gene_ia_norm_CMP),1:2])
gene_ia_norm_CMP$expr_CLP =  rowMeans(expr[rownames(gene_ia_norm_CMP),3:4])
gene_ia_norm_CMP$interaction = (gene_ia_norm_CMP$ST) - (gene_ia_norm_CMP$MEP)
gene_ia_norm_CMP[gene_ia_norm_CMP$interaction > 0.5, "interaction"] <- 0.5
gene_ia_norm_CMP[gene_ia_norm_CMP$interaction < -0.5, "interaction"] <- -0.5
gene_ia_norm_CMP$Gene = rownames(gene_ia_norm_CMP)
g=ggplot(gene_ia_norm_CMP, aes(log(expr_CMP+1), log(expr_CLP+1), name=Gene))+
  geom_point(aes(col=interaction), shape=3, stroke=1.5) + xlim(0,6)+ylim(0,6)+
  xlab("log(ST expr.)")+ylab("log(MEP expr.)")+
  geom_abline(slope = 1, color="red", linetype="dashed", size=1)+
  theme_classic(base_size=15)+
  scale_colour_gradient2()
ggplotly(g) 

g + geom_text_repel(
  data = gene_ia_norm_CMP[c("Slc25c21","Sptb", "Ank1", "Myh10","Atp8a1",
                          "Fli1","Gcnt2","Stat4","Meis1", "Cdk19"),],
  aes(label = Gene),
  size = 5
)




### GAD for human cell


dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  d
} 

plot_GAD_define_hg19 <- function(file){
  gene_ia_hg19 =  read.table(file, stringsAsFactors = F)
  gene_ia_hg19 = gene_ia_hg19[!duplicated(gene_ia_hg19$V1),]
  rownames(gene_ia_hg19) = gene_ia_hg19$V1
  gene_ia_hg19 = gene_ia_hg19[rowMeans(gene_ia_hg19[,2:4]) !=0,]
  gene_ia_hg19$GAD = 2*gene_ia_hg19[, 2] / (gene_ia_hg19[, 3]+gene_ia_hg19[, 4])
  gene_ia_hg19 = gene_ia_hg19[!rownames(gene_ia_hg19) %in% badGenes, ]
  gene_ia_hg19 = gene_ia_hg19[gene_ia_hg19$GAD >=1.5, 1]
  return(gene_ia_hg19)
}


HiC101_GAD = plot_GAD_define_hg19("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/HIC010/SRR1658581_allValidPairs.gene_IA")
GM_GAD = plot_GAD_define_hg19("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_bulk/GM12878.allValidPairs.gene_IA")
GM_lu_GAD = plot_GAD_define_hg19("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_bulk_lumeng/bulk_allValidPairs.gene_IA")
Hela_GAD = plot_GAD_define_hg19("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/Hela_bulk/Hela.allValidPairs.gene_IA")
GM1k_GAD = plot_GAD_define_hg19("/lustre/user/liclab/zhangc/proj/blood3D/HiC/snHiC/hic/GM12878_1k_rep1.hic/GM12878_11.gene_IA")


