library(gplots) 
library(RColorBrewer)

rm(list = ls(all.names = TRUE)) # clear workspace
dev.off() # clear all plots
con = odbcDriverConnect('driver={SQL Server};server=localhost\\SQLEXPRESS;database=MouseDB;trusted_connection=true')
setwd("C:/Users/hadadi/Dropbox/UniGen/manuscripts/All_tissue/")

#########################FIGURE1#################
### 6 BAT studies
#Connecting to the SQL database to retrieve the data
data <- sqlQuery(con,  "  select * from omics_results_FC_subset 
                     where tissue like '%BAT%' and name_for_PCA like '%cold%' and 
                      PValue < 0.06 and
                  omics_name not like '%Fajas%'  ")


logFC= dcast (data, gene_id_name ~ comparison, value.var='logFC' )
Pvalue= dcast (data, gene_id_name ~ comparison, value.var='PValue' )

#LogFC for all 6 studies
write.xlsx(logFC, "BAT_6studies_LogFC_sig.xlsx", row.names = FALSE)

####Plot jitter (FIGURE 1 B & C)
setwd("C:/Users/hadadi/Dropbox/UniGen/manuscripts/All_tissue/")
data_markers = read_excel("C:/Users/hadadi/Dropbox/UniGen/manuscripts/All_tissue/BAT_6studies_LogFC_sig.xlsx", sheet=2)

markers_melt =  melt (data_markers)

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),plot.title = element_text(size=20),
             axis.title.x = element_blank(),axis.title.y = element_text(size=12),axis.text.x=element_text(colour="black",size=10,angle = 90,hjust = 1),
             axis.text.y=element_text(colour="black",size=11),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

ggplot(markers_melt, aes(x = reorder(Symbol, -value, na.rm = TRUE),y=value)) + geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = variable),size = 2.8,alpha = 0.8, shape = 16,position=position_jitter(0.2)) + scale_color_manual(values = viridis(7)) +
  theme_bw() + theme(legend.title = element_blank(),legend.position="top") + theme

############################FIGURE2#########################
#Connecting to the SQL database to retrieve the data
data_raw <- sqlQuery(con,  "   select * from 
	(select gene_id_name, sample_name , raw_value
	from test_for_PCA
	 ) as s
pivot (
	avg(raw_value) for sample_name in (
	Spleen_RT_rep2,	BMC2,	Hypothalamus_CE_rep2,	Hypothalamus_RT_rep1,	BMR5,	Spleen_RT_rep3,	BAT_CE_rep2,	Brain_RT_rep1,	
	Liver_CE_rep1,	Hypothalamus_CE_rep3,	SAT_RT_rep2,	Brain_CE_rep2,	BAT_RT_rep2,	BAT_RT_rep3,	BAT_CE_rep3,	BMC6,	
	SAT_RT_rep1,	VAT_RT_rep1,	VAT_CE_rep2,	Hypothalamus_RT_rep2,	Hypothalamus_CE_rep1,	Liver_RT_rep2,	Ileum_CE_rep3,	
	SAT_CE_rep1,	BMC1,	SAT_CE_rep3,	Hypothalamus_RT_rep3,	Spleen_CE_rep1,	BMC4,	Ileum_CE_rep2,	Ileum_RT_rep3,	Liver_RT_rep1,	
	Ileum_RT_rep2,	SAT_RT_rep3,	BMR1,	BMR3,	VAT_RT_rep3,	Brain_CE_rep1,	Spleen_CE_rep2,	BMC5,	VAT_CE_rep3,	BMR6,	VAT_RT_rep2,	
	Liver_RT_rep3,	BAT_CE_rep1,	Spleen_CE_rep3,	Liver_CE_rep2,	SAT_CE_rep2,	Brain_CE_rep3,	Spleen_RT_rep1,	Ileum_CE_rep1,	BAT_RT_rep1,	
	Brain_RT_rep3,	Liver_CE_rep3,	BMC3,	Ileum_RT_rep1,	VAT_CE_rep1,	BMR2,	BMR4, D0RT6,
D0CE3,D0RT3,D0CE2,D0RT1,D0CE5,D0RT2,D0RT4,D0CE4,D0CE1
)
	) as t
order by 1;
 ")

x=as.matrix(data_raw[,-1])
rownames(x)=data_raw$gene_id_name

colnames(x) = gsub(colnames(x), pattern = "D0CE", replacement = "Spinalcord_CE_Rep")
colnames(x) = gsub(colnames(x), pattern = "D0RT", replacement = "Spinalcord_RT_Rep")
colnames(x) = gsub(colnames(x), pattern = "BMC", replacement = "Bonemarrow_CE_Rep")
colnames(x) = gsub(colnames(x), pattern = "BMR", replacement = "Bonemarrow_RT_Rep")


x2 = x[,grep("Spl", colnames(x))]
x2 = na.omit(x)
x2[is.na(x2)] <- 0
x2 = x2[!(apply(x2, 1, function(y) any(y == 0))),]

cts=x2[rowSums(x2[, ] > 20) != 0, ]

cond= factor(c(gsub(colnames(cts), pattern = "[[:digit:]]", replacement = "")))
y <- DGEList(cts, group= cond,  genes= row.names(cts) )
y <- calcNormFactors(y)
options(digits=3)
y= estimateCommonDisp(y,robust=TRUE)
CPM <- as.data.frame(cpm(y, normalized.lib.sizes = T, log = T))


pca.normcount <- prcomp(t(CPM), center=TRUE, scale=FALSE)
a = factor(c(rownames(pca.normcount$x)))

Cold_RT <- as.data.frame (do.call(rbind, strsplit(as.character(a), '_')))
Cold_RT$V2 = factor(c(gsub(Cold_RT$V2, pattern = "CE", replacement = "Cold")))

pca_data<- data.frame(Tissue = str_sub(a, end=-9), Condition= Cold_RT$V2, PC1 = pca.normcount$x[,1], PC2=pca.normcount$x[,2])
summary(pca.normcount)

percentage <- round((pca.normcount$sdev^2) / (sum(pca.normcount$sdev^2))*100)
percentage <- paste( colnames(pca_data)[3:4], "(", paste( as.character(percentage), "%", ")", sep="") )


#PCA for each tissue
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),legend.title = element_blank(),legend.position = "none",
             axis.text.x=element_text(colour="black",size=2),axis.text.y=element_text(colour="black",size=2),axis.title.x = element_text(size=3),
             axis.title.y = element_text(size=3))

ggplot(pca_data,aes(x=PC1, y=PC2, shape= Condition)) + geom_point(size=1) +theme  +    scale_shape_manual(values=c(0,2)) +
  xlab(percentage[1]) + ylab(percentage[2]) 


ggsave("PCA_new_Spleen.eps", width = 3.5, height = 3.5, units = "cm",scale = 0.9)
ggsave("PCA_new_Spleen.pdf", width = 3.5, height = 3.5, units = "cm",scale = 0.9)

#PCA all tissues
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),plot.title = element_text(size=12),axis.title.x = element_text(size=14),
             axis.title.y = element_text(size=14),axis.text.x=element_text(colour="black",size=12),axis.text.y=element_text(colour="black",size=12),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),legend.direction = "vertical", legend.box = "vertical", 
             legend.background=element_blank(),legend.key=element_blank(),
             legend.title = element_text(colour="Black", size=14, face="bold"))


ggplot(pca_data,aes(x=PC1, y=PC2, col=Tissue, shape= Condition)) + geom_point(size=4) +theme+  xlab(percentage[1]) + ylab(percentage[2]) + 
  scale_shape_manual(values=c(0,2)) + scale_color_manual(values = c(rainbow(11))) +
  geom_text_repel(aes(label = Tissue), hjust=0, vjust=0,size=3,parse=TRUE)
############################FIGURE3#########################
#Connecting to the SQL database to retrieve the data

data  <- sqlQuery(con, "select * from 
	 FC_original_wo_Brain_RTrep2 where  logCPM > 0 ")

data2= dcast (data, gene_id_name ~ comparison, value.var='logFC' )
data3= dcast (data, gene_id_name ~ comparison, value.var='PValue' )

data2[is.na(data2)] <- 0
data3[is.na(data3)] <- 0

colnames(data3) = gsub(colnames(data3), pattern = "D0CE", replacement = "Spinalcord_CE_Rep")
colnames(data3) = gsub(colnames(data3), pattern = "D0RT", replacement = "Spinalcord_RT_Rep")
colnames(data3) = gsub(colnames(data3), pattern = "BMC", replacement = "Bonemarrow_CE_Rep")
colnames(data3) = gsub(colnames(data3), pattern = "BMR", replacement = "Bonemarrow_RT_Rep")
colnames(data3) = gsub(colnames(data3), pattern = " VS ", replacement = "_VS_")

colnames(data2) = gsub(colnames(data2), pattern = "D0CE", replacement = "Spinalcord_CE_Rep")
colnames(data2) = gsub(colnames(data2), pattern = "D0RT", replacement = "Spinalcord_RT_Rep")
colnames(data2) = gsub(colnames(data2), pattern = "BMC", replacement = "Bonemarrow_CE_Rep")
colnames(data2) = gsub(colnames(data2), pattern = "BMR", replacement = "Bonemarrow_RT_Rep")
colnames(data2) = gsub(colnames(data2), pattern = " VS ", replacement = "_VS_")

data_FC= data2[c( "Brain_CE_rep_VS_Brain_RT_rep","Spinalcord_CE_Rep_VS_Spinalcord_RT_Rep", "Hypothalamus_CE_rep_VS_Hypothalamus_RT_rep",
                  "Bonemarrow_CE_Rep_VS_Bonemarrow_RT_Rep","Spleen_CE_rep_VS_Spleen_RT_rep", "Ileum_CE_rep_VS_Ileum_RT_rep","Liver_CE_rep_VS_Liver_RT_rep", 
                  "BAT_CE_rep_VS_BAT_RT_rep",  
                  "VAT_CE_rep_VS_VAT_RT_rep", "SAT_CE_rep_VS_SAT_RT_rep")]


data_Pvalue= data3[c( "Brain_CE_rep_VS_Brain_RT_rep","Spinalcord_CE_Rep_VS_Spinalcord_RT_Rep", "Hypothalamus_CE_rep_VS_Hypothalamus_RT_rep",
                      "Bonemarrow_CE_Rep_VS_Bonemarrow_RT_Rep","Spleen_CE_rep_VS_Spleen_RT_rep", "Ileum_CE_rep_VS_Ileum_RT_rep","Liver_CE_rep_VS_Liver_RT_rep", 
                      "BAT_CE_rep_VS_BAT_RT_rep",  
                      "VAT_CE_rep_VS_VAT_RT_rep", "SAT_CE_rep_VS_SAT_RT_rep")]

data_final=melt(data_FC)
Pvalue = melt(data_Pvalue)
data_final$pval = Pvalue$value

data_final= data_final[!apply(data_final[,2:3] == 0, 1, FUN = any, na.rm = TRUE),]
rownames(data_final) <- NULL

data_final$labels=round(data_final$value,3)

down_label = data %>%
  group_by(comparison) %>%
  top_n(n = -20, wt = logFC)

top_label = data %>%
  group_by(comparison) %>%
  top_n(n = 20, wt = logFC)

write.csv(down_label, file ="down_label.csv")
write.csv(top_label, file ="top_label.csv")


theme<-theme(panel.background = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
             strip.background=element_blank(),axis.text.x=element_text(colour="black",size=16),axis.text.y=element_text(colour="black",size=16),
             axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none")


grey= unique(subset(data_final,data_final$pval>=0.05)[c("variable", "labels")])
color= unique(subset(data_final,data_final$pval<0.05)[c("variable", "labels")])
grey$direction[grey$labels<0]="down"
grey$direction[grey$labels>=0]="up"

color$direction[color$labels<0]="down"
color$direction[color$labels>=0]="up"

cbPalette = c("#BF4343", "#4E9571")

ggplot(data = grey , aes(x = labels, y = variable)) +
  geom_jitter(color = "darkgrey", size = 1.25,height = 0.3) +
  geom_jitter(
    data = color,
    aes(x = labels, y = variable,color=direction), size = 1.3,height = 0.3) + scale_color_manual(values = cbPalette) +
  
  scale_y_discrete(labels=c("VAT_CE_rep_VS_VAT_RT_rep" = "VAT", "Spleen_CE_rep_VS_Spleen_RT_rep" = "Spleen","Spinalcord_CE_Rep_VS_Spinalcord_RT_Rep" = "Spinal cord",
                            "SAT_CE_rep_VS_SAT_RT_rep" = "SAT", "Liver_CE_rep_VS_Liver_RT_rep" = "Liver","Ileum_CE_rep_VS_Ileum_RT_rep" = "leum",
                            "Hypothalamus_CE_rep_VS_Hypothalamus_RT_rep" = "Hypothalamus","Brain_CE_rep_VS_Brain_RT_rep" = "Brain",
                            "Bonemarrow_CE_Rep_VS_Bonemarrow_RT_Rep" = "Bonemarrow","BAT_CE_rep_VS_BAT_RT_rep" = "BAT")) + 
  
  scale_x_continuous(breaks = pretty(data_final$value, n = 10)) +  labs(title= "", x= "\nLog2FC", y= "") + coord_cartesian(xlim = c(-10, 10)) +
  theme(
    axis.title.x = element_text( size=16 ),
    axis.title.y = element_text(size=16)
  ) +theme + 
  geom_vline(
    xintercept = 0.6,
    col = "red",
    linetype = "dotted",
    size = 1.2) +
  geom_vline(
    xintercept = -0.6,
    col = "red",
    linetype = "dotted",
    size = 1.2
  )

############################FIGURE4#########################
DE_pathway <- enrichPathway(unique(bitr(data_scatter$SYMBOL,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")$ENTREZID),organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 1)

A=bitr(data_scatter$SYMBOL,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
names(A)[names(A) == "ENTREZID"] <- "geneID"

B = separate_rows(as.data.frame(DE_pathway) ,geneID, convert = TRUE)
B$geneID=as.character(B$geneID)

with_symbol= B %>%
  left_join(A, by="geneID") %>%
  select(Description,geneID,SYMBOL,pvalue)

####Jitter Plot 
for_jitter= with_symbol %>%
  left_join(data_scatter, by="SYMBOL") %>%
  select(Description,geneID,SYMBOL,pvalue,VAT,BAT,SAT)

for_jitter$Description[for_jitter$Description == "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)"	] <- "Regulation of Insulin-like Growth Factor (IGF)"
for_jitter<-for_jitter[!(for_jitter$Description== "Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.") ,]



for_jitter_melt=melt(for_jitter,measure.vars = c("VAT", "BAT", "SAT"))
names(for_jitter_melt)[names(for_jitter_melt) == "value"] <- "log2FC"

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),plot.title = element_text(size=20),
             axis.title.x = element_blank(),axis.title.y = element_text(size=12),axis.text.x=element_text(colour="black",size=7,angle = 90,hjust = 1),
             axis.text.y=element_text(colour="black",size=11),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

#cbPalette=c("forestgreen", "darkred")

ggplot(for_jitter_melt, aes(x= reorder(Description,pvalue), y=log2FC, fill=variable)) + 
  geom_jitter(aes(colour = variable),size = 1.5,alpha = 0.5, position=position_dodge(0.8),stat="identity") + scale_color_manual(values = viridis(3)) +
  theme_bw() + theme(legend.title = element_blank(),legend.position="top")+ theme 


############################FIGURE7#########################

data_raw <- sqlQuery(con,  "   select * from 
	(select gene_id_name, sample_name , raw_value
	from test_for_PCA
	 ) as s
pivot (
	avg(raw_value) for sample_name in (
	Spleen_RT_rep2,	BMC2,	Hypothalamus_CE_rep2,	Hypothalamus_RT_rep1,	BMR5,	Spleen_RT_rep3,	BAT_CE_rep2,	Brain_RT_rep1,	
	Liver_CE_rep1,	Hypothalamus_CE_rep3,	SAT_RT_rep2,	Brain_CE_rep2,	BAT_RT_rep2,	BAT_RT_rep3,	BAT_CE_rep3,	BMC6,	
	SAT_RT_rep1,	VAT_RT_rep1,	VAT_CE_rep2,	Hypothalamus_RT_rep2,	Hypothalamus_CE_rep1,	Liver_RT_rep2,	Ileum_CE_rep3,	
	SAT_CE_rep1,	BMC1,	SAT_CE_rep3,	Hypothalamus_RT_rep3,	Spleen_CE_rep1,	BMC4,	Ileum_CE_rep2,	Ileum_RT_rep3,	Liver_RT_rep1,	
	Ileum_RT_rep2,	SAT_RT_rep3,	BMR1,	BMR3,	VAT_RT_rep3,	Brain_CE_rep1,	Spleen_CE_rep2,	BMC5,	VAT_CE_rep3,	BMR6,	VAT_RT_rep2,	
	Liver_RT_rep3,	BAT_CE_rep1,	Spleen_CE_rep3,	Liver_CE_rep2,	SAT_CE_rep2,	Brain_CE_rep3,	Spleen_RT_rep1,	Ileum_CE_rep1,	BAT_RT_rep1,	
	Brain_RT_rep3,	Liver_CE_rep3,	BMC3,	Ileum_RT_rep1,	VAT_CE_rep1,	BMR2,	BMR4, D0RT6,
D0CE3,D0RT3,D0CE2,D0RT1,D0CE5,D0RT2,D0RT4,D0CE4,D0CE1
)
	) as t
order by 1;
 ")

x=as.data.frame(data_raw[,-1])
rownames(x)=data_raw$gene_id_name

colnames(x) = gsub(colnames(x), pattern = "D0CE", replacement = "Spinalcord_CE_Rep")
colnames(x) = gsub(colnames(x), pattern = "D0RT", replacement = "Spinalcord_RT_Rep")
colnames(x) = gsub(colnames(x), pattern = "BMC", replacement = "Bonemarrow_CE_Rep")
colnames(x) = gsub(colnames(x), pattern = "BMR", replacement = "Bonemarrow_RT_Rep")

colnames(x) = gsub(colnames(x), pattern = "_rep", replacement = "")
colnames(x) = gsub(colnames(x), pattern = "_Rep", replacement = "")

colnames= gsub(colnames(x), pattern = "[[:digit:]]", replacement = "")


x[is.na(x)] <- 0
x = x[!(apply(x, 1, function(x) all(x == 0))),]

MAT_4ave = x[rowSums(x[, ] > 10) != 0, ]

colnames(MAT_4ave) = colnames
MAT_4ave=as.data.frame(cpm(MAT_4ave))

nms <- unique(colnames(MAT_4ave))
MAT_averaged = as.data.frame(sapply(nms, function(m)  rowMeans(MAT_4ave[names(MAT_4ave) %in% m])))

#####Ploting subsets
subset = MAT_averaged[grep("Slc25a", rownames(MAT_averaged)),]
y=as.matrix((subset))

#gene family Heatmap
setwd("C:/Users/hadadi/Dropbox/UniGen/Martina/All_tissue/Original_wo_Brain_RTrep2/Ave_Expre_Heatmaps/")

for_HM = melt(y)

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),
             axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x=element_text(colour="black",size=10,angle = 90,hjust = 1),
             axis.text.y=element_text(colour="black",size=11),axis.ticks=element_line(colour="black"))


ggplot(for_HM, aes(x = reorder(Var2), y= Var1)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 1)),size = 2.5) +
  scale_fill_gradient(low = "white", high = "red",guide = guide_colourbar(direction = "horizontal")) + theme

ggsave("Abca.pdf", width = 20, height = 15, units = "cm")


#heatmap of all gens with heatmap.2 (Supp Figure)

expressed_all = x[!(apply(x, 1, function(y) any(y <= 5))),]
colnames(expressed_all) = colnames

nms <- unique(colnames(expressed_all))
expressed_all_ave = as.data.frame(sapply(nms, function(m)  rowMeans(expressed_all[names(expressed_all) %in% m])))

expressed_all_ave = expressed_all_ave[c("Brain_RT", "Brain_CE", "Hypothalamus_CE", "Hypothalamus_RT","Bonemarrow_RT", "Bonemarrow_CE","Spinalcord_RT","Spinalcord_CE","Spleen_CE","Spleen_RT",
                                        "Ileum_RT" , "Ileum_CE"  ,"Liver_RT", "Liver_CE","VAT_RT","VAT_CE","SAT_RT","SAT_CE","BAT_RT" ,"BAT_CE")]


MAT_averaged = MAT_averaged[c("Brain_RT", "Brain_CE", "Hypothalamus_CE", "Hypothalamus_RT","Bonemarrow_RT", "Bonemarrow_CE","Spinalcord_RT","Spinalcord_CE","Spleen_CE","Spleen_RT",
                              "Ileum_RT" , "Ileum_CE"  ,"Liver_RT", "Liver_CE","VAT_RT","VAT_CE","SAT_RT","SAT_CE","BAT_RT" ,"BAT_CE")]


MAT_averaged [MAT_averaged < 5] <- 0
MAT_averaged = MAT_averaged[!(apply(MAT_averaged, 1, function(x) all(x == 0))),]

write.xlsx(MAT_averaged, "averaged_expression_10Tissues_18479Genes.xlsx",row.names= TRUE)
write.xlsx(expressed_all_ave, "averaged_expression_10Tissues_10040Genes.xlsx",row.names= TRUE)


y=as.matrix((expressed_all_ave))

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycol <- sample(rainbow(256))
mycol <- mycol[as.vector(mycl)]


## Plot heatmap 
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(20)

heatmap.2(y, Rowv=as.dendrogram(hr), Colv=FALSE, dendrogram="row", col=Colors, scale="row", labRow=FALSE,
          density.info="none", trace="none", RowSideColors=mycol, cexCol = 0.8, margins=c(6,10),srtCol=90, adjCol = c(0.8,0.5),  
          lhei=c(2,16), lwid=c(2,10), keysize=0.55, key.par = list(cex=0.5))


############################FIGURE7#########################
cond= gsub(colnames(x), pattern = "[[:digit:]]", replacement = "")

tissue= gsub(cond, pattern = "_RT", replacement = "")
tissue= gsub(tissue, pattern = "_CE", replacement = "")

tissues = unique(tissue)

i0 <- NULL

for ( i in 1:length(tissues)) 
{
  
  set =  x[,grep(tissues[[i]], colnames(x))]
  keep <- filterByExpr(set,group = gsub(colnames(set), pattern = "[[:digit:]]", replacement = ""), min.count = 5, min.total.count = 10)
  res <- set[keep,]
  
  colnames(res) = gsub(colnames(res), pattern = "[[:digit:]]", replacement = "") 
  MAT_4ave=as.data.frame(cpm(res))
  nms <- unique(colnames(MAT_4ave))
  MAT_averaged = as.data.frame(sapply(nms, function(m)  rowMeans(MAT_4ave[names(MAT_4ave) %in% m])))
  colnames(MAT_averaged) = str_sub(colnames(MAT_averaged),-2,-1)
  towrite = cbind(tissue = tissues[[i]],symbol = rownames(res),MAT_averaged)
  i0 <- rbind(i0,towrite)
}
write.xlsx(i0, "counts.xlsx")


#### PieChart
for_pie = read_excel("C:/Users/hadadi/Dropbox/UniGen/manuscripts/All_tissue/GlobalExpression_Analysis.xlsx", sheet=2)

pie <-readClipboard()
pie = as.matrix(pie)
mycols <- c(viridis(10))

ggplot(for_pie, aes(x = "", y = gene, fill = tissue)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  
  
  coord_polar("y", start = 0)+ scale_fill_manual(values = mycols) +
  theme_void()

pie(for_pie$gene , border="white", col=mycols )
pie(pie , border="white", col=mycols )

ggplot(for_pie, aes(x = "", y = gene, fill = tissue)) +
  geom_bar(width = 1, stat = "identity", color = "white") + coord_polar("y", start = 0)+ 
  scale_fill_brewer(palette="Blues") +
  theme_void()


############################FIGURE8#########################
data  <- sqlQuery(con, "select * from 
	(select gene_id_name, comparison , logFC
	from FC_original_wo_Brain_RTrep2 where  logCPM > 0 and PValue <= 0.05 and logFC  NOT BETWEEN -0.599 AND 0.599  
	 ) as s
pivot (
	avg(logFC) for comparison in (BMC_VS_BMR,
Hypothalamus_CE_rep_VS_Hypothalamus_RT_rep,
Spleen_CE_rep_VS_Spleen_RT_rep,
Liver_CE_rep_VS_Liver_RT_rep,
VAT_CE_rep_VS_VAT_RT_rep,
Ileum_CE_rep_VS_Ileum_RT_rep,
D0CE_VS_D0RT,
BAT_CE_rep_VS_BAT_RT_rep,
SAT_CE_rep_VS_SAT_RT_rep,
Brain_CE_rep_VS_Brain_RT_rep)
	) as t
order by 1; ")
data= dcast (data, gene_id_name ~ comparison, value.var='logFC' )
write.xlsx(data, "5951-sig_changes_allTissues.xlsx", row.names = FALSE)


ave_expression_heatmap  = read.xlsx("averaged_expression_10Tissues_18479Genes.xlsx",1) 


for_heatmap = ave_expression_heatmap[grep("GO",ave_expression_heatmap$marker) ,]
rownames(for_heatmap) = for_heatmap$symbol
for_heatmap = for_heatmap[, -c(1:2)]

y=as.matrix(for_heatmap)

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") 
## Tree cutting
mycl <- cutree(hr, h=max(hr$height)/1.5)
mycol <- sample(rainbow(256))
mycol <- mycol[as.vector(mycl)]


## Plot heatmap 
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(20)

heatmap.2(y, Rowv=as.dendrogram(hr), Colv=FALSE, dendrogram="row", col=Colors, scale="row", labRow=TRUE,
          density.info="none", trace="none", RowSideColors=mycol, cexCol = 0.8, cexRow = 0.6,
          
          margins=c(6,10),srtCol=90, adjCol = c(0.8,0.5),  
          lhei=c(2,16), lwid=c(2,10), keysize=0.55, key.par = list(cex=0.5))


##Boxplot_FC
FC_boxplot = read.xlsx("5951-sig_changes_allTissues.xlsx",1) 

for_boxplot = FC_boxplot[grep("GO",FC_boxplot$marker) ,]
for_boxplot = for_boxplot[, -c(2)]


for_boxplot = for_boxplot[c("symbol", "Count" ,"SAT_CE_VS_RT","BAT_CE_VS_RT","VAT_CE_VS_RT","Liver_CE_VS_RT"
                            ,"Ileum_CE_VS_RT","Spleen_CE_VS_RT","Bonemarrow_CE_VS_RT", "SpinalCord_CE_VS_RT","Brain_CE_VS_RT","Hypothalamus_CE_VS_RT")]
colnames(for_boxplot) = gsub(colnames(for_boxplot), pattern = "_CE_VS_RT", replacement = "")

for_boxplot = melt(for_boxplot,id.vars = c("symbol", "Count") )

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),strip.background=element_blank(),plot.title = element_text(size=20),
             axis.title.x = element_blank(),axis.title.y = element_text(size=12),axis.text.x=element_text(colour="black",size=10,angle = 90,hjust = 1),
             axis.text.y=element_text(colour="black",size=11),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

ggplot(for_boxplot, aes(x=variable ,y=value)) + geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(colour = Count),size = 2.8,alpha = 0.8, shape = 16,position=position_jitter(0.2)) + 
  scale_color_manual(values = viridis(5)) +
  theme_bw() + theme(legend.title = element_blank(),legend.position="top") + theme




