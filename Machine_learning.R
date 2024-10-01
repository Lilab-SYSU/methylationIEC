
library(data.table)
library(e1071)
library(pROC)
library(ggplot2)

dat <- read.csv('Select_DMP_promoter_beta_value.csv',row.names = 1)
CG <- fread('CG_promoter_methylation_and_expression_correlation.csv')

loc=CG$correlation <= -0.27
select <- CG[loc,]
dim(select)
loc = match(select$CG,rownames(dat))
dat=dat[loc,]
DMP <- fread('All_DMP_ExhVsRest.csv')
loc <- match(rownames(dat),DMP$CG)
DMPSelect <- DMP[loc,]
fwrite(DMPSelect,file = "DMP_select.csv")
##########Prediction T cell exhaustion score SVM
score <- read.table('CD8_Tcell_exhausted_v1.xlsimmue_score_NMF.xls',sep = "\t",row.names = 1)
loc <- match(colnames(dat),colnames(score))
train_dat <- t(dat)
score1= score[,loc]
loc <- match(colnames(dat),colnames(score1))

score1T = t(score1)

model=svm(train_dat,score1T[,3],kernel = "polynomial")

LGG=predict(model)
saveRDS(model,file = "MethylationIEC.rds")
####correlation analysis
library(ggpubr)
library(ggsci)
cor.test(LGG,score1T[,3])
df = data.frame(Predicit=LGG,LGG_TEX=score1T[,3])
write.table(df,file = "Predict_LGG_TEX.xls",row.names = T,sep = "\t")
p <- ggplot(data = df, aes(x = Predicit, y = LGG_TEX)) +
  geom_point(colour="#DF8F44FF",size=2)+geom_smooth(method = lm,colour="#374E55FF",size=2) +
  scale_color_jama() +
  labs(x = 'Predict score',y="LGG_TEX") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'pearson', aes(x = Predicit, y = LGG_TEX))
p
ggsave(p,filename = paste0('Predict',"_LGG_TEX_correlation.pdf"),width = 6,height = 6)

mutation = fread('LGG_Mutations_number_stat.xls')
head(mutation)
predict = fread('Predict_LGG_TEX.xls')

###Age correlation
head(predict)
Clinical <- fread('clinical.xls')
loc <- match(predict$V1,gsub('-','.',Clinical$case_submitter_id))
predictage=cbind(predict,Age=Clinical$age_at_index[loc])
cor.test(predictage$Predicit,predictage$Age)
p <- ggplot(data = predictage, aes(x = Predicit, y = Age)) +
  geom_point(colour="#DF8F44FF",size=2)+geom_smooth(method = lm,colour="#374E55FF",size=2) +
  scale_color_jama() +
  labs(x = 'Predict score',y="Age") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'pearson', aes(x = Predicit, y = Age))
p
ggsave(p,filename = paste0('LGG',"_MethylationIEC_Age_correlation.pdf"),width = 6,height = 6)

###TMB correlation
loc <- match(gsub('\\.','-',predict$V1),mutation$samples)
predictMatch <- predict[!is.na(loc),]
loc <- match(gsub('\\.','-',predictMatch$V1),mutation$samples)
predictMatch[,'TMB']= mutation$Mutations[loc]
predictMatch=predictMatch[predictMatch$TMB < 1000,]
cor.test(predictMatch$Predicit,predictMatch$TMB,method = "spearman")

p <- ggplot(data = predictMatch, aes(x = Predicit, y = TMB)) +
  geom_point(colour="#DF8F44FF",size=2)+geom_smooth(method = lm,colour="#374E55FF",size=2) +
  scale_color_jama() +
  labs(x = 'Predict score',y="TMB") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'spearman', aes(x = Predicit, y = TMB))
p
ggsave(p,filename = paste0('LGG',"_MethylationIEC_TMB_correlation.pdf"),width = 6,height = 6)

#####LGG trainning ROC 
group <- fread('NMFconsushctree_cluster_group.xls')
loc <- match(group$samples,gsub('\\.','-',names(LGG)))
group[,'SVMscore']=LGG
group = as.data.frame(group)
group= group[order(group$SVMscore),]
modelroc=roc(group$group,group$SVMscore)
pdf('Trainning_ROC_LGG_v2.pdf',width = 5,height = 5)
plot(modelroc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()

####LGG survival analysis
dat <- fread('LGG_survival.txt')
head(dat)
nmf <- fread('NMFconsushctree_cluster_group.xls')

loc <- match(nmf$samples,dat$`_PATIENT`)
datPrimary <- dat[loc,]
Score <- fread('Predict_LGG_TEX.xls')
Score[Score$Predicit >= quantile(Score$Predicit,0.75),'group']="High IEC Score"
Score[Score$Predicit < quantile(Score$Predicit,0.75),'group']="Low"
loc <- match(gsub('\\.','-',Score$V1),dat$`_PATIENT`)
datPrimary <- dat[loc,]
matchednmfsurclass <- cbind(datPrimary,class=Score$group)
matchednmfsurclass$OS.time=matchednmfsurclass$OS.time/30
matchednmfsurclass$PFI.time=matchednmfsurclass$PFI.time/30
matchednmfsurclass$DSS.time=matchednmfsurclass$DSS.time/30
library(survival)
library(survminer)
fit<- survfit(Surv(DSS.time, DSS) ~ class, data = matchednmfsurclass)



p=ggsurvplot(fit, data = matchednmfsurclass,
             pval = T, # 在图上添加log rank检验的p值
             # conf.int = T,# 添加置信区间
             risk.table = T, # 在图下方添加风险表
             #legend.labs=c("Immune Exhaustion class","Rest class"), #表头标签注释男女
             #legend.title="strata",#表头标签
             title="Disease-specific survival",#改一下整体名称
             ylab="Disease-specific survival",xlab = "Follow-up(months)",#修改X轴Y轴名称
             risk.table.col = "strata", # 风险表加颜色
             #linetype = "strata", # 生存曲线的线型
             surv.median.line = "hv", # 标注出中位生存时间
             ggtheme = theme_bw(), #背景布局
             palette = c("red", "blue")) # 图形颜色风格

pdf("LGG_DSS_plot_predict.pdf",width = 5,height = 6)
p
dev.off()



fit<- survfit(Surv(PFI.time, PFI) ~ class, data = matchednmfsurclass)



p=ggsurvplot(fit, data = matchednmfsurclass,
             pval = T, # 在图上添加log rank检验的p值
             # conf.int = T,# 添加置信区间
             risk.table = T, # 在图下方添加风险表
             #legend.labs=c("Immune Exhaustion class","Rest class"), #表头标签注释男女
             #legend.title="strata",#表头标签
             title="Progress-free survival",#改一下整体名称
             ylab="Progress-free survival",xlab = "Follow-up(months)",#修改X轴Y轴名称
             risk.table.col = "strata", # 风险表加颜色
             #linetype = "strata", # 生存曲线的线型
             surv.median.line = "hv", # 标注出中位生存时间
             ggtheme = theme_bw(), #背景布局
             palette = c("red", "blue")) # 图形颜色风格
p
pdf("LGG_PFS_plot_predict.pdf",width = 5,height = 6)
p
dev.off()

fit<- survfit(Surv(OS.time, OS) ~ class, data = matchednmfsurclass)



p=ggsurvplot(fit, data = matchednmfsurclass,
             pval = T, # 在图上添加log rank检验的p值
             # conf.int = T,# 添加置信区间
             risk.table = T, # 在图下方添加风险表
             #legend.labs=c("Immune Exhaustion class","Rest class"), #表头标签注释男女
             #legend.title="strata",#表头标签
             title="Overall survival",#改一下整体名称
             ylab="Overall survival",xlab = "Follow-up(months)",#修改X轴Y轴名称
             risk.table.col = "strata", # 风险表加颜色
             #linetype = "strata", # 生存曲线的线型
             surv.median.line = "hv", # 标注出中位生存时间
             ggtheme = theme_bw(), #背景布局
             palette = c("red", "blue")) # 图形颜色风格
p
pdf("LGG_OS_plot_predict.pdf",width = 5,height = 6)
p
dev.off()



####GBM Testing
GBM_Dat <- fread('TCGA.GBM.sampleMap_HumanMethylation450.gz')
GBM_Test <- as.data.frame(GBM_Dat[,-1])
GBM_TestT <- t(GBM_Test)
colnames(GBM_TestT) <- GBM_Dat$sample
loc <- match(colnames(train_dat),colnames(GBM_TestT))
GBM_TestT1 <- GBM_TestT[,loc]
GBM_Predict= predict(model,GBM_TestT1)
as.data.frame(GBM_Predict)
write.table(as.data.frame(GBM_Predict),file = "GBM_MethylationIEC_predict.xls",sep = "\t")

GBM_TMB <- fread('GBM_Mutations_number_stat.xls')
head(GBM_Predict)
loc <- match(GBM_TMB$samples,GBM_Predict$V1)
GBM_TMBMatch <- GBM_TMB[!is.na(loc),]
loc <- match(GBM_TMBMatch$samples,GBM_Predict$V1)
GBM_TMBMatch[,'TMB']= GBM_Predict$GBM_Predict[loc]
cor.test(GBM_TMBMatch$Mutations,GBM_TMBMatch$TMB,method = "spearman")
GBM_TMBMatch=GBM_TMBMatch[GBM_TMBMatch$Mutations <2000,]

p <- ggplot(data = GBM_TMBMatch, aes(x = TMB, y = Mutations)) +
  geom_point(colour="#DF8F44FF",size=2)+geom_smooth(method = lm,colour="#374E55FF",size=2) +
  scale_color_jama() +
  labs(x = 'Predict score',y="TMB") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'spearman', aes(x = TMB, y = Mutations))
p
ggsave(p,filename = paste0('GBM',"_MethylationIEC_TMB_correlation.pdf"),width = 6,height = 6)

#####GBM prediction survival analysis
surDat <- fread('GBM_survival.txt')
head(surDat)

names(GBM_Predict)
loc <- match(names(GBM_Predict),surDat$sample)

GBM_Predict1 = GBM_Predict[!is.na(loc)]

loc <- match(names(GBM_Predict1),surDat$sample)
surDatPredict = surDat[loc,]
surDatPredict[,'Predict']= GBM_Predict1
surDatPredict[surDatPredict$Predict >= quantile(surDatPredict$Predict,0.75),'group'] = "High"
surDatPredict[surDatPredict$Predict < quantile(surDatPredict$Predict,0.75),'group'] = "Low"
dim(surDatPredict)
fwrite(surDatPredict,file = "Survival_data_predict_group.csv")

library(survival)
library(survminer)
surDatPredict$OS.time = surDatPredict$OS.time/30
fit<- survfit(Surv(OS.time, OS) ~ group, data = surDatPredict)
p=ggsurvplot(fit, data = surDatPredict,
             pval = T, # 在图上添加log rank检验的p值
             # conf.int = T,# 添加置信区间
             risk.table = T, # 在图下方添加风险表
             legend.labs=c("High IEC score","Low"), #表头标签注释男女
             #legend.title="strata",#表头标签
             title="Overall survival",#改一下整体名称
             ylab="Overall survival",xlab = "Follow-up(months)",#修改X轴Y轴名称
             risk.table.col = "strata", # 风险表加颜色
             #linetype = "strata", # 生存曲线的线型
             surv.median.line = "hv", # 标注出中位生存时间
             ggtheme = theme_bw(), #背景布局
             palette = c("red", "blue")) # 图形颜色风格
pdf('Survival_GBM_SVM_predict_LGG.pdf',width = 5,height = 6)
p
dev.off()

surDatPredict[surDatPredict$Predict >= quantile(surDatPredict$Predict,0.7),'group'] = "High"
surDatPredict[surDatPredict$Predict < quantile(surDatPredict$Predict,0.7),'group'] = "Low"

surDatPredict$PFI.time = surDatPredict$PFI.time/30
fit<- survfit(Surv(PFI.time, PFI) ~ group, data = surDatPredict)
p=ggsurvplot(fit, data = surDatPredict,
             pval = T, # 在图上添加log rank检验的p值
             # conf.int = T,# 添加置信区间
             risk.table = T, # 在图下方添加风险表
             legend.labs=c("High IEC score","Low"), #表头标签注释男女
             #legend.title="strata",#表头标签
             title="Overall survival",#改一下整体名称
             ylab="Progress-free survival",xlab = "Follow-up(months)",#修改X轴Y轴名称
             risk.table.col = "strata", # 风险表加颜色
             #linetype = "strata", # 生存曲线的线型
             surv.median.line = "hv", # 标注出中位生存时间
             ggtheme = theme_bw(), #背景布局
             palette = c("red", "blue")) # 图形颜色风格
p
pdf('PFS_GBM_SVM_predict_LGG.pdf',width = 5,height = 6)
p
dev.off()

surDatPredict$DSS.time = surDatPredict$DSS.time/30
fit<- survfit(Surv(DSS.time, DSS) ~ group, data = surDatPredict)
p=ggsurvplot(fit, data = surDatPredict,
             pval = T, # 在图上添加log rank检验的p值
             # conf.int = T,# 添加置信区间
             risk.table = T, # 在图下方添加风险表
             legend.labs=c("High IEC score","Rest class"), #表头标签注释男女
             #legend.title="strata",#表头标签
             title="DSS",#改一下整体名称
             ylab="Disease-specific survival",xlab = "Follow-up(months)",#修改X轴Y轴名称
             risk.table.col = "strata", # 风险表加颜色
             #linetype = "strata", # 生存曲线的线型
             surv.median.line = "hv", # 标注出中位生存时间
             ggtheme = theme_bw(), #背景布局
             palette = c("red", "blue")) # 图形颜色风格
p
pdf('DSS_GBM_SVM_predict_LGG.pdf',width = 5,height = 6)
p
dev.off()


es.dif <- read.table('GBM_Exhaustion_TEX_immune_cells_cytokine.xlsimmue_score_NMF.xls',row.names = 1,sep="\t")
es.dif[1:5,1:5]
GBM_Predict <- fread('GBM_MethylationIEC_predict.csv')
GBM_Clinical <- fread('TCGA-GBM.GDC_phenotype.tsv.gz')
samples=gsub("-01$","",GBM_Predict$V1)
samples = gsub("-02$","",samples)
samples = gsub("-11$","",samples)
loc <- match(samples,GBM_Clinical$submitter_id)
GBM_Predict_Age <- cbind(GBM_Predict,Age=GBM_Clinical$age_at_initial_pathologic_diagnosis[loc])
cor.test(GBM_Predict_Age$GBM_Predict,GBM_Predict_Age$Age,method = "spearman")
p <- ggplot(data = GBM_Predict_Age, aes(x = GBM_Predict, y = Age)) +
  geom_point(colour="#DF8F44FF",size=2)+geom_smooth(method = lm,colour="#374E55FF",size=2) +
  scale_color_jama() +
  labs(x = 'Predict score',y="Age") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = 'pearson', aes(x = GBM_Predict, y = Age))
p
ggsave(p,filename = paste0('GBM',"_MethylationIEC_Age_correlation.pdf"),width = 6,height = 6)




samples=gsub("-01$","",GBM_Predict$V1)
samples = gsub("-02$","",samples)
samples = gsub("-11$","",samples)
TEX_score=gsub("\\.","-",colnames(es.dif))

GBM_Predict$V1=samples
loc = match(samples,TEX_score)
samplesMatch = GBM_Predict[!is.na(loc),]

loc = match(samplesMatch$V1,TEX_score)

es.difDF= as.data.frame(t(es.dif))
TEX_scoreMatch= es.difDF$LGG_TEX[loc]
plotdf=data.frame(TEX_LGG=TEX_scoreMatch,GBM_Predict_IEC_Score=samplesMatch$GBM_Predict)
p <- ggplot(data = plotdf, aes(x = TEX_LGG, y = GBM_Predict_IEC_Score)) +
  geom_point(colour="#DF8F44FF",size=2) + geom_smooth(method = lm) +
  scale_color_aaas()  +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                     panel.border = element_rect(colour="black", fill=NA, size = 1.3)) +
  stat_cor(method = 'pearson', aes(x = TEX_LGG, y = GBM_Predict_IEC_Score))
p
ggsave(p,filename = "GBM_Predict_TEX_correlation.pdf",width = 6,height = 6)

nmfcluster <- fread('GBM_NMFconsushctree_cluster.xls')
loc=match(samplesMatch$V1,nmfcluster$samples1)
samplesMatch[,'group'] = nmfcluster$group[loc]

testROC=roc(samplesMatch$group,samplesMatch$GBM_Predict)
pdf('Testing_ROC_GBM_v2.pdf',width = 5,height = 5)
plot(testROC, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()