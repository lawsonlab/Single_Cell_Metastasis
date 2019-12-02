##Single Cell Classifier

allGenes <- read.table("~/AllPatients_MetvTum_Markers_tobit.txt", row.names=1,header=T)
data <- t(read.table("~/Combined_Expression_Matrix.txt", row.names=1,header=T))
metaData <- read.table("~/Combined_Expression_Matrix_MetaData.txt", row.names=1,header=T)


dim(data)
head(metaData)
data[1:5,1:15]
summary(allGenes)

length(which(colnames(data) %in% as.character(allGenes$gene)))
length(which(colnames(data) %in% rownames(allGenes)))

datRed = cbind(metaData,data[,which(colnames(data) %in% allGenes$gene)])
dim(datRed); dim(data)
datRed[,"burden"] = 1*(datRed[,"burden"] != "T")

head(datRed)
#cor(datRed[,-(1:8)])
datRed[,9:dim(datRed)[2]] = scale(datRed[,9:dim(datRed)[2]])
summary(datRed)

#######################
## Formulate Groups
#####################

samples = datRed	##this is all the data

########
detach(posGenes)
attach(posGenes)

zVals = data.frame()
accs = data.frame()
multiRunBest = data.frame()
multiRunPval = numeric()
multiRunAuc = numeric()
multiRunAic = numeric()
multiRunMccr = numeric()
for(k in 1:10){
  set.seed(k)
  table(datRed$mouse,datRed$burden)
  #max_rep = max(table(mouse,burden)); table(mouse,burden)
  min_rep = min(table(datRed$mouse,datRed$burden))
  
  samples = data.frame()
  for(i in levels(as.factor(datRed$mouse))){
    for(j in levels(as.factor(datRed$burden))){
      hits = (datRed$burden == j & i == datRed$mouse)
      nHits = sum(hits*1)
      samples = rbind(samples, datRed[hits,][sample(1:nHits,min_rep,replace=FALSE),]) #min_rep
    }
  }
  dim(samples)
  
  train = sample(1:dim(samples)[1],1*dim(samples)[1])		#.8*dim(samples)[1])
  #test = (1:dim(samples)[1])[-train]
  
  table(samples[train,6])
  ############
  
  ################################
  i = 1
  nGenes = ncol(samples)-8
  bestGenes = numeric()
  pvalBestGenes = numeric()
  aicBestGenes = numeric()
  aucs = numeric()
  aics = numeric()
  mccrs = numeric()
  
  for(j in 1:30){
    pvals = numeric()
    indAic = numeric()
    
    remaining = (1:nGenes)[-bestGenes]
    for(i in if(length(bestGenes)==0){1:nGenes}else{remaining}){#(dim(datRed)[2]-8)){
      
      modelDat = as.data.frame(cbind(samples[train,6],samples[train,8+c(bestGenes,i)]))
      #modelDat = as.data.frame(samples[,c(6,8+c(bestGenes,i))])
      colnames(modelDat)[1] = c("burden") #,colnames(datRed)[i+8])
      head(modelDat)
      model = glm(burden ~ .,family=binomial(link='logit'),data=modelDat)
      pvals = c(pvals,summary(model)$coeff[nrow(summary(model)$coeff),4])
      indAic = c(indAic,AIC(model))
      
    }
    pvalBestGenes = c(pvalBestGenes, min(pvals))
    aicBestGenes = c(aicBestGenes, min(indAic))
    
    if(length(bestGenes)==0){bestGenes = which(indAic %in% min(indAic))}else{ 
      bestGenes = c(bestGenes, remaining[which(indAic %in% min(indAic))]) }
    aucs = c(aucs, roc(samples$burden[train], 1*(predict(model)>0))$auc)
    #aucs = c(aucs, roc(samples$burden, 1*(predict(model)>0))$auc)
    aics = c(aics, AIC(model))
    mccrs = c(mccrs, mccr(samples$burden[train], 1*(predict(model)>0)))
    
  }
  #[pvalBestGenes < (.05/(nGenes:(nGenes-j+1)))]
  multiRunBest = rbind(multiRunBest, bestGenes)
  multiRunPval = rbind(multiRunPval, pvalBestGenes*(nGenes:(nGenes-j+1)))
  multiRunAic = rbind(multiRunAic, aics)
  multiRunAuc = rbind(multiRunAuc, aucs)
  multiRunMccr = rbind(multiRunAuc, mccrs)
  
} #end

meltedMccr = melt(multiRunMccr)
ggplot(meltedMccr, aes(Var2,value)) +  stat_summary(fun.data=mean_se, geom="errorbar") + theme_classic()
meltedAuc = melt(multiRunAuc)
ggplot(meltedAuc, aes(Var2,value)) +  stat_summary(fun.data=mean_se, geom="errorbar") + theme_classic()
meltedAic = melt(multiRunAic)
ggplot(meltedAic, aes(Var2,value)) +  stat_summary(fun.data=mean_se, geom="errorbar") + theme_classic()
meltedPval = melt(-log(multiRunPval,10))
ggplot(meltedPval, aes(Var2,value)) +  stat_summary(fun.data=mean_se, geom="errorbar") + theme_classic()+geom_hline(yintercept=-log(0.05,10))

sort(table(melt(multiRunBest[,1:5])$value),decreasing=T)

par(mfrow=c(2,2))
#plot(-log(pvalBestGenes*(nGenes:1)[1:100],10))
plot(-log(pvalBestGenes*(nGenes:1)[1:j],10))
abline(h = -log(0.05,10))

plot(aucs)
plot(aics)
plot(mccrs)

selection = 1:10
modelDat = as.data.frame(cbind(samples[train,6],samples[train,8+bestGenes[selection]]))
colnames(modelDat) = c("burden",colnames(datRed)[8+bestGenes[selection]])
model = glm(burden ~ .,family=binomial(link='logit'),data=modelDat)
summary(model)

trainPreds = predict(model,modelDat)
sortedPreds = sort(trainPreds)
par(mfrow=c(1,1))
plot(1/(1+exp(-sortedPreds)), col = 1+1*(samples[ train[order(trainPreds)],-c(1:5,7:8)]$burden == 1), 
     main = "Ordered Met Likelihood", ylab="Likelihood of being a MET")


## check for directional consistancy with original selection

head(allGenes)
directionDF = allGenes[,c("gene","cluster")]
directionDF$cluster = 1*(directionDF$cluster=="Metastatic")

directionDF[directionDF$gene %in% colnames(modelDat),]


##################################

multiGeneNames = character()
for(i in 1:nrow(multiRunBest)){
  multiGeneNames = c(multiGeneNames,colnames(datRed)[8+as.numeric(multiRunBest[i,1:5])])
}
geneTable = table(multiGeneNames)
sort(geneTable,decreasing=T)

modelDat = as.data.frame(cbind(samples[,6],samples[,c("LDHA", "PHLDA2","BHLHE40")]))
colnames(modelDat) = c("burden","LDHA", "PHLDA2","BHLHE40")
model = glm(burden ~ .,family=binomial(link='logit'),data=modelDat)
summary(model)

roc(samples$burden, 1*(predict(model)>0))$auc
AIC(model)
mccr(samples$burden, 1*(predict(model)>0))
sum(1*(predict(model)>0)== modelDat[,"burden"])/length(modelDat[,"burden"])

##

tumorAcc = sum((trainPreds<0) * (samples[train,-c(1:5,7:8)]$burden == 0))
metAcc = sum((trainPreds>0) * (samples[train,-c(1:5,7:8)]$burden == 1))
#tumorAcc = sum((testPreds<0) * (samples[test,-c(1:5,7:8)]$burden == 0))
#metAcc = sum((testPreds>0) * (samples[test,-c(1:5,7:8)]$burden == 1))

zVals = rbind(zVals, summary(model)$coeff[,"z value"])
accs = rbind(accs,c(tumorAcc,metAcc))

#}

boxplot(zVals, axes = F, range = 3,main="Z-score ranges for 10 models")
abline(h=0,col = 2)
axis(2)
axis(1, at=seq_along(1:31),labels=c("Intercept",as.character(usedTFs)), las=2)
