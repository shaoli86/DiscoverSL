#' findInteractor <- function()
#'
#' It shortlists genes having mutations in at-least 5 samples in cancers. This gene list is used in the later functions to predict SL partners
#' @return A character vector of gene symbols, of the shortlisted genes
#' @export
#' @examples
#' gene <- findInteractor()


findInteractor <- function(){
gene <- DiscoverSL:::gene
gene
}
#' Calculate differential expression
#'
#' Given the primary gene of interest and the cancer type, it calculates the differential expression of the gene list (output from findIteractor) in samples with or without mutation in the primary gene for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param cancer character. cancer type code. e.g. "BRCA"
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), logFC(log fold change of expression of gene2 in presence of mutation in gene1), PValue(p-value for differential expression of gene2 in presence of mutation in gene1)
#' @export
#' @examples
#' dg <- calculateDiffExp("BRCA2","BRCA")

calculateDiffExp <- function(gene1,cancer){
dfgene <- data.frame(Gene2="g",Gene1="g",logFC=0,logCPM=0,PValue=0)
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = paste(tolower(cancer),"_tcga_rna_seq_v2_mrna",sep="")
datax <- cgdsr::getProfileData(mycgds,gene1,mygeneticprofile,mycaselist)
},
error=function(cond){
datax <- data.frame(g=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
datax <- data.frame(g=0)
}
)
if(nrow(datax)>2){
colnames(datax) <- "RNASeqExp"
allcells <- row.names(datax)
for(i in 1:length(allcells)){
alls <- unlist(strsplit(allcells[i],"[.]"))
allcells[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
row.names(datax) <- allcells
tryCatch(
{mymutprofile = paste(tolower(cancer),"_tcga_mutations",sep="")
datam1 <- subset(cgdsr::getMutationData(mycgds,mycaselist,mymutprofile,gene1), select=c(gene_symbol,case_id,mutation_type))
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
dataexpr <- DiscoverSL:::dataexpr
genemut <- subset(datam1, mutation_type==1, select=case_id)
genemut <- as.character(genemut[,1])
genenonmut <- allcells
if(length(genemut)>=1&&length(genenonmut)>=1){
genemut <- unique(genemut)
genenonmut <- unique(genenonmut)
genenonmut <- setdiff(genenonmut, genemut)
if(length(genemut)>=4 && length(genenonmut)>=4){
mut <- genemut
nonmut <- genenonmut
for(i in 1:length(mut)){
saml <- unlist(strsplit(mut[i],"[-]"))
mut[i] <- paste(saml[1],saml[2],saml[3],sep="-")
}
for(i in 1:length(nonmut)){
saml <- unlist(strsplit(nonmut[i],"[-]"))
nonmut[i] <- paste(saml[1],saml[2],saml[3],sep="-")
}
mutexpr <- subset(dataexpr, select=names(dataexpr) %in% mut)
nonmutexpr <- subset(dataexpr, select=names(dataexpr) %in% nonmut)
if(ncol(mutexpr)>=3&&ncol(nonmutexpr)>=3){
dataexp <- cbind(mutexpr, nonmutexpr)
ms <- numeric(ncol(mutexpr))
nms <- numeric(ncol(nonmutexpr))
for(i in 1:ncol(mutexpr))
ms[i] <- 1
for(i in 1:ncol(nonmutexpr))
nms[i] <- 2
group <- factor(c(ms,nms))
y <- edgeR::DGEList(counts=dataexp,group=group)
y <- edgeR::calcNormFactors(y)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
de.tgw <- edgeR::exactTest( y )
dt <- de.tgw$table
genesig <- subset(dt,PValue<0.01)
dfgene <- data.frame(Gene1=rep(gene1, nrow(dt)), dt)
gene2 <- row.names(dfgene)
dfgene <- data.frame(Gene2=gene2,dfgene)
row.names(dfgene) <- NULL
}
}
}
},
error=function(cond){
message(paste("Error while accessing gene mutation data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
}
)
}
dfgene <- dfgene[-1,]
dfgene
}

#' Calculate Expression correlation
#'
#' Given the primary gene of interest and the cancer type, it calculates the expression correlation of the primary gene with the interactor genes from the supplied gene list (output from findIteractor) for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("CHEK1","ATM","NF1","EGFR")
#' @param cancer character. cancer type code. e.g. "BRCA"
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), exp.correlation(Pearson's coefficient for expression correlation of gene2 and gene1), correlation.pvalue(p-value for expression correlation of gene2 and gene1)
#' @export
#' @examples
#' dfgene <- calculateDiffExp("BRCA2","BRCA")
#' gene <- findInteractor()
#' cg <- calculateExpCorrelation("BRCA2",gene,"BRCA")

calculateExpCorrelation <- function(gene1,gene2,cancer){
slg <- gene2
slcor <- numeric(length(slg))
slcorp <- numeric(length(slg))
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = paste(tolower(cancer),"_tcga_rna_seq_v2_mrna",sep="")
data1 <- cgdsr::getProfileData(mycgds,gene1,mygeneticprofile,mycaselist)
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "RNASeqExp"
data1 <- subset(data1, RNASeqExp!="NaN")
data1 <- as.numeric(data1[,1])
for(i in 1:length(slg)){
gene2 <- slg[i]
data2 <- data.frame(m=0)
tryCatch(
{data2 <- cgdsr::getProfileData(mycgds,gene2,mygeneticprofile,mycaselist)
},
error=function(cond){
data2 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data2 <- data.frame(m=0)
}
)
if(nrow(data2)>2){
colnames(data2) <- "RNASeqExp"
data2 <- subset(data2, RNASeqExp!="NaN")
data2 <- as.numeric(data2[,1])
if(length(data1)>=4 && length(data2)>=4){
genemat <- as.matrix(c(data1,data2),nrow=length(data1),ncol=2, bycol=T)
cv <- cor(data1,data2,method="pearson")
if(!is.na(cv)){
slcor[i] <- as.numeric(cor.test(data1,data2,method="pearson")$estimate)
slcorp[i] <- cor.test(data1,data2,method="pearson")$p.value
}
}
else{
slcor[i] <- 0
slcorp[i] <- 0
}
}
}
}
cordata <- data.frame(Gene1=rep(gene1,length(slg)),Gene2=slg,exp.correlation=slcor,correlation.pvalue=slcorp)
cordata
}

#' Calculate Mutual Exclusivity
#'
#' Given the primary gene of interest and the cancer type, it calculates the mutual exclusivity of amplification, deletion and mutation events of the primary gene and the interactor genes from the supplied gene list (output from findIteractor) for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @param cancer character. cancer type code. e.g. "BRCA"
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), MutexAmp(p-value for mutual exclusive amplification of gene1 and gene2), MutexDel(p-value for mutual exclusive deletion of gene1 and gene2), MutexMut(p-value for mutual exclusive mutation of gene1 and gene2)
#' @export
#' @examples
#' dfgene <- calculateDiffExp("BRCA2","BRCA")
#' gene <- findInteractor()
#' cg <- calculateExpCorrelation("BRCA2",gene,"BRCA")
#' gene <- as.character(cg$Gene2)
#' mutdata <- calculateMutex("BRCA2",gene,"BRCA")

calculateMutex <- function(gene1,gene2,cancer){
slg <- gene2
mutexmut <- numeric(length(slg))
mutexamp <- numeric(length(slg))
mutexdel <- numeric(length(slg))
data1 <- data.frame(m=0)
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = paste(tolower(cancer),"_tcga_gistic",sep="")
data1 <- cgdsr::getProfileData(mycgds,gene1,mygeneticprofile,mycaselist)
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "CNA"
allcells <- row.names(data1)
for(i in 1:length(allcells)){
alls <- unlist(strsplit(allcells[i],"[.]"))
allcells[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
datam1 <- data.frame(m=0)
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
mymutprofile = paste(tolower(cancer),"_tcga_mutations",sep="")
datam1 <- subset(cgdsr::getMutationData(mycgds,mycaselist,mymutprofile,gene1), select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
datam1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
datam1 <- data.frame(m=0)
}
)
if(nrow(datam1)>2){
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
cellmut <- unique(as.character(datam1$case_id))
allcells <- unique(c(allcells,cellmut))
gene1mut <- subset(datam1, gene_symbol==gene1 & mutation_type==1, select="case_id")
gene1mut <- unique(as.character(gene1mut$case_id))
gene1amp <- row.names(subset(data1, CNA>=2))
gene1amp <- unique(as.character(gene1amp))
if(length(gene1amp)>=1){
for(i in 1:length(gene1amp)){
alls <- unlist(strsplit(gene1amp[i],"[.]"))
gene1amp[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
}
gene1del <- row.names(subset(data1, CNA<= -2))
gene1del <- unique(as.character(gene1del))
if(length(gene1del)>=1){
for(i in 1:length(gene1del)){
alls <- unlist(strsplit(gene1del[i],"[.]"))
gene1del[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
}
gene1alt <- unique(c(gene1mut,gene1amp,gene1del))
for(j in 1:length(slg)){
gene2 <- slg[j]
data2 <- data.frame(m=0)
datam2 <- data.frame(m=0)
tryCatch(
{data2 <- cgdsr::getProfileData(mycgds,gene2,mygeneticprofile,mycaselist)
datam2 <- subset(cgdsr::getMutationData(mycgds,mycaselist,mymutprofile,gene2), select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
data2 <- data.frame(m=0)
datam2 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data2 <- data.frame(m=0)
datam2 <- data.frame(m=0)
}
)
if(nrow(datam2)>2 && nrow(data2)>2){
colnames(data2) <- "CNA"
datam2 <- within(datam2, mutation_type[mutation_type != 'Silent'] <- '1')
datam2 <- within(datam2, mutation_type[mutation_type == 'Silent'] <- '0')
gene2mut <- subset(datam2, gene_symbol==gene2 & mutation_type==1, select="case_id")
gene2mut <- unique(as.character(gene2mut$case_id))
gene2amp <- row.names(subset(data2, CNA>=2))
gene2amp <- unique(as.character(gene2amp))
if(length(gene2amp)>=1){
for(i in 1:length(gene2amp)){
alls <- unlist(strsplit(gene2amp[i],"[.]"))
gene2amp[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
}
gene2del <- row.names(subset(data2, CNA<= -2))
gene2del <- unique(as.character(gene2del))
if(length(gene2del)>=1){
for(i in 1:length(gene2del)){
alls <- unlist(strsplit(gene2del[i],"[.]"))
gene2del[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
}
gene2alt <- unique(c(gene2mut,gene2amp,gene2del))
gene1gene2mut <- intersect(gene1mut,gene2mut)
gene1gene2amp <- intersect(gene1amp,gene2amp)
gene1gene2del <- intersect(gene1del,gene2del)
gene1gene2alt <- unique(c(gene1gene2amp,gene1gene2mut,gene1gene2del))
mut1 <- length(gene1amp)
mut2 <- length(gene2amp)
mutboth <- length(gene1gene2amp)-1
allc <- length(allcells)
nonmut1 <- allc-mut1
mutexamp[j] <- round(phyper(mutboth,mut1,nonmut1,mut2,lower.tail=FALSE,log.p=FALSE),2)
mut1 <- length(gene1mut)
mut2 <- length(gene2mut)
mutboth <- length(gene1gene2mut)-1
allc <- length(allcells)
nonmut1 <- allc-mut1
mutexmut[j] <- round(phyper(mutboth,mut1,nonmut1,mut2,lower.tail=FALSE,log.p=FALSE),2)
mut1 <- length(gene1del)
mut2 <- length(gene2del)
mutboth <- length(gene1gene2del)-1
allc <- length(allcells)
nonmut1 <- allc-mut1
mutexdel[j] <- round(phyper(mutboth,mut1,nonmut1,mut2,lower.tail=FALSE,log.p=FALSE),2)
}
}
}
}
mutdata <- data.frame(Gene1=rep(gene1,length(slg)),Gene2=slg,MutexAmp=mutexamp,MutexDel=mutexdel,MutexMut=mutexmut)
mutdata
}

#' Calculate probability of sharing pathways
#'
#' Given the primary gene of interest and the interactor gene list, it calculates the probability of the primary gene and the interactor gene to share common pathways
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), PvalPathway(p-value for gene1 and gene2 sharing common pathways)
#' @export
#' @examples
#' gene <- findInteractor()
#' genepath <- calculatePathway("BRCA2",gene)

calculatePathway <- function(gene1,gene){
genesets <- DiscoverSL:::genesets
genesets <- subset(genesets, 1==startsWith(genesets[,1],"REACTOME_") | 1==startsWith(genesets[,1],"KEGG_"))
slg <- gene
pathpval <- numeric(length(slg))
numboth <- numeric(length(slg))
numgene1 <- numeric(length(slg))
numgene2 <- numeric(length(slg))
for(i in 1:length(slg)){
numboth[i] <- 0
numgene1[i] <- 0
numgene2[i] <- 0
 for(j in 1:nrow(genesets)){
 genes <- unlist(strsplit(genesets[j,3],","))
 if(gene1 %in% genes)
 numgene1[i]=numgene1[i]+1
 if(slg[i] %in% genes)
 numgene2[i]=numgene2[i]+1
 if(gene1 %in% genes && slg[i] %in% genes)
 numboth[i]=numboth[i]+1
 }
allc <- nrow(genesets)
nongene1 <- allc-numgene1[i]
pathpval[i] <- phyper((numboth[i]-1),numgene1[i],nongene1,numgene2[i],lower.tail=FALSE,log.p=FALSE)
}
genepath <- data.frame(Gene1=rep(gene1,length(slg)), Gene2=slg, PvalPathway=pathpval)
genepath
}

#' Predict synthetic lethality
#'
#' Given the primary gene of interest and the cancer type, it calculates the Synthetic lethality score for all interactor genes from the supplied gene list (output from findIteractor) for the cancer type, using the four parameters calculated from functions calculateDiffExp, calculateExpCorrelation, calculateMutex and calculatePathway
#' @import randomForest
#' @param dg data frame. Output from calculateDiffExp
#' @param cg data frame. Output from calculateExpCorrelation
#' @param mutdata data frame. Output from calculateMutex
#' @param genepath data frame. Output from calculatePathway
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), Diffexp.Pvalue(p-value for differential expression of gene2 in presence of mutation in gene1), correlation.pvalue(p-value for expression correlation of gene2 and gene1), MutexMut(p-value for mutual exclusive mutation of gene1 and gene2), PvalPathway(p-value for gene1 and gene2 sharing common pathways), MetaFDR (FDR corrected combined P-value from four parameters), SyntheticLethalScore(Predicted synthetic lethality score between gene1 and gene2)
#' @export
#' @examples
#' dg <- calculateDiffExp("BRCA2","BRCA")
#' gene <- findInteractor()
#' cg <- calculateExpCorrelation("BRCA2",gene,"BRCA")
#' gene <- as.character(cg$Gene2)
#' mutdata <- calculateMutex("BRCA2",gene,"BRCA")
#' genepath <- calculatePathway("BRCA2",gene)
#' sl <- predictSL(dg,cg,mutdata,genepath)

predictSL <- function(dg,cg,mutdata,genepath){
synleth <- dplyr::inner_join(as.data.frame(dg),cg)
synleth <- dplyr::inner_join(synleth,mutdata)
synleth <- dplyr::inner_join(synleth,genepath)
synleth$MutexAmp <- 1 - synleth$MutexAmp
synleth$MutexDel <- 1 - synleth$MutexDel
synleth$MutexMut <- 1 - synleth$MutexMut
synleth <- na.omit(synleth)
pvalm <- numeric(nrow(synleth))
for(i in 1:nrow(synleth)){
allp <- c(synleth$MutexMut[i],synleth$MutexAmp[i],synleth$MutexDel[i])
tryCatch(
{pvalm[i]<- metap::sumlog(allp)$p
},
error=function(cond){
pvalm[i] <- 0
},
warning=function(cond){
pvalm[i] <- 0
}
)
}
synleth <- data.frame(synleth, Mutex=pvalm)
pval <- numeric(nrow(synleth))
for(i in 1:nrow(synleth)){
allp <- c(synleth$PValue[i],synleth$correlation.pvalue[i],synleth$Mutex[i],synleth$PvalPathway[i])
tryCatch(
{pval[i]<- metap::sumlog(allp)$p
},
error=function(cond){
pval[i] <- 0
},
warning=function(cond){
pval[i] <- 0
}
)
}
synleth <- data.frame(synleth, MetaPval=pval)
mf <- p.adjust(synleth$MetaPval, method="fdr")
synleth <- data.frame(synleth, MetaFDR=mf)
msl <- synleth
msl$PvalPathway <- 1 - msl$PvalPathway
testdata <- matrix(c(as.numeric(msl$PValue),as.numeric(msl$Mutex),as.numeric(msl$correlation.pvalue),as.numeric(msl$PvalPathway)), nrow=nrow(msl), ncol=4)
model <- DiscoverSL:::model
pred <- predict(model, newdata = testdata, type="response")
synleth <- data.frame(synleth, SyntheticLethalScore=pred)
synleth <- unique(synleth)
synleth <- subset(synleth, select=c("Gene1","Gene2","PValue","correlation.pvalue","MutexMut","PvalPathway","MetaFDR","SyntheticLethalScore"))
colnames(synleth)[3] <- "Diffexp.Pvalue"
if(nrow(synleth)>100)
synleth <- subset(synleth, SyntheticLethalScore>0.5)
synleth
}

#' Calculate conditional gene essentiality from RNAi screening
#'
#' Given the primary gene of interest and the interactor gene list, it calculates the difference in RNAi score for targeting the interactor gene in cases where the primary gene is mutated vs not mutated
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), PvalRNAi(p-value for gene2 having less RNAi score in cells where gene1 is mutated compared to cells where gene1 is not mutated)
#' @export
#' @examples
#' dg <- calculateDiffExp("BRCA2","BRCA")
#' gene <- as.character(dfgene$Gene2)
#' cg <- calculateExpCorrelation("BRCA2",gene,"BRCA")
#' gene <- as.character(cg$Gene2)
#' mutdata <- calculateMutex("BRCA2",gene,"BRCA")
#' genepath <- calculatePathway("BRCA2",gene)
#' sl <- predictSL(dg,cg,mutdata,genepath)
#' gene <- as.character(sl$Gene2)
#' rnaiscore <- calculateRNAiScore("BRCA2",gene)

calculateRNAiScore <- function(gene1,gene2){
mcle <- DiscoverSL:::mcle
mrnai <- DiscoverSL:::mrnai
slgene <- unique(gene2)
pval <- numeric(length(slgene))
prgene <- character(length(slgene))
meanmut <- numeric(length(slgene))
meannonmut <- numeric(length(slgene))
meanfc <- numeric(length(slgene))
i=1;
shalt <- subset(mcle, Name==paste(gene1,"_MUT",sep="") & value==1, select=variable)
shalt <- unique(shalt)
shunalt <- subset(mcle, Name==paste(gene1,"_MUT",sep="") & value==0, select=variable)
shunalt <- unique(shunalt)
rnaimut <- subset(mrnai, variable %in% shalt[,1])
rnainon <- subset(mrnai, variable %in% shunalt[,1])
for(i in 1:length(slgene)){
emutalt <- subset(rnaimut, Description==slgene[i], select=value)
enmutunalt <- subset(rnainon, Description==slgene[i], select=value)
emut <- as.numeric(emutalt[,1])
enmut <- as.numeric(enmutunalt[,1])
prgene[i] <- gene1
if(length(emut)>4 && length(enmut)>4 && emut[1]!="NaN" && enmut[1]!="NaN"){
tp <- t.test(emut,y=enmut,alternative="less")
pval[i] <- tp$p.value
meanmut[i] <- mean(emut)
meannonmut[i] <- mean(enmut)
}
else{
pval[i] <- 0
meanmut[i] <- 0
meannonmut[i] <- 0
}
}
d <- data.frame(Gene1=prgene,Gene2=slgene,PvalRNAi=pval)
d
}

#' Calculate relative targetability by checking amplification difference
#'
#' Given the primary gene of interest and the cancer type, it calculates the difference in amplification status of the interactor genes from the supplied gene list (output from findIteractor) in cases where the primary gene is mutated vs not mutated for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @param cancer character. cancer type code. e.g. "BRCA"
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), TTestMutAmp(p-value for gene2 having higher amplification of gene2 in cases where gene1 is mutated compared to cases where gene1 is not mutated)
#' @export
#' @examples
#' gene <- findInteractor()
#' ampdata <- calculateTTestAmp("BRCA2",gene,"BRCA")

calculateTTestAmp <- function(gene1,gene2,cancer){
slg <- gene2
ttestamp <- numeric(length(slg))
datam1 <- data.frame(m=0)
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = paste(tolower(cancer),"_tcga_gistic",sep="")
mymutprofile = paste(tolower(cancer),"_tcga_mutations",sep="")
datam1 <- subset(cgdsr::getMutationData(mycgds,mycaselist,mymutprofile,gene1), select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
datam1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
datam1 <- data.frame(m=0)
}
)
if(nrow(datam1)>2){
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
cellmut <- unique(as.character(datam1$case_id))
colnames(datam1) <- c("Name","Cell_line","Mutflag")
gene1mut <- subset(datam1, Name==gene1 & Mutflag==1, select="Cell_line") #Gene1 mutated samples
gene1mut <- unique(as.character(gene1mut$Cell_line))
for(j in 1:length(slg)){
gene2 <- slg[j]
data1 <- data.frame(m=0)
tryCatch(
{data1 <- cgdsr::getProfileData(mycgds,gene2,mygeneticprofile,mycaselist)
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "CNA"
allcells <- row.names(data1)
for(i in 1:length(allcells)){
alls <- unlist(strsplit(allcells[i],"[.]"))
allcells[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
cmutl <- data.frame(Name=rep(gene2,nrow(data1)),Cell_line=allcells,Mutflag=data1$CNA)
gene2ampgene1mut <- subset(cmutl, Name==gene2 & Cell_line %in% gene1mut, select="Mutflag")
emut <- as.numeric(gene2ampgene1mut[,1])
gene2ampgene1nonmut <- subset(cmutl, Name==gene2 & !(Cell_line %in% gene1mut), select="Mutflag")
enmut <- as.numeric(gene2ampgene1nonmut[,1])
ttestamp[j] <- 1
if(length(emut)>2 && length(enmut)>2){
try(tp <- t.test(emut,y=enmut,alternative="greater"))
try(ttestamp[j] <- tp$p.value)
}
}
}
}
df <- data.frame(Gene1=rep(gene1,length(slg)),Gene2=slg,TTestMutAmp=ttestamp)
df
}

#' Calculate conditional drug sensitivity
#'
#' Given the primary gene of interest and the interactor gene list, it calculates the difference in Drug activity score (IC50) for targeting the interactor gene in cases where the primary gene is mutated vs not mutated
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), Drug (Drug name), PvalDrug(p-value for gene2 having less RNAi score in cells where gene1 is mutated compared to cells where gene1 is not mutated)
#' @export
#' @examples
#' gene <- findInteractor()
#' drugscore <- calculateDrugSensitivity("BRCA2",gene)

calculateDrugSensitivity <- function(gene,gene2){
ccle <- DiscoverSL:::ccle
rnaidrug <- DiscoverSL:::rnaidrug
shalt <- subset(ccle, Gene1==gene & Classification!="silent", select=COSMIC_ID)
shalt <- unique(shalt)
shunalt <- subset(ccle, Gene1!=gene, select=COSMIC_ID)
shunalt <- unique(shunalt)
rnaimut <- dplyr::inner_join(shalt,rnaidrug)
rnainon <- dplyr::inner_join(shunalt,rnaidrug)
prgene <- gene
slgene <- unique(gene2)
sldrug <- unique(as.character(rnaidrug$DrugName))
slen <- length(slgene) + length(sldrug)
pval <- numeric(slen)
gene1 <- character(slen)
gene2 <- character(slen)
meanmut <- numeric(slen)
meannonmut <- numeric(slen)
meanfc <- numeric(slen)
l=1;
drug <- character(slen)
for(i in 1:length(slgene)){
emutalt <- subset(rnaimut, Gene2==slgene[i])
enmutunalt <- subset(rnainon, Gene2==slgene[i])
unidrug <- intersect(unique(as.character(emutalt$DrugName)), unique(as.character(enmutunalt$DrugName)))
unidrug <- unique(unidrug)
for(k in 1:length(unidrug)){
mutsub <- subset(emutalt, DrugName==unidrug[k])
nonmutsub <- subset(enmutunalt, DrugName==unidrug[k])
emut <- as.numeric(mutsub$LN_IC50)
enmut <- as.numeric(nonmutsub$LN_IC50)
if(length(emut)>2 && length(enmut)>2){
gene1[l] <- prgene
gene2[l] <- slgene[i]
try(tp <- t.test(emut,y=enmut,alternative="greater"))
try(pval[l] <- tp$p.value)
meanmut[l] <- mean(emut)
meannonmut[l] <- mean(enmut)
drug[l] <- unidrug[k]
l = l + 1
}
}
}
d <- data.frame(Gene1=gene1,Gene2=gene2,Drug=drug,PvalDrug=pval)
d <- subset(d, Gene1!="" & Gene2!="" & Drug!="")
d
}

#' Summarize potential synthetic lethals by all parameters
#'
#' Given the prior calculated parameters SyntheticLethalScore, RNAiscore,DrugSensitivity and TargetAmplification, it summarizes all these parameters in a single table for all interactor genes
#' @param sl data frame. Output from predictSL
#' @param rnaiscore data frame. Output from calculateRNAiScore
#' @param drugscore data frame. Output from calculateDrugSensitivity
#' @param ampdata data frame. Output from calculateTTestAmp
#' @return A data frame with following fields: Gene1(primary gene), Gene2(interactor gene), Diffexp.Pvalue(p-value for differential expression of gene2 in presence of mutation in gene1), correlation.pvalue(p-value for expression correlation of gene2 and gene1), MutexMut(p-value for mutual exclusive mutation of gene1 and gene2), PvalPathway(p-value for gene1 and gene2 sharing common pathways), MetaFDR (FDR corrected combined P-value from four parameters), SyntheticLethalScore(Predicted synthetic lethality score between gene1 and gene2), PvalRNAi (P-value for conditional essentiality from RNAi screens), Drug (Drug name), PvalDrug (P-value for conditional drugsensitivity), TTestMutAmp (p-value for conditional targetability)
#' @export
#' @examples
#' dg <- calculateDiffExp("BRCA2","BRCA")
#' gene <- findInteractor()
#' cg <- calculateExpCorrelation("BRCA2",gene,"BRCA")
#' gene <- as.character(cg$Gene2)
#' mutdata <- calculateMutex("BRCA2",gene,"BRCA")
#' genepath <- calculatePathway("BRCA2",gene)
#' sl <- predictSL(dg,cg,mutdata,genepath)
#' rnaiscore <- calculateRNAiScore(gene1,gene)
#' drugscore <- calculateDrugSensitivity(gene1,gene)
#' ampdata <- calculateTTestAmp(gene1,gene,cancer)
#' finalist <- summarizeSLParams(sl,rnaiscore,drugscore,ampdata)

summarizeSLParams <- function(sl,rnaiscore,drugscore,ampdata){
slfinal <- dplyr::left_join(sl,rnaiscore)
slfinal <- dplyr::left_join(slfinal,drugscore)
slfinal <- dplyr::left_join(slfinal,ampdata)
if(nrow(slfinal)>100)
slfinal <- subset(slfinal, (PvalRNAi<0.1 & TTestMutAmp<0.1) | PvalDrug<0.1)
slfinal
}

#' Plot conditional gene essentiality from Achilles RNAi screening
#'
#' Given the primary gene of interest and the interactor gene, it calculates the difference in RNAi score for targeting the interactor gene in cases where the primary gene is mutated vs not mutated and shows the difference in a box plot
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character. gene symbol for the interactor gene. e.g. "PARP1"
#' @param mcle data frame. output from loadMutationDataCCLE
#' @param mrnai data frame. output from loadRNAidataAchilles
#' @export
#' @examples
#' plotRNAi("BRCA1","PARP1")

plotRNAi <- function(gene1,gene2){
mcle <- DiscoverSL:::mcle
mrnai <- DiscoverSL:::mrnai
shalt <- subset(mcle, Name==paste(gene1,"_MUT",sep="") & value==1, select=variable)
shalt <- unique(shalt)
shunalt <- subset(mcle, Name==paste(gene1,"_MUT",sep="") & value==0, select=variable)
shunalt <- unique(shunalt)
rnaimut <- subset(mrnai, variable %in% shalt[,1])
rnainon <- subset(mrnai, variable %in% shunalt[,1])
emutalt <- subset(rnaimut, Description==gene2, select=value)
enmutunalt <- subset(rnainon, Description==gene2, select=value)
emut <- as.numeric(emutalt[,1])
enmut <- as.numeric(enmutunalt[,1])
tryCatch(
{tp <- t.test(emut,y=enmut,alternative="less")
pval <- tp$p.value
meanmut <- mean(emut)
meannonmut <- mean(enmut)
mutrnai <- data.frame(RNAiScore=emut, Status=rep(paste(gene1,"Mutated\nSamples:",length(emut)),length(emut)))
nonmutrnai <- data.frame(RNAiScore=enmut, Status=rep(paste(gene1,"Not Mutated\nSamples:",length(enmut)),length(enmut)))
mdata <- rbind(mutrnai,nonmutrnai)
p <- ggplot2::ggplot(mdata,ggplot2::aes(Status,RNAiScore))
p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill=Status)) +
ggplot2::ggtitle(paste("Cancer Cell lines Response to targeting",gene2,"in cases",gene1,"mutated vs not mutated", sep=" ")) +
ggplot2::annotate("text", x = 1.5, y = 6, label = paste("P-value:",formatC(pval,format="e",digits=2)))
p
},
error=function(cond){
message(paste("Not enough data in Achilles for targeting",gene2,"in cells carrying",gene1,"mutation for performing t-test for conditional essentiality"))
},
warning=function(cond){
}
)
}

#' Plot conditional sensitivity for a given drug
#'
#' Given the primary gene of interest, the interactor gene and the drug name targeting the interactor gene, it calculates the difference in Drug activity score (IC50) for targeting the interactor gene in cases where the primary gene is mutated vs not mutated and shows the difference in a box plot
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character. gene symbol for the interactor gene. e.g. "PARP1"
#' @param drug character. drug name. e.g. "Olaparib"
#' @export
#' @examples
#' plotSensitivitybyDrug("BRCA1","PARP1","Olaparib")

plotSensitivitybyDrug <- function(gene1,gene2,drug){
ccle <- DiscoverSL:::ccle
rnaidrug <- DiscoverSL:::rnaidrug
shalt <- subset(ccle, Gene1==gene1 & Classification!="silent", select=COSMIC_ID)
shalt <- unique(shalt)
shunalt <- subset(ccle, Gene1!=gene1, select=COSMIC_ID)
shunalt <- unique(shunalt)
rnaimut <- dplyr::inner_join(shalt,rnaidrug)
rnainon <- dplyr::inner_join(shunalt,rnaidrug)
emutalt <- subset(rnaimut, Gene2==gene2)
enmutunalt <- subset(rnainon, Gene2==gene2)
mutsub <- subset(emutalt, DrugName==drug)
nonmutsub <- subset(enmutunalt, DrugName==drug)
emut <- as.numeric(mutsub$LN_IC50)
enmut <- as.numeric(nonmutsub$LN_IC50)
tryCatch(
{tp <- t.test(emut,y=enmut,alternative="greater")
pval <- tp$p.value
meanmut <- mean(emut)
meannonmut <- mean(enmut)
mutrnai <- data.frame(DrugScore=emut, Status=rep(paste(gene1,"Mutated\nSamples:",length(emut)),length(emut)))
nonmutrnai <- data.frame(DrugScore=enmut, Status=rep(paste(gene1,"Not Mutated\nSamples:",length(enmut)),length(enmut)))
mdata <- rbind(mutrnai,nonmutrnai)
p <- ggplot2::ggplot(mdata,ggplot2::aes(Status,DrugScore))
p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill=Status)) +
ggplot2::ggtitle(paste("Cancer Cell lines Response to targeting",gene2,"in cases",gene1,"mutated vs not mutated", sep=" ")) +
ggplot2::annotate("text", x = 1.5, y = 6, label = paste("P-value:",formatC(pval,format="e",digits=2)))
p
},
error=function(cond){
message(paste("Not enough data in GDSC for drug",drug,"targeting",gene2,"in cells carrying",gene1,"mutation for performing t-test for conditional drug sensitivity"))
},
warning=function(cond){
}
)
}

#' Plot conditional sensitivity to any drugs targeting SL gene
#'
#' Given the primary gene of interest and the interactor gene, it calculates the difference in Drug activity score (IC50) for targeting the interactor gene in cases where the primary gene is mutated vs not mutated and shows the difference in a box plot
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character. gene symbol for the interactor gene. e.g. "PARP1"
#' @export
#' @examples
#' plotDrugSensitivity("BRCA1","PARP1")

plotDrugSensitivity <- function(gene1,gene2){
ccle <- DiscoverSL:::ccle
rnaidrug <- DiscoverSL:::rnaidrug
shalt <- subset(ccle, Gene1==gene1 & Classification!="silent", select=COSMIC_ID)
shalt <- unique(shalt)
shunalt <- subset(ccle, Gene1!=gene1, select=COSMIC_ID)
shunalt <- unique(shunalt)
rnaimut <- dplyr::inner_join(shalt,rnaidrug)
rnainon <- dplyr::inner_join(shunalt,rnaidrug)
emutalt <- subset(rnaimut, Gene2==gene2)
enmutunalt <- subset(rnainon, Gene2==gene2)
mutsub <- emutalt
nonmutsub <- enmutunalt
emut <- as.numeric(mutsub$LN_IC50)
enmut <- as.numeric(nonmutsub$LN_IC50)
pval <- 1
tryCatch(
{tp <- t.test(emut,y=enmut,alternative="greater")
pval <- tp$p.value
mutrnai <- data.frame(DrugScore=emut, Status=rep(paste(gene1,"Mutated\nSamples:",length(emut)),length(emut)))
nonmutrnai <- data.frame(DrugScore=enmut, Status=rep(paste(gene1,"Not Mutated\nSamples:",length(enmut)),length(enmut)))
mdata <- rbind(mutrnai,nonmutrnai)
p <- ggplot2::ggplot(mdata,ggplot2::aes(Status,DrugScore))
p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill=Status)) +
ggplot2::ggtitle(paste("Cancer Cell lines Response to targeting",gene2,"in cases",gene1,"mutated vs not mutated", sep=" ")) +
ggplot2::annotate("text", x = 1.5, y = 6, label = paste("P-value:",formatC(pval,format="e",digits=2)))
p
},
error=function(cond){
message(paste("Not enough data in GDSC for drugs targeting",gene2,"in cells carrying",gene1,"mutation for performing t-test for conditional drug sensitivity"))
},
warning=function(cond){
}
)
}

#' Plot relative amplification difference for SL gene
#'
#' Given the primary gene of interest and the cancer type, it calculates the difference in amplification status of the interactor gene in cases where the primary gene is mutated vs not mutated for the cancer type and shows the difference in a box plot
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character. gene symbol for the primary gene. e.g. "PARP1"
#' @param cancer character. cancer type code. e.g. "BRCA"
#' @export
#' @examples
#' plotAmplificationDiff("BRCA1","PARP1","BRCA")

plotAmplificationDiff <- function(gene1,gene2,cancer){
data1 <- data.frame(m=0)
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = paste(tolower(cancer),"_tcga_gistic",sep="")
data1 <- cgdsr::getProfileData(mycgds,gene2,mygeneticprofile,mycaselist)
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "CNA"
allcells <- row.names(data1)
for(i in 1:length(allcells)){
alls <- unlist(strsplit(allcells[i],"[.]"))
allcells[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep="-")
}
cmutl <- data.frame(Name=rep(gene2,nrow(data1)),Cell_line=allcells,Mutflag=data1$CNA)
mymutprofile = paste(tolower(cancer),"_tcga_mutations",sep="")
datam1 <- data.frame(m=0)
tryCatch(
{datam1 <- subset(cgdsr::getMutationData(mycgds,mycaselist,mymutprofile,gene1), select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
datam1 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
datam1 <- data.frame(m=0)
}
)
if(nrow(datam1)>2){
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
cellmut <- unique(as.character(datam1$case_id))
allcells <- unique(c(allcells,cellmut))
colnames(datam1) <- c("Name","Cell_line","Mutflag")
gene1mut <- subset(datam1, Name==gene1 & Mutflag==1, select="Cell_line")
gene1mut <- unique(as.character(gene1mut$Cell_line))
gene2ampgene1mut <- subset(cmutl, Name==gene2 & Cell_line %in% gene1mut, select="Mutflag")
emut <- as.numeric(gene2ampgene1mut[,1])
gene2ampgene1nonmut <- subset(cmutl, Name==gene2 & !(Cell_line %in% gene1mut), select="Mutflag")
enmut <- as.numeric(gene2ampgene1nonmut[,1])
tryCatch(
{tp <- t.test(emut,y=enmut,alternative="greater")
ttestamp <- tp$p.value
mutrnai <- data.frame(AmpDiff=emut, Status=rep(paste(gene1,"Mutated\nSamples:",length(emut)),length(emut)))
nonmutrnai <- data.frame(AmpDiff=enmut, Status=rep(paste(gene1," Not Mutated\nSamples:",length(enmut)),length(enmut)))
mdata <- rbind(mutrnai,nonmutrnai)
p <- ggplot2::ggplot(mdata,ggplot2::aes(Status,AmpDiff))
p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill=Status)) +
ggplot2::ggtitle(paste("Difference in amplification of",gene2,"in cases",gene1,"mutated vs not mutated", sep=" ")) +
ggplot2::annotate("text", x = 1.5, y = 6, label = paste("P-value:",formatC(ttestamp,format="e",digits=2)))
p
},
error=function(cond){
message(paste("Not enough data in",cancer,"for amplification of",gene2,"in samples carrying",gene1,"mutation for performing t-test for conditional targetability"))
},
warning=function(cond){
}
)
}
else{
message(paste("Not enough data in",cancer,"for amplification of",gene2,"in samples carrying",gene1,"mutation for performing t-test for conditional targetability"))
}
}
else{
message(paste("Not enough data in",cancer,"for amplification of",gene2,"in samples carrying",gene1,"mutation for performing t-test for conditional targetability"))
}
}

#' Plot Survival difference for SL co-inactivation
#'
#' Given the primary gene of interest and the cancer type, it calculates the difference in patient disease-free survival depending on the co-inactivation of the primary gene and the interactor gene for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character. gene symbol for the interactor gene. e.g. "PARP1"
#' @param cancer character. cancer type code. e.g. "BRCA"
#' @export
#' @examples
#' plotSurvivalCurveSL("BRCA1","PARP1","BRCA")

plotSurvivalCurveSL <- function(gene1,gene2,cancer){
tryCatch(
{mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
mycancerstudy <- paste(tolower(cancer),"tcga",sep="_")
mycaselist = cgdsr::getCaseLists(mycgds,mycancerstudy)[1,1]
datac <- cgdsr::getClinicalData(mycgds,mycaselist)
mygeneticprofile = paste(tolower(cancer),"_tcga_rna_seq_v2_mrna_median_Zscores",sep="")
data2 <- cgdsr::getProfileData(mycgds,gene2,mygeneticprofile,mycaselist)
data2 <- na.omit(data2)
mymutprofile = paste(tolower(cancer),"_tcga_mutations",sep="")
datam1 <- subset(cgdsr::getMutationData(mycgds,mycaselist,mymutprofile,gene1), select=c(gene_symbol,case_id,mutation_type))
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
},
error=function(cond){
message(paste("Error retrieving data from CBioportal for",gene1,"and",gene2,"in",cancer, "\n"))
},
warning=function(cond){
}
)
if(nrow(datam1)>2 && nrow(data2)>2){
genemut1 <- as.character(subset(datam1, mutation_type==1, select=case_id)[,1])
for(i in 1:length(genemut1)){
alls <- unlist(strsplit(genemut1[i],"-"))
genemut1[i] <- paste(alls[1],alls[2],alls[3],alls[4],sep=".")
}
gene2mutexp <- subset(data2, row.names(data2) %in% genemut1)
geneup2 <- row.names(subset(gene2mutexp, gene2mutexp[,1]>=median(gene2mutexp[,1])))
genedown2 <- row.names(subset(gene2mutexp, gene2mutexp[,1]<median(gene2mutexp[,1])))
if("DFS_MONTHS" %in% colnames(datac) && length(geneup2)>2 && length(genedown2)>2){
clinical <- subset(datac, SAMPLE_TYPE!="Blood Derived Normal" & SAMPLE_TYPE!="Solid Tissue Normal", select=c(DFS_STATUS,DFS_MONTHS,AGE))
condition <- character(nrow(clinical))
mutcondition <- character(nrow(clinical))
for(i in 1:nrow(clinical)){
if(row.names(clinical)[i] %in% genemut1)
condition[i] <- paste(gene1,"_mut",sep="")
else if(!row.names(clinical)[i] %in% genemut1)
condition[i] <- paste(gene1,"_nomut",sep="")
else
conditionn[i] <- "NA"
}
for(i in 1:nrow(clinical)){
if(row.names(clinical)[i] %in% genemut1 && row.names(clinical)[i] %in% geneup2)
mutcondition[i] <- paste(gene2,"_up",sep="")
else if(row.names(clinical)[i] %in% genemut1 && row.names(clinical)[i] %in% genedown2)
mutcondition[i] <- paste(gene2,"_down",sep="")
else
mutcondition[i] <- "NA"
}
clinical <- data.frame(clinical,Condition1=condition, Condition2=mutcondition, stringsAsFactors=FALSE)
clinical <- subset(clinical, DFS_STATUS!="")
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == '[Not Available]'] <- '0')
clinical <- within(clinical, DFS_MONTHS[DFS_MONTHS == 'NA'] <- '0')
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == '[Not Available]'] <- 'NA')
clinical <- subset(clinical, DFS_STATUS!="NA")
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'DiseaseFree'] <- 1)
clinical <- within(clinical, DFS_STATUS[DFS_STATUS == 'Recurred/Progressed'] <- 0)
clinical <- subset(clinical, Condition1!="NA")
clin1 <- subset(clinical, Condition2==paste(gene2,"_down",sep="")|Condition2==paste(gene2,"_up",sep=""), select=c("DFS_MONTHS","DFS_STATUS","AGE","Condition2"))
cl1 <- subset(clinical, Condition1!="NA",select=c("DFS_MONTHS","DFS_STATUS","AGE","Condition1"))
if(nrow(clin1)>4 && nrow(cl1)>4){
pb <- list()
blank <- grid::rectGrob(gp=grid::gpar(col="white"))
tryCatch(
{mutsurv1 <- survival::survfit(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ Condition1, data=cl1)
pb[[1]] <- survminer::ggsurvplot(mutsurv1, data=cl1, pval=TRUE, risk.table=TRUE, title=paste(gene1,"_Mutated+Not_Mutated", sep=""), legend="none")
mutsurv2 <- survival::survfit(survival::Surv(as.numeric(DFS_MONTHS), as.numeric(DFS_STATUS)) ~ Condition2, data=clin1)
pb[[2]] <- survminer::ggsurvplot(mutsurv2, data=clin1, pval=TRUE, risk.table=TRUE, title=paste(gene1,"_Mutated", sep=""), legend="none")
print(gridExtra::grid.arrange(pb[[1]]$plot,blank,pb[[2]]$plot,pb[[1]]$table,blank,pb[[2]]$table,ncol=3,nrow=2, heights=c(1,0.4), widths=c(1,0.1,1)))
},
error=function(cond){
message(paste("Cannot plot Survival curves for",gene1,"and",gene2,"in",cancer,"\n",cond,"\n"))
},
warning=function(cond){
message(paste("Warning: for",gene1,"and",gene2,"in",cancer,"\n",cond,"\n"))
}
)
}
else{
message(paste("Cannot plot Survival curves for",gene1,"and",gene2,"in",cancer,"\n","Not enough clinical data in mutation subgroups"))
}
}
else{
message(paste("Cannot plot Survival curves for",gene1,"and",gene2,"in",cancer,"\n","Not enough clinical data in mutation subgroups"))
}
}
else{
message(paste("Cannot plot Survival curves for",gene1,"and",gene2,"in",cancer,"\n","Not enough data in mutation subgroups"))
}
}

#' Calculate Mutation associated differential expression for user provided mutation and gene expression data
#'
#' Given the primary gene of interest and the cancer type, it calculates the differential expression of the gene list (output from findIteractor) in samples with or without mutation in the primary gene for the cancer type
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param mymutdata data frame. gene mutation data in tumor samples in the three column format: (gene_symbol, tumor_sample_id, variant_classification).
#' @param myexpression data frame. gene expression data in tumor samples in a matrix format. The rows in the matrix will correspond to genes and the columns will correspond to tumor samples.
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), logFC(log fold change of expression of gene2 in presence of mutation in gene1), PValue(p-value for differential expression of gene2 in presence of mutation in gene1)
#' @export
#' @examples
#' mymutdata
#'   Hugo_Symbol Tumor_Sample_Barcode Variant_Classification
#' 1       TKTL1      TCGA.AB.2988.03                 Silent
#' 2      GUCA2A      TCGA.AB.2988.03      Missense_Mutation
#' 3      DNMT3A      TCGA.AB.2988.03      Missense_Mutation
#' 4     SLC17A3      TCGA.AB.2988.03      Missense_Mutation
#' 5       LRWD1      TCGA.AB.2988.03      Missense_Mutation
#' 6        BAAT      TCGA.AB.2988.03      Missense_Mutation
#' myexpression
#'   Hugo_Symbol TCGA.AB.2803.03 TCGA.AB.2805.03 TCGA.AB.2806.03 TCGA.AB.2807.03 TCGA.AB.2808.03 TCGA.AB.2810.03
#' 1      TKTL1        563.7860        480.5890        415.6588        407.1500        638.6489        350.2327
#' 2      GUCA2A       232.5103        0.6693          103.7556        204.5680        27.4941         36.8666
#' 3      DNMT3A       3169.7531       5251.6734       3499.6817       4340.6157       3394.3441       2516.1456
#' 4      SLC17A3      0.0000          0.0000          0.0000          0.0000          0.0000          0.0000
#' 5      LRWD1        341.5638        167.3360        309.9936        55.6107         197.1720        21.0666
#' 6      BAAT         279.4136        106.4726        218.0204        378.0139        164.0927        155.0899
#' dg <- customDataDiffExp("BRCA2",mymutdata,myexpression)


customDataDiffExp <- function(gene1,mymutdata,myexpression){
dfgene <- data.frame(Gene2="g",Gene1="g",logFC=0,logCPM=0,PValue=0)
tryCatch(
{names(mymutdata) <- c("gene_symbol","case_id","mutation_type")
datam1 <- subset(mymutdata, select=c(gene_symbol,case_id,mutation_type))
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
dataexpr <- myexpression
row.names(dataexpr) <- dataexpr[,1]
dataexpr <- round(dataexpr[,-1])
genemut <- subset(datam1, gene_symbol==gene1 & mutation_type==1, select=case_id)
genemut <- as.character(genemut[,1])
genenonmut <- unique(names(dataexpr))
if(length(genemut)>=1&&length(genenonmut)>=1){
genemut <- unique(genemut)
genenonmut <- unique(genenonmut)
genenonmut <- setdiff(genenonmut, genemut)
if(length(genemut)>=1 && length(genenonmut)>=1){
mut <- genemut
nonmut <- genenonmut
mutexpr <- subset(dataexpr, select=names(dataexpr) %in% mut)
nonmutexpr <- subset(dataexpr, select=names(dataexpr) %in% nonmut)
if(ncol(mutexpr)>=1&&ncol(nonmutexpr)>=1){
dataexp <- cbind(mutexpr, nonmutexpr)
ms <- numeric(ncol(mutexpr))
nms <- numeric(ncol(nonmutexpr))
for(i in 1:ncol(mutexpr))
ms[i] <- 1
for(i in 1:ncol(nonmutexpr))
nms[i] <- 2
group <- factor(c(ms,nms))
y <- edgeR::DGEList(counts=dataexp,group=group)
y <- edgeR::calcNormFactors(y)
design <- model.matrix(~group)
y <- edgeR::estimateDisp(y, design)
de.tgw <- edgeR::exactTest( y )
dt <- de.tgw$table
genesig <- subset(dt,PValue<0.01)
dfgene <- data.frame(Gene1=rep(gene1, nrow(dt)), dt)
gene2 <- row.names(dfgene)
dfgene <- data.frame(Gene2=gene2,dfgene)
row.names(dfgene) <- NULL
}
}
}
},
error=function(cond){
message(paste("Error while reading custom dataset:",cond,"\n"))
},
warning=function(cond){
}
)
dfgene <- dfgene[-1,]
dfgene
}

#' Calculate Expression correlation for user provided gene expression data
#'
#' Given the primary gene of interest and the cancer type, it calculates the expression correlation of the primary gene with the interactor genes from the supplied gene list (output from findIteractor) for the cancer type
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("CHEK1","ATM","NF1","EGFR")
#' @param myexpression data frame. gene expression data in tumor samples in a matrix format. The rows in the matrix will correspond to genes and the columns will correspond to tumor samples.
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), exp.correlation(Pearson's coefficient for expression correlation of gene2 and gene1), correlation.pvalue(p-value for expression correlation of gene2 and gene1)
#' @export
#' @examples
#' myexpression
#'   Hugo_Symbol TCGA.AB.2803.03 TCGA.AB.2805.03 TCGA.AB.2806.03 TCGA.AB.2807.03 TCGA.AB.2808.03 TCGA.AB.2810.03
#' 1      TKTL1        563.7860        480.5890        415.6588        407.1500        638.6489        350.2327
#' 2      GUCA2A       232.5103        0.6693          103.7556        204.5680        27.4941         36.8666
#' 3      DNMT3A       3169.7531       5251.6734       3499.6817       4340.6157       3394.3441       2516.1456
#' 4      SLC17A3      0.0000          0.0000          0.0000          0.0000          0.0000          0.0000
#' 5      LRWD1        341.5638        167.3360        309.9936        55.6107         197.1720        21.0666
#' 6      BAAT         279.4136        106.4726        218.0204        378.0139        164.0927        155.0899
#' gene2 <- findInteractor()
#' cg <- customDataExpCorrelation("BRCA2",gene2,myexpression)

customDataExpCorrelation <- function(gene1,gene2,myexpression){
slg <- gene2
slcor <- numeric(length(slg))
slcorp <- numeric(length(slg))
tryCatch(
{names(myexpression)[1] <- "Hugo_Symbol"
data1 <- subset(myexpression, Hugo_Symbol==gene1)
data1 <- as.data.frame(t(data1[,-1]))
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while reading gene expression data:",cond,"\n"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "RNASeqExp"
data1 <- subset(data1, RNASeqExp!="NaN")
data1 <- as.numeric(data1[,1])
for(i in 1:length(slg)){
gene2 <- slg[i]
data2 <- data.frame(m=0)
tryCatch(
{data2 <- subset(myexpression, Hugo_Symbol==gene2)
data2 <- as.data.frame(t(data2[,-1]))
},
error=function(cond){
data2 <- data.frame(m=0)
message(paste("Error while reading gene expression data:",cond,"\n"))
},
warning=function(cond){
data2 <- data.frame(m=0)
}
)
if(nrow(data2)>2){
colnames(data2) <- "RNASeqExp"
data2 <- subset(data2, RNASeqExp!="NaN")
data2 <- as.numeric(data2[,1])
if(length(data1)>=4 && length(data2)>=4){
genemat <- as.matrix(c(data1,data2),nrow=length(data1),ncol=2, bycol=T)
cv <- cor(data1,data2,method="pearson")
if(!is.na(cv)){
slcor[i] <- as.numeric(cor.test(data1,data2,method="pearson")$estimate)
slcorp[i] <- cor.test(data1,data2,method="pearson")$p.value
}
}
else{
slcor[i] <- 0
slcorp[i] <- 0
}
}
}
}
cordata <- data.frame(Gene1=rep(gene1,length(slg)),Gene2=slg,exp.correlation=slcor,correlation.pvalue=slcorp)
cordata
}

#' Calculate Mutual exclusivity for user provided mutation and CNA data
#'
#' Given the primary gene of interest and the cancer type, it calculates the mutual exclusivity of amplification, deletion and mutation events of the primary gene and the interactor genes from the supplied gene list (output from findIteractor) for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @param mycnadata data frame. gene-wise copy number alteration data (GISTIC processed) in tumor samples in a matrix format. The rows in the matrix will correspond to genes and the columns will correspond to tumor samples.
#' @param mymutdata data frame. gene mutation data in tumor samples in the three column format: (gene_symbol, tumor_sample_id, variant_classification).
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), MutexAmp(p-value for mutual exclusive amplification of gene1 and gene2), MutexDel(p-value for mutual exclusive deletion of gene1 and gene2), MutexMut(p-value for mutual exclusive mutation of gene1 and gene2)
#' @export
#' @examples
#' gene <- findInteractor()
#' mycnadata
#'   Hugo_Symbol TCGA.3C.AAAU.01 TCGA.3C.AALI.01 TCGA.3C.AALJ.01 TCGA.3C.AALK.01 TCGA.4H.AAAK.01 TCGA.5L.AAT0.01
#' 1 TKTL1               1              -1              -1               0               0               0
#' 2 GUCA2A              1               1               0               0               0               0
#' 3 DNMT3A              1               1               2               0               1               0
#' 4 SLC17A3             1               0               0              -1              -1              -1
#' 5 LRWD1              -1              -1               0               0               0               0
#' 6 BAAT                1               1               0               0               0               0
#' mymutdata
#'   Hugo_Symbol Tumor_Sample_Barcode Variant_Classification
#' 1       TKTL1      TCGA.AB.2988.03                 Silent
#' 2      GUCA2A      TCGA.AB.2988.03      Missense_Mutation
#' 3      DNMT3A      TCGA.AB.2988.03      Missense_Mutation
#' 4     SLC17A3      TCGA.AB.2988.03      Missense_Mutation
#' 5       LRWD1      TCGA.AB.2988.03      Missense_Mutation
#' 6        BAAT      TCGA.AB.2988.03      Missense_Mutation
#' mutdata <- customDataMutex("BRCA2",gene,mycnadata,mymutdata)

customDataMutex <- function(gene1,gene2,mycnadata,mymutdata){
slg <- gene2
mutexmut <- numeric(length(slg))
mutexamp <- numeric(length(slg))
mutexdel <- numeric(length(slg))
data1 <- data.frame(m=0)
tryCatch(
{names(mycnadata)[1] <- "Hugo_Symbol"
data1 <- subset(mycnadata, Hugo_Symbol==gene1)
data1 <- as.data.frame(t(data1[,-1]))
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while reading CNA data:",cond,"\n"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "CNA"
allcells <- row.names(data1)
datam1 <- data.frame(m=0)
tryCatch(
{names(mymutdata) <- c("gene_symbol","case_id","mutation_type")
datam1 <- subset(mymutdata, gene_symbol==gene1, select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
datam1 <- data.frame(m=0)
message(paste("Error while reading mutation data:",cond,"\n"))
},
warning=function(cond){
datam1 <- data.frame(m=0)
}
)
if(nrow(datam1)>2){
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
cellmut <- unique(as.character(datam1$case_id))
allcells <- unique(c(allcells,cellmut))
gene1mut <- subset(datam1, gene_symbol==gene1 & mutation_type==1, select="case_id")
gene1mut <- unique(as.character(gene1mut$case_id))
gene1amp <- row.names(subset(data1, CNA>=2))
gene1amp <- unique(as.character(gene1amp))
gene1del <- row.names(subset(data1, CNA<= -2))
gene1del <- unique(as.character(gene1del))
gene1alt <- unique(c(gene1mut,gene1amp,gene1del))
for(j in 1:length(slg)){
gene2 <- slg[j]
data2 <- data.frame(m=0)
datam2 <- data.frame(m=0)
tryCatch(
{names(mycnadata)[1] <- "Hugo_Symbol"
data2 <- subset(mycnadata, Hugo_Symbol==gene2)
data2 <- as.data.frame(t(data2[,-1]))
names(mymutdata) <- c("gene_symbol","case_id","mutation_type")
datam2 <- subset(mymutdata, gene_symbol==gene2, select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
data2 <- data.frame(m=0)
datam2 <- data.frame(m=0)
message(paste("Error while accessing gene expression data from cBioPortal:",cond,"\n"))
},
warning=function(cond){
data2 <- data.frame(m=0)
datam2 <- data.frame(m=0)
}
)
if(nrow(datam2)>2 && nrow(data2)>2){
colnames(data2) <- "CNA"
datam2 <- within(datam2, mutation_type[mutation_type != 'Silent'] <- '1')
datam2 <- within(datam2, mutation_type[mutation_type == 'Silent'] <- '0')
gene2mut <- subset(datam2, gene_symbol==gene2 & mutation_type==1, select="case_id")
gene2mut <- unique(as.character(gene2mut$case_id))
gene2amp <- row.names(subset(data2, CNA>=2))
gene2amp <- unique(as.character(gene2amp))
gene2del <- row.names(subset(data2, CNA<= -2))
gene2del <- unique(as.character(gene2del))
gene2alt <- unique(c(gene2mut,gene2amp,gene2del))
gene1gene2mut <- intersect(gene1mut,gene2mut)
gene1gene2amp <- intersect(gene1amp,gene2amp)
gene1gene2del <- intersect(gene1del,gene2del)
gene1gene2alt <- unique(c(gene1gene2amp,gene1gene2mut,gene1gene2del))
mut1 <- length(gene1amp)
mut2 <- length(gene2amp)
mutboth <- length(gene1gene2amp)-1
allc <- length(allcells)
nonmut1 <- allc-mut1
mutexamp[j] <- round(phyper(mutboth,mut1,nonmut1,mut2,lower.tail=FALSE,log.p=FALSE),2)
mut1 <- length(gene1mut)
mut2 <- length(gene2mut)
mutboth <- length(gene1gene2mut)-1
allc <- length(allcells)
nonmut1 <- allc-mut1
mutexmut[j] <- round(phyper(mutboth,mut1,nonmut1,mut2,lower.tail=FALSE,log.p=FALSE),2)
mut1 <- length(gene1del)
mut2 <- length(gene2del)
mutboth <- length(gene1gene2del)-1
allc <- length(allcells)
nonmut1 <- allc-mut1
mutexdel[j] <- round(phyper(mutboth,mut1,nonmut1,mut2,lower.tail=FALSE,log.p=FALSE),2)
}
}
}
}
mutdata <- data.frame(Gene1=rep(gene1,length(slg)),Gene2=slg,MutexAmp=mutexamp,MutexDel=mutexdel,MutexMut=mutexmut)
mutdata
}

#' Calculate conditional essentiality (difference in RNAi activity) for user provided mutation and RNAi screening data
#'
#' Given the primary gene of interest and the interactor gene list, it calculates the difference in RNAi score for targeting the interactor gene in cases where the primary gene is mutated vs not mutated
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @param celmutdata data frame. gene mutation data in tumor cell lines in the three column format: (gene_symbol, cell_line_id, variant_classification).
#' @param cellrnai data frame. RNAi activity data in tumor cell lines in a matrix format. The rows in the matrix will correspond to genes and the columns will correspond to tumor cell lines.
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), PvalRNAi(p-value for gene2 having less RNAi score in cells where gene1 is mutated compared to cells where gene1 is not mutated)
#' @export
#' @examples
#' celmutdata
#'        gene_symbol    case_id mutation
#' 1         AKT3 DMS53_LUNG        0
#' 2         ABI1 DMS53_LUNG        0
#' 3         CDH2 DMS53_LUNG        0
#' 4 LOC100130776 DMS53_LUNG        0
#' 5        HDAC6 DMS53_LUNG        0
#' 6      BCL2L11 DMS53_LUNG        0
#' cellrnai
#'                           shRNA target_gene X22RV1_PROSTATE X697_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE X786O_KIDNEY
#' 1 GCCCACCAAAGCAAGAAACTT_61E3.4      61E3.4      -1.0261547                              -2.3936993   -0.9911353
#' 2   CAGTGACAGAAGCAGCCATAT_A1BG        A1BG      -1.6216855                              -1.7048357   -2.8186215
#' 3   CCGCCTGTGCTGATGCACCAT_A1BG        A1BG      -1.4392042                              -1.2772873   -2.1917983
#' 4   CTCTTCGAGCTGCACAACATT_A1BG        A1BG      -1.3452530                              -1.1004493   -2.9640288
#' 5   GAGTCCGAATCACTGCTGAAA_A1BG        A1BG      -0.7926535                              -0.3372055    0.2775138
#' 6   CCATGCTGCAAGGAGAGTATA_A1CF        A1CF      -1.4126040                              -0.7056490    1.3978028
#' gene <- findInteractor()
#' rnaiscore <- customDataRNAiScore("BRCA2",gene,celmutdata,cellrnai)

customDataRNAiScore <- function(gene1,gene2,celmutdata,cellrnai){
tryCatch(
{mrnai <- reshape2::melt(cellrnai)
names(mrnai) <- c("shRNA","gene_symbol","case_id","rnaiscore")
names(celmutdata) <- c("gene_symbol","case_id","mutation")
},
error=function(cond){
mrnai <- data.frame(m=0)
message(paste("Error while reading RNAi data:",cond,"\n","Please type help(customDataRNAiscore) to see the appropriate format of input data"))
},
warning=function(cond){
}
)
slgene <- unique(gene2)
pval <- numeric(length(slgene))
prgene <- character(length(slgene))
meanmut <- numeric(length(slgene))
meannonmut <- numeric(length(slgene))
meanfc <- numeric(length(slgene))
i=1;
shalt <- subset(celmutdata, gene_symbol==gene1 & mutation==1, select=case_id)
shalt <- unique(shalt)
shunalt <- subset(celmutdata, gene_symbol==gene1 & mutation==0, select=case_id)
shunalt <- unique(shunalt)
rnaimut <- subset(mrnai, case_id %in% shalt[,1])
rnainon <- subset(mrnai, case_id %in% shunalt[,1])
for(i in 1:length(slgene)){
emutalt <- subset(rnaimut, gene_symbol==slgene[i], select=rnaiscore)
enmutunalt <- subset(rnainon, gene_symbol==slgene[i], select=rnaiscore)
emut <- as.numeric(emutalt[,1])
enmut <- as.numeric(enmutunalt[,1])
prgene[i] <- gene1
if(length(emut)>4 && length(enmut)>4 && emut[1]!="NaN" && enmut[1]!="NaN"){
tp <- t.test(emut,y=enmut,alternative="less")
pval[i] <- tp$p.mutation
meanmut[i] <- mean(emut)
meannonmut[i] <- mean(enmut)
}
else{
pval[i] <- 0
meanmut[i] <- 0
meannonmut[i] <- 0
}
}
d <- data.frame(Gene1=prgene,Gene2=slgene,PvalRNAi=pval)
d
}


#' Calculate relative targetability (diffrence in amplification w.r.t. presence of mutation in primary gene) for user provided mutation and CNA data
#'
#' Given the primary gene of interest and the cancer type, it calculates the difference in amplification status of the interactor genes from the supplied gene list (output from findIteractor) in cases where the primary gene is mutated vs not mutated for the cancer type
#' @import cgdsr
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @param mymutdata data frame. gene mutation data in tumor samples in the three column format: (gene_symbol, tumor_sample_id, variant_classification).
#' @param mycnadata data frame. gene-wise copy number alteration data (GISTIC processed) in tumor samples in a matrix format. The rows in the matrix will correspond to genes and the columns will correspond to tumor samples.
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), TTestMutAmp(p-value for gene2 having higher amplification of gene2 in cases where gene1 is mutated compared to cases where gene1 is not mutated)
#' @export
#' @examples
#' gene <- findInteractor()
#' mymutdata
#'   Hugo_Symbol Tumor_Sample_Barcode Variant_Classification
#' 1       TKTL1      TCGA.AB.2988.03                 Silent
#' 2      GUCA2A      TCGA.AB.2988.03      Missense_Mutation
#' 3      DNMT3A      TCGA.AB.2988.03      Missense_Mutation
#' 4     SLC17A3      TCGA.AB.2988.03      Missense_Mutation
#' 5       LRWD1      TCGA.AB.2988.03      Missense_Mutation
#' 6        BAAT      TCGA.AB.2988.03      Missense_Mutation
#' mycnadata
#'   Hugo_Symbol TCGA.3C.AAAU.01 TCGA.3C.AALI.01 TCGA.3C.AALJ.01 TCGA.3C.AALK.01 TCGA.4H.AAAK.01 TCGA.5L.AAT0.01
#' 1 TKTL1               1              -1              -1               0               0               0
#' 2 GUCA2A              1               1               0               0               0               0
#' 3 DNMT3A              1               1               2               0               1               0
#' 4 SLC17A3             1               0               0              -1              -1              -1
#' 5 LRWD1              -1              -1               0               0               0               0
#' 6 BAAT                1               1               0               0               0               0
#' ampdata <- customDataTTestAmp("BRCA2",gene,mymutdata,mycnadata)

customDataTTestAmp <- function(gene1,gene2,mymutdata,mycnadata){
slg <- gene2
ttestamp <- numeric(length(slg))
datam1 <- data.frame(m=0)
tryCatch(
{names(mymutdata) <- c("gene_symbol","case_id","mutation_type")
datam1 <- subset(mymutdata, gene_symbol==gene1, select=c(gene_symbol,case_id,mutation_type))
},
error=function(cond){
datam1 <- data.frame(m=0)
message(paste("Error while reading mutation data:",cond,"\n","Please type help(customDataTTestAmp) to see the appropriate format of input data"))
},
warning=function(cond){
datam1 <- data.frame(m=0)
}
)
if(nrow(datam1)>2){
datam1 <- within(datam1, mutation_type[mutation_type != 'Silent'] <- '1')
datam1 <- within(datam1, mutation_type[mutation_type == 'Silent'] <- '0')
cellmut <- unique(as.character(datam1$case_id))
colnames(datam1) <- c("Name","Cell_line","Mutflag")
gene1mut <- subset(datam1, Name==gene1 & Mutflag==1, select="Cell_line") #' Gene1 mutated samples
gene1mut <- unique(as.character(gene1mut$Cell_line))
for(j in 1:length(slg)){
gene2 <- slg[j]
data1 <- data.frame(m=0)
tryCatch(
{names(mycnadata)[1] <- "Hugo_Symbol"
data1 <- subset(mycnadata, Hugo_Symbol==gene2)
data1 <- as.data.frame(t(data1[,-1]))
},
error=function(cond){
data1 <- data.frame(m=0)
message(paste("Error while reading CNA data:",cond,"\n","Please type help(customDataTTestAmp) to see the appropriate format of input data"))
},
warning=function(cond){
data1 <- data.frame(m=0)
}
)
if(nrow(data1)>2){
colnames(data1) <- "CNA"
allcells <- row.names(data1)
cmutl <- data.frame(Name=rep(gene2,nrow(data1)),Cell_line=allcells,Mutflag=data1$CNA)
gene2ampgene1mut <- subset(cmutl, Name==gene2 & Cell_line %in% gene1mut, select="Mutflag")
emut <- as.numeric(gene2ampgene1mut[,1])
gene2ampgene1nonmut <- subset(cmutl, Name==gene2 & !(Cell_line %in% gene1mut), select="Mutflag")
enmut <- as.numeric(gene2ampgene1nonmut[,1])
ttestamp[j] <- 1
if(length(emut)>2 && length(enmut)>2){
try(tp <- t.test(emut,y=enmut,alternative="greater"))
try(ttestamp[j] <- tp$p.value)
}
}
}
}
df <- data.frame(Gene1=rep(gene1,length(slg)),Gene2=slg,TTestMutAmp=ttestamp)
df
}


#' Calculate relative drug sensitivity (difference in drug sensitivity) for user provided mutation and drug screening data
#'
#' Given the primary gene of interest and the interactor gene list, it calculates the difference in Drug activity score (IC50) for targeting the interactor gene in cases where the primary gene is mutated vs not mutated
#' @param gene1 character. gene symbol for the primary gene. e.g. "BRCA2"
#' @param gene2 character vector. It can be either the output gene list from findInteractor, or a user-supplied gene list. list of gene symbols for the interactor gene. e.g. c("PARP1","ATM","NF1","EGFR")
#' @param celmutdata data frame. gene mutation data in tumor cell lines in the three column format: (gene_symbol, cell_line_id, variant_classification).
#' @param celldrug data frame. Drug response data (IC50 values) in tumor cell lines in a matrix format. The rows in the matrix will correspond to genes and the columns will correspond to tumor cell lines.
#' @return A data frame with following fields: Gene2(interactor gene), Gene1(primary gene), Drug (Drug name), PvalDrug(p-value for gene2 having less RNAi score in cells where gene1 is mutated compared to cells where gene1 is not mutated)
#' @export
#' @examples
#' gene <- findInteractor()
#' celmutdata
#'        gene_symbol    case_id mutation
#' 1         AKT3 DMS53_LUNG        0
#' 2         ABI1 DMS53_LUNG        0
#' 3         CDH2 DMS53_LUNG        0
#' 4 LOC100130776 DMS53_LUNG        0
#' 5        HDAC6 DMS53_LUNG        0
#' 6      BCL2L11 DMS53_LUNG        0
#' celldrug
#' 		Drug	target_gene	COSMIC_906800	COSMIC_910921	COSMIC_687562
#' 1	Erlotinib	EGFR	0.6922463	       2.5381818	2.5602048
#' 2	Sunitinib	PDGFRB	-1.6216855	       -1.7048357   -2.8186215
#' 3	Sunitinib	FLT1	-1.4392042		   -1.2772873   -2.1917983
#' 4	Sunitinib	KIT		-1.3452530		   -1.1004493   -2.9640288
#' 5	Sunitinib	FLT3	-0.7926535		   -0.3372055    0.2775138
#' 6	Sunitinib	CSF1R	-1.4126040		   -0.7056490    1.3978028
#' drugscore <- calculateDrugSensitivity("BRCA2",gene,celmutdata,celldrug)

customDataDrugSensitivity <- function(gene,gene2,celmutdata,celldrug){
tryCatch(
{mdrug <- reshape2::melt(celldrug)
names(mdrug) <- c("DrugName","gene_symbol","case_id","LN_IC50")
names(celmutdata) <- c("gene_symbol","case_id","mutation")
},
error=function(cond){
mdrug <- data.frame(m=0)
message(paste("Error while reading Drug screening data:",cond,"\n","Please type help(customDataDrugSensitivity) to see the appropriate format of input data"))
},
warning=function(cond){
}
)
shalt <- subset(celmutdata, gene_symbol==gene1 & mutation==1, select=case_id)
shalt <- unique(shalt)
shunalt <- subset(celmutdata, gene_symbol==gene1 & mutation==0, select=case_id)
shunalt <- unique(shunalt)
rnaimut <- subset(mdrug, case_id %in% shalt[,1])
rnainon <- subset(mdrug, case_id %in% shunalt[,1])
prgene <- gene
slgene <- unique(gene2)
sldrug <- unique(as.character(mdrug$DrugName))
slen <- length(slgene) + length(sldrug)
pval <- numeric(slen)
gene1 <- character(slen)
gene2 <- character(slen)
meanmut <- numeric(slen)
meannonmut <- numeric(slen)
meanfc <- numeric(slen)
l=1;
drug <- character(slen)
for(i in 1:length(slgene)){
emutalt <- subset(rnaimut, Gene2==slgene[i])
enmutunalt <- subset(rnainon, Gene2==slgene[i])
unidrug <- intersect(unique(as.character(emutalt$DrugName)), unique(as.character(enmutunalt$DrugName)))
unidrug <- unique(unidrug)
for(k in 1:length(unidrug)){
mutsub <- subset(emutalt, DrugName==unidrug[k])
nonmutsub <- subset(enmutunalt, DrugName==unidrug[k])
emut <- as.numeric(mutsub$LN_IC50)
enmut <- as.numeric(nonmutsub$LN_IC50)
if(length(emut)>2 && length(enmut)>2){
gene1[l] <- prgene
gene2[l] <- slgene[i]
try(tp <- t.test(emut,y=enmut,alternative="greater"))
try(pval[l] <- tp$p.value)
meanmut[l] <- mean(emut)
meannonmut[l] <- mean(enmut)
drug[l] <- unidrug[k]
l = l + 1
}
}
}
d <- data.frame(Gene1=gene1,Gene2=gene2,Drug=drug,PvalDrug=pval)
d <- subset(d, Gene1!="" & Gene2!="" & Drug!="")
d
}
