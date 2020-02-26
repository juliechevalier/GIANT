# A command-line interface to plot heatmap based on expression or diff. exp. analysis 
# written by Jimmy Vandel
# one of these arguments is required:
#
#
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
source(file.path(script.basename, "utils.R"))
source(file.path(script.basename, "getopt.R"))

#addComment("Welcome R!")

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat(geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
loc <- Sys.setlocale("LC_NUMERIC", "C")

#get starting time
start.time <- Sys.time()


options(stringAsfactors = FALSE, useFancyQuotes = FALSE, OutDec=".")

#get options
args <- commandArgs()

# get options, using the spec as defined by the enclosed list.
# we read the options   from the default: commandArgs(TRUE).
spec <- matrix(c(
  "expressionFile", "x", 1, "character",
  "diffAnalyseFile", "x", 1, "character",
  "factorInfo","x", 1, "character",
  "genericData","x", 0, "logical",
  "comparisonName","x",1,"character",
  "comparisonNameLow","x",1,"character",
  "comparisonNameHigh","x",1,"character",
  "filterInputOutput","x", 1, "character",
  "FCthreshold","x", 1, "double",
  "pvalThreshold","x", 1, "double",
  "geneListFiltering","x",1,"character",
  "clusterNumber","x",1,"integer",
  "maxRows","x",1,"integer",
  "sampleClusterNumber","x",1,"integer",
  "dataTransformation","x",1,"character",
  "distanceMeasure","x",1,"character",
  "aggloMethod","x",1,"character",
  "personalColors","x",1,"character",
  "sideBarColorPalette","x",1,"character",
  "format", "x", 1, "character",
  "quiet", "x", 0, "logical",
  "log", "x", 1, "character",
  "outputFile" , "x", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getoptLong(spec)

# enforce the following required arguments
if (is.null(opt$log)) {
  addComment("[ERROR]'log file' is required")
  q( "no", 1, F )
}
addComment("[INFO]Start of R script",T,opt$log,display=FALSE)
if (is.null(opt$format)) {
  addComment("[ERROR]'output format' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$outputFile)) {
  addComment("[ERROR]'output file' is required",T,opt$log)
  q( "no", 1, F )
}

if(is.null(opt$expressionFile) && !is.null(opt$genericData)){
  addComment("[ERROR]generic data clustering is based on expression clustering",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$clusterNumber) || opt$clusterNumber<2) {
  addComment("[ERROR]valid genes clusters number is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$sampleClusterNumber) || opt$sampleClusterNumber<1) {
  addComment("[ERROR]valid samples clusters number is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$dataTransformation)) {
  addComment("[ERROR]data transformation option is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$distanceMeasure)) {
  addComment("[ERROR]distance measure option is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$aggloMethod)) {
  addComment("[ERROR]agglomeration method option is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$maxRows) || opt$maxRows<2) {
  addComment("[ERROR]valid plotted row number is required",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt[["comparisonName"]]) && nchar(opt[["comparisonName"]])==0){
  addComment("[ERROR]you have to specify comparison",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt$comparisonNameLow) && nchar(opt$comparisonNameLow)==0){
  addComment("[ERROR]you have to specify comparisonLow",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt$comparisonNameHigh) && nchar(opt$comparisonNameHigh)==0){
  addComment("[ERROR]you have to specify comparisonHigh",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$genericData) && (!is.null(opt$comparisonNameLow) || !is.null(opt$comparisonNameHigh))){
  addComment("[ERROR]comparisonLow and comparisonHigh can be specified only with generic data",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt$genericData) && !is.null(opt[["comparisonName"]])){
  addComment("[ERROR]basic comparison cannot be specified for generic data",T,opt$log)
  q( "no", 1, F )
}

if ((!is.null(opt[["comparisonName"]]) || !is.null(opt$comparisonNameLow) || !is.null(opt$comparisonNameHigh)) && is.null(opt$diffAnalyseFile)) {
  addComment("[ERROR]'diff. exp. analysis file' is required",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt$genericData) && !is.null(opt$diffAnalyseFile) && is.null(opt$comparisonNameLow) && is.null(opt$comparisonNameHigh)){
  addComment("[ERROR]Missing comparison information for filtering",T,opt$log)
  q( "no", 1, F )
}

if ((!is.null(opt$FCthreshold) || !is.null(opt$pvalThreshold)) && (is.null(opt[["comparisonName"]]) && is.null(opt$comparisonNameLow) && is.null(opt$comparisonNameHigh))) {
  addComment("[ERROR]'comparisons' are missing for filtering",T,opt$log)
  q( "no", 1, F )
}

if ((!is.null(opt$FCthreshold) || !is.null(opt$pvalThreshold)) && !is.null(opt$geneListFiltering)) {
  addComment("[ERROR]Cannot have two filtering strategies",T,opt$log)
  q( "no", 1, F )
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
}else{
  FALSE}

addComment("[INFO]Parameters checked!",T,opt$log,display=FALSE)

addComment(c("[INFO]Working directory: ",getwd()),TRUE,opt$log,display=FALSE)
addComment(c("[INFO]Command line: ",args),TRUE,opt$log,display=FALSE)

#directory for plots and HTML
dir.create(file.path(getwd(), "plotDir"))
dir.create(file.path(getwd(), "plotLyDir"))

#silent package loading
suppressPackageStartupMessages({
  library("plotly")
  library("dendextend")
  #library("ggdendro")
  #library("plyr")
  library("ggplot2")
  library("heatmaply")
  library("circlize")
  #library("RColorBrewer")
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("ComplexHeatmap")
  library("ComplexHeatmap")
  #library("processx")
})

expressionToCluster=!is.null(opt$expressionFile)

#load input data files
if(expressionToCluster){
  #first expression data
  expressionMatrix=read.csv(file=opt$expressionFile,header=F,sep="\t",colClasses="character")
  #remove first row to convert it as colnames (to avoid X before colnames with header=T)
  colNamesData=expressionMatrix[1,-1]
  expressionMatrix=expressionMatrix[-1,]
  #remove first colum to convert it as rownames
  rowNamesData=expressionMatrix[,1]
  expressionMatrix=expressionMatrix[,-1]
  if(is.data.frame(expressionMatrix)){
    expressionMatrix=data.matrix(expressionMatrix)
  }else{
    expressionMatrix=data.matrix(as.numeric(expressionMatrix))
  }
  dimnames(expressionMatrix)=list(rowNamesData,colNamesData)
  
  #check input files
  if (!is.numeric(expressionMatrix)) {
    addComment("[ERROR]Expression data is not fully numeric!",T,opt$log,display=FALSE)
    q( "no", 1, F )
  }
  
  addComment("[INFO]Expression data loaded and checked")
  addComment(c("[INFO]Dim of expression matrix:",dim(expressionMatrix)),T,opt$log,display=FALSE)
}

nbComparisons=0
nbColPerContrast=5
comparisonMatrix=NULL
comparisonMatrixInfoGene=NULL
#if available comparisons
if(!is.null(opt[["comparisonName"]])){
    #load results from differential expression analysis
    #consider first row contains column names
    comparisonMatrix=read.csv(file=opt$diffAnalyseFile,header=F,sep="\t")
    colnames(comparisonMatrix)=as.character(unlist(comparisonMatrix[1,]))
    #remove the second line also as it's information line (p-val,FDR.p-val,FC,logFC)
    comparisonMatrix=comparisonMatrix[-c(1,2),]
    #remove first and second colums, convert the first one as rownames
    rownames(comparisonMatrix)=as.character(unlist(comparisonMatrix[,1]))
    #and save second column content that contain geneInfo
    comparisonMatrixInfoGene=as.character(unlist(comparisonMatrix[,2]))
    names(comparisonMatrixInfoGene)=as.character(unlist(comparisonMatrix[,1]))
    comparisonMatrix=comparisonMatrix[,-c(1,2)]
    
    comparisonMatrix=matrix(as.numeric(as.matrix(comparisonMatrix)),ncol=ncol(comparisonMatrix),dimnames = dimnames(comparisonMatrix))
    
    if (ncol(comparisonMatrix)%%nbColPerContrast != 0) {
      addComment("[ERROR]Diff. exp. data does not contain good number of columns per contrast, should contains in this order:p-val,FDR.p-val,FC,log2(FC) and t-stat",T,opt$log,display=FALSE)
      q( "no", 1, F )
    }
    
    if(max(comparisonMatrix[,c(seq(1,ncol(comparisonMatrix),nbColPerContrast),seq(2,ncol(comparisonMatrix),nbColPerContrast))])>1 || min(comparisonMatrix[,c(seq(1,ncol(comparisonMatrix),nbColPerContrast),seq(2,ncol(comparisonMatrix),nbColPerContrast))])<0){
      addComment("[ERROR]Seem that diff. exp. data does not contain correct values for p-val and FDR.p-val columns, should be including in [0,1] interval",T,opt$log,display=FALSE)
      q( "no", 1, F )
    }
    
    if (!is.numeric(comparisonMatrix)) {
      addComment("[ERROR]Diff. exp. data is not fully numeric!",T,opt$log,display=FALSE)
      q( "no", 1, F )
    }
    
    if(expressionToCluster && length(setdiff(rownames(comparisonMatrix),rownames(expressionMatrix)))!=0){
      addComment("[WARNING]All genes from diff. exp. file are not included in expression file",T,opt$log,display=FALSE)
    }
    
    if(expressionToCluster && length(setdiff(rownames(expressionMatrix),rownames(comparisonMatrix)))!=0){
      addComment("[WARNING]All genes from expression file are not included in diff. exp. file",T,opt$log,display=FALSE)
    }
    
    addComment("[INFO]Diff. exp. analysis loaded and checked",T,opt$log,display=FALSE)
    addComment(c("[INFO]Dim of original comparison matrix:",dim(comparisonMatrix)),T,opt$log,display=FALSE)
    
    #restrict to user specified comparisons
    restrictedComparisons=unlist(strsplit(opt[["comparisonName"]],","))
    #should be improved to avoid selection of column names starting too similarly  
    colToKeep=which(unlist(lapply(colnames(comparisonMatrix),function(x)any(startsWith(x,restrictedComparisons)))))
    comparisonMatrix=matrix(comparisonMatrix[,colToKeep],ncol=length(colToKeep),dimnames = list(rownames(comparisonMatrix),colnames(comparisonMatrix)[colToKeep]))
    
    #get number of required comparisons
    nbComparisons=ncol(comparisonMatrix)/nbColPerContrast
    
    addComment(c("[INFO]Dim of effective filtering matrix:",dim(comparisonMatrix)),T,opt$log,display=FALSE)
}

#should be only the case with generic data
if(!is.null(opt$comparisonNameLow) || !is.null(opt$comparisonNameHigh)){
    #load generic data used for filtering
    nbColPerContrast=1
    #consider first row contains column names
    comparisonMatrix=read.csv(file=opt$diffAnalyseFile,header=F,sep="\t")
    colnames(comparisonMatrix)=as.character(unlist(comparisonMatrix[1,]))
    #remove first colum, convert the first one as rownames
    rownames(comparisonMatrix)=as.character(unlist(comparisonMatrix[,1]))
    comparisonMatrix=comparisonMatrix[-1,-1]
    
    comparisonMatrix=matrix(as.numeric(as.matrix(comparisonMatrix)),ncol=ncol(comparisonMatrix),dimnames = dimnames(comparisonMatrix))
    
    if (!is.numeric(comparisonMatrix)) {
      addComment("[ERROR]Filtering matrix is not fully numeric!",T,opt$log,display=FALSE)
      q( "no", 1, F )
    }
    
    if(expressionToCluster && length(setdiff(rownames(comparisonMatrix),rownames(expressionMatrix)))!=0){
      addComment("[WARNING]All genes from filtering file are not included in expression file",T,opt$log,display=FALSE)
    }
    
    if(expressionToCluster && length(setdiff(rownames(expressionMatrix),rownames(comparisonMatrix)))!=0){
      addComment("[WARNING]All genes from expression file are not included in filtering file",T,opt$log,display=FALSE)
    }
    
    addComment("[INFO]Filtering file loaded and checked",T,opt$log,display=FALSE)
    addComment(c("[INFO]Dim of original filtering matrix:",dim(comparisonMatrix)),T,opt$log,display=FALSE)
    
    #restrict to user specified comparisons
    restrictedComparisons=c()
    if(!is.null(opt$comparisonNameLow))restrictedComparisons=unique(c(restrictedComparisons,unlist(strsplit(opt$comparisonNameLow,","))))
    if(!is.null(opt$comparisonNameHigh))restrictedComparisons=unique(c(restrictedComparisons,unlist(strsplit(opt$comparisonNameHigh,","))))
    
    if (!all(restrictedComparisons%in%colnames(comparisonMatrix))){
      addComment("[ERROR]Selected columns in filtering file are not present in filtering matrix!",T,opt$log,display=FALSE)
      q( "no", 1, F )
    }
    comparisonMatrix=matrix(comparisonMatrix[,restrictedComparisons],ncol=length(restrictedComparisons),dimnames = list(rownames(comparisonMatrix),restrictedComparisons))
    
    #get number of required comparisons
    nbComparisons=ncol(comparisonMatrix)
    
    addComment(c("[INFO]Dim of effective filtering matrix:",dim(comparisonMatrix)),T,opt$log,display=FALSE)
}



factorInfoMatrix=NULL
if(!is.null(opt$factorInfo)){
  #get group information
  #load factors file
  factorInfoMatrix=read.csv(file=opt$factorInfo,header=F,sep="\t",colClasses="character")
  #remove first row to convert it as colnames
  colnames(factorInfoMatrix)=factorInfoMatrix[1,]
  factorInfoMatrix=factorInfoMatrix[-1,]
  #use first colum to convert it as rownames but not removing it to avoid conversion as vector in unique factor case
  rownames(factorInfoMatrix)=factorInfoMatrix[,1]
  
  factorBarColor=colnames(factorInfoMatrix)[2]
  
  if(ncol(factorInfoMatrix)>2){
    addComment("[ERROR]Factors file should not contain more than 2 columns",T,opt$log,display=FALSE)
    q( "no", 1, F )
  }
  
  #factor file is used for color band on heatmap, so all expression matrix column should be in the factor file
  if(expressionToCluster && length(setdiff(colnames(expressionMatrix),rownames(factorInfoMatrix)))!=0){
    addComment("[ERROR]Missing samples in factor file",T,opt$log,display=FALSE)
    q( "no", 1, F )
  }
  
  #factor file is used for color band on heatmap, so all comparison matrix column should be in the factor file
  if(!expressionToCluster && length(setdiff(colnames(comparisonMatrix),rownames(factorInfoMatrix)))!=0){
    addComment("[ERROR]Missing differential contrasts in factor file",T,opt$log,display=FALSE)
    q( "no", 1, F )
  }
  
  addComment("[INFO]Factors OK",T,opt$log,display=FALSE)
  addComment(c("[INFO]Dim of factorInfo matrix:",dim(factorInfoMatrix)),T,opt$log,display=FALSE)
}

if(!is.null(opt$personalColors)){
 ##parse personal colors
  personalColors=unlist(strsplit(opt$personalColors,","))
  if(length(personalColors)==2){
    ##add medium color between two to get three colors
    personalColors=c(personalColors[1],paste(c("#",as.character(as.hexmode(floor(apply(col2rgb(personalColors),1,mean))))),collapse=""),personalColors[2])
  }
  if(length(personalColors)!=3){
    addComment("[ERROR]Personalized colors doesn't contain enough colors",T,opt$log,display=FALSE)
    q( "no", 1, F )
  }
    
}


if(!is.null(opt$filterInputOutput) && opt$filterInputOutput=="input"){
  #filter input data
  
    if(is.null(opt$geneListFiltering)){
      #filtering using stat thresholds
      #rowToKeep=intersect(which(comparisonMatrix[,seq(2,ncol(comparisonMatrix),4)]<=opt$pvalThreshold),which(abs(comparisonMatrix[,seq(4,ncol(comparisonMatrix),4)])>=log2(opt$FCthreshold)))
      if(is.null(opt$genericData)){
        #diff. expression matrix
        rowToKeep=names(which(unlist(apply(comparisonMatrix,1,function(x)length(intersect(which(x[seq(2,length(x),nbColPerContrast)]<opt$pvalThreshold),which(abs(x[seq(4,length(x),nbColPerContrast)])>log2(opt$FCthreshold))))!=0))))
      }else{
        #generic filtering matrix
        rowToKeep=rownames(comparisonMatrix)
        if(!is.null(opt$comparisonNameLow)){
          restrictedLowComparisons=unlist(strsplit(opt$comparisonNameLow,","))
          rowToKeep=intersect(rowToKeep,names(which(unlist(apply(comparisonMatrix,1,function(x)length(which(x[restrictedLowComparisons]>opt$FCthreshold))!=0)))))
        }
        if(!is.null(opt$comparisonNameHigh)){
          restrictedHighComparisons=unlist(strsplit(opt$comparisonNameHigh,","))
          rowToKeep=intersect(rowToKeep,names(which(unlist(apply(comparisonMatrix,1,function(x)length(which(x[restrictedHighComparisons]<opt$pvalThreshold))!=0)))))
        }
      }
    }else{
      #filtering using user gene list
      geneListFiltering=read.csv(opt$geneListFiltering,as.is = 1,header=F)
      rowToKeep=unlist(c(geneListFiltering))
    }
    
    if(!is.null(comparisonMatrix) && !all(rowToKeep%in%rownames(comparisonMatrix))){
      #should arrive only with user gene list filtering with diff.exp. results clustering
      addComment("[WARNING] some genes of the user defined list are not in the diff. exp. input file",T,opt$log)
      rowToKeep=intersect(rowToKeep,rownames(comparisonMatrix))
    }
  
    if(expressionToCluster && !all(rowToKeep%in%rownames(expressionMatrix))){
      addComment("[WARNING] some genes selected by the input filter are not in the expression file",T,opt$log)
      rowToKeep=intersect(rowToKeep,rownames(expressionMatrix))
    }
  
    if(length(rowToKeep)==0){
      addComment("[ERROR]No gene survived to the input filtering thresholds, execution will be aborted.
                 Please consider to change threshold values and re-run the tool.",T,opt$log)
      q( "no", 1, F )
    }

    #filter comparison matrix 
    if(!is.null(comparisonMatrix)){
      comparisonMatrix=matrix(comparisonMatrix[rowToKeep,],ncol=ncol(comparisonMatrix),dimnames = list(rowToKeep,colnames(comparisonMatrix)))
      if(!is.null(comparisonMatrixInfoGene))comparisonMatrixInfoGene=comparisonMatrixInfoGene[rowToKeep]
    }
    #then expression matrix
    if(expressionToCluster)expressionMatrix=matrix(expressionMatrix[rowToKeep,],ncol=ncol(expressionMatrix),dimnames = list(rowToKeep,colnames(expressionMatrix)))

    if(!is.null(comparisonMatrix) && expressionToCluster && nrow(comparisonMatrix)!=nrow(expressionMatrix)){
      addComment("[ERROR]Problem during input filtering, please check code",T,opt$log,display=FALSE)
      q( "no", 1, F )
    }
    
    addComment("[INFO]Filtering step done",T,opt$log,display=FALSE)
    addComment(c("[INFO]Input filtering step:",length(rowToKeep),"remaining rows"),T,opt$log,display=FALSE)
}


addComment("[INFO]Ready to plot",T,opt$log,display=FALSE)

##---------------------

#plot heatmap
if(expressionToCluster){
  #will make clustering based on expression value or generic value
  dataToHeatMap=expressionMatrix
  valueMeaning="Intensity"
  if(!is.null(opt$genericData))valueMeaning="Value"
}else{
  #will make clustering on log2(FC) values
  dataToHeatMap=matrix(comparisonMatrix[,seq(4,ncol(comparisonMatrix),nbColPerContrast)],ncol=nbComparisons,dimnames = list(rownames(comparisonMatrix),colnames(comparisonMatrix)[seq(1,ncol(comparisonMatrix),nbColPerContrast)]))
  valueMeaning="Log2(FC)"
}
addComment(c("[INFO]Dim of heatmap matrix:",dim(dataToHeatMap)),T,opt$log,display=FALSE)

if(nrow(dataToHeatMap)==1 && ncol(dataToHeatMap)==1){
  addComment("[ERROR]Cannot make clustering with unique cell tab",T,opt$log,display=FALSE)
  q( "no", 1, F )
}


#apply data transformation if needed
if(opt$dataTransformation=="log"){
  dataToHeatMap=log(dataToHeatMap)
  valueMeaning=paste(c("log(",valueMeaning,")"),collapse="")
  addComment("[INFO]Data to cluster and to display in the heatmap are log transformed",T,opt$log,display=FALSE)
}
if(opt$dataTransformation=="log2"){
  dataToHeatMap=log2(dataToHeatMap)
  valueMeaning=paste(c("log2(",valueMeaning,")"),collapse="")
  addComment("[INFO]Data to cluster and to display in the heatmap are log2 transformed",T,opt$log,display=FALSE)
}

maxRowsToDisplay=opt$maxRows

nbClusters=opt$clusterNumber
if(nbClusters>nrow(dataToHeatMap)){
  #correct number of clusters if needed
  nbClusters=nrow(dataToHeatMap)
  addComment(c("[WARNING]Not enough rows to reach required clusters number, it is reduced to number of rows:",nbClusters),T,opt$log,display=FALSE)
}

nbSampleClusters=opt$sampleClusterNumber
if(nbSampleClusters>ncol(dataToHeatMap)){
  #correct number of clusters if needed
  nbSampleClusters=ncol(dataToHeatMap)
  addComment(c("[WARNING]Not enough columns to reach required conditions clusters number, it is reduced to number of columns:",nbSampleClusters),T,opt$log,display=FALSE)
}

colClust=FALSE
rowClust=FALSE
effectiveRowClust=FALSE

#make appropriate clustering if needed
if(nrow(dataToHeatMap)>1 && nbClusters>1)rowClust=hclust(distExtended(dataToHeatMap,method = opt$distanceMeasure),method = opt$aggloMethod)
if(ncol(dataToHeatMap)>1 && nbSampleClusters>1)colClust=hclust(distExtended(t(dataToHeatMap),method = opt$distanceMeasure),method = opt$aggloMethod)

if(nrow(dataToHeatMap)>maxRowsToDisplay){
  #make subsampling based on preliminary global clustering
  #clusteringResults=cutree(rowClust,nbClusters)
  #heatMapGenesToKeep=unlist(lapply(seq(1,nbClusters),function(x)sample(which(clusteringResults==x),min(length(which(clusteringResults==x)),round(maxRowsToDisplay/nbClusters)))))
  ##OR
  #basic subsampling
  heatMapGenesToKeep=sample(rownames(dataToHeatMap),maxRowsToDisplay)
  effectiveDataToHeatMap=matrix(dataToHeatMap[heatMapGenesToKeep,],ncol=ncol(dataToHeatMap),dimnames=list(heatMapGenesToKeep,colnames(dataToHeatMap)))
  effectiveNbClusters=min(nbClusters,maxRowsToDisplay)
  if(nrow(effectiveDataToHeatMap)>1 && effectiveNbClusters>1)effectiveRowClust=hclust(distExtended(effectiveDataToHeatMap, method = opt$distanceMeasure),method = opt$aggloMethod)
  addComment(c("[WARNING]Too many rows for efficient heatmap drawing",maxRowsToDisplay,"subsampling is done for vizualization only"),T,opt$log,display=FALSE)
  rm(heatMapGenesToKeep)
}else{
  effectiveDataToHeatMap=dataToHeatMap
  effectiveRowClust=rowClust 
  effectiveNbClusters=nbClusters
}

addComment(c("[INFO]Dim of plotted heatmap matrix:",dim(effectiveDataToHeatMap)),T,opt$log,display=FALSE)

personalized_hoverinfo=matrix("",ncol = ncol(effectiveDataToHeatMap),nrow = nrow(effectiveDataToHeatMap),dimnames = dimnames(effectiveDataToHeatMap))
if(expressionToCluster){
  for(iCol in colnames(effectiveDataToHeatMap)){for(iRow in rownames(effectiveDataToHeatMap)){personalized_hoverinfo[iRow,iCol]=paste(c("Probe: ",iRow,"\nCondition: ",iCol,"\n",valueMeaning,": ",effectiveDataToHeatMap[iRow,iCol]),collapse="")}}
}else{
  for(iCol in colnames(effectiveDataToHeatMap)){for(iRow in rownames(effectiveDataToHeatMap)){personalized_hoverinfo[iRow,iCol]=paste(c("Probe: ",iRow,"\nCondition: ",iCol,"\nFC: ",round(2^effectiveDataToHeatMap[iRow,iCol],2)),collapse="")}}
}

#trying to overcome limitation of heatmaply package to modify xtick and ytick label, using directly plotly functions, but for now plotly do not permit to have personalized color for each x/y tick separately
test=FALSE
if(test==TRUE){
  
  #define dendogram shapes
  dd.row <- as.dendrogram(effectiveRowClust)
  dd.col <- as.dendrogram(colClust)
  
  #and color them
  dd.row=color_branches(dd.row, k = effectiveNbClusters, groupLabels = T)
  dd.col=color_branches(dd.col, k = nbSampleClusters, groupLabels = T)
  
  #generating function for dendogram from segment list
  ggdend <- function(df) {
    ggplot() +
      geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
      labs(x = "", y = "") + theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
  }
  
  # generate x/y dendogram plots
  px <- ggdend(dendro_data(dd.col)$segments)
  py <- ggdend(dendro_data(dd.row)$segments) + coord_flip()
  
  # reshape data matrix
  col.ord <- order.dendrogram(dd.col)
  row.ord <- order.dendrogram(dd.row)
  xx <- effectiveDataToHeatMap[row.ord, col.ord]
  # and also personalized_hoverinfo
  personalized_hoverinfo=personalized_hoverinfo[row.ord, col.ord]
  
  # hide axis ticks and grid lines
  eaxis <- list(
    showticklabels = FALSE,
    showgrid = FALSE,
    zeroline = FALSE
  )
  
  #make the empty plot
  p_empty <- plot_ly() %>%
    layout(margin = list(l = 200),
           xaxis = eaxis,
           yaxis = eaxis)
  
  heatmap.plotly <- plot_ly(
    z = xx, x = 1:ncol(xx), y = 1:nrow(xx), colors = viridis(n = 101, alpha = 1, begin = 0, end = 1, option = "inferno"),
    type = "heatmap", showlegend = FALSE, text = personalized_hoverinfo, hoverinfo = "text",
    colorbar = list(
      # Capitalise first letter
      title = valueMeaning,
      tickmode = "array",
      len = 0.3
    )
  ) %>%
    layout(
      xaxis = list(
        tickfont = list(size = 10,color=get_leaves_branches_col(dd.row)),
        tickangle = 45,
        tickvals = 1:ncol(xx), ticktext = colnames(xx),
        linecolor = "#ffffff",
        range = c(0.5, ncol(xx) + 0.5),
        showticklabels = TRUE
      ),
      yaxis = list(
        tickfont = list(size = 10, color=get_leaves_branches_col(dd.col)),
        tickangle = 0,
        tickvals = 1:nrow(xx), ticktext = rownames(xx),
        linecolor = "#ffffff",
        range = c(0.5, nrow(xx) + 0.5),
        showticklabels = TRUE
      )
    )
  
  #generate plotly 
  pp <- subplot(px, p_empty, heatmap.plotly, py, nrows = 2, margin = 0,widths = c(0.8,0.2),heights = c(0.2,0.8), shareX = TRUE, 
                shareY = TRUE)
  
  #save image file
  export(pp, file =  paste(c(file.path(getwd(), "plotDir"),"/Heatmap.",opt$format),collapse=""))
  #rise a bug due to token stuf
  #orca(pp, file =  paste(c(file.path(getwd(), "plotDir"),"/Heatmap.",opt$format),collapse=""))
  
  
  #save plotLy file
  htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/Heatmap.html"),collapse=""),selfcontained = F)
  
  #htmlwidgets::saveWidget(as_widget(pp),"~/Bureau/test.html",selfcontained = F)
  
}else{ #test
  label_names=c("Probe","Condition",valueMeaning)
  
  # #color hclust objects
  # dd.row=color_branches(effectiveRowClust, k = effectiveNbClusters)
  # #rowColors=get_leaves_branches_col(dd.row)
  # #rowColors[order.dendrogram(dd.row)]=rowColors
  # rowGroup=cutree(effectiveRowClust, k = effectiveNbClusters)
  # 
  # #get order of class as they will be displayed on the dendogram
  # rowGroupRenamed=data.frame(cluster=mapvalues(rowGroup, unique(rowGroup[order.dendrogram(dd.row)[nleaves(dd.row):1]]), 1:effectiveNbClusters))
  #
  #  dd.col=color_branches(colClust, k = nbSampleClusters)
  #  #colColors=get_leaves_branches_col(dd.col)
  #  #colColors[order.dendrogram(dd.col)]=colColors
  #  colGroup=cutree(colClust, k = nbSampleClusters)
  #  
  # # #get order of class as they will be displayed on the dendogram
  #  colGroupRenamed=data.frame(sampleCluster=mapvalues(colGroup, unique(colGroup[order.dendrogram(dd.col)[nleaves(dd.col):1]]), 1:nbSampleClusters))
  
  
  #while option is not correctly managed by heatmap apply, put personalized_hoverinfo to NULL
  personalized_hoverinfo=NULL
  
  if(is.null(opt$personalColors)){
    heatmapColors=viridis(n = 101, alpha = 1, begin = 0, end = 1, option = "inferno")
  }else{
    heatmapColors=personalColors
  }
  
  colGroupRenamed=NULL
  if(!is.null(factorInfoMatrix)){
    colGroupRenamed=eval(parse(text=(paste("data.frame(",factorBarColor,"=factorInfoMatrix[colnames(effectiveDataToHeatMap),2])",sep=""))))
    sideBarGroupNb=length(table(factorInfoMatrix[colnames(effectiveDataToHeatMap),2]))
    sideBarColorPaletteName="Spectral"
    if(!is.null(opt$sideBarColorPalette) && opt$sideBarColorPalette%in%rownames(RColorBrewer::brewer.pal.info)){
      sideBarColorPaletteName=opt$sideBarColorPalette
    }
    sideBarColorPalette=setNames(colorRampPalette(RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[sideBarColorPaletteName,"maxcolors"], sideBarColorPaletteName))(sideBarGroupNb),unique(factorInfoMatrix[colnames(effectiveDataToHeatMap),2]))
  }
  
  if(!is.null(colGroupRenamed)){
    pp <- heatmaply(effectiveDataToHeatMap,key.title = valueMeaning,k_row=effectiveNbClusters,k_col=nbSampleClusters,col_side_colors=colGroupRenamed,col_side_palette=sideBarColorPalette,Rowv=effectiveRowClust,Colv=colClust,label_names=label_names,custom_hovertext=personalized_hoverinfo,plot_method = "plotly",colors = heatmapColors)
  }else{
    pp <- heatmaply(effectiveDataToHeatMap,key.title = valueMeaning,k_row=effectiveNbClusters,k_col=nbSampleClusters,Rowv=effectiveRowClust,Colv=colClust,label_names=label_names,custom_hovertext=personalized_hoverinfo,plot_method = "plotly",colors = heatmapColors)
  }
  
  
  #save image file
  export(pp, file =  paste(c(file.path(getwd(), "plotDir"),"/Heatmap.",opt$format),collapse=""))
  #rise a bug due to token stuf
  #orca(pp, file =  paste(c(file.path(getwd(), "plotDir"),"/Heatmap.",opt$format),collapse=""))
  
  
  #save plotLy file
  htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/Heatmap.html"),collapse=""),selfcontained = F)
  
}
addComment("[INFO]Heatmap drawn",T,opt$log,display=FALSE)  


#plot circular heatmap
if(!class(effectiveRowClust)=="logical"){
  dendo=as.dendrogram(effectiveRowClust)
  
  if(is.null(opt$personalColors)){
    col_fun = colorRamp2(quantile(effectiveDataToHeatMap,probs = seq(0,1,0.01)), viridis(101,option = "inferno"))
  }else{
    col_fun = colorRamp2(quantile(effectiveDataToHeatMap,probs = seq(0,1,0.5)), personalColors)
  }
  
  if(opt$format=="pdf"){
    pdf(paste(c("./plotDir/circularPlot.pdf"),collapse=""))}else{
      png(paste(c("./plotDir/circularPlot.png"),collapse=""))
    }
  
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5)
  circos.initialize(c(rep("a",nrow(effectiveDataToHeatMap)),"b"),xlim=cbind(c(0,0),c(nrow(effectiveDataToHeatMap),5)))
  circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
    if(CELL_META$sector.index=="a"){
      nr = ncol(effectiveDataToHeatMap)
      nc = nrow(effectiveDataToHeatMap)
      circos.text(1:nc- 0.5, rep(0,nc), adj = c(0, 0), 
                  rownames(effectiveDataToHeatMap)[order.dendrogram(dendo)], facing = "clockwise", niceFacing = TRUE, cex = 0.3)
    }
  })
  
  circos.track(ylim = c(0, ncol(effectiveDataToHeatMap)), bg.border = NA, panel.fun = function(x, y) {
    
    m = t(matrix(effectiveDataToHeatMap[order.dendrogram(dendo),],ncol=ncol(effectiveDataToHeatMap)))
    col_mat = col_fun(m)
    nr = nrow(m)
    nc = ncol(m)
    if(CELL_META$sector.index=="a"){
      for(i in 1:nr) {
        circos.rect(1:nc - 1, rep(nr - i, nc), 
                    1:nc, rep(nr - i + 1, nc), 
                    border = col_mat[i, ], col = col_mat[i, ])
      }
    }else{
      circos.text(rep(1,nr), seq(nr,1,-1) , colnames(effectiveDataToHeatMap),cex = 0.3)
    }
  })
  
  #dendo = color_branches(dendo, k = effectiveNbClusters, col = colorRampPalette(brewer.pal(12,"Set3"))(effectiveNbClusters))
  dendo = color_branches(dendo, k = effectiveNbClusters, col = rev(colorspace::rainbow_hcl(effectiveNbClusters)))
  
  
  circos.track(ylim = c(0, attributes(dendo)$height), bg.border = NA, track.height = 0.25, 
               panel.fun = function(x, y) {
                 if(CELL_META$sector.index=="a")circos.dendrogram(dendo)} )
  
  circos.clear()
  ##add legend
  lgd_links = Legend(at = seq(ceiling(min(effectiveDataToHeatMap)),floor(max(effectiveDataToHeatMap)),ceiling((floor(max(effectiveDataToHeatMap))-ceiling(min(effectiveDataToHeatMap)))/4)), col_fun = col_fun, 
                     title_position = "topleft", grid_width = unit(5, "mm") ,title = valueMeaning)
  
  pushViewport(viewport(x = 0.85, y = 0.80, 
                        width = 0.1, 
                        height = 0.1, 
                        just = c("left", "bottom")))
  grid.draw(lgd_links)
  upViewport()
  
  
  dev.off()
  
  addComment("[INFO]Circular heatmap drawn",T,opt$log,display=FALSE)  
  loc <- Sys.setlocale("LC_NUMERIC","C")
}else{
  addComment(c("[WARNING]Circular plot will not be plotted considering row or cluster number < 2"),T,opt$log,display=FALSE)
}
rm(effectiveDataToHeatMap,effectiveRowClust,effectiveNbClusters)

#plot screeplot 
if(class(rowClust)!="logical" && nrow(dataToHeatMap)>2){
  screePlotData=c()
  for(iNbClusters in 2:(nbClusters+min(10,max(0,nrow(dataToHeatMap)-nbClusters)))){
    clusteringResults=cutree(rowClust,iNbClusters)
    #clusteringResults=kmeans(dataToHeatMap,iNbClusters)$cluster
    
    #compute variance between each intra-class points amongst themselves (need at least 3 points by cluster)
    #screePlotData=c(screePlotData,sum(unlist(lapply(seq(1,iNbClusters),function(x){temp=which(clusteringResults==x);if(length(temp)>2){var(dist(dataToHeatMap[temp,]))}else{0}}))) )
    #compute variance between each intra-class points and fictive mean point (need at least 2 points by cluster)
    #screePlotData=c(screePlotData,sum(unlist(lapply(seq(1,iNbClusters),function(x){temp=which(clusteringResults==x);if(length(temp)>1){   var(dist(rbind(apply(dataToHeatMap[temp,],2,mean),dataToHeatMap[temp,]))[1:length(temp)]) }else{0}}))) )
    if(ncol(dataToHeatMap)>1)screePlotData=c(screePlotData,sum(unlist(lapply(seq(1,iNbClusters),function(x){temp=which(clusteringResults==x);if(length(temp)>1){   sum((distExtended(rbind(apply(dataToHeatMap[temp,],2,mean),dataToHeatMap[temp,]),method = opt$distanceMeasure)[1:length(temp)])^2) }else{0}}))) )
    else screePlotData=c(screePlotData,sum(unlist(lapply(seq(1,iNbClusters),function(x){temp=which(clusteringResults==x);if(length(temp)>1){   sum((dataToHeatMap[temp,]-mean(dataToHeatMap[temp,]))^2) }else{0}}))) )
  }
  
  dataToPlot=data.frame(clusterNb=seq(2,length(screePlotData)+1),wcss=screePlotData)
  p <- ggplot(data=dataToPlot, aes(clusterNb,wcss)) + geom_point(colour="#EE4444") + geom_line(colour="#DD9999") +
    ggtitle("Scree plot") + theme_bw() + xlab(label="Cluster number") + ylab(label="Within cluster sum of squares") + 
    theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none") +
    scale_x_continuous(breaks=seq(min(dataToPlot$clusterNb), max(dataToPlot$clusterNb), 1))
  
  #save plotly files   
  pp <- ggplotly(p)
  
  if(opt$format=="pdf"){
    pdf(paste(c("./plotDir/screePlot.pdf"),collapse=""))}else{
      png(paste(c("./plotDir/screePlot.png"),collapse=""))
    }
  plot(p)
  dev.off()
  
  #save plotly files 
  htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/screePlot.html"),collapse=""),selfcontained = F)
  
  addComment("[INFO]Scree plot drawn",T,opt$log,display=FALSE)  
}else{
  addComment(c("[WARNING]Scree plot will not be plotted considering row number <= 2"),T,opt$log,display=FALSE)
}

##----------------------
  
#filter output based on parameters

rowToKeep=rownames(dataToHeatMap)
if(!is.null(opt$filterInputOutput) && opt$filterInputOutput=="output"){
  #rowToKeep=intersect(which(comparisonMatrix[,seq(2,ncol(comparisonMatrix),4)]<=opt$pvalThreshold),which(abs(comparisonMatrix[,seq(4,ncol(comparisonMatrix),4)])>=log2(opt$FCthreshold)))
  if(is.null(opt$geneListFiltering)){
    if(is.null(opt$genericData)){
      #diff. expression matrix
      rowToKeep=names(which(unlist(apply(comparisonMatrix,1,function(x)length(intersect(which(x[seq(2,length(x),nbColPerContrast)]<=opt$pvalThreshold),which(abs(x[seq(4,length(x),nbColPerContrast)])>=log2(opt$FCthreshold))))!=0))))
    }else{
      #generic filtering matrix
      rowToKeep=rownames(comparisonMatrix)
      if(!is.null(opt$comparisonNameLow)){
        restrictedLowComparisons=unlist(strsplit(opt$comparisonNameLow,","))
        rowToKeep=intersect(rowToKeep,names(which(unlist(apply(comparisonMatrix,1,function(x)length(which(x[restrictedLowComparisons]>opt$FCthreshold))!=0)))))
      }
      if(!is.null(opt$comparisonNameHigh)){
        restrictedHighComparisons=unlist(strsplit(opt$comparisonNameHigh,","))
        rowToKeep=intersect(rowToKeep,names(which(unlist(apply(comparisonMatrix,1,function(x)length(which(x[restrictedHighComparisons]<opt$pvalThreshold))!=0)))))
      }
    }
  }else{
    geneListFiltering=read.csv(opt$geneListFiltering,as.is = 1,header=F)
    rowToKeep=unlist(c(geneListFiltering))
  }
  if(!is.null(comparisonMatrix) && !all(rowToKeep%in%rownames(comparisonMatrix))){
    #should arrive only with user gene list filtering with diff.exp. results clustering
    addComment("[WARNING] some genes of the user defined list are not in the diff. exp. input file",T,opt$log)
    rowToKeep=intersect(rowToKeep,rownames(comparisonMatrix))
  }
  
  if(expressionToCluster && !all(rowToKeep%in%rownames(expressionMatrix))){
    addComment("[WARNING] some genes selected by the output filter are not in the expression file",T,opt$log)
    rowToKeep=intersect(rowToKeep,rownames(expressionMatrix))
  }
  addComment(c("[INFO]Output filtering step:",length(rowToKeep),"remaining rows"),T,opt$log,display=FALSE) 
}

#we add differential analysis info in output if it was directly used for clustering or when it was used for filtering with expression

#in case of expression or generic data clustering without filtering based on external stats
if(expressionToCluster && is.null(comparisonMatrix)){
  if(length(rowToKeep)==0){
    addComment("[WARNING]No more gene after output filtering step, tabular output will be empty",T,opt$log,display=FALSE)
    outputData=matrix(c("Gene","Cluster","noGene","noClustering"),ncol=2,nrow=2,byrow = TRUE)
  }else{
      outputData=matrix(0,ncol=2,nrow=length(rowToKeep)+1)
      outputData[1,]=c("Gene","Cluster")
      outputData[2:(length(rowToKeep)+1),1]=rowToKeep
      if(class(rowClust)!="logical" ){
        outputData[2:(length(rowToKeep)+1),2]=cutree(rowClust,nbClusters)[rowToKeep]
      }else{
        outputData[2:(length(rowToKeep)+1),2]=0
      }
  }
}

#in case of generic data clustering with filtering based on generic external data
if(!is.null(opt$genericData) && !is.null(comparisonMatrix)){
  if(length(rowToKeep)==0){
    addComment("[WARNING]No more gene after output filtering step, tabular output will be empty",T,opt$log,display=FALSE)
    outputData=matrix(c("Gene","Cluster","noGene","noClustering"),ncol=2,nrow=2,byrow = TRUE)
  }else{
    outputData=matrix(0,ncol=2+nbComparisons,nrow=length(rowToKeep)+1)
    outputData[1,]=c("Gene","Cluster",colnames(comparisonMatrix))
    outputData[2:(length(rowToKeep)+1),1]=rowToKeep
    if(class(rowClust)!="logical" ){
      outputData[2:(length(rowToKeep)+1),2]=cutree(rowClust,nbClusters)[rowToKeep]
    }else{
      outputData[2:(length(rowToKeep)+1),2]=0
    }
    outputData[2:(length(rowToKeep)+1),3:(ncol(comparisonMatrix)+2)]=prettyNum(comparisonMatrix[rowToKeep,],digits=4)
  }
}

#in case of expression data clustering with filtering based on diff. exp. results or diff. exp. results clustering
if(is.null(opt$genericData) && !is.null(comparisonMatrix)){
  if(length(rowToKeep)==0){
    addComment("[WARNING]No more gene after output filtering step, tabular output will be empty",T,opt$log,display=FALSE)
    outputData=matrix(0,ncol=3,nrow=3)
    outputData[1,]=c("","","Comparison")
    outputData[2,]=c("Gene","Info","Cluster")
    outputData[3,]=c("noGene","noInfo","noClustering")
  }else{
      outputData=matrix(0,ncol=3+nbComparisons*nbColPerContrast,nrow=length(rowToKeep)+2)
      outputData[1,]=c("","","Comparison",rep(colnames(comparisonMatrix)[seq(1,ncol(comparisonMatrix),nbColPerContrast)],each=nbColPerContrast))
      outputData[2,]=c("Gene","Info","Cluster",rep(c("p-val","FDR.p-val","FC","log2(FC)","t-stat"),nbComparisons))
      outputData[3:(length(rowToKeep)+2),1]=rowToKeep
      outputData[3:(length(rowToKeep)+2),2]=comparisonMatrixInfoGene[rowToKeep]
      if(class(rowClust)!="logical" ){
        outputData[3:(length(rowToKeep)+2),3]=cutree(rowClust,nbClusters)[rowToKeep]
      }else{
        outputData[3:(length(rowToKeep)+2),3]=0
      }
      outputData[3:(length(rowToKeep)+2),4:(ncol(comparisonMatrix)+3)]=prettyNum(comparisonMatrix[rowToKeep,],digits=4)
  }
}

addComment("[INFO]Formated output",T,opt$log,display=FALSE) 
write.table(outputData,file=opt$outputFile,quote=FALSE,sep="\t",col.names = F,row.names = F)
  
##----------------------

end.time <- Sys.time()
addComment(c("[INFO]Total execution time for R script:",as.numeric(end.time - start.time,units="mins"),"mins"),T,opt$log,display=FALSE)


addComment("[INFO]End of R script",T,opt$log,display=FALSE)

printSessionInfo(opt$log)

#sessionInfo()



