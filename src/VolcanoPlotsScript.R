# R script to plot volcanos through Galaxy based GIANT tool 
# written by Jimmy Vandel
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

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs()

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "statisticsFile", "i", 1, "character",
  "volcanoName" , "n", 1, "character",
  "pvalColumnName" , "p", 1, "character",
  "fdrColumnName" , "m", 1, "character",
  "fcColumnName" , "c", 1, "character",
  "fcKind","d", 1, "character",
  "fdrThreshold","s", 1, "double",
  "fcThreshold","e", 1, "double",
  "organismID","x",1,"character",
  "rowNameType","y",1,"character",
  "log", "l", 1, "character",
  "outputFile" , "o", 1, "character",
  "format", "f", 1, "character",
  "quiet", "q", 0, "logical"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# enforce the following required arguments
if (is.null(opt$log)) {
  addComment("[ERROR]'log file' is required\n")
  q( "no", 1, F )
}
addComment("[INFO]Start of R script",T,opt$log,display=FALSE)
if (is.null(opt$statisticsFile)) {
  addComment("[ERROR]'statisticsFile' is required",T,opt$log)
  q( "no", 1, F )
}
if (length(opt$pvalColumnName)==0 || length(opt$fdrColumnName)==0  || length(opt$fcColumnName)==0) {
  addComment("[ERROR]no selected columns",T,opt$log)
  q( "no", 1, F )
}
if (length(opt$pvalColumnName)!=length(opt$fcColumnName) || length(opt$pvalColumnName)!=length(opt$fdrColumnName)) {
  addComment("[ERROR]different number of selected columns between p.val, adj-p.val and FC ",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$fcKind)) {
  addComment("[ERROR]'fcKind' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$fdrThreshold)) {
  addComment("[ERROR]'FDR threshold' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$fcThreshold)) {
  addComment("[ERROR]'FC threshold' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$outputFile)) {
  addComment("[ERROR]'output file' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$format)) {
  addComment("[ERROR]'output format' is required",T,opt$log)
  q( "no", 1, F )
}

#demande si le script sera bavard
verbose <- if (is.null(opt$quiet)) {
  TRUE
}else{
  FALSE
}

#paramètres internes
addComment("[INFO]Parameters checked test mode !",T,opt$log,display=FALSE)

addComment(c("[INFO]Working directory: ",getwd()),TRUE,opt$log,display=FALSE)
addComment(c("[INFO]Command line: ",args),TRUE,opt$log,display=FALSE)

#directory for plots
dir.create(file.path(getwd(), "plotDir"))
dir.create(file.path(getwd(), "plotLyDir"))

#charge des packages silencieusement
suppressPackageStartupMessages({
  library("methods")
  library("biomaRt")
  library("ggplot2")
  library("plotly")
  library("stringr")
})

#define some usefull variable
nbVolcanosToPlot=length(opt$pvalColumnName)

#load input file
statDataMatrix=read.csv(file=file.path(getwd(), opt$statisticsFile),header=F,sep="\t",colClasses="character")
#remove first colum to convert it as rownames
rownames(statDataMatrix)=statDataMatrix[,1]
statDataMatrix=statDataMatrix[,-1]

#identify lines without adjusted p-value info (should contain the same content as rownames) and replace them with NA values
FDRinfo=rep(TRUE,nbVolcanosToPlot)
for(iVolcano in 1:nbVolcanosToPlot){
  #input parameter should be None when adjusted p-val are not available
  if(opt$fdrColumnName[iVolcano]=="None"){
    #content of the corresponding column should also be the same as rownames
    if(!all(statDataMatrix[,(iVolcano-1)*3+2]==rownames(statDataMatrix))){
      addComment(c("[ERROR]It seems that input stat matrix contains adjusted p-values for volcano",iVolcano,"whereas input parameter indicates that not."),T,opt$log)
      q( "no", 1, F )
    }
    FDRinfo[iVolcano]=FALSE
    statDataMatrix[,(iVolcano-1)*3+2]=NA
  }
}

if(is.data.frame(statDataMatrix)){
  statDataMatrix=data.matrix(statDataMatrix)
}else{
  statDataMatrix=data.matrix(as.numeric(statDataMatrix))
}

#check if available column number match with volcano requested number
if(ncol(statDataMatrix)!=3*nbVolcanosToPlot){
  addComment("[ERROR]Input file column number is different from requested volcano number",T,opt$log)
  q( "no", 1, F )
}

#build global dataFrame with data and fill with p.val and log2(FC) and FDR
dataFrame=data.frame(row.names = rownames(statDataMatrix))
#start with p-value
dataFrame$p.value=statDataMatrix[,seq(1,nbVolcanosToPlot*3,3),drop=FALSE]
#compute FDR if needed or just get available info
dataFrame$adj_p.value=dataFrame$p.value
for(iVolcano in 1:nbVolcanosToPlot){
  #adjusted p-value are already computed
  if(FDRinfo[iVolcano]){
    dataFrame$adj_p.value[,iVolcano]=statDataMatrix[,(iVolcano-1)*3+2,drop=FALSE]
  }else{
    #adjusted p-value should be computed based on p-val using FDR
    dataFrame$adj_p.value[,iVolcano]=p.adjust(dataFrame$p.value[,iVolcano,drop=FALSE],"fdr")
    addComment(c("[INFO]Adjusted p-values are not available in input for volcano",iVolcano,", FDR approach will be used on available raw p-values"),T,opt$log)
  } 
}
if(opt$fcKind=="FC"){
  #we should transform as Log2FC
  dataFrame$coefficients=log2(statDataMatrix[,seq(3,nbVolcanosToPlot*3,3),drop=FALSE])
  addComment(c("[INFO]FC are converted in log2(FC) for plotting"),T,opt$log)
}else{
  dataFrame$coefficients=statDataMatrix[,seq(3,nbVolcanosToPlot*3,3),drop=FALSE]
}

addComment(c("[INFO]Input data available for",nbVolcanosToPlot,"volcano(s) with",nrow(statDataMatrix),"rows"),T,opt$log)


#plot VOLCANOs
volcanoPerPage=1
logFCthreshold=log2(opt$fcThreshold)
iToPlot=1
plotVector=list()
volcanoNameList=c()
for (iVolcano in 1:nbVolcanosToPlot){
  
  if(nchar(opt$volcanoName[iVolcano])>0){
    curentVolcanoName=opt$volcanoName[iVolcano]
  }else{
    curentVolcanoName=paste(iVolcano,opt$pvalColumnName[iVolcano],sep="_")
  }
  
  #keep only rows without NA for p-val, adjusted p-val and coeff
  pValToPlot=dataFrame$p.value[,iVolcano]
  fdrToPlot=dataFrame$adj_p.value[,iVolcano]
  coeffToPlot=dataFrame$coefficients[,iVolcano]

  rowToRemove=unique(c(which(is.na(pValToPlot)),which(is.na(fdrToPlot)),which(is.na(coeffToPlot))))
  if(length(rowToRemove)>0){
    pValToPlot=pValToPlot[-rowToRemove]
    fdrToPlot=fdrToPlot[-rowToRemove]
    coeffToPlot=coeffToPlot[-rowToRemove]
  }
  addComment(c("[INFO]For",curentVolcanoName,"volcano,",length(rowToRemove),"rows are discarded due to NA values,",length(pValToPlot),"remaining rows."),T,opt$log)
  
  #save volcano name
  volcanoNameList=c(volcanoNameList,curentVolcanoName)
  
  #remove characters possibly troubling
  volcanoFileName=iVolcano
  
  #define the log10(p-val) threshold corresponding to FDR threshold fixed by user
  probeWithLowFDR=-log10(pValToPlot[which(fdrToPlot<=opt$fdrThreshold)])
  pvalThresholdFDR=NULL
  if(length(probeWithLowFDR)>0)pvalThresholdFDR=min(probeWithLowFDR)
  
  #get significant points over FC and FDR thresholds
  significativePoints=intersect(which(abs(coeffToPlot)>=logFCthreshold),which(fdrToPlot<=opt$fdrThreshold))
  
  #to reduce size of html plot, we keep 20000 points maximum sampled amongst genes with pval>=33%(pval) and abs(log2(FC))<=66%(abs(log2(FC)))
  htmlPointsToRemove=intersect(which(abs(coeffToPlot)<=quantile(abs(coeffToPlot),c(0.66))),which(pValToPlot>=quantile(abs(pValToPlot),c(0.33))))
  if(length(htmlPointsToRemove)>20000){
    htmlPointsToRemove=setdiff(htmlPointsToRemove,sample(htmlPointsToRemove,20000))
  }else{
    htmlPointsToRemove=c() 
  }

  xMinLimPlot=min(coeffToPlot)-0.2
  xMaxLimPlot=max(coeffToPlot)+0.2
  yMaxLimPlot= max(-log10(pValToPlot))+0.2
  
  if(length(significativePoints)>0){
    dataSignifToPlot=data.frame(pval=-log10(pValToPlot[significativePoints]),FC=coeffToPlot[significativePoints],description=paste(names(coeffToPlot[significativePoints]),"\n","FC: " , round(2^coeffToPlot[significativePoints],2) , " | Adjusted p-val: ",prettyNum(fdrToPlot[significativePoints],digits=4), sep=""))
    #to test if remains any normal points to draw
    if(length(significativePoints)<length(pValToPlot)){
      dataToPlot=data.frame(pval=-log10(pValToPlot[-significativePoints]),FC=coeffToPlot[-significativePoints],description=paste("FC: " , round(2^coeffToPlot[-significativePoints],2) , " | Adjusted p-val: ",prettyNum(fdrToPlot[-significativePoints],digits=4), sep=""))
    }else{
      dataToPlot=data.frame(pval=0,FC=0,description="null")
    }
  }else{
    dataToPlot=data.frame(pval=-log10(pValToPlot),FC=coeffToPlot,description=paste("FC: " , round(2^coeffToPlot,2) , " | Adjusted p-val: ",prettyNum(fdrToPlot,digits=4), sep=""))
  }
  
  ##traditional plot
  
  p <- ggplot(data=dataToPlot, aes(x=FC, y=pval)) + geom_point() + 
    theme_bw() + ggtitle(curentVolcanoName) + ylab(label="-Log10(p-val)") + xlab(label="Log2 Fold Change") +
    theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none")
  if(logFCthreshold!=0) p <- p + geom_vline(xintercept=-logFCthreshold, color="salmon",linetype="dotted", size=1) +  geom_vline(xintercept=logFCthreshold, color="salmon",linetype="dotted", size=1) + geom_text(data.frame(text=c(paste(c("log2(1/FC=",opt$fcThreshold,")"),collapse=""),paste(c("log2(FC=",opt$fcThreshold,")"),collapse="")),x=c(-logFCthreshold,logFCthreshold),y=c(0,0)),mapping=aes(x=x, y=y, label=text), size=4, angle=90, vjust=-0.4, hjust=0, color="salmon")
  if(!is.null(pvalThresholdFDR)) p <- p + geom_hline(yintercept=pvalThresholdFDR, color="skyblue1",linetype="dotted", size=0.5) + geom_text(data.frame(text=c(paste(c("Adjusted pval limit(",opt$fdrThreshold,")"),collapse="")),x=c(xMinLimPlot),y=c(pvalThresholdFDR)),mapping=aes(x=x, y=y, label=text), size=4, vjust=0, hjust=0, color="skyblue3")
  if(length(significativePoints)>0)p <- p + geom_point(data=dataSignifToPlot,aes(colour=description))
  
  ##interactive plot
  
  if(length(htmlPointsToRemove)>0){
    pointToRemove=union(htmlPointsToRemove,significativePoints)
    #to test if it remains any normal points to draw
    if(length(pointToRemove)<length(pValToPlot)){
      dataToPlot=data.frame(pval=-log10(pValToPlot[-pointToRemove]),FC=coeffToPlot[-pointToRemove],description=paste("FC: " , round(2^coeffToPlot[-pointToRemove],2) , " | Adjusted p-val: ", prettyNum(fdrToPlot[-pointToRemove],digits=4), sep=""))
    }else{
      dataToPlot=data.frame(pval=0,FC=0,description="null")
    }
  }
  
  if((nrow(dataToPlot)+length(significativePoints))>40000)addComment(c("[WARNING]For",curentVolcanoName,"volcano, numerous points to plot(",nrow(dataToPlot)+nrow(dataSignifToPlot),"), resulting volcano could be heavy, using more stringent thresholds could be helpful."),T,opt$log)
  
  phtml <- plot_ly(data=dataToPlot, x=~FC, y=~pval,type="scatter", mode="markers",showlegend = FALSE, marker = list(color="gray",opacity=0.5), text=~description, hoverinfo="text") %>%
    layout(title = curentVolcanoName[iVolcano],xaxis=list(title="Log2 Fold Change",showgrid=TRUE, zeroline=FALSE),yaxis=list(title="-Log10(p-val)", showgrid=TRUE, zeroline=FALSE))
  if(length(significativePoints)>0) phtml=add_markers(phtml,data=dataSignifToPlot, x=~FC, y=~pval, mode="markers" , marker=list( color=log10(abs(dataSignifToPlot$FC)*dataSignifToPlot$pval),colorscale='Rainbow'), text=~description, hoverinfo="text", inherit = FALSE) %>% hide_colorbar()
  if(logFCthreshold!=0){
    phtml=add_trace(phtml,x=c(-logFCthreshold,-logFCthreshold), y=c(0,yMaxLimPlot), type="scatter", mode = "lines", line=list(color="coral",dash="dash"), hoverinfo='none', showlegend = FALSE,inherit = FALSE)
    phtml=add_annotations(phtml,x=-logFCthreshold,y=0,xref = "x",yref = "y",text = paste(c("log2(1/FC=",opt$fcThreshold,")"),collapse=""),xanchor = 'right',showarrow = F,textangle=270,font=list(color="coral"))
    phtml=add_trace(phtml,x=c(logFCthreshold,logFCthreshold), y=c(0, yMaxLimPlot), type="scatter",  mode = "lines", line=list(color="coral",dash="dash"), hoverinfo='none', showlegend = FALSE,inherit = FALSE)
    phtml=add_annotations(phtml,x=logFCthreshold,y=0,xref = "x",yref = "y",text = paste(c("log2(FC=",opt$fcThreshold,")"),collapse=""),xanchor = 'right',showarrow = F,textangle=270,font=list(color="coral"))
  }
  if(!is.null(pvalThresholdFDR)){
    phtml=add_trace(phtml,x=c(xMinLimPlot,xMaxLimPlot), y=c(pvalThresholdFDR,pvalThresholdFDR), type="scatter",  mode = "lines", line=list(color="cornflowerblue",dash="dash"), hoverinfo='none', showlegend = FALSE,inherit = FALSE)
    phtml=add_annotations(phtml,x=xMinLimPlot,y=pvalThresholdFDR+0.1,xref = "x",yref = "y",text = paste(c("Adjusted pval limit(",opt$fdrThreshold,")"),collapse=""),xanchor = 'left',showarrow = F,font=list(color="cornflowerblue"))
  }
  plotVector[[length(plotVector)+1]]=p
  
  #save plotly files
  pp <- ggplotly(phtml)
  htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/Volcanos_",volcanoFileName,".html"),collapse=""),selfcontained = F)
  
  
  if(iVolcano==nbVolcanosToPlot || length(plotVector)==volcanoPerPage){
    #plot and close the actual plot
    if(opt$format=="pdf"){
      pdf(paste(c("./plotDir/Volcanos_",volcanoFileName,".pdf"),collapse=""))}else{
        png(paste(c("./plotDir/Volcanos_",volcanoFileName,".png"),collapse=""))
      }
    multiplot(plotlist=plotVector,cols=1)
    dev.off()
    if(iVolcano<nbVolcanosToPlot){
      #prepare for a new ploting file if necessary
      plotVector=list()
      iToPlot=iToPlot+1
    }
  }
}
remove(dataToPlot,dataSignifToPlot)
addComment("[INFO]Volcanos drawn",T,opt$log,T,display=FALSE)


#now add anotation infos about genes

rowItemInfo=NULL
if(!is.null(opt$rowNameType) && !is.null(opt$organismID)){
  ##get gene information from BioMart
  #if(!require("biomaRt")){
  #    source("https://bioconductor.org/biocLite.R")
  #    biocLite("biomaRt")
  #}
  
  ensembl_hs_mart <- useMart(biomart="ensembl", dataset=opt$organismID)
  ensembl_df <- getBM(attributes=c(opt$rowNameType,"description"),mart=ensembl_hs_mart)
  rowItemInfo=ensembl_df[which(ensembl_df[,1]!=""),2]
  rowItemInfo=unlist(lapply(rowItemInfo,function(x)substr(unlist(strsplit(x," \\[Source"))[1],1,30)))
  names(rowItemInfo)=ensembl_df[which(ensembl_df[,1]!=""),1]
}

#filter out genes with higher p-values for all comparisons
genesToKeep=names(which(apply(dataFrame$adj_p.value,1,function(x)length(which(x<=opt$fdrThreshold))>0)))
#filter out genes with lower FC for all comparisons
genesToKeep=intersect(genesToKeep,names(which(apply(dataFrame$coefficients,1,function(x)length(which(abs(x)>=logFCthreshold))>0))))

if(length(genesToKeep)>0){
  dataFrameNew=data.frame(row.names=genesToKeep)
  
  dataFrameNew$adj_p.value=matrix(dataFrame$adj_p.value[genesToKeep,,drop=FALSE],ncol=ncol(dataFrame$adj_p.value))
  rownames(dataFrameNew$adj_p.value)=genesToKeep
  colnames(dataFrameNew$adj_p.value)=colnames(dataFrame$p.value)
  
  dataFrameNew$p.value=matrix(dataFrame$p.value[genesToKeep,,drop=FALSE],ncol=ncol(dataFrame$p.value))
  rownames(dataFrameNew$p.value)=genesToKeep
  colnames(dataFrameNew$p.value)=colnames(dataFrame$adj_p.value)
  
  dataFrameNew$coefficients=matrix(dataFrame$coefficients[genesToKeep,,drop=FALSE],ncol=ncol(dataFrame$coefficients))
  rownames(dataFrameNew$coefficients)=genesToKeep
  colnames(dataFrameNew$coefficients)=colnames(dataFrame$adj_p.value)
  
  dataFrame=dataFrameNew
  rm(dataFrameNew)
}else{
  addComment("[WARNING]No significative genes",T,opt$log,display=FALSE)
}

addComment("[INFO]Significant genes filtering done",T,opt$log,T,display=FALSE)


#plot VennDiagramm for genes below threshold between comparisons
#t=apply(dataFrame$adj_p.value[,1:4],2,function(x)names(which(x<=opt$threshold)))
#get.venn.partitions(t)
#vennCounts(dataFrame$adj_p.value[,1:4]<=opt$threshold)

#make a simple sort genes based only on the first comparison
#newOrder=order(dataFrame$adj_p.value[,1])
#dataFrame$adj_p.value=dataFrame$adj_p.value[newOrder,]

#alternative sorting strategy based on the mean gene rank over all comparisons
if(length(genesToKeep)>1){
  currentRank=rep(0,nrow(dataFrame$adj_p.value))
  for(iVolcano in 1:ncol(dataFrame$adj_p.value)){
    currentRank=currentRank+rank(dataFrame$adj_p.value[,iVolcano])
  }
  currentRank=currentRank/ncol(dataFrame$adj_p.value)
  newOrder=order(currentRank)
  rownames(dataFrame)=rownames(dataFrame)[newOrder]
  
  dataFrame$adj_p.value=matrix(dataFrame$adj_p.value[newOrder,],ncol=ncol(dataFrame$adj_p.value))
  rownames(dataFrame$adj_p.value)=rownames(dataFrame$p.value)[newOrder]
  colnames(dataFrame$adj_p.value)=colnames(dataFrame$p.value)
  
  dataFrame$p.value=matrix(dataFrame$p.value[newOrder,],ncol=ncol(dataFrame$p.value))
  rownames(dataFrame$p.value)=rownames(dataFrame$adj_p.value)
  colnames(dataFrame$p.value)=colnames(dataFrame$adj_p.value)
  
  dataFrame$coefficients=matrix(dataFrame$coefficients[newOrder,],ncol=ncol(dataFrame$coefficients))
  rownames(dataFrame$coefficients)=rownames(dataFrame$adj_p.value)
  colnames(dataFrame$coefficients)=colnames(dataFrame$adj_p.value)
}

#formating output matrix depending on genes to keep
if(length(genesToKeep)==0){
  outputData=matrix(0,ncol=ncol(dataFrame$adj_p.value)*4+2,nrow=3)
  outputData[1,]=c("X","X",rep(volcanoNameList,each=4))
  outputData[2,]=c("X","X",rep(c("p-val","Adjusted.p-val","FC","log2(FC)"),ncol(dataFrame$adj_p.value)))
  outputData[,1]=c("Volcano","Gene","noGene")
  outputData[,2]=c("Comparison","Info","noInfo")
}else{
  if(length(genesToKeep)==1){
    outputData=matrix(0,ncol=ncol(dataFrame$adj_p.value)*4+2,nrow=3)
    outputData[1,]=c("X","X",rep(volcanoNameList,each=4))
    outputData[2,]=c("X","X",rep(c("p-val","Adjusted.p-val","FC","log2(FC)"),ncol(dataFrame$adj_p.value)))
    outputData[,1]=c("Volcano","Gene",genesToKeep)
    outputData[,2]=c("Comparison","Info","na")
    if(!is.null(rowItemInfo))outputData[3,2]=rowItemInfo[genesToKeep]
    outputData[3,seq(3,ncol(outputData),4)]=prettyNum(dataFrame$p.value,digits=4)
    outputData[3,seq(4,ncol(outputData),4)]=prettyNum(dataFrame$adj_p.value,digits=4)
    outputData[3,seq(5,ncol(outputData),4)]=prettyNum(2^dataFrame$coefficients,digits=4)
    outputData[3,seq(6,ncol(outputData),4)]=prettyNum(dataFrame$coefficients,digits=4)
  }else{
    #format matrix to be correctly read by galaxy (move headers in first column and row)
    outputData=matrix(0,ncol=ncol(dataFrame$adj_p.value)*4+2,nrow=nrow(dataFrame$adj_p.value)+2)
    outputData[1,]=c("X","X",rep(volcanoNameList,each=4))
    outputData[2,]=c("X","X",rep(c("p-val","Adjusted.p-val","FC","log2(FC)"),ncol(dataFrame$adj_p.value)))
    outputData[,1]=c("Volcano","Gene",rownames(dataFrame$adj_p.value))
    outputData[,2]=c("Comparison","Info",rep("na",nrow(dataFrame$adj_p.value)))
    if(!is.null(rowItemInfo))outputData[3:nrow(outputData),2]=rowItemInfo[rownames(dataFrame$adj_p.value)]
    outputData[3:nrow(outputData),seq(3,ncol(outputData),4)]=prettyNum(dataFrame$p.value,digits=4)
    outputData[3:nrow(outputData),seq(4,ncol(outputData),4)]=prettyNum(dataFrame$adj_p.value,digits=4)
    outputData[3:nrow(outputData),seq(5,ncol(outputData),4)]=prettyNum(2^dataFrame$coefficients,digits=4)
    outputData[3:nrow(outputData),seq(6,ncol(outputData),4)]=prettyNum(dataFrame$coefficients,digits=4)
  }
}
addComment("[INFO]Formated output",T,opt$log,display=FALSE) 

#write output results
write.table(outputData,file=opt$outputFile,quote=FALSE,sep="\t",col.names = F,row.names = F)


end.time <- Sys.time()
addComment(c("[INFO]Total execution time for R script:",as.numeric(end.time - start.time,units="mins"),"mins"),T,opt$log,display=FALSE)

addComment("[INFO]End of R script",T,opt$log,display=FALSE)

printSessionInfo(opt$log)

#sessionInfo()