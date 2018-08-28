# A command-line interface to plot heatmap based on expression or diff. exp. analysis 
# written by Jimmy Vandel
# one of these arguments is required:
#
#
# the output file has columns:
#
#
#
source("./utils.R")
source("./getopt.R")

addComment("Welcome R!")

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat(geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

#get starting time
start.time <- Sys.time()

#get options
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs()


# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "expressionFile", "i", 1, "character",
  "diffAnalyseFile", "d", 1, "character",
  "factorInfo","g", 1, "character",
  "comparisonName","c",1,"character",
  "filterInputOutput","r", 1, "character",
  "FCthreshold","p", 1, "double",
  "pvalThreshold","t", 1, "double",
  "geneListFiltering","e",1,"character",
  "clusterNumber","n",1,"integer",
  "maxRows","m",1,"integer",
  "sampleClustering","s",0,"logical",
  "personalColors","k",1,"character",
  "format", "f", 1, "character",
  "quiet", "q", 0, "logical",
  "log", "l", 1, "character",
  "outputFile" , "o", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# enforce the following required arguments
if (is.null(opt$log)) {
  addComment("'log file' is required")
  q( "no", 1, F )
}
if (is.null(opt$format)) {
  addComment("'output format' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$outputFile)) {
  addComment("'output file' is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$clusterNumber) || opt$clusterNumber<2) {
  addComment("'valid cluster number' is required",T,opt$log)
  q( "no", 1, F )
}

if (is.null(opt$maxRows) || opt$maxRows<2) {
  addComment("'valid plotted row number' is required",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt$comparisonName) && nchar(opt$comparisonName)==0){
  addComment("'you have to specify comparison",T,opt$log)
  q( "no", 1, F )
}

if (!is.null(opt$comparisonName) && is.null(opt$diffAnalyseFile)) {
  addComment("'diff. exp. analysis file' is required",T,opt$log)
  q( "no", 1, F )
}

if ((!is.null(opt$FCthreshold) || !is.null(opt$pvalThreshold)) && is.null(opt$comparisonName)) {
  addComment("'comparisons' are missing for filtering",T,opt$log)
  q( "no", 1, F )
}

if ((!is.null(opt$FCthreshold) || !is.null(opt$pvalThreshold)) && !is.null(opt$geneListFiltering)) {
  addComment("cannot have two filtering strategies",T,opt$log)
  q( "no", 1, F )
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
}else{
  FALSE}

addComment("Parameters checked!")

addComment(getwd(),TRUE,opt$log,display=FALSE)
addComment(args,TRUE,opt$log,display=FALSE)

#directory for plots and HTML
dir.create(file.path(getwd(), "plotDir"))
dir.create(file.path(getwd(), "plotLyDir"))

#silent package loading
suppressPackageStartupMessages({
  library("plotly")
  library("dendextend")
  library("ggplot2")
  library("heatmaply")
  library("circlize")
  #library("RColorBrewer")
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("ComplexHeatmap")
  library("ComplexHeatmap")
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
    addComment("expression data is not fully numeric!",T,opt$log)
    q( "no", 1, F )
  }
  
  addComment("Expression data loaded and checked")
  addComment(c("dim of expression matrix:",dim(expressionMatrix)),T,opt$log,display=F)
}


if(!is.null(opt$factorInfo)){
  #get group information
  #chargement du fichier des facteurs
  factorInfoMatrix=read.csv(file=opt$factorInfo,header=F,sep="\t",colClasses="character")
  #remove first row to convert it as colnames
  colnames(factorInfoMatrix)=factorInfoMatrix[1,]
  factorInfoMatrix=factorInfoMatrix[-1,]
  #use first colum to convert it as rownames but not removing it to avoid conversion as vector in unique factor case
  rownames(factorInfoMatrix)=factorInfoMatrix[,1]
  
  
  if(expressionToCluster && length(setdiff(colnames(expressionMatrix),rownames(factorInfoMatrix)))!=0){
    addComment("Missing samples in factor file",T,opt$log)
    q( "no", 1, F )
  }
  
  #order sample as in expression matrix and remove spurious sample
  if(expressionToCluster)factorInfoMatrix=factorInfoMatrix[colnames(expressionMatrix),]
  
  addComment("Factors OK")
  addComment(c("dim of factorInfo matrix:",dim(factorInfoMatrix)),T,opt$log,display=F)
}

nbComparisons=0
#if available comparisons
if(!is.null(opt$comparisonName)){
  
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
  
  if (!is.numeric(comparisonMatrix)) {
    addComment("diff. exp. data is not fully numeric!",T,opt$log)
    q( "no", 1, F )
  }
  
  if(expressionToCluster && length(setdiff(rownames(comparisonMatrix),rownames(expressionMatrix)))!=0){
    addComment("geneSet from diff. exp. file should be included in expression file",T,opt$log)
    q( "no", 1, F )
  }
  
  addComment("Diff. exp. analysis loaded and checked")
  addComment(c("dim of original comparison matrix:",dim(comparisonMatrix)),T,opt$log,display=F)
  
  #restrict to user specified comparisons
  restrictedComparisons=unlist(strsplit(opt$comparisonName,","))
  colToKeep=which(unlist(lapply(colnames(comparisonMatrix),function(x)any(startsWith(x,restrictedComparisons)))))
  comparisonMatrix=matrix(comparisonMatrix[,colToKeep],ncol=length(colToKeep),dimnames = list(rownames(comparisonMatrix),colnames(comparisonMatrix)[colToKeep]))
  
  addComment(c("dim of effective comparison matrix:",dim(comparisonMatrix)),T,opt$log,display=F)
  
  #get number of required comparisons
  nbComparisons=ncol(comparisonMatrix)/4 
}

if(!is.null(opt$personalColors)){
 ##parse personal colors
  personalColors=unlist(strsplit(opt$personalColors,","))
  if(length(personalColors)==2){
    ##add medium color between two to get three colors
    personalColors=c(personalColors[1],paste(c("#",as.character(as.hexmode(floor(apply(col2rgb(personalColors),1,mean))))),collapse=""),personalColors[2])
  }
  if(length(personalColors)!=3){
    addComment("personalized colors doesn't contain enough colors",T,opt$log)
    q( "no", 1, F )
  }
    
}


if(!is.null(opt$filterInputOutput) && opt$filterInputOutput=="input"){
  #filter input data
  
    if(is.null(opt$geneListFiltering)){
      #filtering using stat thresholds
      #rowToKeep=intersect(which(comparisonMatrix[,seq(2,ncol(comparisonMatrix),4)]<=opt$pvalThreshold),which(abs(comparisonMatrix[,seq(4,ncol(comparisonMatrix),4)])>=log2(opt$FCthreshold)))
      rowToKeep=names(which(unlist(apply(comparisonMatrix,1,function(x)length(intersect(which(x[seq(2,length(x),4)]<=opt$pvalThreshold),which(abs(x[seq(4,length(x),4)])>=log2(opt$FCthreshold))))!=0))))
    }else{
      #filtering using gene list
      geneListFiltering=read.csv(opt$geneListFiltering,as.is = 1,header=F)
      rowToKeep=unlist(c(geneListFiltering))
      if(!is.null(opt$comparisonName)){
        rowToKeep=intersect(rowToKeep,rownames(comparisonMatrix))
      }else{
        rowToKeep=intersect(rowToKeep,rownames(expressionMatrix))
      }
    }
    
    if(length(rowToKeep)==0){
      addComment("No gene survive from input filtering, execution will be aborted.",T,opt$log)
      q( "no", 1, F )
    }
    
    #filter comparison matrix 
    if(!is.null(opt$comparisonName)){
      comparisonMatrix=matrix(comparisonMatrix[rowToKeep,],ncol=ncol(comparisonMatrix),dimnames = list(rowToKeep,colnames(comparisonMatrix)))
      comparisonMatrixInfoGene=comparisonMatrixInfoGene[rowToKeep]
    }
    #then expression matrix
    if(expressionToCluster)expressionMatrix=matrix(expressionMatrix[rowToKeep,],ncol=ncol(expressionMatrix),dimnames = list(rowToKeep,colnames(expressionMatrix)))

    if(!is.null(opt$comparisonName) && expressionToCluster && nrow(comparisonMatrix)!=nrow(expressionMatrix)){
      addComment("Problem during input filtering, please check code",T,opt$log)
      q( "no", 1, F )
    }
    
    addComment("Filtering step done")
    addComment(c("input filtering step:",length(rowToKeep),"remaining rows"),T,opt$log,display=F)
}

addComment("Ready to plot")

##---------------------

#plot heatmap
if(expressionToCluster){
    #will make clustering based on expression value
    dataToHeatMap=expressionMatrix
    valueMeaning="Intensities"
  }else{
    #will make clustering on log2(FC) values
    dataToHeatMap=matrix(comparisonMatrix[,seq(4,ncol(comparisonMatrix),4)],ncol=nbComparisons,dimnames = list(rownames(comparisonMatrix),colnames(comparisonMatrix)[seq(1,ncol(comparisonMatrix),4)]))
    valueMeaning="Log2(FC)"
  }
  addComment(c("dim of heatmap matrix:",dim(dataToHeatMap)),T,opt$log,display=F)
  
  if(nrow(dataToHeatMap)==1 && ncol(dataToHeatMap)==1){
    addComment("cannot make clustering with unique cell tab",T,opt$log)
    q( "no", 1, F )
  }

  #heatMapGenesToKeep=genesToKeep
  #if(length(heatMapGenesToKeep)>5000){
  #  heatMapGenesToKeep=sample(genesToKeep,5000)
  #  addComment("WARNING limit signifiaddCommentive genes to 5000 for heatmap!")
  #}
  

  
  #define dendogram shapes
  #dd.col <- as.dendrogram(hclust(dist(dataToHeatMap)))
  #dd.row <- as.dendrogram(hclust(dist(t(dataToHeatMap))))
  #dx <- dendro_data(dd.row)
  #dy <- dendro_data(dd.col)
  
  #ggdend <- function(df) {
  #  ggplot() +
  #    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
  #    labs(x = "", y = "") + theme_minimal() +
  #    theme(axis.text = element_blank(), axis.ticks = element_blank(),
  #          panel.grid = element_blank())
  #}
  
  # x/y dendograms
  #px <- ggdend(dx$segments)
  #py <- ggdend(dy$segments) + coord_flip()
  
  # reshape data matrix
  #col.ord <- order.dendrogram(dd.col)
  #row.ord <- order.dendrogram(dd.row)
  #xx <- scale(dataToHeatMap)[col.ord, row.ord]
  #xx_names <- attr(xx, "dimnames")
  #df <- as.data.frame(xx)
  #colnames(df) <- xx_names[[2]]
  #df$gene <- xx_names[[1]]
  #df$gene <- with(df, factor(gene, levels=gene, ordered=TRUE))
  #mdf <- reshape2::melt(df, id.vars="gene")
  
  #plot heatmap
  #p <- ggplot(mdf, aes(x = variable, y = gene)) + geom_tile(aes(fill = value))
  
  #mat <- matrix(unlist(dplyr::select(df,-gene)),nrow=nrow(df))
  #colnames(mat) <- colnames(df)[1:ncol(df)-1]
  
  # hide axis ticks and grid lines
  #eaxis <- list(
  #  showticklabels = FALSE,
  #  showgrid = FALSE,
  #  zeroline = FALSE
  #)
  
  #make the empty plot
  # p_empty <- plot_ly() %>%
  # note that margin applies to entire plot, so we can
  # add it here to make tick labels more readable
  #  layout(margin = list(l = 200),
  #          xaxis = eaxis,
  #          yaxis = eaxis)
  
  #make the heatmap through plotLy
  #heatmap.plotly <- plot_ly() %>% add_heatmap(z=~xx,x=factor(colnames(xx),lev=colnames(xx)),y=factor(rownames(xx),lev=rownames(xx)))
  
  ###pp <- subplot(px, p_empty, p, py, nrows = 2, margin = 0.01)
  #pp <- subplot(px, p_empty, heatmap.plotly, py, nrows = 2, margin = 0.01)
  

  #heatmaply seems to make all the job by its own
  maxRowsToDisplay=opt$maxRows
  
  nbClusters=opt$clusterNumber
  if(nbClusters>nrow(dataToHeatMap)){
    #correct number of cluster if needed
    nbClusters=nrow(dataToHeatMap)
    addComment(c("WARNING not enough rows to reach required cluster number, it is reduced to number of rows:",nbClusters),T,opt$log)
  }
  
  colClust=FALSE
  rowClust=FALSE
  
  #make appropriate clustering if needed
  if(nrow(dataToHeatMap)>1)rowClust=hclust(dist(dataToHeatMap))
  if(ncol(dataToHeatMap)>1 && !is.null(opt$sampleClustering))colClust=hclust(dist(t(dataToHeatMap)))
  
  if(nrow(dataToHeatMap)>maxRowsToDisplay){
    #make subsampling based on preliminary global clustering
    #clusteringResults=cutree(rowClust,nbClusters)
    #heatMapGenesToKeep=unlist(lapply(seq(1,nbClusters),function(x)sample(which(clusteringResults==x),min(length(which(clusteringResults==x)),round(maxRowsToDisplay/nbClusters)))))
    ##OR
    #basic subsampling
    heatMapGenesToKeep=sample(rownames(dataToHeatMap),maxRowsToDisplay)
    effectiveDataToHeatMap=matrix(dataToHeatMap[heatMapGenesToKeep,],ncol=ncol(dataToHeatMap),dimnames=list(heatMapGenesToKeep,colnames(dataToHeatMap)))
    effectiveRowClust=hclust(dist(effectiveDataToHeatMap))
    effectiveNbClusters=min(nbClusters,maxRowsToDisplay)
    addComment(c("WARNING too many rows for efficient heatmap drawing",maxRowsToDisplay,"subsampling is done for vizualization only"),T,opt$log)
    rm(heatMapGenesToKeep)
  }else{
    effectiveDataToHeatMap=dataToHeatMap
    effectiveRowClust=rowClust 
    effectiveNbClusters=nbClusters
  }
  
  addComment(c("dim of plotted heatmap matrix:",dim(effectiveDataToHeatMap)),T,opt$log,display=F)
  
  personalized_hoverinfo=matrix("",ncol = ncol(effectiveDataToHeatMap),nrow = nrow(effectiveDataToHeatMap),dimnames = dimnames(effectiveDataToHeatMap))
  if(expressionToCluster){
    for(iCol in colnames(effectiveDataToHeatMap)){for(iRow in rownames(effectiveDataToHeatMap)){personalized_hoverinfo[iRow,iCol]=paste(c("Probe: ",iRow," | Condition: ",iCol,"\n Intensity: ",effectiveDataToHeatMap[iRow,iCol]),collapse="")}}
  }else{
    for(iCol in colnames(effectiveDataToHeatMap)){for(iRow in rownames(effectiveDataToHeatMap)){personalized_hoverinfo[iRow,iCol]=paste(c("Probe: ",iRow," | Condition: ",iCol,"\n FC: ",round(2^effectiveDataToHeatMap[iRow,iCol],2)),collapse="")}}
  }
  
  if(is.null(opt$personalColors)){
    pp <- heatmaply(effectiveDataToHeatMap,key.title = valueMeaning,k_row=effectiveNbClusters,Rowv=effectiveRowClust,Colv=colClust,custom_hovertext=personalized_hoverinfo,plot_method = "plotly",colors = viridis(n = 101, alpha = 1, begin = 0, end = 1, option = "inferno"),margins = c(9*max(unlist(lapply(colnames(effectiveDataToHeatMap),nchar))),9*max(unlist(lapply(rownames(effectiveDataToHeatMap),nchar))),0,0))
  }else{
    pp <- heatmaply(effectiveDataToHeatMap,key.title = valueMeaning,k_row=effectiveNbClusters,Rowv=effectiveRowClust,Colv=colClust,custom_hovertext=personalized_hoverinfo,plot_method = "plotly",colors = personalColors,margins = c(9*max(unlist(lapply(colnames(effectiveDataToHeatMap),nchar))),9*max(unlist(lapply(rownames(effectiveDataToHeatMap),nchar))),0,0))
  }
  
  
  #save image file
  export(pp, file =  paste(c(file.path(getwd(), "plotDir"),"/Heatmap.",opt$format),collapse=""))
  
  #save plotLy file
  htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/Heatmap.html"),collapse=""),selfcontained = F)
  
  addComment("Heatmap drawn")  
  
  
  
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
    
    addComment("Circular heatmap drawn")  
    
    rm(effectiveDataToHeatMap,effectiveRowClust,effectiveNbClusters)
  }
  
  
  
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
      if(ncol(dataToHeatMap)>1)screePlotData=c(screePlotData,sum(unlist(lapply(seq(1,iNbClusters),function(x){temp=which(clusteringResults==x);if(length(temp)>1){   sum((dist(rbind(apply(dataToHeatMap[temp,],2,mean),dataToHeatMap[temp,]))[1:length(temp)])^2) }else{0}}))) )
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
    
    addComment("Scree plot drawn")  
  }else{
    addComment(c("Scree plot will not be plotted considering row number <= 2"),T,opt$log)
  }

##----------------------
  
#filter output based on parameters

rowToKeep=rownames(dataToHeatMap)
if(!is.null(opt$filterInputOutput) && opt$filterInputOutput=="output"){
  #rowToKeep=intersect(which(comparisonMatrix[,seq(2,ncol(comparisonMatrix),4)]<=opt$pvalThreshold),which(abs(comparisonMatrix[,seq(4,ncol(comparisonMatrix),4)])>=log2(opt$FCthreshold)))
  
  if(is.null(opt$geneListFiltering)){
    rowToKeep=rownames(comparisonMatrix)[which(unlist(apply(comparisonMatrix,1,function(x)length(intersect(which(x[seq(2,length(x),4)]<=opt$pvalThreshold),which(abs(x[seq(4,length(x),4)])>=log2(opt$FCthreshold))))!=0)))]
  }else{
    geneListFiltering=read.csv(opt$geneListFiltering,as.is = 1)
    rowToKeep=unlist(c(geneListFiltering))
    if(!is.null(opt$comparisonName)){
      rowToKeep=intersect(rowToKeep,rownames(comparisonMatrix))
    }else{
      rowToKeep=intersect(rowToKeep,rownames(expressionMatrix))
    }
  }
  addComment(c("output filtering step:",length(rowToKeep),"remaining rows"),T,opt$log,display=F) 
}
  

#we add differential analysis info in output if it was directly used for clustering or when it was used for filtering with expression
addingDiffExpStat=nbComparisons>=1
  
#formating output matrix depending on genes to keep
if(length(rowToKeep)==0){
  addComment("No more gene after output filtering step, tabular output will be empty",T,opt$log)
  outputData=matrix(0,ncol=nbComparisons*4+3,nrow=3-!addingDiffExpStat)
  if(addingDiffExpStat)outputData[1,]=c("","","Comparison",rep(colnames(comparisonMatrix)[seq(1,ncol(comparisonMatrix),4)],each=4))
  outputData[2-!addingDiffExpStat,]=c("Gene","Info","Cluster",rep(c("p-val","FDR.p-val","FC","log2(FC)"),nbComparisons))
  outputData[3-!addingDiffExpStat,1:3]=c("noGene","noInfo","noClustering")
}else{
  if(addingDiffExpStat){
    outputData=matrix(0,ncol=nbComparisons*4+3,nrow=length(rowToKeep)+2)
    outputData[1,]=c("","","Comparison",rep(colnames(comparisonMatrix)[seq(1,ncol(comparisonMatrix),4)],each=4))
    outputData[2,]=c("Gene","Info","Cluster",rep(c("p-val","FDR.p-val","FC","log2(FC)"),nbComparisons))
    outputData[3:(length(rowToKeep)+2),1]=rowToKeep
    outputData[3:(length(rowToKeep)+2),2]=comparisonMatrixInfoGene[rowToKeep]
    if(class(rowClust)!="logical" ){outputData[3:(length(rowToKeep)+2),3]=cutree(rowClust,nbClusters)[rowToKeep]
    }else{outputData[3:(length(rowToKeep)+2),3]=0}
    outputData[3:(length(rowToKeep)+2),4:(ncol(comparisonMatrix)+3)]=prettyNum(comparisonMatrix[rowToKeep,],digits=4)
  }else{
    outputData=matrix(0,ncol=2,nrow=length(rowToKeep)+1)
    outputData[1,]=c("Gene","Cluster")
    outputData[2:(length(rowToKeep)+1),1]=rowToKeep
    if(class(rowClust)!="logical" ){outputData[2:(length(rowToKeep)+1),2]=cutree(rowClust,nbClusters)[rowToKeep]
    }else{outputData[2:(length(rowToKeep)+1),2]=0}
  }
}
addComment("Formated output") 
  
write.table(outputData,file=opt$outputFile,quote=FALSE,sep="\t",col.names = F,row.names = F)
  
##----------------------

end.time <- Sys.time()
addComment(c("Total execution time:",as.numeric(end.time - start.time,units="mins"),"mins"),T,opt$log,display=F)


addComment("End of R script",T,opt$log)

#sessionInfo()




