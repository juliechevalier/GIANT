# A command-line interface to basic plots for use with Galaxy
# written by Jimmy Vandel
# one of these arguments is required:
#
#
# the output file has columns:
#
#
#
source("/data/galaxy-dist/tools/jimmytools/utils.R")
source("/data/galaxy-dist/tools/jimmytools/getopt.R")

addComment("Welcome R!")

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

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
  "dataFile", "i", 1, "character",
  "factorInfo","t", 1, "character",
  "dataFileFormat","j",1,"character",
  "conditionNames","c",1,"character",
  "format", "f", 1, "character",
  "quiet", "q", 0, "logical",
  "log", "l", 1, "character",
  "histo" , "h", 1, "character",
  "maPlot" , "a", 1, "character",
  "boxplot" , "b", 1, "character",
  "microarray" , "m", 1, "character",
  "acp" , "p" , 1, "character",
  "screePlot" , "s" , 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# enforce the following required arguments
if (is.null(opt$log)) {
  addComment("'log file' is required")
  q( "no", 1, F )
}
if (is.null(opt$dataFile) || is.null(opt$dataFileFormat)) {
  addComment("'dataFile' and it format are required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$format)) {
  addComment("'output format' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$histo) & is.null(opt$maPlot) & is.null(opt$boxplot) & is.null(opt$microarray) & is.null(opt$acp)){
  addComment("select at least one plot to draw",T,opt$log)
  q( "no", 1, F )
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
}else{
  FALSE}

if(verbose)addComment("parameters checked")

addComment(getwd(),TRUE,opt$log,display=FALSE)
addComment(args,TRUE,opt$log,display=FALSE)

#directory for plots
dir.create(file.path(getwd(), "plotDir"))
dir.create(file.path(getwd(), "plotLyDir"))

#silent package loading
suppressPackageStartupMessages({
  library("oligo")
  library("ff")
  library("ggplot2")
  library("plotly")
})


#chargement des fichiers en entrée
#fichier de type CEL
dataAreFromCel=FALSE
if(toupper(opt$dataFileFormat)=="CEL"){
  dataAreFromCel=TRUE
  celData=read.celfiles(unlist(strsplit(opt$dataFile,",")))
  #load all expressions
  dataMatrix=exprs(celData)
  #select "pm" probes
  probeInfo=getProbeInfo(celData,probeType = c("pm"),target="probeset")
  #reduce dataMatrix to log expression matrix for a randomly probe selection
  dataMatrix=log2(dataMatrix[sample(unique(probeInfo[,1]),min(100000,length(unique(probeInfo[,1])))),])
  remove(probeInfo)
}else{
  #fichier deja tabule
  dataMatrix=read.csv(file=opt$dataFile,header=F,sep="\t",colClasses="character")
  #remove first row to convert it as colnames (to avoid X before colnames with header=T)
  colNamesData=dataMatrix[1,-1]
  dataMatrix=dataMatrix[-1,]
  #remove first colum to convert it as rownames
  rowNamesData=dataMatrix[,1]
  dataMatrix=dataMatrix[,-1]
  if(is.data.frame(dataMatrix)){
    dataMatrix=data.matrix(dataMatrix)
  }else{
    dataMatrix=data.matrix(as.numeric(dataMatrix))
  }
  dimnames(dataMatrix)=list(rowNamesData,colNamesData)
}

if(verbose)addComment("Input data loaded")

#get number of conditions
nbConditions=ncol(dataMatrix)

#get condition names if they are specified
if(!is.null(opt$conditionNames) && length(opt$conditionNames)==nbConditions){
  nameConditions=opt$conditionNames
  colnames(dataMatrix)=nameConditions
  #rownames(phenoData(celData)@data)=nameConditions
  #rownames(protocolData(celData)@data)=nameConditions
}else{
  nameConditions=colnames(dataMatrix)
}

#create a correspondance table between plot file names and name displayed in figure legend and html items 
correspondanceNameTable=matrix("",ncol=2,nrow=nbConditions)
correspondanceNameTable[,1]=paste("Condition",1:nbConditions,sep="")
correspondanceNameTable[,2]=nameConditions
rownames(correspondanceNameTable)=correspondanceNameTable[,2]

if(verbose)addComment("Conditions names OK")

if(!is.null(opt$factorInfo)){
  #chargement du fichier des facteurs
  factorInfoMatrix=read.csv(file=file.path(getwd(), opt$factorInfo),header=F,sep="\t",colClasses="character")
  #remove first row to convert it as colnames
  colnames(factorInfoMatrix)=factorInfoMatrix[1,]
  factorInfoMatrix=factorInfoMatrix[-1,]
  #use first colum to convert it as rownames but not removing it to avoid conversion as vector in unique factor case
  rownames(factorInfoMatrix)=factorInfoMatrix[,1]
  
  
  if(length(setdiff(colnames(dataMatrix),rownames(factorInfoMatrix)))!=0){
    addComment("missing samples in factor file",T,opt$log)
    q( "no", 1, F )
  }
  
  #order sample as in expression matrix and remove spurious sample
  factorInfoMatrix=factorInfoMatrix[colnames(dataMatrix),]
  
  if(verbose)addComment("Factors OK")
  
}

if(verbose)addComment("Ready to plot")


##----------------------

###plot histograms###
histogramPerFigure=50
if (!is.null(opt$histo)) {
  for(iToPlot in 1:(((nbConditions-1)%/%histogramPerFigure)+1)){
    firstPlot=1+histogramPerFigure*(iToPlot-1)
    lastPlot=min(nbConditions,histogramPerFigure*iToPlot)
    dataToPlot=data.frame(x=c(dataMatrix[,firstPlot:lastPlot]),Experiment=rep(colnames(dataMatrix)[firstPlot:lastPlot],each=nrow(dataMatrix)))
    p <- ggplot(data=dataToPlot, aes(x = x, color=Experiment)) + stat_density(geom="line", size=1, position="identity") +
      ggtitle("Intensity densities") + theme_bw() + ylab(label="Density") + 
      theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5))
    if(dataAreFromCel){
      #original ploting function
      #hist(celData[,firstPlot:lastPlot],lty=rep(1,nbConditions)[firstPlot:lastPlot],lwd=2,which='pm',target="probeset",transfo=log2,col=rainbow(nbConditions)[firstPlot:lastPlot])
      p <- p + xlab(label="Log2 intensities") 
    }else{
      p <- p + xlab(label="Intensities") 
    }
    if(opt$format=="pdf"){
      pdf(paste(c("./plotDir/",opt$histo,iToPlot,".pdf"),collapse=""))}else{
        png(paste(c("./plotDir/",opt$histo,iToPlot,".png"),collapse=""))
      }
    print(p)
    dev.off()
    #save plotly files
    pp <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$histo,iToPlot,".html"),collapse=""),selfcontained = F)
  }
  remove(p,dataToPlot)
  if(verbose)addComment("Histograms drawn")
}

##----------------------

###plot MAplots###
MAplotPerPage=4
if (!is.null(opt$maPlot)) {
  iToPlot=1
  plotVector=list()
  toTake=sample(nrow(dataMatrix),min(200000,nrow(dataMatrix)))
  refMedianColumn=rowMedians(as.matrix(dataMatrix[toTake,]))
  for (iCondition in 1:nbConditions){
    #MAplot(celData,which=i,what=pm,transfo=log2)
    #smoothScatter(x=xToPlot,y=yToPlot,main=nameConditions[iCondition])
    dataA=dataMatrix[toTake,iCondition]
    dataB=refMedianColumn####ATTENTION PAR DEFAUT
    xToPlot=0.5*(dataA+dataB)
    yToPlot=dataA-dataB
    tempX=seq(min(xToPlot),max(xToPlot),0.1)
    tempY=unlist(lapply(tempX,function(x){median(yToPlot[intersect(which(xToPlot>=(x-0.1/2)),which(xToPlot<(x+0.1/2)))])}))
    
    dataToPlot=data.frame(x=xToPlot,y=yToPlot)
    dataMedianToPlot=data.frame(x=tempX,y=tempY)
    p <- ggplot(data=dataToPlot, aes(x,y)) + stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 100) +
      scale_fill_continuous(low = "white", high = "dodgerblue4") + geom_smooth(data=dataMedianToPlot,colour="red", size=0.5, se=FALSE) +
      ggtitle(correspondanceNameTable[iCondition,2]) + theme_bw() + xlab(label="") + ylab(label="") + 
      theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5),legend.position = "none")
    plotVector[[length(plotVector)+1]]=p
  
    #save plotly files   
    pp <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$maPlot,"_",correspondanceNameTable[iCondition,1],".html"),collapse=""),selfcontained = F)

    if(iCondition==nbConditions || length(plotVector)==MAplotPerPage){
      #define a new plotting file
      if(opt$format=="pdf"){
        pdf(paste(c("./plotDir/",opt$maPlot,iToPlot,".pdf"),collapse=""))}else{
          png(paste(c("./plotDir/",opt$maPlot,iToPlot,".png"),collapse=""))
        }
      multiplot(plotlist=plotVector,cols=2)
      dev.off()
      if(iCondition<nbConditions){
        #prepare for a new plotting file if necessary
        plotVector=list()
        iToPlot=iToPlot+1
      }
    }
  }
  remove(p,dataToPlot,dataA,dataB,toTake,xToPlot,yToPlot)
  if(verbose)addComment("MAplots drawn")
}

##----------------------

###plot boxplots###
boxplotPerFigure=50
if (!is.null(opt$boxplot)) {
  for(iToPlot in 1:(((nbConditions-1)%/%boxplotPerFigure)+1)){
    firstPlot=1+boxplotPerFigure*(iToPlot-1)
    lastPlot=min(nbConditions,boxplotPerFigure*iToPlot)
    dataToPlot=data.frame(intensities=c(dataMatrix[,firstPlot:lastPlot]),Experiment=rep(colnames(dataMatrix)[firstPlot:lastPlot],each=nrow(dataMatrix)))
    #to make HTML file lighter, sampling will be done amongst outliers
    #get outliers for each boxplot
    boxplotsOutliers=apply(dataMatrix[,firstPlot:lastPlot],2,function(x)boxplot.stats(x)$out)
    #sample amongst them to keep at maximum of 1000 points and include both min and max outliers values
    boxplotsOutliers=lapply(boxplotsOutliers,function(x)if(length(x)>0)c(sample(c(x),min(length(x),1000)),max(c(x)),min(c(x))))
    dataOutliers=data.frame(yVal=unlist(boxplotsOutliers),xVal=unlist(lapply(seq_along(boxplotsOutliers),function(x)rep(names(boxplotsOutliers)[x],length(boxplotsOutliers[[x]])))))
    #plot boxplots without outliers
    p <- ggplot(data=dataToPlot, aes(y = intensities, x=Experiment ,color=Experiment)) + geom_boxplot(outlier.colour=NA,outlier.shape =NA) +
      ggtitle("Intensities") + theme_bw() + xlab(label="") + 
      theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))
    #add to plot sampled outliers
    p <- p + geom_point(data=dataOutliers,aes(x=xVal,y=yVal,color=xVal),inherit.aes = F)
    if(dataAreFromCel){
      #original plotting function
  	  #boxplot(celData[,firstPlot:lastPlot],which='pm',col=rainbow(nbConditions)[firstPlot:lastPlot],target="probeset",transfo=log2,names=nameConditions[firstPlot:lastPlot],main="Intensities") 
      p <- p + ylab(label="Log2 intensities")
    }else{
      p <- p + ylab(label="Intensities")
    }
    if(opt$format=="pdf"){
      pdf(paste(c("./plotDir/",opt$boxplot,iToPlot,".pdf"),collapse=""))}else{
        png(paste(c("./plotDir/",opt$boxplot,iToPlot,".png"),collapse=""))
      }
    print(p)
    dev.off()
    #save plotly files    
    pp <- ggplotly(p)
    
    #modify plotly object to get HTML file not too heavy for loading
    for(iData in 1:length(pp$x$data)){
      ##get kept outliers y values
      #yPointsToKeep=dataOutliers$yVal[which(dataOutliers$xVal==pp$x$data[[iData]]$name)]
      if(pp$x$data[[iData]]$type=="scatter"){
        ##scatter plot represent outliers points added to boxplot through geom_point
        ##nothing to do as outliers have been sampled allready, just have to modify hover text
        #if(length(yPointsToKeep)>0){
          #pointsToKeep=which(pp$x$data[[iData]]$y %in% yPointsToKeep)
          #pp$x$data[[iData]]$x=pp$x$data[[iData]]$x[pointsToKeep]
          #pp$x$data[[iData]]$y=pp$x$data[[iData]]$y[pointsToKeep]
          #pp$x$data[[iData]]$text=pp$x$data[[iData]]$text[pointsToKeep]
        #}else{
          #pp$x$data[[iData]]$x=NULL
          #pp$x$data[[iData]]$y=NULL
          #pp$x$data[[iData]]$marker$opacity=0
          #pp$x$data[[iData]]$hoverinfo=NULL
          #pp$x$data[[iData]]$text=NULL
        #}
        #modify text to display
        if(dataAreFromCel){
          pp$x$data[[iData]]$text=unlist(lapply(seq_along(pp$x$data[[iData]]$y),function(x)return(paste(c("log2(intensity) ",prettyNum(pp$x$data[[iData]]$y[x],digits=4)),collapse = ""))))
        }else{
          pp$x$data[[iData]]$text=unlist(lapply(seq_along(pp$x$data[[iData]]$y),function(x)return(paste(c("intensity ",prettyNum(pp$x$data[[iData]]$y[x],digits=4)),collapse = ""))))
        }
      }else{
        ##disable marker plotting to keep only box and whiskers plot (outliers are displayed through scatter plot)
        pp$x$data[[iData]]$marker$opacity=0
        
        #sample 50000 points amongst all data to get a lighter html file, sampling size should not be too low to avoid modifying limit of boxplots
        pp$x$data[[iData]]$y=c(sample(dataMatrix[,pp$x$data[[iData]]$name],min(length(dataMatrix[,pp$x$data[[iData]]$name]),50000)),min(dataMatrix[,pp$x$data[[iData]]$name]),max(dataMatrix[,pp$x$data[[iData]]$name]))
        pp$x$data[[iData]]$x=rep(pp$x$data[[iData]]$x[1],length(pp$x$data[[iData]]$y))
        
        ##first remove outliers info
        #downUpValues=boxplot.stats(dataMatrix[,pp$x$data[[iData]]$name])$stats
        #if(verbose)addComment(c("filter values for boxplot",pp$x$data[[iData]]$name,"between",min(downUpValues),"and",max(downUpValues)),T,opt$log)
        #pointsToRemove=which(pp$x$data[[iData]]$y<min(downUpValues))
        #if(length(pointsToRemove)>0)pp$x$data[[iData]]$y=pp$x$data[[iData]]$y[-pointsToRemove]
        #pointsToRemove=which(pp$x$data[[iData]]$y>max(downUpValues))
        #if(length(pointsToRemove)>0)pp$x$data[[iData]]$y=pp$x$data[[iData]]$y[-pointsToRemove]
        #then add sampled outliers info
        #pp$x$data[[iData]]$y=c(yPointsToKeep,pp$x$data[[iData]]$y)
        #pp$x$data[[iData]]$x=rep(pp$x$data[[iData]]$x[1],length(pp$x$data[[iData]]$y))
      }
    }
    
    htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$boxplot,iToPlot,".html"),collapse=""),selfcontained = F)
  }
  remove(p,dataToPlot)
  if(verbose)addComment("Boxplots drawn")
  
}

##----------------------

###plot microarrays (only for .CEL files)###
if (!is.null(opt$microarray) && dataAreFromCel) {
  for (iCondition in 1:nbConditions){
    if(opt$format=="pdf"){
      pdf(paste(c("./plotDir/",opt$microarray,"_",correspondanceNameTable[iCondition,1],".pdf"),collapse=""),onefile = F,width = 5,height = 5)}else{
        png(paste(c("./plotDir/",opt$microarray,"_",correspondanceNameTable[iCondition,1],".png"),collapse=""))
      }
    image(celData[,iCondition],main=correspondanceNameTable[iCondition,2])
    dev.off()
  }
  if(verbose)addComment("Microarray drawn")
}

##----------------------

###plot PCA plot###
if (!is.null(opt$acp)){
  ##to avoid error when nrow is too large, results quite stable with 200k random selected rows
  randomSelection=sample(nrow(dataMatrix),min(200000,nrow(dataMatrix)))
  #remove constant variables
  
  dataFiltered=dataMatrix[randomSelection,]
  toRemove=which(unlist(apply(dataFiltered,1,var))==0)
  if(length(toRemove)>0){
    dataFiltered=dataFiltered[-toRemove,]
  }
  ##geom_text(aes(label=Experiments,hjust=1, vjust=1.3), y = PC2+0.01)
  PACres = prcomp(t(dataFiltered),scale.=TRUE)
  
  if(!is.null(opt$screePlot)){
    #scree plot
    #p <- fviz_eig(PACres)
    dataToPlot=data.frame(compo=seq(1,length(PACres$sdev)),var=(PACres$sdev^2/sum(PACres$sdev^2))*100)
    p<-ggplot(data=dataToPlot, aes(x=compo, y=var)) + geom_bar(stat="identity", fill="steelblue") + geom_line() + geom_point() +
      ggtitle("Scree plot") + theme_bw() + theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5)) +
      xlab(label="Dimensions") + ylab(label="% explained variances") + scale_x_discrete(limits=dataToPlot$compo)
    pp <- ggplotly(p)
    
    if(opt$format=="pdf"){
        pdf(paste(c("./plotDir/",opt$screePlot,".pdf"),collapse=""))}else{
        png(paste(c("./plotDir/",opt$screePlot,".png"),collapse=""))
      }
    plot(p)
    dev.off()
    htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$screePlot,".html"),collapse=""),selfcontained = F)
  }
  
  #now plot pca plots
  
  if(!is.null(opt$factorInfo)){
    fileIdent=""
    symbolset = c("circle","cross","square","diamond","circle-open","square-open","diamond-open","x")
    
    #save equivalence between real factor names and generic ones in correspondanceNameTable
    correspondanceNameTable=rbind(correspondanceNameTable,matrix(c(paste("Factor",1:(ncol(factorInfoMatrix)-1),sep=""),colnames(factorInfoMatrix)[-1]),ncol=2,nrow=ncol(factorInfoMatrix)-1))
    rownames(correspondanceNameTable)=correspondanceNameTable[,2]
    
    #first order factors from decreasing groups number
    orderedFactors=colnames(factorInfoMatrix)[-1][order(unlist(lapply(colnames(factorInfoMatrix)[-1],function(x)length(table(factorInfoMatrix[,x])))),decreasing = T)]
    allFactorsBigger=length(table(factorInfoMatrix[,orderedFactors[length(orderedFactors)]]))>length(symbolset)
    if(allFactorsBigger)addComment("All factors are composed of too many groups to display two factors at same time, each PCA plot will display only one factor groups",T,opt$log) 
    for(iFactor in 1:length(orderedFactors)){
      #if it is the last factor of the list or if all factor 
      if(iFactor==length(orderedFactors) || allFactorsBigger){
       if(length(orderedFactors)==1 || allFactorsBigger){ 
          dataToPlot=data.frame(PC1=PACres$x[,1],PC2=PACres$x[,2],PC3=PACres$x[,3],Experiments=rownames(PACres$x), Attribute1=factorInfoMatrix[rownames(PACres$x),orderedFactors[iFactor]], hoverLabel=unlist(lapply(rownames(PACres$x),function(x)paste(factorInfoMatrix[x,-1],collapse=","))))
          p <- plot_ly(dataToPlot,x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode="markers", color=~Attribute1,colors=rainbow(length(levels(dataToPlot$Attribute1))+2),hoverinfo = 'text', text = ~paste(Experiments,"\n",hoverLabel),marker=list(size=5))%>%
            layout(title = "Principal Component Analysis", scene = list(xaxis = list(title = "Component 1"),yaxis = list(title = "Component 2"),zaxis = list(title = "Component 3")),
                   legend=list(font = list(family = "sans-serif",size = 15,color = "#000")))
          fileIdent=correspondanceNameTable[orderedFactors[iFactor],1]
          #add text label to plot
          ##p <- add_text(p,x = dataToPlot$PC1, y = dataToPlot$PC2 + (max(PACres$x[,2])-min(PACres$x[,2]))*0.02, z = dataToPlot$PC3, mode = 'text', inherit = F, text=rownames(PACres$x), hoverinfo='skip', showlegend = FALSE, color=I('black'))
          #save the plotly plot
          htmlwidgets::saveWidget(as_widget(p), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$acp,"_",fileIdent,".html"),collapse=""),selfcontained = F)
        }
      }else{
        for(iFactorBis in (iFactor+1):length(orderedFactors)){
          if(length(table(factorInfoMatrix[,orderedFactors[iFactorBis]]))<=length(symbolset)){
            dataToPlot=data.frame(PC1=PACres$x[,1],PC2=PACres$x[,2],PC3=PACres$x[,3],Experiments=rownames(PACres$x), Attribute1=factorInfoMatrix[rownames(PACres$x),orderedFactors[iFactor]], Attribute2=factorInfoMatrix[rownames(PACres$x),orderedFactors[iFactorBis]], hoverLabel=unlist(lapply(rownames(PACres$x),function(x)paste(factorInfoMatrix[x,-1],collapse=","))))
            p <- plot_ly(dataToPlot,x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode="markers", color=~Attribute1,colors=rainbow(length(levels(dataToPlot$Attribute1))+2),symbol=~Attribute2,symbols = symbolset,hoverinfo = 'text', text = ~paste(Experiments,"\n",hoverLabel),marker=list(size=5))%>%
              layout(title = "Principal Component Analysis", scene = list(xaxis = list(title = "Component 1"),yaxis = list(title = "Component 2"),zaxis = list(title = "Component 3")),
                     legend=list(font = list(family = "sans-serif",size = 15,color = "#000")))
            fileIdent=paste(correspondanceNameTable[orderedFactors[c(iFactor,iFactorBis)],1],collapse="_AND_")
            #add text label to plot
            ##p <- add_text(p,x = dataToPlot$PC1, y = dataToPlot$PC2 + (max(PACres$x[,2])-min(PACres$x[,2]))*0.02, z = dataToPlot$PC3, mode = 'text', inherit = F, text=rownames(PACres$x), hoverinfo='skip', showlegend = FALSE, color=I('black'))
            #save the plotly plot
            htmlwidgets::saveWidget(as_widget(p), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$acp,"_",fileIdent,".html"),collapse=""),selfcontained = F)
          }else{
            addComment(c("PCA with",orderedFactors[iFactor],"and",orderedFactors[iFactorBis],"groups cannot be displayed, too many groups (max",length(symbolset),")"),T,opt$log) 
          }
        }
      }
    }
  }else{
    dataToPlot=data.frame(PC1=PACres$x[,1],PC2=PACres$x[,2],PC3=PACres$x[,3],Experiments=rownames(PACres$x))
    p <- plot_ly(dataToPlot,x = ~PC1, y = ~PC2, z = ~PC3, type = 'scatter3d', mode="markers",marker=list(size=5,color="salmon"),hoverinfo = 'text',text = ~paste(Experiments))%>%
      layout(title = "Principal Component Analysis", scene = list(xaxis = list(title = "Component 1"),yaxis = list(title = "Component 2"),zaxis = list(title = "Component 3")),
             legend=list(font = list(family = "sans-serif",size = 15,color = "#000")))
    ##p <- add_text(p,x = dataToPlot$PC1, y = dataToPlot$PC2 + (max(PACres$x[,2])-min(PACres$x[,2]))*0.02, z = dataToPlot$PC3, mode = 'text', inherit = F, text=rownames(PACres$x), hoverinfo='skip', showlegend = FALSE, color=I('black'))
    
    #save plotly files 
    htmlwidgets::saveWidget(as_widget(p), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$acp,"_plot.html"),collapse=""),selfcontained = F)
  }
  remove(p,dataToPlot,dataFiltered)
  if(verbose)addComment("ACP plot drawn")
}

#write correspondances between plot file names and displayed names in figure legends, usefull to define html items in xml file
write.table(correspondanceNameTable,file=file.path(getwd(), "correspondanceFileNames.csv"),quote=FALSE,sep="\t",col.names = F,row.names = F)


end.time <- Sys.time()
addComment(c("Total execution time:",as.numeric(end.time - start.time,units="mins"),"mins"),T,opt$log,display=F)

addComment("End of R script",T,opt$log)
#sessionInfo()
