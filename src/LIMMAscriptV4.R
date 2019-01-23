# A command-line interface for LIMMA to use with Galaxy
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

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs()

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "dataFile", "i", 1, "character",
  "factorInfo","a", 1, "character",
  "blockingInfo","b", 1, "character",
  "dicoRenaming","g",1,"character",
  "blockingPolicy","u", 1, "character",
  "thresholdPval","t", 1, "double",
  "thresholdFC","d", 1, "double",
  "format", "f", 1, "character",
  "histo","h", 1, "character",
  "volcano","v", 1, "character",
  "factorsContrast","r", 1, "character",
  "contrastNames","p", 1, "character",
  "firstGroupContrast","m", 1, "character",
  "secondGroupContrast","n", 1, "character",
  "controlGroups","c", 1, "character",
  "fratioFile","s",1,"character",
  "organismID","x",1,"character",
  "rowNameType","y",1,"character",
  "quiet", "q", 0, "logical",
  "log", "l", 1, "character",
  "outputFile" , "o", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# enforce the following required arguments
if (is.null(opt$log)) {
  addComment("[ERROR]'log file' is required\n")
  q( "no", 1, F )
}
addComment("[INFO]Start of R script",T,opt$log,display=FALSE)
if (is.null(opt$dataFile)) {
  addComment("[ERROR]'dataFile' is required",T,opt$log)
  q( "no", 1, F )
}
if (!is.null(opt$blockingInfo) && is.null(opt$blockingPolicy) ) {
  addComment("[ERROR]blocking policy is missing",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$dicoRenaming)) {
  addComment("[ERROR]renaming dictionnary is missing",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$factorsContrast) && is.null(opt$firstGroupContrast) && is.null(opt$secondGroupContrast)) {
  addComment("[ERROR]factor and contrast informations are missing",T,opt$log)
  q( "no", 1, F )
}
if (length(opt$firstGroupContrast)!=length(opt$secondGroupContrast)) {
  addComment("[ERROR]some contrast groups seems to be empty",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$factorInfo)) {
  addComment("[ERROR]factors info is missing",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$format)) {
  addComment("[ERROR]'output format' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$thresholdPval)) {
  addComment("[ERROR]'p-val threshold' is required",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$outputFile)) {
  addComment("[ERROR]'output file' is required",T,opt$log)
  q( "no", 1, F )
}
if (!is.null(opt$volcano) && is.null(opt$thresholdFC)){
  addComment("[ERROR]'FC threshold' is required for volcanos",T,opt$log)
  q( "no", 1, F )
}
if (is.null(opt$fratioFile)) {
  addComment("[ERROR]F-ratio parameter is missing",T,opt$log)
  q( "no", 1, F )
}

#demande si le script sera bavard
verbose <- if (is.null(opt$quiet)) {
  TRUE
}else{
  FALSE
}

#paramètres internes
#pour savoir si on remplace les FC calculés par LIMMA par un calcul du LS-MEAN (ie moyenne de moyennes de chaque groupe dans chaque terme du contraste plutôt qu'une moyenne globale dans chaque terme)
useLSmean=FALSE

addComment("[INFO]Parameters checked!",T,opt$log,display=FALSE)

addComment(c("[INFO]Working directory: ",getwd()),TRUE,opt$log,display=FALSE)
addComment(c("[INFO]Command line: ",args),TRUE,opt$log,display=FALSE)

#directory for plots
dir.create(file.path(getwd(), "plotDir"))
dir.create(file.path(getwd(), "plotLyDir"))

#charge des packages silencieusement
suppressPackageStartupMessages({
  library("methods")
  library("limma")
  library("biomaRt")
  library("ggplot2")
  library("plotly")
  library("stringr")
})


#chargement du fichier dictionnaire de renommage
renamingDico=read.csv(file=file.path(getwd(), opt$dicoRenaming),header=F,sep="\t",colClasses="character")
rownames(renamingDico)=renamingDico[,2]


#chargement des fichiers en entrée
expDataMatrix=read.csv(file=file.path(getwd(), opt$dataFile),header=F,sep="\t",colClasses="character")
#remove first row to convert it as colnames (to avoid X before colnames with header=T)
colNamesData=expDataMatrix[1,-1]
expDataMatrix=expDataMatrix[-1,]
#remove first colum to convert it as rownames
rowNamesData=expDataMatrix[,1]
expDataMatrix=expDataMatrix[,-1]
if(is.data.frame(expDataMatrix)){
  expDataMatrix=data.matrix(expDataMatrix)
}else{
  expDataMatrix=data.matrix(as.numeric(expDataMatrix))
}
dimnames(expDataMatrix)=list(rowNamesData,colNamesData)

#test the number of rows that are constant in dataMatrix
nbConstantRows=length(which(unlist(apply(expDataMatrix,1,var))==0))
if(nbConstantRows>0){
  addComment(c("[WARNING]",nbConstantRows,"rows are constant across conditions in input data file"),T,opt$log,display=FALSE)
}

#test if all condition names are present in dico
if(!all(colnames(expDataMatrix) %in% rownames(renamingDico))){
  addComment("[ERROR]Missing condition names in renaming dictionary",T,opt$log)
  q( "no", 1, F )
}

addComment("[INFO]Expression data loaded and checked",T,opt$log,display=FALSE)
addComment(c("[INFO]Dim of expression matrix:",dim(expDataMatrix)),T,opt$log,display=FALSE)

#chargement du fichier des facteurs
factorInfoMatrix=read.csv(file=file.path(getwd(), opt$factorInfo),header=F,sep="\t",colClasses="character")
#remove first row to convert it as colnames
colnames(factorInfoMatrix)=factorInfoMatrix[1,]
factorInfoMatrix=factorInfoMatrix[-1,]
#use first colum to convert it as rownames but not removing it to avoid conversion as vector in unique factor case
rownames(factorInfoMatrix)=factorInfoMatrix[,1]

if(length(setdiff(colnames(expDataMatrix),rownames(factorInfoMatrix)))!=0){
  addComment("[ERROR]Missing samples in factor file",T,opt$log)
  q( "no", 1, F )
}

#order sample as in expression matrix and remove spurious sample
factorInfoMatrix=factorInfoMatrix[colnames(expDataMatrix),]

#test if all values names are present in dico
if(!all(unlist(factorInfoMatrix) %in% rownames(renamingDico))){
  addComment("[ERROR]Missing factor names in renaming dictionary",T,opt$log)
  q( "no", 1, F )
}

addComment("[INFO]Factors OK",T,opt$log,display=FALSE)
addComment(c("[INFO]Dim of factorInfo matrix:",dim(factorInfoMatrix)),T,opt$log,display=FALSE)
  
##manage blocking factor
blockingFactor=NULL
blockinFactorsList=NULL
if(!is.null(opt$blockingInfo)){
  
  #chargement du fichier des blocking factors
  blockingInfoMatrix=read.csv(file=file.path(getwd(), opt$blockingInfo),header=F,sep="\t",colClasses="character")
  #remove first row to convert it as colnames
  colnames(blockingInfoMatrix)=blockingInfoMatrix[1,]
  blockingInfoMatrix=blockingInfoMatrix[-1,]
  #use first colum to convert it as rownames but not removing it to avoid conversion as vector in unique factor case
  rownames(blockingInfoMatrix)=blockingInfoMatrix[,1]
  
  
  if(length(setdiff(colnames(expDataMatrix),rownames(blockingInfoMatrix)))!=0){
    addComment("[ERROR]Missing samples in blocking factor file",T,opt$log)
    q( "no", 1, F )
  }
  
  #order sample as in expression matrix
  blockingInfoMatrix=blockingInfoMatrix[colnames(expDataMatrix),]
  
  #test if all blocking names are present in dico
  if(!all(unlist(blockingInfoMatrix) %in% rownames(renamingDico))){
    addComment("[ERROR]Missing blocking names in renaming dictionary",T,opt$log)
    q( "no", 1, F )
  }
  
  #remove blocking factors allready present as real factors
  blockingNotInMainFactors=setdiff(colnames(blockingInfoMatrix)[-1],colnames(factorInfoMatrix)[-1])
  
  if(length(blockingNotInMainFactors)<(ncol(blockingInfoMatrix)-1))addComment("[WARNING]Blocking factors cannot be principal factors",T,opt$log,display=FALSE)
  
  if(length(blockingNotInMainFactors)>0){
    
    blockingInfoMatrix=blockingInfoMatrix[,c(colnames(blockingInfoMatrix)[1],blockingNotInMainFactors)]
    
    groupBlocking=rep("c",ncol(expDataMatrix))
    #for each blocking factor
    for(blockingFact in blockingNotInMainFactors){
      if(opt$blockingPolicy=="correlated"){
        indNewFact=as.numeric(factor(blockingInfoMatrix[,blockingFact]))
        groupBlocking=paste(groupBlocking,indNewFact,sep="_")
      }else{
        if(is.null(blockinFactorsList))blockinFactorsList=list()
        blockinFactorsList[[blockingFact]]=factor(unlist(lapply(blockingInfoMatrix[,blockingFact],function(x)paste(c(blockingFact,"_",x),collapse=""))))
      }
    }
    if(opt$blockingPolicy=="correlated"){
      blockingFactor=factor(groupBlocking)
      if(length(levels(blockingFactor))==1){
        addComment("[ERROR]Selected blocking factors seems to be constant",T,opt$log)
        q( "no", 1, F )
      }
    }
    addComment("[INFO]Blocking info OK",T,opt$log,display=FALSE)
  }else{
    addComment("[WARNING]No blocking factors will be considered",T,opt$log,display=FALSE)
  }
}


##rename different input parameters using renamingDictionary
opt$factorsContrast=renamingDico[unlist(lapply(unlist(strsplit(opt$factorsContrast,",")),function(x)which(renamingDico[,1]==x))),2]

for(iContrast in 1:length(opt$firstGroupContrast)){
  opt$firstGroupContrast[iContrast]=paste(unlist(lapply(unlist(strsplit(opt$firstGroupContrast[iContrast],",")),function(x)paste(renamingDico[unlist(lapply(unlist(strsplit(x,"\\*")),function(x)which(renamingDico[,1]==x))),2],collapse="*"))),collapse=",")
  opt$secondGroupContrast[iContrast]=paste(unlist(lapply(unlist(strsplit(opt$secondGroupContrast[iContrast],",")),function(x)paste(renamingDico[unlist(lapply(unlist(strsplit(x,"\\*")),function(x)which(renamingDico[,1]==x))),2],collapse="*"))),collapse=",")
  }

if(!is.null(opt$controlGroups)){
  renamedGroups=c()
  for(iGroup in unlist(strsplit(opt$controlGroups,","))){
    renamedGroups=c(renamedGroups,paste(renamingDico[unlist(lapply(unlist(strsplit(iGroup,":")),function(x)which(renamingDico[,1]==x))),2],collapse=":"))
  }
  opt$controlGroups=renamedGroups
}
addComment("[INFO]Contrast variables are renamed to avoid confusion",T,opt$log,display=FALSE)
##renaming done

#to convert factor as numeric value --> useless now ?
#expDataMatrix=apply(expDataMatrix,c(1,2),function(x)as.numeric(paste(x)))

#get factors info for LIMMA
factorsList=list()
for(iFactor in opt$factorsContrast){
  if(!(iFactor %in% colnames(factorInfoMatrix))){
    addComment("[ERROR]Required factors are missing in input file",T,opt$log)
    q( "no", 1, F )
  }
  factorsList[[iFactor]]=factor(unlist(lapply(factorInfoMatrix[,iFactor],function(x)paste(c(iFactor,"_",x),collapse=""))))
  if(length(levels(factorsList[[iFactor]]))==1){
    addComment("[ERROR]One selected factor seems to be constant",T,opt$log)
    q( "no", 1, F )
  }
}

#check if there is at least 2 factors to allow interaction computation
if(!is.null(opt$controlGroups) && length(factorsList)<2){
  addComment("[ERROR]You cannot ask for interaction with less than 2 factors",T,opt$log)
  q( "no", 1, F )
}

#merge all factors as a single one
factorsMerged=as.character(factorsList[[opt$factorsContrast[1]]])
for(iFactor in opt$factorsContrast[-1]){
  factorsMerged=paste(factorsMerged,as.character(factorsList[[iFactor]]),sep=".")
}
factorsMerged=factor(factorsMerged)

#checked that coefficient number (ie. factorsMerged levels) is strictly smaller than sample size
if(length(levels(factorsMerged))>=length(factorsMerged)){
  addComment(c("[ERROR]No enough samples (",length(factorsMerged),") to estimate ",length(levels(factorsMerged))," coefficients"),T,opt$log)
  q( "no", 1, F )
}

#get the sample size of each factor values
sampleSizeFactor=table(factorsMerged)


if(!is.null(blockinFactorsList)){
  factorString=c("blockinFactorsList[['", names(blockinFactorsList)[1],"']]")
  for(blockingFact in names(blockinFactorsList)[-1]){
    factorString=c(factorString," + blockinFactorsList[['",blockingFact,"']]")
  }
  design = model.matrix(as.formula(paste(c("~ factorsMerged +",factorString," + 0"),collapse="")))
  
  #rename design columns
  coeffMeaning = levels(factorsMerged)
  for(blockingFact in blockinFactorsList){
    coeffMeaning=c(coeffMeaning,levels(blockingFact)[-1])
  }
  colnames(design) = coeffMeaning
}else{
  design = model.matrix(as.formula( ~ factorsMerged + 0))
  
  #rename degin columns
  coeffMeaning = levels(factorsMerged)
  colnames(design) = coeffMeaning
}
  
addComment(c("[INFO]Available coefficients: ",coeffMeaning),T,opt$log,display=F)

estimableCoeff=which(colSums(design)!=0)
  
addComment("[INFO]Design done",T,opt$log,display=F)
  
  #use blocking factor if exists
if(!is.null(blockingFactor)){
  corfit <- duplicateCorrelation(expDataMatrix, design, block=blockingFactor)
    
  #run linear model  fit
  data.fit = lmFit(expDataMatrix,design,block = blockingFactor, correlation=corfit$consensus.correlation)
}else{
  #run linear model  fit  
  data.fit = lmFit(expDataMatrix,design)
}
  
estimatedCoeff=which(!is.na(data.fit$coefficients[1,]))
  
addComment("[INFO]Lmfit done",T,opt$log,display=F)

#catch situation where some coefficients cannot be estimated, probably due to dependances between design columns 
#if(length(setdiff(estimableCoeff,estimatedCoeff))>0){
#  addComment("[ERROR]Error in design matrix, check your group definitions",T,opt$log)
#  q( "no", 1, F )
#}
#to strong condition, should return ERROR only when coefficients relative to principal factors cannot be estimated, otherwise, return a simple WARNING
  
#define requested contrasts 
requiredContrasts=c()
humanReadingContrasts=c()
persoContrastName=c()
  for(iContrast in 1:length(opt$firstGroupContrast)){
    posGroup=unlist(lapply(unlist(strsplit(opt$firstGroupContrast[iContrast],",")),function(x)paste(paste(opt$factorsContrast,unlist(strsplit(x,"\\*")),sep="_"),collapse=".")))
    negGroup=unlist(lapply(unlist(strsplit(opt$secondGroupContrast[iContrast],",")),function(x)paste(paste(opt$factorsContrast,unlist(strsplit(x,"\\*")),sep="_"),collapse=".")))
    #clear posGroup and negGroup from empty groups
    emptyPosGroups=which(!(posGroup%in%coeffMeaning))
    if(length(emptyPosGroups)>0){
      addComment(c("[WARNING]The group(s)",posGroup[emptyPosGroups],"is/are removed from contrast as it/they is/are empty"),T,opt$log,display=FALSE)
      posGroup=posGroup[-emptyPosGroups]
      currentHumanContrast=paste(unlist(strsplit(opt$firstGroupContrast[iContrast],","))[-emptyPosGroups],collapse="+") 
    }else{
      currentHumanContrast=paste(unlist(strsplit(opt$firstGroupContrast[iContrast],",")),collapse="+")  
    }
    emptyNegGroups=which(!(negGroup%in%coeffMeaning))
    if(length(emptyNegGroups)>0){
      addComment(c("[WARNING]The group(s)",negGroup[emptyNegGroups],"is/are removed from contrast as it/they is/are empty"),T,opt$log,display=FALSE)
      negGroup=negGroup[-emptyNegGroups]
      currentHumanContrast=paste(c(currentHumanContrast,unlist(strsplit(opt$secondGroupContrast[iContrast],","))[-emptyNegGroups]),collapse="-")
    }else{
      currentHumanContrast=paste(c(currentHumanContrast,unlist(strsplit(opt$secondGroupContrast[iContrast],","))),collapse="-")
    }
    if(length(posGroup)==0 || length(negGroup)==0 ){
      addComment(c("[WARNING]Contrast",posGroup,"-",negGroup,"cannot be estimated due to empty group"),T,opt$log,display=FALSE)
    }else{
      if(all(posGroup%in%negGroup) && all(negGroup%in%posGroup)){
        addComment(c("[WARNING]Contrast",posGroup,"-",negGroup,"cannot be estimated due to null contrast"),T,opt$log,display=FALSE)
      }else{
        #get coefficients required for first group added as positive
        positiveCoeffWeights=sampleSizeFactor[posGroup]/sum(sampleSizeFactor[posGroup])
        #positiveCoeffWeights=rep(1,length(posGroup))
        #names(positiveCoeffWeights)=names(sampleSizeFactor[posGroup])
        #get coefficients required for second group added as negative
        negativeCoeffWeights=sampleSizeFactor[negGroup]/sum(sampleSizeFactor[negGroup])
        #negativeCoeffWeights=rep(1,length(negGroup))
        #names(negativeCoeffWeights)=names(sampleSizeFactor[negGroup])
        #build the resulting contrast
        currentContrast=paste(paste(positiveCoeffWeights[posGroup],posGroup,sep="*"),collapse="+")
        currentContrast=paste(c(currentContrast,paste(paste(negativeCoeffWeights[negGroup],negGroup,sep="*"),collapse="-")),collapse="-")
        requiredContrasts=c(requiredContrasts,currentContrast)
        
        #build the human reading contrast
        humanReadingContrasts=c(humanReadingContrasts,currentHumanContrast)
        if(!is.null(opt$contrastNames) && nchar(opt$contrastNames[iContrast])>0){
          persoContrastName=c(persoContrastName,opt$contrastNames[iContrast])
        }else{
          persoContrastName=c(persoContrastName,"")
        }
        
        addComment(c("[INFO]Contrast added : ",currentHumanContrast),T,opt$log,display=F)
        addComment(c("with complete formula ",currentContrast),T,opt$log,display=F)
      }
    }
  }
  
  
  #define the true formula with interactions to get interaction coefficients
  factorString=c("factorsList[['", names(factorsList)[1],"']]")
  for(iFactor in names(factorsList)[-1]){
    factorString=c(factorString," * factorsList[['",iFactor,"']]")
  }
  
  if(!is.null(blockinFactorsList)){
    for(blockingFact in names(blockinFactorsList)){
      factorString=c(factorString," + blockinFactorsList[['",blockingFact,"']]")
    }
  }
  
  #should not be null at the end 
  allFtestMeanSquare=NULL
  #to get the F-test values
  estimatedInteractions=rownames(anova(lm(as.formula(paste(c("expDataMatrix[1,] ~ ",factorString),collapse="")))))
  estimatedInteractions=c(unlist(lapply(estimatedInteractions[-length(estimatedInteractions)],function(x){temp=unlist(strsplit(x,"[ \" | : ]"));paste(temp[seq(2,length(temp),3)],collapse=":")})),estimatedInteractions[length(estimatedInteractions)])
  #rename estimated interaction terms using renamingDico
  estimatedInteractions=c(unlist(lapply(estimatedInteractions[-length(estimatedInteractions)],function(x)paste(renamingDico[unlist(strsplit(x,":")),1],collapse=":"))),estimatedInteractions[length(estimatedInteractions)])
  t <- unlist(apply(expDataMatrix,1,function(x){temp=anova(lm(as.formula(paste(c("x ~ ",factorString),collapse=""))))$`Mean Sq`;temp/temp[length(temp)]}))
  allFtestMeanSquare <- t(matrix(t,nrow=length(estimatedInteractions)))
  #remove from allFtest rows containing NA
  if(length(which(is.na(allFtestMeanSquare[,1])))>0)allFtestMeanSquare=allFtestMeanSquare[-(which(is.na(allFtestMeanSquare[,1]))),]
  colnames(allFtestMeanSquare)=estimatedInteractions
  
  #add contrasts corresponding to interaction terms
  if(!is.null(opt$controlGroups)){
    #first load user defined control group for each factor
    controlGroup=rep(NA,length(factorsList))
    names(controlGroup)=names(factorsList)
    for(iGroup in opt$controlGroups){
      splitGroup=unlist(strsplit(iGroup,":"))
      splitGroup[2]=paste(splitGroup[1],splitGroup[2],sep = "_")
      #check if defined control group is really a level of the corresponding factor
      if(!splitGroup[1]%in%names(controlGroup) || !splitGroup[2]%in%factorsList[[splitGroup[1]]]){
        addComment(c("[ERROR]The factor name",splitGroup[1],"does not exist or group name",splitGroup[2]),T,opt$log)
        q( "no", 1, F )
      }
      if(!is.na(controlGroup[splitGroup[1]])){
        addComment("[ERROR]Several control groups are defined for the same factor",T,opt$log)
        q( "no", 1, F )
      }
      controlGroup[splitGroup[1]]=splitGroup[2]
    }
    
    #check if all factor have a defined control group
    if(any(is.na(controlGroup))){
      addComment("[ERROR]Missing control group for some factors",T,opt$log)
      q( "no", 1, F )
    }
    
    nbFactors=length(factorsList)
    interactionContrasts=c()
    contrastClass=c()
    #initialize list for the first level
    newPreviousLoopContrast=list()
    for(iFactorA in 1:(nbFactors-1)){
      nameFactorA=names(factorsList)[iFactorA]
      compA=c()
      for(levelA in setdiff(levels(factorsList[[iFactorA]]),controlGroup[nameFactorA])){
        compA=c(compA,paste(levelA,controlGroup[nameFactorA],sep="-"))
      }
      newPreviousLoopContrast[[as.character(iFactorA)]]=compA
    }
    #make a loop for growing interaction set
    for(globalIfactor in 1:(nbFactors-1)){
      previousLoopContrast=newPreviousLoopContrast
      newPreviousLoopContrast=list()
      #factor A reuse contrasts made at previsous loop
      for(iFactorA in names(previousLoopContrast)){
        compA=previousLoopContrast[[iFactorA]]
  
        if(max(as.integer(unlist(strsplit(iFactorA,"\\."))))<nbFactors){
          #factor B is the new factor to include in intreraction set
          for(iFactorB in (max(as.integer(unlist(strsplit(iFactorA,"\\."))))+1):nbFactors){
            nameFactorB=names(factorsList)[iFactorB]
            #keep contrasts identified trough interacting factors set
            newPreviousLoopContrast[[paste(iFactorA,iFactorB,sep=".")]]=c()
              for(iCompA in compA){
                for(levelB in setdiff(levels(factorsList[[iFactorB]]),controlGroup[nameFactorB])){
                  #decompose the contrast compA to apply the new level of factor B on each term
                  temp=unlist(strsplit(iCompA,"[ + ]"))
                  splitCompA=temp[1]
                  for(iTemp in temp[-1])splitCompA=c(splitCompA,"+",iTemp)
                  splitCompA=unlist(lapply(splitCompA,function(x){temp=unlist(strsplit(x,"-"));splitCompB=temp[1];for(iTemp in temp[-1])splitCompB=c(splitCompB,"-",iTemp);splitCompB}))
                  #apply on each contrast term the new level of factor B
                  firstTerm=paste(unlist(lapply(splitCompA,function(x)if(x!="+" && x!="-"){paste(x,levelB,sep=".")}else{x})),collapse="")
                  secondTerm=negativeExpression(paste(unlist(lapply(splitCompA,function(x)if(x!="+" && x!="-"){paste(x,controlGroup[nameFactorB],sep=".")}else{x})),collapse=""))
                  currentContrast=paste(c(firstTerm,secondTerm),collapse="")
                  
                  newPreviousLoopContrast[[paste(iFactorA,iFactorB,sep=".")]]=c(newPreviousLoopContrast[[paste(iFactorA,iFactorB,sep=".")]],currentContrast)
                }
              }
            }
        }
      }
      for(iContrast in names(newPreviousLoopContrast)){
        contrastClass=c(contrastClass,rep(iContrast,length(newPreviousLoopContrast[[iContrast]])))
      }
      interactionContrasts=c(interactionContrasts,unlist(newPreviousLoopContrast))
    }
    #make human title for interactions
    names(interactionContrasts)=contrastClass
    humanReadingInteraction=unlist(lapply(interactionContrasts,function(x)gsub("\\.",":",unlist(strsplit(x,"[+-]"))[1])))
    
    contrastToIgnore=c()
    
    #complete with control groups and order to match to coeffs
    for(iContrast in 1:length(interactionContrasts)){
      missingFactor=setdiff(1:nbFactors,as.integer(unlist(strsplit(names(interactionContrasts[iContrast]),"\\."))))
      #decompose the contrast
      temp=unlist(strsplit(interactionContrasts[iContrast],"[ + ]"))
      splitContrast=temp[1]
      for(iTemp in temp[-1])splitContrast=c(splitContrast,"+",iTemp)
      splitContrast=unlist(lapply(splitContrast,function(x){temp=unlist(strsplit(x,"-"));splitCompB=temp[1];for(iTemp in temp[-1])splitCompB=c(splitCompB,"-",iTemp);splitCompB}))
      for(iFactor in missingFactor){
        for(iTerm in 1:length(splitContrast)){
          if(splitContrast[iTerm]!="+" && splitContrast[iTerm]!="-"){
            splitTerm=unlist(strsplit(splitContrast[iTerm],"\\."))
            if(iFactor==1)splitContrast[iTerm]=paste(c(controlGroup[names(factorsList)[iFactor]],splitTerm),collapse=".")
            if(iFactor==nbFactors)splitContrast[iTerm]=paste(c(splitTerm,controlGroup[names(factorsList)[iFactor]]),collapse=".")
            if(iFactor>1 && iFactor<nbFactors)splitContrast[iTerm]=paste(c(splitTerm[1:(iFactor-1)],controlGroup[names(factorsList)[iFactor]],splitTerm[iFactor:length(splitTerm)]),collapse=".")
          }
        }
      }
      interactionContrasts[iContrast]=paste(splitContrast,collapse="")
      if(all(splitContrast[seq(1,length(splitContrast),2)]%in%coeffMeaning)){
        addComment(c("[INFO]Interaction contrast added : ",humanReadingInteraction[iContrast]),T,opt$log,display=F)
        addComment(c("with complete formula ",interactionContrasts[iContrast]),T,opt$log,display=F)
      }else{
        contrastToIgnore=c(contrastToIgnore,iContrast)
        addComment(c("[WARNING]Interaction contrast",humanReadingInteraction[iContrast],"is removed due to empty group"),T,opt$log,display=F)
      }
    }
    
    #add interaction contrasts to global contrast list
    if(length(contrastToIgnore)>0){
      requiredContrasts=c(requiredContrasts,interactionContrasts[-contrastToIgnore])
      humanReadingContrasts=c(humanReadingContrasts,humanReadingInteraction[-contrastToIgnore])
      persoContrastName=c(persoContrastName,rep("",length(humanReadingInteraction[-contrastToIgnore])))
    }else{
      requiredContrasts=c(requiredContrasts,interactionContrasts)
      humanReadingContrasts=c(humanReadingContrasts,humanReadingInteraction)
      persoContrastName=c(persoContrastName,rep("",length(humanReadingInteraction)))
    }
  }
  
  #remove from requiredContrasts contrasts that cannot be estimated
  toRemove=unique(unlist(lapply(setdiff(coeffMeaning,names(estimatedCoeff)),function(x)grep(x,requiredContrasts))))
  if(length(toRemove)>0){
    addComment(c("[WARNING]",length(toRemove)," contrasts are removed, due to missing coefficients"),T,opt$log,display=FALSE)
    requiredContrasts=requiredContrasts[-toRemove]
    humanReadingContrasts=humanReadingContrasts[-toRemove]
    persoContrastName=persoContrastName[-toRemove]
  }
  
  #compute for each contrast mean of coefficients in posGroup and negGroup for FC computation of log(FC) with LSmean as in Partek
  meanPosGroup=list()
  meanNegGroup=list()
  for(iContrast in 1:length(requiredContrasts)){
    #define posGroup and negGroup
    #first split contrast 
    temp=unlist(strsplit(requiredContrasts[iContrast],"[ + ]"))
    splitContrast=temp[1]
    for(iTemp in temp[-1])splitContrast=c(splitContrast,"+",iTemp)
    splitContrast=unlist(lapply(splitContrast,function(x){temp=unlist(strsplit(x,"-"));splitCompB=temp[1];for(iTemp in temp[-1])splitCompB=c(splitCompB,"-",iTemp);splitCompB}))
    #and then put each term in good group
    posGroup=c()
    negGroup=c()
    nextIsPos=TRUE
    for(iSplit in splitContrast){
      if(iSplit=="+")nextIsPos=TRUE
      if(iSplit=="-")nextIsPos=FALSE
      if(iSplit!="-" && iSplit!="+"){
        #remove weights of contrast terms
        iSplitBis=unlist(strsplit(iSplit,"[*]"))
        iSplitBis=iSplitBis[length(iSplitBis)]
        if(nextIsPos)posGroup=c(posGroup,iSplitBis)
        else negGroup=c(negGroup,iSplitBis)
      }
    }
    #compute means for each group
    meanPosGroup[[iContrast]]=apply(as.matrix(data.fit$coefficients[,posGroup],ncol=length(posGroup)),1,mean)
    meanNegGroup[[iContrast]]=apply(as.matrix(data.fit$coefficients[,negGroup],ncol=length(negGroup)),1,mean)
  }

  
  contrast.matrix = makeContrasts(contrasts=requiredContrasts,levels=design)
  data.fit.con = contrasts.fit(data.fit,contrast.matrix)
  
  addComment("[INFO]Contrast definition done",T,opt$log,T,display=FALSE)
  
  #compute LIMMA statistics
  data.fit.eb = eBayes(data.fit.con)
  
  addComment("[INFO]Estimation done",T,opt$log,T,display=FALSE)
  
  #adjust p.value through FDR
  data.fit.eb$adj_p.value=data.fit.eb$p.value
  for(iComparison in 1:ncol(data.fit.eb$adj_p.value)){
    data.fit.eb$adj_p.value[,iComparison]=p.adjust(data.fit.eb$p.value[,iComparison],"fdr")
  }

  #add a new field based on LS-means for each contrast instead of global mean like they were calculated in coefficients field
  data.fit.eb$coefficientsLS=data.fit.eb$coefficients
  if(ncol(data.fit.eb$coefficients)!=length(meanPosGroup)){
    addComment("[ERROR]Estimated contrasts number unexpected",T,opt$log)
    q( "no", 1, F )
  }
  for(iContrast in 1:length(meanPosGroup)){
    data.fit.eb$coefficientsLS[,iContrast]=meanPosGroup[[iContrast]][rownames(data.fit.eb$coefficientsLS)]-meanNegGroup[[iContrast]][rownames(data.fit.eb$coefficientsLS)]
  }
  
  #if requested replace coefficient computed as global mean by LS-means values
  if(useLSmean)data.fit.eb$coefficients=data.fit.eb$coefficientsLS

addComment("[INFO]Core treatment done",T,opt$log,T,display=FALSE)
  
  
##convert humanReadingContrasts with namingDictionary to create humanReadingContrastsRenamed and keep original humanReadingContrasts names for file names 
humanReadingContrastsRenamed=rep("",length(humanReadingContrasts))
for(iContrast in 1:length(humanReadingContrasts)){
  if(persoContrastName[iContrast]==""){
    #if(verbose)addComment(humanReadingContrasts[iContrast])
    specialCharacters=str_extract_all(humanReadingContrasts[iContrast],"[+|*|_|:|-]")[[1]]
    #if(verbose)addComment(specialCharacters)
    nameConverted=unlist(lapply(strsplit(humanReadingContrasts[iContrast],"[+|*|_|:|-]")[[1]],function(x)renamingDico[x,1]))
    #if(verbose)addComment(nameConverted)
    humanReadingContrastsRenamed[iContrast]=paste(nameConverted,specialCharacters,collapse="",sep="")
    #if(verbose)addComment(humanReadingContrastsRenamed[iContrast])
    humanReadingContrastsRenamed[iContrast]=substr(humanReadingContrastsRenamed[iContrast],1,nchar(humanReadingContrastsRenamed[iContrast])-1)
  }else{
    humanReadingContrastsRenamed[iContrast]=persoContrastName[iContrast]
  }
}

#write correspondances between plot file names (humanReadingContrasts) and displayed names in figure legends (humanReadingContrastsRenamed), usefull to define html items in xml file
correspondanceTable=matrix("",ncol=2,nrow=ncol(data.fit.eb$p.value))
correspondanceTable[,1]=unlist(lapply(humanReadingContrasts,function(x)gsub(":","_INT_",gsub("\\+","_PLUS_",gsub("\\*","_AND_",x)))))
correspondanceTable[,2]=humanReadingContrastsRenamed
rownames(correspondanceTable)=correspondanceTable[,2]
write.table(correspondanceTable,file=file.path(getwd(), "correspondanceFileNames.csv"),quote=FALSE,sep="\t",col.names = F,row.names = F)

#plot nominal p-val histograms for selected comparisons
histogramPerPage=6
if (!is.null(opt$histo)) {
  iToPlot=1
  plotVector=list()
  nbComparisons=ncol(data.fit.eb$p.value)
  for (iComparison in 1:nbComparisons){
    dataToPlot=data.frame(pval=data.fit.eb$p.value[,iComparison],id=rownames(data.fit.eb$p.value))
    p <- ggplot(data=dataToPlot, aes(x=pval)) + geom_histogram(colour="red", fill="salmon") +
      theme_bw() + ggtitle(humanReadingContrastsRenamed[iComparison]) + ylab(label="Frequencies") + xlab(label="Nominal p-val") +
      theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5))
    plotVector[[length(plotVector)+1]]=p
    
    pp <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$histo,"_",correspondanceTable[humanReadingContrastsRenamed[iComparison],1],".html"),collapse=""),selfcontained = F)
    
    if(iComparison==nbComparisons || length(plotVector)==histogramPerPage){
      #plot and close the actual plot
      if(opt$format=="pdf"){
        pdf(paste(c("./plotDir/",opt$histo,iToPlot,".pdf"),collapse=""))}else{
          png(paste(c("./plotDir/",opt$histo,iToPlot,".png"),collapse=""))
        }
      multiplot(plotlist=plotVector,cols=2)
      dev.off()
      if(iComparison<nbComparisons){
        #prepare for a new plotting file if necessary
        plotVector=list()
        iToPlot=iToPlot+1
      }
    }
  }
  addComment("[INFO]Histograms drawn",T,opt$log,T,display=FALSE)
  
}

#plot F-test sum square barplot
if(!is.null(allFtestMeanSquare)){
  dataToPlot=data.frame(Fratio=apply(allFtestMeanSquare,2,mean),Factors=factor(colnames(allFtestMeanSquare),levels = colnames(allFtestMeanSquare)))

  p <- ggplot(data=dataToPlot, aes(x=Factors, y=Fratio, fill=Factors)) +
    geom_bar(stat="identity") + scale_fill_brewer(palette="Set1") + ylab(label="mean F-ratio") +
    theme_bw() + theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5)) + ggtitle("Source of variation")
  
   if(opt$format=="pdf"){
    pdf(paste(c("./plotDir/",opt$fratioFile,".pdf"),collapse=""))}else{
      png(paste(c("./plotDir/",opt$fratioFile,".png"),collapse=""))
    }
  plot(p)
  dev.off()
  
  pp <- ggplotly(p)
  htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$fratioFile,".html"),collapse=""),selfcontained = F)
  
  addComment("[INFO]SumSquareTest drawn",T,opt$log,T,display=FALSE)
}

#plot VOLCANO plot
#volcanoplot(data.fit.eb,coef=1,highlight=10)
volcanoPerPage=1
FCthreshold=log2(opt$thresholdFC)
if (!is.null(opt$volcano)) {
  iToPlot=1
  plotVector=list()
  nbComparisons=ncol(data.fit.eb$adj_p.value)
  for (iComparison in 1:nbComparisons){
    
    #define the log10(p-val) threshold corresponding to FDR threshold fixed by user
    probeWithLowFDR=-log10(data.fit.eb$p.value[which(data.fit.eb$adj_p.value[,iComparison]<=opt$thresholdPval),iComparison])
    pvalThresholdFDR=NULL
    if(length(probeWithLowFDR)>0)pvalThresholdFDR=min(probeWithLowFDR)
    
    #get significative points over FC thresholds
    significativePoints=which(abs(data.fit.eb$coefficients[,iComparison])>=FCthreshold)
    
    #to reduce size of html plot, we keep 20000 points maximum sampled amongst genes with pval>=33%(pval) and abs(log2(FC))<=66%(abs(log2(FC)))
    htmlPointsToRemove=intersect(which(abs(data.fit.eb$coefficients[,iComparison])<=quantile(abs(data.fit.eb$coefficients[,iComparison]),c(0.66))),which(data.fit.eb$p.value[,iComparison]>=quantile(abs(data.fit.eb$p.value[,iComparison]),c(0.33))))
    if(length(htmlPointsToRemove)>20000){
      htmlPointsToRemove=setdiff(htmlPointsToRemove,sample(htmlPointsToRemove,20000))
    }else{
      htmlPointsToRemove=c() 
    }
      
      xMinLimPlot=min(data.fit.eb$coefficients[,iComparison])-0.2
      xMaxLimPlot=max(data.fit.eb$coefficients[,iComparison])+0.2
      yMaxLimPlot= max(-log10(data.fit.eb$p.value[,iComparison]))+0.2
      
      if(length(significativePoints)>0){
        dataSignifToPlot=data.frame(pval=-log10(data.fit.eb$p.value[significativePoints,iComparison]),FC=data.fit.eb$coefficients[significativePoints,iComparison],description=paste(names(data.fit.eb$coefficients[significativePoints,iComparison]),"\n","FC: " , round(2^data.fit.eb$coefficients[significativePoints,iComparison],2) , " | FDR p-val: ",prettyNum(data.fit.eb$adj_p.value[significativePoints,iComparison],digits=4), sep=""))
        #to test if remains any normal points to draw
        if(length(significativePoints)<nrow(data.fit.eb$p.value)){
          dataToPlot=data.frame(pval=-log10(data.fit.eb$p.value[-significativePoints,iComparison]),FC=data.fit.eb$coefficients[-significativePoints,iComparison],description=paste("FC: " , round(2^data.fit.eb$coefficients[-significativePoints,iComparison],2) , " | FDR p-val: ",prettyNum(data.fit.eb$adj_p.value[-significativePoints,iComparison],digits=4), sep=""))
        }else{
          dataToPlot=data.frame(pval=0,FC=0,description="null")
        }
      }else{
        dataToPlot=data.frame(pval=-log10(data.fit.eb$p.value[,iComparison]),FC=data.fit.eb$coefficients[,iComparison],description=paste("FC: " , round(2^data.fit.eb$coefficients[,iComparison],2) , " | FDR p-val: ",prettyNum(data.fit.eb$adj_p.value[,iComparison],digits=4), sep=""))
      }
        
      p <- ggplot(data=dataToPlot, aes(x=FC, y=pval)) + geom_point() + geom_vline(xintercept=-FCthreshold, color="salmon",linetype="dotted", size=1) +
        geom_vline(xintercept=FCthreshold, color="salmon",linetype="dotted", size=1) +
        geom_text(data.frame(text=c(paste(c("log2(1/FC=",opt$thresholdFC,")"),collapse=""),paste(c("log2(FC=",opt$thresholdFC,")"),collapse="")),x=c(-FCthreshold,FCthreshold),y=c(0,0)),mapping=aes(x=x, y=y, label=text), size=4, angle=90, vjust=-0.4, hjust=0, color="salmon") +
        theme_bw() + ggtitle(humanReadingContrastsRenamed[iComparison]) + ylab(label="-log10(p-val)") + xlab(label="Log2 Fold Change") +
        theme(panel.border=element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none") 
      if(!is.null(pvalThresholdFDR)) p <- p + geom_hline(yintercept=pvalThresholdFDR, color="skyblue1",linetype="dotted", size=0.5) + geom_text(data.frame(text=c(paste(c("FDR pval limit(",opt$thresholdPval,")"),collapse="")),x=c(xMinLimPlot),y=c(pvalThresholdFDR)),mapping=aes(x=x, y=y, label=text), size=4, vjust=0, hjust=0, color="skyblue3")
      if(length(significativePoints)>0)p <- p + geom_point(data=dataSignifToPlot,aes(colour=description))
      
      if(length(htmlPointsToRemove)>0){
        pointToRemove=union(htmlPointsToRemove,significativePoints)
        #to test if remains any normal points to draw
        if(length(pointToRemove)<nrow(data.fit.eb$p.value)){
          dataToPlot=data.frame(pval=-log10(data.fit.eb$p.value[-pointToRemove,iComparison]),FC=data.fit.eb$coefficients[-pointToRemove,iComparison],description=paste("FC: " , round(2^data.fit.eb$coefficients[-pointToRemove,iComparison],2) , " | FDR p-val: ", prettyNum(data.fit.eb$adj_p.value[-pointToRemove,iComparison],digits=4), sep=""))
        }else{
          dataToPlot=data.frame(pval=0,FC=0,description="null")
        }
      }
      
      phtml <- plot_ly(data=dataToPlot, x=~FC, y=~pval,type="scatter", mode="markers",showlegend = FALSE, marker = list(color="gray",opacity=0.5), text=~description, hoverinfo="text") %>%
        layout(title = humanReadingContrastsRenamed[iComparison],xaxis=list(title="Log2 Fold Change",showgrid=TRUE, zeroline=FALSE),yaxis=list(title="-log10(p-val)", showgrid=TRUE, zeroline=FALSE))
      if(length(significativePoints)>0) phtml=add_markers(phtml,data=dataSignifToPlot, x=~FC, y=~pval, mode="markers" , marker=list( color=log10(abs(dataSignifToPlot$FC)*dataSignifToPlot$pval),colorscale='Rainbow'), text=~description, hoverinfo="text", inherit = FALSE) %>% hide_colorbar()
      phtml=add_trace(phtml,x=c(-FCthreshold,-FCthreshold), y=c(0,yMaxLimPlot), type="scatter", mode = "lines", line=list(color="coral",dash="dash"), hoverinfo='none', showlegend = FALSE,inherit = FALSE)
      phtml=add_annotations(phtml,x=-FCthreshold,y=0,xref = "x",yref = "y",text = paste(c("log2(1/FC=",opt$thresholdFC,")"),collapse=""),xanchor = 'right',showarrow = F,textangle=270,font=list(color="coral"))
      phtml=add_trace(phtml,x=c(FCthreshold,FCthreshold), y=c(0, yMaxLimPlot), type="scatter",  mode = "lines", line=list(color="coral",dash="dash"), hoverinfo='none', showlegend = FALSE,inherit = FALSE)
      phtml=add_annotations(phtml,x=FCthreshold,y=0,xref = "x",yref = "y",text = paste(c("log2(FC=",opt$thresholdFC,")"),collapse=""),xanchor = 'right',showarrow = F,textangle=270,font=list(color="coral"))
      if(!is.null(pvalThresholdFDR)){
        phtml=add_trace(phtml,x=c(xMinLimPlot,xMaxLimPlot), y=c(pvalThresholdFDR,pvalThresholdFDR), type="scatter",  mode = "lines", line=list(color="cornflowerblue",dash="dash"), hoverinfo='none', showlegend = FALSE,inherit = FALSE)
        phtml=add_annotations(phtml,x=xMinLimPlot,y=pvalThresholdFDR+0.1,xref = "x",yref = "y",text = paste(c("FDR pval limit(",opt$thresholdPval,")"),collapse=""),xanchor = 'left',showarrow = F,font=list(color="cornflowerblue"))
      }
    plotVector[[length(plotVector)+1]]=p
    
    #save plotly files
    pp <- ggplotly(phtml)
    htmlwidgets::saveWidget(as_widget(pp), paste(c(file.path(getwd(), "plotLyDir"),"/",opt$volcano,"_",correspondanceTable[humanReadingContrastsRenamed[iComparison],1],".html"),collapse=""),selfcontained = F)
    
    
    if(iComparison==nbComparisons || length(plotVector)==volcanoPerPage){
      #plot and close the actual plot
      if(opt$format=="pdf"){
        pdf(paste(c("./plotDir/",opt$volcano,"_",correspondanceTable[humanReadingContrastsRenamed[iComparison],1],".pdf"),collapse=""))}else{
          png(paste(c("./plotDir/",opt$volcano,"_",correspondanceTable[humanReadingContrastsRenamed[iComparison],1],".png"),collapse=""))
        }
      multiplot(plotlist=plotVector,cols=1)
      dev.off()
      if(iComparison<nbComparisons){
        #prepare for a new ploting file if necessary
        plotVector=list()
        iToPlot=iToPlot+1
      }
    }
  }
  remove(dataToPlot,dataSignifToPlot)
  addComment("[INFO]Volcanos drawn",T,opt$log,T,display=FALSE)
}

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
  
#write(unlist(dimnames(data.fit.eb$adj_p.value)),opt$log,append = T)

#filter out genes with higher p-values for all comparisons
genesToKeep=names(which(apply(data.fit.eb$adj_p.value,1,function(x)length(which(x<=opt$thresholdPval))>0)))
if(length(genesToKeep)>0){
  data.fit.eb$adj_p.value=matrix(data.fit.eb$adj_p.value[genesToKeep,],ncol=ncol(data.fit.eb$adj_p.value))
  rownames(data.fit.eb$adj_p.value)=genesToKeep
  colnames(data.fit.eb$adj_p.value)=colnames(data.fit.eb$p.value)
  
  data.fit.eb$p.value=matrix(data.fit.eb$p.value[genesToKeep,],ncol=ncol(data.fit.eb$p.value))
  rownames(data.fit.eb$p.value)=genesToKeep
  colnames(data.fit.eb$p.value)=colnames(data.fit.eb$adj_p.value)
  
  data.fit.eb$coefficients=matrix(data.fit.eb$coefficients[genesToKeep,],ncol=ncol(data.fit.eb$coefficients))
  rownames(data.fit.eb$coefficients)=genesToKeep
  colnames(data.fit.eb$coefficients)=colnames(data.fit.eb$adj_p.value)
}else{
  addComment("[WARNING]No significative genes",T,opt$log,display=FALSE)
}

addComment("[INFO]Significant genes filtering done",T,opt$log,T,display=FALSE)


#plot VennDiagramm for genes below threshold between comparisons
#t=apply(data.fit.eb$adj_p.value[,1:4],2,function(x)names(which(x<=opt$threshold)))
#get.venn.partitions(t)
#vennCounts(data.fit.eb$adj_p.value[,1:4]<=opt$threshold)

#make a simple sort genes based only on the first comparison
#newOrder=order(data.fit.eb$adj_p.value[,1])
#data.fit.eb$adj_p.value=data.fit.eb$adj_p.value[newOrder,]

#alternative sorting strategy based on the mean gene rank over all comparisons
if(length(genesToKeep)>1){
  currentRank=rep(0,nrow(data.fit.eb$adj_p.value))
  for(iComparison in 1:ncol(data.fit.eb$adj_p.value)){
    currentRank=currentRank+rank(data.fit.eb$adj_p.value[,iComparison])
  }
  currentRank=currentRank/ncol(data.fit.eb$adj_p.value)
  newOrder=order(currentRank)
  
  data.fit.eb$adj_p.value=matrix(data.fit.eb$adj_p.value[newOrder,],ncol=ncol(data.fit.eb$adj_p.value))
  rownames(data.fit.eb$adj_p.value)=rownames(data.fit.eb$p.value)[newOrder]
  colnames(data.fit.eb$adj_p.value)=colnames(data.fit.eb$p.value)
  
  data.fit.eb$p.value=matrix(data.fit.eb$p.value[newOrder,],ncol=ncol(data.fit.eb$p.value))
  rownames(data.fit.eb$p.value)=rownames(data.fit.eb$adj_p.value)
  colnames(data.fit.eb$p.value)=colnames(data.fit.eb$adj_p.value)
  
  data.fit.eb$coefficients=matrix(data.fit.eb$coefficients[newOrder,],ncol=ncol(data.fit.eb$coefficients))
  rownames(data.fit.eb$coefficients)=rownames(data.fit.eb$adj_p.value)
  colnames(data.fit.eb$coefficients)=colnames(data.fit.eb$adj_p.value)
}


#formating output matrix depending on genes to keep
if(length(genesToKeep)==0){
  outputData=matrix(0,ncol=ncol(data.fit.eb$adj_p.value)*4+2,nrow=3)
  outputData[1,]=c("X","X",rep(humanReadingContrastsRenamed,each=4))
  outputData[2,]=c("X","X",rep(c("p-val","FDR.p-val","FC","log2(FC)"),ncol(data.fit.eb$adj_p.value)))
  outputData[,1]=c("LIMMA","Gene","noGene")
  outputData[,2]=c("Comparison","Info","noInfo")
}else{
  if(length(genesToKeep)==1){
    outputData=matrix(0,ncol=ncol(data.fit.eb$adj_p.value)*4+2,nrow=3)
    outputData[1,]=c("X","X",rep(humanReadingContrastsRenamed,each=4))
    outputData[2,]=c("X","X",rep(c("p-val","FDR.p-val","FC","log2(FC)"),ncol(data.fit.eb$adj_p.value)))
    outputData[,1]=c("LIMMA","Gene",genesToKeep)
    outputData[,2]=c("Comparison","Info","na")
    if(!is.null(rowItemInfo))outputData[3,2]=rowItemInfo[genesToKeep]
    outputData[3,seq(3,ncol(outputData),4)]=prettyNum(data.fit.eb$p.value,digits=4)
    outputData[3,seq(4,ncol(outputData),4)]=prettyNum(data.fit.eb$adj_p.value,digits=4)
    outputData[3,seq(5,ncol(outputData),4)]=prettyNum(2^data.fit.eb$coefficients,digits=4)
    outputData[3,seq(6,ncol(outputData),4)]=prettyNum(data.fit.eb$coefficients,digits=4)
  }else{
    #format matrix to be correctly read by galaxy (move headers in first column and row)
    outputData=matrix(0,ncol=ncol(data.fit.eb$adj_p.value)*4+2,nrow=nrow(data.fit.eb$adj_p.value)+2)
    outputData[1,]=c("X","X",rep(humanReadingContrastsRenamed,each=4))
    outputData[2,]=c("X","X",rep(c("p-val","FDR.p-val","FC","log2(FC)"),ncol(data.fit.eb$adj_p.value)))
    outputData[,1]=c("LIMMA","Gene",rownames(data.fit.eb$adj_p.value))
    outputData[,2]=c("Comparison","Info",rep("na",nrow(data.fit.eb$adj_p.value)))
    if(!is.null(rowItemInfo))outputData[3:nrow(outputData),2]=rowItemInfo[rownames(data.fit.eb$adj_p.value)]
    outputData[3:nrow(outputData),seq(3,ncol(outputData),4)]=prettyNum(data.fit.eb$p.value,digits=4)
    outputData[3:nrow(outputData),seq(4,ncol(outputData),4)]=prettyNum(data.fit.eb$adj_p.value,digits=4)
    outputData[3:nrow(outputData),seq(5,ncol(outputData),4)]=prettyNum(2^data.fit.eb$coefficients,digits=4)
    outputData[3:nrow(outputData),seq(6,ncol(outputData),4)]=prettyNum(data.fit.eb$coefficients,digits=4)
  }
}
addComment("[INFO]Formated output",T,opt$log,display=FALSE) 

#write output results
write.table(outputData,file=opt$outputFile,quote=FALSE,sep="\t",col.names = F,row.names = F)

end.time <- Sys.time()
addComment(c("[INFO]Total execution time for R script:",as.numeric(end.time - start.time,units="mins"),"mins"),T,opt$log,display=FALSE)

addComment("[INFO]End of R script",T,opt$log,display=FALSE)

#sessionInfo()



