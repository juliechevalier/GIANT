import re
import numpy as np

def get_column_names( file_path, toNotConsider=-1, each=1):
	options=[]
	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	cpt=0
	for i, field_component in enumerate( firstLine ):
		if i!=toNotConsider:#to squeeze the first column
			if cpt==0:
				options.append( ( field_component, field_component, False ) )
			cpt+=1
			if cpt==each:
				cpt=0
	inputfile.close()
	return options

def get_column_names_filteredList( file_path, toNotConsider=[], each=1):
	options=[]
	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	cpt=0
	for i, field_component in enumerate( firstLine ):
		if i not in toNotConsider:#to squeeze the first columns
			if cpt==0:
				options.append( ( field_component, field_component, False ) )
			cpt+=1
			if cpt==each:
				cpt=0
	inputfile.close()
	return options

def get_column_names_mergeNumber(file_path, numberToMerge=1, toNotConsider=[]):
	options=[]
	inputfile = open(file_path)
	if int(numberToMerge)>0:
		iHeader=0
		for iCurrentLine in inputfile:
			iHeader=iHeader+1
			if iHeader>int(numberToMerge):
				break
			currentLine=iCurrentLine.strip().split("\t")
			iOption=-1
			for i, field_component in enumerate( currentLine ):
				if i not in toNotConsider:#to squeeze specified columns
					iOption=iOption+1
					if iHeader==1:
						options.append( ( str(field_component), str(field_component), False ) )
					else:
						options[iOption]=(options[iOption][0]+"_"+str(field_component),options[iOption][1]+"_"+str(field_component),False)
	else:
		currentLine = next(inputfile).strip().split("\t")
		for i, field_component in enumerate( currentLine ):
			if i not in toNotConsider:#to squeeze specified columns
				options.append( ( "Column_"+str(i), "Column_"+str(i), False ) )
	inputfile.close()
	return options

def get_row_names( file_path, factorName ):
	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	iColumn=-1
	for i, field_component in enumerate( firstLine ):
		if field_component==factorName:#to test
			iColumn=i
	options=[]
	if iColumn!=-1:
		for nextLine in inputfile:
			nextLine=nextLine.strip().split("\t")
			if len(nextLine)>1:
				if (nextLine[iColumn], nextLine[iColumn], False) not in options:
					options.append( (nextLine[iColumn], nextLine[iColumn], False) )
	inputfile.close()
	return options

def get_condition_file_names( file_list, toNotConsider=-1, each=1):
	options=[]
	if not isinstance(file_list,list):#if input file is a tabular file, act as get_column_names
		inputfile = open(file_list.file_name)
		firstLine = next(inputfile).strip().split("\t")
		cpt=0
		for i, field_component in enumerate( firstLine ):
			if i!=toNotConsider:#to squeeze the first column
				if cpt==0:
					options.append( ( field_component, field_component, False ) )
				cpt+=1
				if cpt==each:
					cpt=0
		inputfile.close()
	else:#if input file is a .cel file list
		for i, field_component in enumerate( file_list ):
			options.append( ( field_component.name, field_component.name, False ) )
	return options

def generateFactorFile( file_list, factor_list, outputFileName, logFile):
	forbidenCharacters={"*",":",",","|"}
	outputfile = open(outputFileName, 'w')
	outputLog = open(logFile, 'w')
	sampleList=[]
	if not isinstance(file_list,list):
		conditionNames=get_condition_file_names(file_list,0)
	else :
		conditionNames=get_condition_file_names(file_list)	
	for iSample, sample_component in enumerate (conditionNames):
		sampleList.append(str(sample_component[1]))
	outputLog.write("[INFO] "+str(len(sampleList))+" sample are detected as input\n")
	globalDict=dict()
	factorNameList=[]
	firstLine="Conditions"
	if len(factor_list)==0:#check if there is at least one factor available
		outputLog.write("[ERROR] no factor was defined !\n")
		return 1
	else:
		for iFactor, factor_component in enumerate( factor_list ):
			currentSampleList=list(sampleList)
			currentFactor=str(factor_component['factorName'])
			#check if factor name contains forbidden characters
			for specialCharacter in forbidenCharacters:
				if currentFactor.find(specialCharacter)!=-1:
					outputLog.write("[ERROR] '"+specialCharacter+"' character is forbidden in factor name : '"+currentFactor+"'\n")	
					return 4
			#check if factor allready named like that
			if not globalDict.get(currentFactor) is None:
				outputLog.write("[ERROR] '"+currentFactor+"' is used several times as factor name\n")	
				return 3
			globalDict[currentFactor]=dict()
			firstLine=firstLine+"\t"+currentFactor
			factorNameList.append(currentFactor)
			if len(factor_component['valueList'])<=1:#check if there is at least two value available
				outputLog.write("[ERROR] at least two different values are necessary for '"+currentFactor+"' factor\n")
				return 1
			else:
				for iValue, value_component in enumerate( factor_component['valueList'] ):
					currentValue=str(value_component['valueName'])
					#check if factor name contains forbidden characters
					for specialCharacter in forbidenCharacters:
						if currentValue.find(specialCharacter)!=-1:
							outputLog.write("[ERROR] '"+specialCharacter+"' character is forbidden in value name : '"+currentValue+"'\n")	
							return 4
					currentSample=str(value_component['valueConditions']).split(",")
					for iSample, sample_component in enumerate (currentSample):
						if not sample_component in currentSampleList:
							outputLog.write("[ERROR] sample "+sample_component+" was assigned several times for factor '"+currentFactor+"'\n")
							return 2
						currentSampleList.remove(sample_component)
						globalDict[currentFactor][sample_component]=currentValue
			if(len(currentSampleList)>0):
				outputLog.write("[ERROR] for factor '"+currentFactor+"'' sample "+str(currentSampleList)+" are not assigned to any value\n")
				return 2
	outputLog.write("[INFO] "+str(len(globalDict))+" factors are detected\n")
	#start writing the factor file
	outputfile.write(firstLine+"\n") 
	for iSample, sample_component in enumerate(sampleList):
		newLine=sample_component
		for iFactor, factor_component in enumerate(factorNameList):
			newLine=newLine+"\t"+globalDict[factor_component][sample_component]
		outputfile.write(newLine+"\n") 
	outputfile.close()
	outputLog.close()
	return 0

def selectSubSetTable(file_path,headerLine_number,columnsToAdd,columnNamesToKeep,outputFileName,logFile):
	outputLog = open(logFile, 'w')
	outputLog.write("[INFO] header line number : "+ headerLine_number+" lines\n")	
	availableColumnsTuple=get_column_names_mergeNumber(file_path, headerLine_number)
	#convert tuple list as a simple array
	availableColumns=[]
	for iTuple, tuple_content in enumerate (availableColumnsTuple): 
		availableColumns.append(str(tuple_content[0]))
	if len(availableColumns)==0:
		outputLog.write("[ERROR] No detected columns in input file\n")	
		return 1
	selectedColumns=list(columnsToAdd)
	for iVolcano, volcano_content in enumerate(columnNamesToKeep):
		selectedColumns.append(availableColumns.index(volcano_content['pvalColumn']))
		selectedColumns.append(availableColumns.index(volcano_content['fcColumn']))
	if len(selectedColumns)!=(2*len(columnNamesToKeep)+len(columnsToAdd)):
		outputLog.write("[ERROR] matching between input file colnames and requested column names failed\n")	
		return 1
	outputLog.write("[INFO] columns kept : "+str(selectedColumns)+"\n")	
	#start writting formatted file
	inputfile = open(file_path)
	outputfile = open(outputFileName, 'w')
	iLineCpt=-1
	for iCurrentLine in inputfile:
		iLineCpt=iLineCpt+1
		if iLineCpt>=int(headerLine_number):
			currentLineFields=np.array(iCurrentLine.strip().split("\t"))
			newLine="\t".join(currentLineFields[selectedColumns])
			outputfile.write(newLine+"\n")
	if iLineCpt<int(headerLine_number):
		outputLog.write("[ERROR] not enough lines in input files ("+(iLineCpt+1)+" lines)\n")	
		return 1
	inputfile.close()
	outputfile.close()
	outputLog.close()
	return 0