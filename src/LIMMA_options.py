import re

def get_column_names( file_path, toNotConsider=None, toNotConsiderBis=None):
	options=[]
	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	for i, field_component in enumerate( firstLine ):
		if i!=0 and field_component!=toNotConsider and field_component!=toNotConsiderBis:#to squeeze the first column
			options.append( ( field_component, field_component, False ) )
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

def get_row_names_interaction( file_path, factorNameA, factorNameB ):
	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	iColumnA=-1
	iColumnB=-1
	for i, field_component in enumerate( firstLine ):
		if field_component==factorNameA:#to test
			iColumnA=i
		if field_component==factorNameB:#to test
			iColumnB=i
	possibleValuesA=[]
	possibleValuesB=[]
	if iColumnA!=-1 and iColumnB!=-1:
		for nextLine in inputfile:
			nextLine=nextLine.strip().split("\t")
			if len(nextLine)>1:
				if nextLine[iColumnA] not in possibleValuesA:
					possibleValuesA.append(nextLine[iColumnA])
				if nextLine[iColumnB] not in possibleValuesB:
					possibleValuesB.append(nextLine[iColumnB])
	inputfile.close()	
	options=[]
	if len(possibleValuesA)>=1 and len(possibleValuesB)>=1 and possibleValuesA[0]!="None" and possibleValuesB[0]!="None":
		for counterA in range(len(possibleValuesA)):
			for counterB in range(len(possibleValuesB)):
				options.append( (possibleValuesA[counterA]+"*"+possibleValuesB[counterB], possibleValuesA[counterA]+"*"+possibleValuesB[counterB], False) )	
	return options

def get_comparisonsA( factorA, valuesA ):
	options=[]
	formatValuesA=re.sub("(^\[u')|('\]$)","", str(valuesA))
	possibleValues=formatValuesA.split("', u'")
	if len(possibleValues)>=2:
		for counter in range(len(possibleValues)-1):
			for innerCounter in range(counter+1,len(possibleValues)):
				options.append( (possibleValues[counter]+" - "+possibleValues[innerCounter], possibleValues[counter]+" - "+possibleValues[innerCounter], False) )
				options.append( (possibleValues[innerCounter]+" - "+possibleValues[counter], possibleValues[innerCounter]+" - "+possibleValues[counter], False) )
	return options

def get_comparisonsAB(factorA, valuesA, factorB, valuesB, interaction):
	options=[]
	formatValuesA=re.sub("(^\[u')|('\]$)","", str(valuesA))
	possibleValuesA=formatValuesA.split("', u'")
	formatValuesB=re.sub("(^\[u')|('\]$)","", str(valuesB))
	possibleValuesB=formatValuesB.split("', u'")
	if str(interaction)=="False":
		if len(possibleValuesA)>=2:
			for counter in range(len(possibleValuesA)-1):
				for innerCounter in range(counter+1,len(possibleValuesA)):
					options.append( (possibleValuesA[counter]+" - "+possibleValuesA[innerCounter], possibleValuesA[counter]+" - "+possibleValuesA[innerCounter], False) )
					options.append( (possibleValuesA[innerCounter]+" - "+possibleValuesA[counter], possibleValuesA[innerCounter]+" - "+possibleValuesA[counter], False) )
		if len(possibleValuesB)>=2:
			for counter in range(len(possibleValuesB)-1):
				for innerCounter in range(counter+1,len(possibleValuesB)):
					options.append( (possibleValuesB[counter]+" - "+possibleValuesB[innerCounter], possibleValuesB[counter]+" - "+possibleValuesB[innerCounter], False) )
					options.append( (possibleValuesB[innerCounter]+" - "+possibleValuesB[counter], possibleValuesB[innerCounter]+" - "+possibleValuesB[counter], False) )
	else:
		if len(possibleValuesA)>=1 and len(possibleValuesB)>=1 and possibleValuesA[0]!="None" and possibleValuesB[0]!="None":
			for counterA in range(len(possibleValuesA)):
				for innerCounterA in range(len(possibleValuesA)):
					for counterB in range(len(possibleValuesB)):
						for innerCounterB in range(len(possibleValuesB)):
							if not(counterA==innerCounterA and counterB==innerCounterB):
								options.append( ("("+possibleValuesA[counterA]+" * "+possibleValuesB[counterB]+") - ("+possibleValuesA[innerCounterA]+" * "+possibleValuesB[innerCounterB]+")","("+possibleValuesA[counterA]+" * "+possibleValuesB[counterB]+") - ("+possibleValuesA[innerCounterA]+" * "+possibleValuesB[innerCounterB]+")", False) )
	return options

def get_row_names_allInteractions( file_path, factorSelected):
	formatFactors=re.sub("(^\[u')|('\]$)","", str(factorSelected))
	factorsList=formatFactors.split("', u'")
	iColumn=[None] * len(factorsList)
	valuesList=[None] * len(factorsList)

	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	for iField, fieldComponent in enumerate( firstLine ):
		for iFactor, factorComponent in enumerate(factorsList):
			if fieldComponent==factorComponent:
				iColumn[iFactor]=iField
				valuesList[iFactor]=[]

	for nextLine in inputfile:
		nextLine=nextLine.strip().split("\t")
		if len(nextLine)>1:
			for iFactor, factorComponent in enumerate(factorsList):
				if nextLine[iColumn[iFactor]] not in valuesList[iFactor]:
					valuesList[iFactor].append(nextLine[iColumn[iFactor]])
	inputfile.close()

	allCombinations=[]
	for iFactor, factorComponent in enumerate(factorsList):
		if iFactor==0:
			allCombinations=valuesList[iFactor]
		else:
			currentCombinations=allCombinations
			allCombinations=[]	
			for iValue, valueComponent in enumerate(valuesList[iFactor]):
				for iCombination, combination in enumerate(currentCombinations):
					allCombinations.append(combination+"*"+valueComponent)	

	options=[]
	for iCombination, combination in enumerate(allCombinations):
		options.append((combination,combination,False))

	return options

def get_allrow_names( file_path, factorSelected ):
	formatFactors=re.sub("(^\[u')|('\]$)","", str(factorSelected))
	factorsList=formatFactors.split("', u'")
	iColumn=[None] * len(factorsList)
	valuesList=[None] * len(factorsList)

	inputfile = open(file_path)
	firstLine = next(inputfile).strip().split("\t")
	for iField, fieldComponent in enumerate( firstLine ):
		for iFactor, factorComponent in enumerate(factorsList):
			if fieldComponent==factorComponent:
				iColumn[iFactor]=iField
				valuesList[iFactor]=[]

	for nextLine in inputfile:
		nextLine=nextLine.strip().split("\t")
		if len(nextLine)>1:
			for iFactor, factorComponent in enumerate(factorsList):
				if nextLine[iColumn[iFactor]] not in valuesList[iFactor]:
					valuesList[iFactor].append(nextLine[iColumn[iFactor]])
	inputfile.close()

	allValues=[]
	for iFactor, factorComponent in enumerate(factorsList):
		for iValue, valueComponent in enumerate(valuesList[iFactor]):
			allValues.append(factorComponent+":"+valueComponent)	

	options=[]
	for iValue, valueComponent in enumerate(allValues):
		options.append((valueComponent,valueComponent,False))

	return options

def replaceNamesInFiles(expressionFile_name,conditionFile_name,outputExpressionFile,outputConditionFile,ouputDictionnary):
	dico={}
	forbidenCharacters={"*",":",",","|"}
	##start with expression file, read only the first line
	inputfile = open(expressionFile_name)
	outputfile = open(outputExpressionFile, 'w')
	firstLine = next(inputfile).rstrip().split("\t")
	iCondition=1
	newFirstLine=""
	for i, field_component in enumerate( firstLine ):
		if (i>0):
			#conditions names should not be redundant with other conditions
			if(field_component not in dico):
				dico[field_component]="Condition"+str(iCondition)
				newFirstLine+="\t"+"Condition"+str(iCondition)
				iCondition+=1
			else:
				raise NameError('condition name allready exists!')
		else:
			newFirstLine+=field_component
	outputfile.write(newFirstLine+"\n")
	for line in inputfile:
		outputfile.write(line)
	outputfile.close()
	inputfile.close()
	#then parse condition file, read all lines in this case
	inputfile = open(conditionFile_name)
	outputfile = open(outputConditionFile, 'w')
	firstLine=1
	iFactor=1
	iValue=1
	for line in inputfile:
		currentLine = line.rstrip().split("\t")
		newCurrentLine=""
		for i, field_component in enumerate( currentLine ):
			#special treatment for the first line
			if (firstLine==1):
				if (i==0):
					newCurrentLine=field_component
				else:
					#factor names should not be redundant with other factors or conditions
					if(field_component not in dico):
						dico[field_component]="Factor"+str(iFactor)
						newCurrentLine+="\t"+"Factor"+str(iFactor)
						iFactor+=1
					else:
						raise NameError('factor name allready exists!')
			else:	
				if (i==0):
					#check if condition name allready exist and used it if it is, or create a new one if not
					if(field_component not in dico):
						dico[field_component]="Condition"+str(iCondition)
						newCurrentLine="Condition"+str(iCondition)
						iCondition+=1
					else:
						newCurrentLine=dico[field_component]
				else:
					if(field_component not in dico):
						dico[field_component]="Value"+str(iValue)
						newCurrentLine+="\tValue"+str(iValue)
						iValue+=1
					else:
						newCurrentLine+="\t"+dico[field_component]
		outputfile.write(newCurrentLine+"\n")
		firstLine=0
	outputfile.close()
	inputfile.close()
	##check if any entries in dictionnary contains forbiden character
	for key, value in dico.iteritems():
		for specialCharacter in forbidenCharacters:
			if value.startswith("Condition")==False and key.find(specialCharacter)!=-1:
				return 1
	##then write dictionnary in a additional file
	outputfile = open(ouputDictionnary, 'w')
	for key, value in dico.iteritems():
		outputfile.write(key+"\t"+value+"\n")
	outputfile.close()
	return 0


def replaceNamesBlockInFiles(expressionFile_name,conditionFile_name,blockingFile_name,outputExpressionFile,outputConditionFile,outputBlockingFile,ouputDictionnary):
	dico={}
	forbidenCharacters={"*",":",",","|"}
	##start with expression file, read only the first line
	inputfile = open(expressionFile_name)
	outputfile = open(outputExpressionFile, 'w')
	firstLine = next(inputfile).rstrip().split("\t")
	iCondition=1
	newFirstLine=""
	for i, field_component in enumerate( firstLine ):
		if (i>0):
			#conditions names should not be redundant with other conditions
			if(field_component not in dico):
				dico[field_component]="Condition"+str(iCondition)
				newFirstLine+="\t"+"Condition"+str(iCondition)
				iCondition+=1
			else:
				raise NameError('condition name allready exists!')
		else:
			newFirstLine+=field_component
	outputfile.write(newFirstLine+"\n")
	for line in inputfile:
		outputfile.write(line)
	outputfile.close()
	inputfile.close()
	#then parse condition file, read all lines in this case
	iFactor=1
	iValue=1
	for fileNum in range(2):
		if fileNum==0:
			inputfile = open(conditionFile_name)
			outputfile = open(outputConditionFile, 'w')
		else:
			inputfile = open(blockingFile_name)
			outputfile = open(outputBlockingFile, 'w')
		firstLine=1
		for line in inputfile:
			currentLine = line.rstrip().split("\t")
			newCurrentLine=""
			for i, field_component in enumerate( currentLine ):
				#special treatment for the first line
				if (firstLine==1):
					if (i==0):
						newCurrentLine=field_component
					else:
						#factor names should not be redundant with other factors or conditions
						if(field_component not in dico):
							dico[field_component]="Factor"+str(iFactor)
							newCurrentLine+="\t"+"Factor"+str(iFactor)
							iFactor+=1
						else:
							raise NameError('factor name allready exists!')
				else:	
					if (i==0):
						#check if condition name allready exist and used it if it is, or create a new one if not
						if(field_component not in dico):
							dico[field_component]="Condition"+str(iCondition)
							newCurrentLine="Condition"+str(iCondition)
							iCondition+=1
						else:
							newCurrentLine=dico[field_component]
					else:
						if(field_component not in dico):
							dico[field_component]="Value"+str(iValue)
							newCurrentLine+="\tValue"+str(iValue)
							iValue+=1
						else:
							newCurrentLine+="\t"+dico[field_component]
			outputfile.write(newCurrentLine+"\n")
			firstLine=0
		outputfile.close()
		inputfile.close()
	##check if any entries in dictionnary contains forbiden character
	for key, value in dico.iteritems():
		for specialCharacter in forbidenCharacters:
			if value.startswith("Condition")==False and key.find(specialCharacter)!=-1:
				return 1
	##then write dictionnary in a additional file
	outputfile = open(ouputDictionnary, 'w')
	for key, value in dico.iteritems():
		outputfile.write(key+"\t"+value+"\n")
	outputfile.close()
	return 0
