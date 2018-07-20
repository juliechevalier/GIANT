import re

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
