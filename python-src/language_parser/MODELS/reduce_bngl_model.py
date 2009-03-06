#!/usr/bin/python
import sys, re, pdb


molLines = ["1	Cell", "12	Ste11(Ste5_site,MAPK_site,S302_S306_T307~none~pS~pSpS~pSpSpT,Feedback_PO4~none~PO4)"]
seedLines = ["15	Ste20(Ste4_site,Ste11_site) Ste20_tot_conc", "21 	Ste4(Gpa1_site!1,Ste5_site,Ste20_site).Gpa1(Ste2_site,Sst2_site,Ste4_site!1,nucleotide~GDP) Ste4_tot_conc"]
rxnLines = ["3	Dig1(MAPK_site,PO4_site~none) + Fus3(target_site,T180~PO4,Y182~PO4) <-> Dig1(MAPK_site!1,PO4_site~none).Fus3(target_site!1,T180~PO4,Y182~PO4) kon_Fus3pTpY_Dig1	koff_Fus3pTpY_Dig1", "28	Ste12(Dig1_site,Dig2_site!+,MAPK_site) + Dig1(Ste12_site,PO4_site~PO4) <-> Ste12(Dig1_site!1,Dig2_site!+,MAPK_site).Dig1(Ste12_site!1,PO4_site~PO4) kon_Ste12Dig2_Dig1PO4 koff_Ste12Dig2_Dig1PO4"]

def uniq(array):
    if len(array) <= 1:
        return array
    
    uniqArray = [array[0]]
    lastEl = array[0]
    
    for ndx in range(1, len(array)):
        if array[ndx] != lastEl:
            lastEl = array[ndx]
            uniqArray.append(lastEl)
    return uniqArray

def removeNumberingOnLine(line):
    return re.sub("^\d*\s*", "", line)

def getBitsFromMoleculeTypeLine(line):
    line = removeNumberingOnLine(line)
    try:
        parenNdx = line.index("(")
        line = line[:parenNdx]
    except:
        pass

    line = line.strip()
    return [line]

def getBitsFromSeedSpecies(line):
    # Input: "12 Ste11(Ste5_site,MAPK_site,S302_S306_T307~none,Feedback_PO4~none) Ste11_tot_conc"

    line = removeNumberingOnLine(line)
    # "Ste11(Ste5_site,MAPK_site,S302_S306_T307~none,Feedback_PO4~none) Ste11_tot_conc"
    
    
    line = line.split()
    
    # "[Ste11(Ste5_site,MAPK_site,S302_S306_T307~none,Feedback_PO4~none), Ste11_tot_conc]"
    # "[Ste11(Ste5_site,MAPK_site,S302_S306_T307~none,Feedback_PO4~none)]"
    line = line[0]
    individualMols = line.split('.')
        
    theBits = []
    for ndx in indexes(individualMols):
        a = getBitsFromMoleculeTypeLine(individualMols[ndx])
        theBits.extend(a)

    return theBits

def combineLinesOnBlock(block):
    pdb.set_trace()
    
    block = [l.strip() for l in block]
    string = ""
    for line in block:
        string += (line + '\n')
    newLine = re.sub(r"\\\s*\n\s*", " ", string)
    return newLine.split('\n')

def getBitsFromRxnRule(line):
    line = removeNumberingOnLine(line)
    pieces = line.split()
    
    if "<->" in pieces:
        pieces = pieces[:-2]
    else:
        pieces = pieces[:-1]
    

    pieces = [x for x in pieces if "->" not in x]

    pieces = [x for x in pieces if x != "+"]

    theBits = []
    for plex in pieces:
        theBits.extend( getBitsFromSeedSpecies(plex))
    

    theBits.sort()
    return uniq(theBits)
    
    
def filterBlock(originalBlock, bitsGetter, excludeTypes, includeTypes):
    firstLine = originalBlock[0]
    opBlock = originalBlock[1:-1]
    lastLine = originalBlock[-1]
    
    for excluded in excludeTypes:
        opBlock = [line for line in opBlock if excluded not in bitsGetter(line)]

    newBlock = [firstLine]
    newBlock.extend(opBlock)
    newBlock.append(lastLine)

    return newBlock

def removeNewLine(aLine):
    while len(aLine) > 0 and aLine[-1] == "\n":
        aLine = aLine[:-1]
    return aLine
    
def writeArrayToFile(anArray, aFile):
    anArray = [removeNewLine(line) for line in anArray]
    for line in anArray:
        aFile.write(line + '\n')
    aFile.write('\n')

def indexes(x):
    return range(len(x))

def print_array(anArray):
    for l in anArray: print l

def barf(message):
    if message: print message
    print "Crashing...."
    sys.exit(1)

def main():
    if not "-f" in sys.argv:
        barf("No file specified (must be in form ./reducemodel -f filename)")

    filename = sys.argv[ sys.argv.index("-f") + 1]
    outputFilename = 'reduced_' + filename
    outputFile = open(outputFilename, 'w')
    linesinfile = open(filename).read().split("\n")
    linesinfile = [line.strip() for line in linesinfile]

    try:
        beginIndex = linesinfile.index("begin molecule types")
        endIndex = linesinfile.index("end molecule types")
    except:
        barf("Strings 'begin molecule types' and 'end molecule types' were not found")

    moleculeTypes = linesinfile[beginIndex + 1 : endIndex]

    regexp1 = re.compile("^\d*")
    regexp2 = re.compile(r"\(.*\Z")

    moleculeTypes = [ regexp1.sub( "", l) for l in moleculeTypes]
    moleculeTypes = [l.strip() for l in moleculeTypes]
    moleculeTypes = [ regexp2.sub( "", l) for l in moleculeTypes]

    moleculeTypes = [ [x, False] for x in moleculeTypes]

    moleculeTypes.sort(lambda x, y: cmp(x[0], y[0]))
    num = 0
    for mol in moleculeTypes:
        num += 1
        sys.stdout.write(mol[0] + ", ")
        if num % 5 == 0:
            sys.stdout.write('\n')

    sys.stdout.write('\n')

    answerArray = [0,0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,1,0]

    for ndx in indexes(moleculeTypes):
        molInput = raw_input(moleculeTypes[ndx][0] + ": ")
        yesno = molInput.lower()[0]
        # yesno = answerArray[ndx]

        if yesno == "t" or yesno == "y" or yesno == 1:
            include = True
        else:
            include = False

        moleculeTypes[ndx][1] = include


    excludedMolTypes = [x[0] for x in moleculeTypes if not x[1]]
    molTypes = [x[0] for x in moleculeTypes if x[1]]
    outputlines = []

    parameters, moleculeTypes, seedSpecies, reactionRules = macroParseFileLines(linesinfile)


    moleculeTypes = filterBlock(moleculeTypes, getBitsFromMoleculeTypeLine, excludedMolTypes, molTypes)
    seedSpecies = filterBlock(seedSpecies, getBitsFromSeedSpecies, excludedMolTypes, molTypes)

    reactionRules = [x for x in reactionRules if not x == '']
    reactionRules = [x.strip() for x in reactionRules]
    reactionRules = [x for x in reactionRules if not x.startswith("#")]
    pdb.set_trace()
    reactionRules = combineLinesOnBlock(reactionRules)
    reactionRules = filterBlock(reactionRules, getBitsFromRxnRule, excludedMolTypes, molTypes)

    writeArrayToFile(parameters, outputFile)
    writeArrayToFile(moleculeTypes, outputFile)    
    writeArrayToFile(seedSpecies, outputFile)
    writeArrayToFile(reactionRules, outputFile)

    outputFile.close()

    
def macroParseFileLines(bngFile):
    parameters = []
    moleculeTypes = []
    seedSpecies = []
    reactionRules = []

    beginParamNdx = bngFile.index('begin parameters')
    endParamNdx = bngFile.index('end parameters') + 1
    beginMolNdx = bngFile.index('begin molecule types')
    endMolNdx = bngFile.index('end molecule types') + 1
    beginSeedNdx = bngFile.index('begin seed species')
    endSeedNdx = bngFile.index('end seed species') + 1
    beginRxnNdx = bngFile.index('begin reaction rules')
    endRxnNdx = bngFile.index('end reaction rules') + 1

    parameters = bngFile[beginParamNdx : endParamNdx]
    moleculeTypes = bngFile[beginMolNdx : endMolNdx]
    seedSpecies = bngFile[beginSeedNdx : endSeedNdx]
    reactionRules = bngFile[beginRxnNdx : endRxnNdx]

    return parameters, moleculeTypes, seedSpecies, reactionRules

if __name__ == "__main__":
    main()

