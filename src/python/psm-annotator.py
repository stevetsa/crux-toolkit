"""
Daniel Chee
Noble Lab
June 22, 2009 - 
PSM Annotator:
This program takes a PSM file as input, and is meant ot annotate the different ions 
that are present 
in the data
 
"""

from optparse import OptionParser 
import subprocess
import sys
import math
import numpy as np
import os


def calculateMass(sequence, fix, type):
    """
    Input: the peptide sequence, prefix or suffrix(B or Y), and the type refering to monoisotopic
        or average mass. 
    Output: This function outputs an array of values that correspond to the different fragments
        that occur from the input peptide. 
    Function: This method takes the input peptide of each PSM and returns an array ot values
        corresponding to each ion produced. Used for annotation.
    """
    global averageMasses
    global monoisotopicMasses
    invalidChar = False
    massArray = []
    mass = 0
    charCounter = 0
    for char in sequence:
        if charCounter < (len(sequence) - 1): 
            if averageMasses.get(char) != None:     
                if type == 'mono':
                    mass = mass + monoisotopicMasses.get(char)
                else:
                    mass = mass + averageMasses.get(char)
                massArray.append(mass)
            else:
                invalidChar = True
        charCounter = charCounter + 1
    if fix == 'suffix':
        massArray.reverse()
    return massArray
   
def transform(m,charge,sym):
    """
    Input: The mass value of the ion array. 
    Output: The mass transormed into the mass of the whichever ion is specified. 
    Function: This method takes the input mass from the ion fragment array, and changes the mass 
        in accordance with the different types of ions that occur.  
    """
    if sym == 'b':
        m = (m + (charge*massH))/charge
        return m
    if sym == 'b-h2o':
        m = (m + (charge * massH) - massh2o)/charge
        return m 
    if sym == 'b-nh3':
        m = (m + ( charge * massH ) - massnh3 )/charge
        return m
    if sym == 'b-h2o-h2o':
        m = (m + ( charge * massH ) - massh2o - massh2o )/charge
        return m 
    if sym == 'b-h2o-nh3':
        m = (m + ( charge * massH ) -  massh2o - massnh3 )/charge
        return m 
    if sym == 'a':
        m = (m + (charge*massH) - massCO)/charge
        return m
    if sym == 'a-h2o':
        m = (m + (charge*massH) - massCO - massh2o)/charge
        return m
    if sym == 'a-nh3':
        m = (m + (charge*massH) - massCO - massnh3)/charge
        return m
    if sym == 'y':
        m = (m +(charge*massH) + massh2o)/charge 
        return m
    if sym == 'y-h2o':
        m = (m + (charge*massH))/charge
        return m 
    if sym == 'y-nh3':
        m = (m + (charge*massH) + massh2o - massnh3)/charge
        return m
    if sym == 'y-h2o-h2o':
        m = (m + (charge*massH) - massh2o)/charge
        return m 
    if sym == 'y-h2o-nh3':
        m = (m + (charge*massH) - massnh3)/charge
        return m 

def annotateIon( type , loss , ending  , ID ):
    """
    Input: This method takes an input a tag that is used for annotating purposes. 
    Output: N/A
    Function: This method searches the spectrum for the calculated ions and annotates the spectrum
        based on the ions that it finds. 
    """
    global amIHere
    global wasFound
    global line
    global charges
    global threshold
    global floatError
    global arrayIndex
    global sym2
    global sequenceLength
    if options.test == 'yes':
        global annotatorFinds
        global f2value
    if options.sistats > 1:
        global plotNumber
    if options.pistats > 0:
        global iontag1
        global iontag2
    sym = symGen( amIHere )
    peptideCounter = 1
    if type =='y':
        fixArray = suffixArrays[arrayIndex]
    else:
        fixArray = prefixArrays[arrayIndex]
    for things in fixArray:
        if type =='y':
            var1 = (sequenceLength-1) - peptideCounter + 1
            var2 = (sequenceLength - ((sequenceLength-1) - peptideCounter + 1)) -1
            var3 = (sequenceLength) - peptideCounter - 1
        else:
            var1 = peptideCounter
            var2 = peptideCounter-1
            var3 = peptideCounter - 1  
        answer =  transform(things,charges, (type+loss)) - float(lines[0])
        if (abs(answer ) < (threshold + floatError)) and\
                options.test == "yes" and ending == '@':
            annotatorFinds.append(transform(things , charges , (type+loss)))
        if abs(answer ) < threshold:
            amIHere = True
            wasFound = True
            wasFound = True
            line = line.rstrip()
            if options.annotate == 'yes':
                annotatorfile.write(sym + sym2 + type +str(var1) + "-ch" +\
                                        str(charges) + loss + '('+\
                                        str(transform(things,charges,(type + loss))) +\
                                        ')' )
            if options.test == "yes":
                f2value = f2value + (sym + sym2 + type +str(var1) + "-ch" +\
                                         str(charges) + loss + ending)
            if options.sistats > 1:
                if ID == generateStatsID(plotNumber):
                    signal.append(float(lines[1]))
            if options.pistats != 0:
                if (type + loss) == iontag1:
                    
                    if float(lines[1]) > ionMatrix1[peptideCounter - 1][charges - 1]:
                        ionMatrix1[peptideCounter - 1][charges - 1] = float(lines[1])
                if (type + loss) == iontag2:
                    if float(lines[1]) > ionMatrix2[peptideCounter - 1][charges - 1]:
                        ionMatrix2[peptideCounter - 1][charges - 1] = float(lines[1])
            if options.visualize == 'yes':
                global visualTag
                if abs(answer) < visualTag[0]:
                    if sym2 == '':
                        visualTag = [abs(answer), (type +str(var1) + "-ch" +\
                                            str(charges) + loss)]
                    else:
                        visualTag = [abs(answer), (type +str(var1) + "-ch" +\
                                            str(charges) + loss) + '(' + str(sym2).rstrip('-') +\
                                            ')'] 
                    if sym2 == 'P1-':
                        visualTag.append('red')
                    elif sym2 == 'P2-':
                        visualTag.append('blue')
                    else:
                        if type  == 'b' or type == 'a':
                            visualTag.append('red')
                        else:
                            visualTag.append('blue')
            if stringArrays[arrayIndex][peptideCounter - 1].count(sym2 + type + str(var1) +\
                                                                        "-ch" + str(charges) + loss + ',')\
                                                                        == 0:
                if len(stringArrays[arrayIndex][peptideCounter - 1]) == 0:
                    stringArrays[arrayIndex][peptideCounter - 1] = (stringArrays[arrayIndex])\
                        [peptideCounter - 1] + sym2 + type + str(var1) + "-ch" + str(charges) + loss + ','
                else:
                    stringArrays[arrayIndex][peptideCounter - 1] = (stringArrays[arrayIndex])\
                        [peptideCounter - 1] + sym2 + type + str(var1) + "-ch" + str(charges) + loss + ','
        peptideCounter = peptideCounter + 1
    
def annotateEnding(sequences, annotatorfile, stringArrays):
    """
    Input: A string array of the peptides, the file currently being annotated, the array of strings
        that consist of the differetn ions that occured for each cut. 
    Output: N/A
    Function: This methods produces the print-out of different family ions at the end of each PSM
    """
    seqIndex = 0
    for eachThing in sequences:
        pepIndex = 0 
        annotatorfile.write(">" + eachThing + "\n")
        for eachChar in eachThing:
            if aminoAcids.find(eachChar) >= 0:
                if len(stringArrays[seqIndex][pepIndex]) != 0: 
                    annotatorfile.write(str(pepIndex +1 ) + " " +\
                                       sequences[seqIndex][pepIndex] +\
                                       " ")
                annotatorfile.write(stringArrays[seqIndex][pepIndex].rstrip(','))
                if len(stringArrays[seqIndex][pepIndex]) != 0: 
                    annotatorfile.write("\n")
                pepIndex = pepIndex + 1
        seqIndex = seqIndex +1

def symGen( amIHere):
    """
    Input: A decider variable that keeps track of whether or not the current ion is the first
        one annotated for the current peak. 
    Output: Outputs a symbol that is either a comma or a space. 
    Function: This method is basically used for organizational purposes, putting a space before
        the first ion for the current peak and a comma in front of the following ions.
    """
    if amIHere == True:
        sym = ","
    else:
        sym = " "
    return sym

def peptideNumberGenerator( arrayIndex , sequences):
    """
    Input: An array index that is used to switch between peptides in multi-peptide PSMs,
        a string array of the sequences
    Output: A symbol that is either blank or that represents which peptide the ion 
        belongs to.
    Function: This mehtod basically puts P1, P2, etc in front of the ions in the (multi-peptide)
       annotated specrtum to make it easier to figure out which peptide the ion came from. 
    """
    if len(sequences) > 1:
        sym2 = "P" + str( arrayIndex +1 ) + "-"
    else:
        sym2 = ""
    return sym2

def provideStatsFile(array ,title , file):
    """
    Input: An array of the intensities or values recieved from the input file, an identifier that 
        describes where the array came from, the output statistics file.
    Output: N/A
    Function: This method is a general method used by the stats options to extract the basic 
        stats from an array and write them to a file. 
    """
    file.write(title)
    file.write('MEAN: ' + str(np.mean(array)) + '\n')
    file.write('MEDIAN: ' + str(np.median(array)) + '\n')
    file.write('STD: ' + str(np.std(array)) + '\n')
    file.write('MAX: ' + str(np.nanmax(array)) + '\n')
    file.write('MIN: ' + str(np.nanmin(array)) + '\n\n')
    file.write('\nData:\n')
    for x in range(0, len(array)):
        file.write(str(array[x]) + '\n')
    file.write('\n')

def cruxConverter( array ):
    """
    Input: An array of output values recieved from Crux
    Output:A string that corresponds to the identified ion.
    Function: This method basically converts the values of ion predictions made by crux
        and converts them into the same format used by the annotator. USed for testing.
    """
    so = ''
    if array[3] == '1':
        part1 = "b"
    elif array[3] == '4':
        part1 = "y"

    part2 = array[4]
    part3 = "-ch" + array[2]
    
    if array[5] == '1':
        part5 = "-nh3"
    else:
        part5 = ""

    if array[6] == '1':
        part4 = "-h2o"
    else:
        part4 = ""
    
    string = part1 + part2 + part3 + part4 + part5

    return string

def plotter( array , title , xlabel , ylabel , filename ):
    """
    Input: An array of data, graph title, x label, y label, output file name
    Output:N/A
    Function: This method is a general plotting function which is used to make a histogram for 
         all of the data extracted by the statistics options.
    """
    if options.overwrite == 'no':
        if os.path.exists(directory + filename + '.png'):
            sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
            sys.exit()
    if options.sistats == 1:
        type ='step'
    else:
        type = 'bar'
    upperBound = percentile(np.sort(array) , 0.90)
    pyplot.hist((array),histtype=type, range=(0, upperBound),bins=100)
    pyplot.title( title.upper() )
    pyplot.xlabel( xlabel )
    pyplot.ylabel( ylabel )
    if title != 'NOISE':
        pyplot.savefig( directory + '/stats/plots/' + filename )
        pyplot.clf()

def generateStatsID( num ):
    """
    Input: A number input by the user that refers to a specific ion. 
    Output: Returns a tag string that represents the ion and corresponds with how the annotator
        identifies eahc ion during annotation. 
    Function: This method is used by the single ion stats option to extract the data for a single 
        specific ion. We organized the tool this way so that more that all the statistics could be 
        extracted in parallel.
    """
    idCharge = str((num / 13) + 1)

    if num%13 == 1:
        ID = 'b-ch' + idCharge
    if num%13 == 2:
        ID = 'b-ch' + idCharge + '-h2o'
    if num%13 == 3:
        ID = 'b-ch' + idCharge + '-nh3'
    if num%13 == 4:
        ID = 'b-ch' + idCharge + '-h2o-h2o'
    if num%13 == 5:
        ID = 'b-ch' + idCharge + '-h2o-nh3'
    if num%13 == 6:
        ID = 'a-ch' + idCharge
    if num%13 == 7:
        ID = 'a-ch' + idCharge + '-h2o'
    if num%13 == 8:
        ID = 'a-ch' + idCharge + '-nh3'
    if num%13 == 9:
        ID = 'y-ch' + idCharge
    if num%13 == 10:
        ID = 'y-ch' + idCharge + '-h2o'
    if num%13 == 11:
        ID = 'b-ch' + idCharge + '-nh3'
    if num%13 == 12:
        ID = 'b-ch' + idCharge + '-h2o-h2o' 
    if num%13 == 0:
        ID = 'b-ch' + str((num / 13)) + '-h2o-nh3'

    return ID

def generatePIStatsTag(num):
    """
    Input: The number recieved as input to the pair ion stats option.
    Output: Two tags that correspond to the ion types that will be compared.
    Function: This method is used to identify which specific comparison will be made this run.
        Once again things were set up this way so that a large data set could be ran in parallel.
    """
    ionTags = ['b','b-h2o','b-nh3','b-h2o-h2o','b-h2o-nh3','a','a-h2o','a-nh3','y',
               'y-h2o','y-nh3','y-h2o-h2o','y-h2o-nh3']
    correlationTags = []
    for tagNum in range(0, (len(ionTags) - 1) ):
        for tagNum2 in range((tagNum +1) , (len(ionTags))):
            correlationTags.append([ionTags[tagNum],ionTags[tagNum2]])
    tag1 = correlationTags[(num%78)-1][0]
    tag2 = correlationTags[(num%78)-1][1]
    return tag1, tag2

def appendIonTypes(iontag1, iontag2, iontype1, iontype2, ionMatrix1, ionMatrix2, charge):
    """
    This method is used to update that matrices that contain the intensities of the ions that are 
    being compared in the pair-ion stats option. 
    """
    global counterX
    global counterY
    if iontag1[0] == 'b' and iontag2[0] == 'b' or\
            iontag1[0] == 'b' and iontag2[0] == 'a' or\
            iontag1[0] == 'a' and iontag2[0] == 'b' or\
            iontag1[0] == 'y' and iontag2[0] == 'y':
        for x in range(0, len(ionMatrix1)):
            for y in range(0, len(ionMatrix1[0])):
                if (ionMatrix1[x][y] + ionMatrix2[x][y]) != 0:
                    iontype1 = np.append(iontype1, ionMatrix1[x][y])
                    iontype2 = np.append(iontype2, ionMatrix2[x][y])
                    if ionMatrix1[x][y] == 0 or ionMatrix2[x][y] == 0:
                        counterY = counterY +1
                    else:
                        counterX = counterX +1
    else: 
        for x in range(0, len(ionMatrix1)):
            for y in range(0, len(ionMatrix1[0])):
                if (x+y) == charge:
                    if (ionMatrix1[x][y] + ionMatrix2[x][y]) != 0:
                        iontype1 = np.append(iontype1, ionMatrix1[x][y])
                        iontype2 = np.append(iontype2, ionMatrix2[x][y])
                        if ionMatrix1[x][y] == 0 or ionMatrix2[x][y] == 0:
                            counterY = counterY +1
                        else:
                            counterX = counterX +1
    return iontype1, iontype2

def getCruxPrediction(cruxValues, sequences, theCharges):
    """
    Input: An empty array that will be filled by this method,an array of peptides for 
        each PSM, an array of the charges for the PSM.
    Output: An array of values that will be used to search each spectrum for ions, 
    Function: This method calls crux and obtains an array of values that will be used to 
        search the spectrum for ions. This search will be used to check whether or not 
        the annotator is working.
    """
    sequenceIndex = 0
    testIndex = 0
    for eachSequence in sequences:
        if options.nl == "none":
            tester = subprocess.Popen(['/net/noble/vol2/home/dchee7/Working-copy/' +\
                                      'ms2denoising/bin/thirdparty/arch/' + arch +\
                                      '/crux-predict-peptide-ions-modified',\
                                      eachSequence, theCharges[testIndex].lstrip('+' and '-')],\
                                      stdout = subprocess.PIPE)      
            testOutput = tester.communicate()
            testStrings = testOutput[0].rsplit('\n')
        if options.nl == 'all':
            tester = subprocess.Popen(['/net/noble/vol2/home/dchee7/Working-copy/' +\
                                      'ms2denoising/bin/thirdparty/arch/' + arch +\
                                      '/crux-predict-peptide-ions-modified',\
                                      '--nh3', '1' , '--h2o' , '1'  , eachSequence ,\
                                      theCharges[testIndex].lstrip('+' and '-')], \
                                      stdout = subprocess.PIPE)
            testOutput = tester.communicate()
            testStrings = testOutput[0].rsplit("\n")
        for strings in testStrings:
            if strings.startswith("#") != True and strings.startswith("m") != True:
                cruxValues.append(strings.rsplit("\t"))
        testIndex = testIndex + 1
    #sequenceIndex = sequenceIndex + 1
    return cruxValues

def testAnnotator():
    """
    Input:N/A
    Output:N/A
    Function: This method actually compares the list of ions found by crux with a list of values 
        found by the annotator. It prints out a report that states whether or not the annotator passed.
        It gives a line by line print out so which makes finding the point of failure very easy. 
    """
    global f2value
    global cruxValue
    global cruxFinds
    global annotatorFinds
    global threshold
    global floatError
    global testFile
    global line 
    global failed        
    f2value  = line + f2value  
    cruxValue = cruxValue.rstrip()
    f2value = f2value.rstrip()
    cruxStrings = cruxValue.rsplit(" ")
    f2Strings = f2value.rsplit(" ")
    if len(cruxStrings) == 2 and len(f2Strings) == 2:
        passAnnotator(line, '', '', '' , '', '', '\n')
    elif len(cruxStrings) == 2 and len(f2Strings) != 2:
        if f2Strings[2].count("@") > 0:
            if (len(cruxFinds) == 0) or (len(annotatorFinds) == 0):
                failAnnotator(line, '\t', '#', 'NONE\t', '*', f2Strings[2], '\n', 1)
            elif len(cruxFinds) != len(annotatorFinds):
                failAnnotator(line, '\t', '#', 'NONE\t', '*', f2Strings[2], '\n', 2)
            else:
                for findsIndex in range(0,len(cruxFinds)):
                    if abs(float(cruxFinds[findsIndex])- float(annotatorFinds[findsIndex])) >\
                            floatError:
                        failAnnotator(line, '\t', '#', 'NONE\t', '*', f2Strings[2], '\n', 3)
                    else:
                        passAnnotator(line, '\t' , '#' , 'NONE\t', '*', f2Strings[2], 
                                      ' Within floating point error threshold.\n')
        else:
            passAnnotator(line, '\t' , '#' , 'NONE\t', '*', f2Strings[2] , '\n')
    elif len(cruxStrings) != 2 and len(f2Strings) == 2:
        if f2Strings[1].count("@") > 0:
            if (len(cruxFinds) == 0) or (len(annotatorFinds) == 0):
                failAnnotator(line, '\t', '#', cruxStrings[2], '*', 'NONE\t', '\n', 4)
            elif len(cruxFinds) != len(annotatorFinds):
                failAnnotator(line, '\t', '#', cruxStrings[2], '*', 'NONE\t', '\n', 5)
            else:
                for findsIndex in range(0,len(cruxFinds)):
                    if abs(float(cruxFinds[findsIndex])- float(annotatorFinds[findsIndex])) >\
                            floatError:
                        failAnnotator(line, '\t', '#', cruxStrings[2], '*', 'NONE\t', '\n', 6)
                    else:
                        passAnnotator(line, '\t', '#', cruxStrings[2], '*', 'NONE\t',
                                      ' Within the floating point error threshold\n') 
    elif len(cruxStrings) > 2 and len(f2Strings) > 2:
        cruxSubStrings = cruxStrings[2].rsplit(",")
        for subStrings in cruxSubStrings:
            if f2Strings[2].find(subStrings) >= 0:
                passAnnotator(line, '\t', '#', (subStrings + '\t'), '*', f2Strings[2], '\n')
            else:
                if len(cruxFinds) != len(annotatorFinds):
                    failAnnotator(line, '\t', '#', cruxStrings[2], '*', f2Strings[2], '\n', 7)
                else:
                    for findsIndex in range(0,len(cruxFinds)):
                        cruxFinds.sort()
                        annotatorFinds.sort()
                        if abs(float(cruxFinds[findsIndex])- float(annotatorFinds[findsIndex])) >\
                                floatError: 
                            failAnnotator(line, '\t', '#', cruxStrings[2], '*', f2Strings[2],\
                                              '\n', 8)
                        else:
                            passAnnotator(line, '\t', '#', cruxStrings[2], '*', f2Strings[2],
                                          ' Within floating point error threshold\n')
    else:
        print "ERROR: INVALID INPUT"
        sys.exit()

def passAnnotator(var1, var2, var3, var4, var5, var6 , var7):
    """
    These methods were added in an effort to reduce repetative code.They basically write pass
    or fail to the test report and then prints out the line that cause each pass or fail.
    """
    global testFile
    testFile.write("PASS\t")
    testFile.write(var1 + var2 + var3 + var4 + var5 + var6 + var7)

def failAnnotator(var1, var2, var3, var4, var5, var6, var7, code):
    """
    These methods were added in an effort to reduce repetative code.They basically write pass
    or fail to the test report and then prints out the line that cause each pass or fail.
    """
    global failed
    global testFile
    testFile.write("FAIL\t")
    testFile.write(var1 + var2 + var3 + var4 + var5 +var6 + var7)
    failed = True
    print 'The Annotator Failed at area ' + str(code)

def percentile(N, percent, key=lambda x:x):
    """
    This method was pulled from an internet posting and was not written originally by Daniel Chee. 

    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

def normalize(array):
    """
    Input: An array of data
    Output: outputs a new normalized array
    Function: This method basically takes in an array of values which are normalized by dividing them 
        by the largest value. Makes every value between 0 and 1. 
    """
    NA = np.array([])
    for items in array:
        NA = np.append(NA, items/(np.max(array)))
    return NA

def rank(array):
    """
    Input: An array of data.
    Output: An array of ranks in the same order as the original array. 
    Function: This method basically ranks each value (1 being the highest) in the input array based on 
        the magnitude of each value compared to reast of the array.
    """
    RA = np.array([])
    array2 = np.sort(array)
    for x in range(0, len(array2)):
        for items in array:
            if items == array2[x]:
                RA = np.append(RA, (len(array2) - x))
    return RA

def relatem2q(array, pepmass):
    """
    Input: An array of data.
    Output: An array of m/q values relative to the total pepmass.
    Function: This method basically takes in an array of m/q values and divides each va;ue by the 
        overall pepmass.
    """
    RM = np.array([])
    for items in array:
        RM = np.append(RM, items/float(pepmass)) # works for single peptides
    return RM

def buildFamily(mass, PSMarray, pepmass):
    familyArray = []
    ##### Has some extra pieces, I don't know yet if the are necessary. 
    #familyArray = [(mass - massh2o), (mass + massh2o), (mass - massnh3), (mass + massnh3), (mass - massh2o - massh2o),
     #              (mass + massh2o + massh2o), (mass - massh2o - massnh3), (mass + massh2o + massnh3),
      #             (mass - massh2o + massnh3), (mass + massh2o - massnh3), (mass - massCO), (mass - massCO - massh2o),
       #            (mass - massCO + massh2o), (mass - massCO - massnh3), (mass - massCO + massnh3), (mass + massCO),
        #           (mass + massCO - massh2o), (mass + massCO + massh2o), (mass + massCO - massnh3),
         #          (mass + massCO + massnh3)]
    familyArray1 = [(mass - massh2o), (mass + massh2o), (mass - massnh3), (mass + massnh3), (mass - massh2o - massh2o),
                   (mass + massh2o + massh2o), (mass - massh2o - massnh3), (mass + massh2o + massnh3),
                   (mass - massCO), (mass - massCO - massh2o),
                   (mass - massCO + massh2o), (mass - massCO - massnh3), (mass - massCO + massnh3), (mass + massCO),
                   (mass + massCO - massh2o), (mass + massCO + massh2o), (mass + massCO - massnh3),
                   (mass + massCO + massnh3)]
    ## This works for the test data but the 2 needs to be replaced with the peptide charge
    mass = (pepmass * 2) - mass
    ##
    familyArray2 = [(mass - massh2o), (mass + massh2o), (mass - massnh3), (mass + massnh3), (mass - massh2o - massh2o),
                   (mass + massh2o + massh2o), (mass - massh2o - massnh3), (mass + massh2o + massnh3),
                   (mass - massCO), (mass - massCO - massh2o),
                   (mass - massCO + massh2o), (mass - massCO - massnh3), (mass - massCO + massnh3), (mass + massCO),
                   (mass + massCO - massh2o), (mass + massCO + massh2o), (mass + massCO - massnh3),
                   (mass + massCO + massnh3)]
    found = []
    for peaks in PSMarray:
        for x in range(0, len(familyArray1)):
            if abs(familyArray1[x] - float(peaks.rsplit(' ')[0])) < .25:
                found.append(float(peaks.rsplit(' ')[0]))
            if abs(familyArray2[x] - float(peaks.rsplit(' ')[0])) < .25:
                found.append(float(peaks.rsplit(' ')[0]))
    while (len(found) < 18):
        found.append(0)
    found.sort(reverse=True)
    print found
     
def extractFeatures(PSMarray, pepmass):
    for peak in PSMarray:
        buildFamily(float(peak.rsplit(' ')[0]), PSMarray, pepmass)
    
def blDenoise(array, denoisedFile):
    """
    Input: An array of peaks recieved from the input file.
    Output: N/A
    Function: This method basically writes a file where the bottom Nth percentage of the 
        intensities have been removed. 
    """
    percent = options.bld * .01
    removed = int (round(len(array) * percent)) - 1
    array2 = []
    for items in array:
        items.reverse()
    array.sort()
    for x in range(removed, len(array)):
        array2.append(array[x])
    for items in array2:
        items.reverse()
    array2.sort()
    for items in array2:
        denoisedFile.write(str(items[0]) + ' ' + str(items[1]) + '\n')





############################################################################################
############################################################################################


####Here I am implementing the option function

########################################
usage = 'psm-annotator.py [options] <input_psm_file> <output_directory>'
parser = OptionParser(usage = usage)

parser.add_option('--annotate' , action = 'store' , type = 'string' , dest = 'annotate' ,\
                      default = 'yes' , help = 'Takes a yes or no answer and determines ' +\
                      'whether or not this tool produces an annotated spectrum. Default=yes')

parser.add_option( '--mz-threshold', action='store' , type = 'float' , dest = 'threshold' ,\
                       default = .25 , help = 'Input a threshold between 0 and 1. Default=.25' )

parser.add_option( '--neutral-losses', action='store' , type = 'string' , dest = 'nl' ,\
                       default ='none', help = 'Input all or none; depending on the ' +\
                       'desired output. None outputs the b and y ions and all outputs' + \
                       ' sevral more. Default=none')

parser.add_option( '--isotopic-mass', action='store' , type = 'string' , dest = 'masstype',\
                       default = 'mono', help = 'input \"mono\" for monoisotopic mass and ' +\
                       'average for \"average\" mass. Default=mono')

parser.add_option( '--general-statistics' , action = 'store' , type = 'string' , dest =\
                       'genstats' , default = 'no', help = 'This option takes the the' +\
                       ' arguments yes or no and determines whether or not a general-statistics' +\
                       ' report is created.' + 'If \"yes\" is selected, then you must' +\
                       ' also provide an output directory. Default=no' )

parser.add_option( '--single-ion-statistics' , action = 'store' , type = 'int' ,\
                       dest = 'sistats' , default = 0, help = 'input a number from 1 to ' +\
                       '13*maxcharge. Corresponds to the particular ion you want ' +\
                       'plotted.1 is for noise vs signal. Default=0') 

parser.add_option( '--pair-ion-statistics', action = 'store', type = 'int',\
                       dest = 'pistats', default = 0, help = 'input a number from 1 to 78*' +\
                       'the max charge(Corresponds to a particular ion relationship).' +\
                       ' Default=0') 

parser.add_option('--feature ' , action = 'store' , type = 'string' , dest = 'feature' ,\
                      default  = 'no', help = 'This option takes a \"all\" or \"no\" ' +\
                      'argument, or an array of selected features ' + \
                      'directory. Default=no')

parser.add_option( '--test' , action = 'store' , type = 'string' , dest = 'test' ,\
                       default  = 'no', help = 'This option tests whether or not the' +\
                       ' annotator is working correctly, It prints out a report expressing' +\
                       ' the details as well as a pass or fail verdict. Default=no')

parser.add_option( '--visualize', action = 'store', type = 'string', dest = 'visualize',\
                       default = 'no', help = 'If yes, this file creates a visual interpretstion' + \
                       'of each annotated PSM(Stores in a directory called \"visualize\")' +\
                       '. Default=no')

parser.add_option('--baseline-denoising', action = 'store', type = 'int', dest = 'bld',\
                      default = 0, help = 'This option when an input value other than zero' +\
                      ' is specified(values can range from 0-99) will remove the peaks with' +\
                      ' intensity in the bottom "n" percent of the spectrum')

parser.add_option('--overwrite', action = 'store', type = 'string', dest = 'overwrite',\
                      default = 'no', help = 'The tool makes subdirectories inside the' +\
                      ' specified directory. If yes is input, if these directories already' +\
                      ' exist their contents will be overwritten. Otherwise the directories' +\
                      ' will be created. If no is input, the tool will only create the' +\
                      ' directories if they do not exist. Under this option it will exit if' +\
                      ' the directories do exist.')

(options, args) = parser.parse_args()

if options.threshold > 1 or options.threshold < 0:
    sys.stderr.write('ERROR: The threshold must be between 1 and 0\n')
    sys.exit()

if options.nl != 'none' and options.nl  != 'all':
    sys.stderr.write('ERROR: The only types of neutral loss covered are none and all.\n')
    sys.exit()

if options.masstype != 'mono' and options.masstype != 'average':
    sys.stderr.write('ERROR: This tool only does monoisotopic and average mass\n')
    sys.exit()

if options.genstats != 'no' and options.annotate != 'yes' or \
        options.sistats != 0 and options.annotate != 'yes' or \
        options.pistats != 0 and options.annotate != 'yes':
    sys.stderr.write('ERROR: In order to do statistics of any kind the annotate' + 
                     ' option must be set to \"yes\"\n')
    sys.exit()

if options.test == 'yes' and options.annotate != 'yes':
    sys.stderr.write('ERROR: In order to test the annotator, the annotate option must be ' +\
                      'set to \"yes\"\n')
    sys.exit()

if len(sys.argv) < 3:
    sys.stderr.write('USAGE: ' + parser.usage + '\n')
    sys.stderr.write('type option "--help" for more info\n\n')
    sys.exit()

if options.bld < 0 or options.bld > 99:
    sys.stderr.write('ERROR: The baseline denoising option anly takes values from 0-99')
    sys.exit()


########
# you still have to add the conditions that will insure that the options will work correctly 
#########

#####################
##### CONSTANTS #####
#####################
MASS_NEUTRON = 1.00866
if options.masstype == "mono":
    massnh3 = 17.02655 
    massh2o = 18.01056
    massH = 1.0078246
    massCO = 27.9949 
else:
    massnh3 = 17.03056 
    massh2o = 18.0153
    massH = 1.00794
    massCO = 28.0101
floatError = .001
threshold = options.threshold
aminoAcids = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
averageMasses={ 'A':71.0788, 'B':114.5962, 'C':103.1388 + 57.0000,
                'D':115.0886, 'E':129.1155, 'F':147.1766, 'G':57.0519,
                'H':137.1411, 'I':113.1594, 'K':128.1741, 'L':113.1594,
                'M':131.1926, 'N':114.1038, 'O':114.1472, 'P':97.1167,
                'Q':128.1307, 'R':156.1875, 'S':87.0782, 'T':101.1051,
                'U':150.0388, 'V':99.1326, 'W':186.2132, 'X':113.1594,
                'Y':163.1760, 'Z':128.6231}
monoisotopicMasses={ 'A':71.03711, 'B':114.53494, 'C':103.00919 + 57.00000,
                     'D':115.02694, 'E':129.04259, 'F':147.06841, 'G':57.02146,
                     'H':137.05891, 'I':113.08406, 'K':128.09496, 'L':113.08406,
                     'M':131.04049, 'N':114.04293, 'O':114.07931, 'P':97.05276,
                     'Q':128.05858, 'R':156.10111, 'S':87.03203, 'T':101.04768,
                     'U':150.04344, 'V':99.06841, 'W':186.07931, 'X':113.08406,
                     'Y':163.06333, 'Z':128.55059}


############################
##### GLOBAL VARIALBES #####
############################
directory = sys.argv[len(sys.argv)-1]
inputFileName = sys.argv[len(sys.argv)-2]
inputfile = open(inputFileName,'r')
isFound = False
amIHere = False
sequences = []
prefixArrays = []
suffixArrays = []
arch = subprocess.Popen(['uname -m'], shell = True, stdout = subprocess.PIPE)\
    .communicate()[0].strip()
if options.feature == 'yes':
    maxCharge = int(subprocess.Popen(['grep CHARGE ' +  inputFileName +\
                                      ' | sort | tail -n 1 | cut -f 2 -d "+"'],\
                                     shell = True, stdout = subprocess.PIPE).communicate()[0].strip())

if options.annotate == 'yes':
    if options.overwrite == 'no':
        if os.path.exists(directory + '/annotated-spectrum.mgf2'):
            sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
            sys.exit()
    annotatorfile = open(directory + '/annotated-spectrum.mgf2','w')

if options.test ==  "yes":
    if options.overwrite == 'no':
        if os.path.exists(directory + '/test_report.txt'):
            sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
            sys.exit()
    testFile = open(directory + '/test_report.txt','w')
    failed = False
    cruxValue = ""
    f2value = ''
    annotatorFinds = []
    cruxFinds = []
    
    
if options.genstats == "yes":
    if not os.path.isdir(directory + '/stats/'):
        os.mkdir(directory + '/stats/')
    if not os.path.isdir(directory + '/stats/plots/'):
        os.mkdir(directory + '/stats/plots/')
    import matplotlib.pyplot as pyplot
    numberOfPSMs = 0
    peptideCounter = np.array([])
    peakNumber = np.array([])
    peakCounter = 0 
    chargeNumber = np.array([])
    peptideLength = np.array([])
    firstPeaks = np.array([])
    firstPeak = 0
    lastPeaks = np.array([])
    lastPeak = 0

if options.sistats != 0:
    if not os.path.isdir(directory + '/stats/'):
        os.mkdir(directory + '/stats/')
    if not os.path.isdir(directory + '/stats/plots/'):
        os.mkdir(directory + '/stats/plots/')
    if not os.path.isdir(directory + '/stats/plots/ion_types/'):
        os.mkdir(directory + '/stats/plots/ion_types/')
    import matplotlib.pyplot as pyplot
    import math
    import functools
    if options.sistats == 1:
        noise  = np.array([])
        signal  = np.array([])
    if options.sistats > 1:
        plotNumber = options.sistats - 1
        signal = []

if options.pistats !=0:
    if not os.path.isdir(directory + '/stats/'):
        os.mkdir(directory + '/stats/')
    if not os.path.isdir(directory + '/stats/plots/'):
        os.mkdir(directory + '/stats/plots/')
    if not os.path.isdir(directory + '/stats/plots/ion_types/'):
        os.mkdir(directory + '/stats/plots/ion_types/')
    if not os.path.isdir(directory + '/stats/plots/ion_types/correlations'):
        os.mkdir(directory + '/stats/plots/ion_types/correlations')
    import matplotlib.pyplot as pyplot
    familyArray = []
    iontype1 = np.array([])
    iontype2 = np.array([])
    plotNumber = options.pistats
    piCharges = np.array([])
    iontag1, iontag2 = generatePIStatsTag(plotNumber)
    counterX = 0
    counterY = 0

if options.feature != 'no':
    if not os.path.isdir(directory + '/learning/'):
        os.mkdir(directory + '/learning/')
        featureFile = open(directory + '/learning/features.txt', 'w')
        if options.annotate == 'yes':
            labelFile = open(directory + '/learning/labels.txt', 'w')
        else:
            labelFile = open(directory + '/learning/labels_null.txt', 'w')
    else:
        if options.overwrite == 'no':
            if os.path.exists(directory + '/learning/features.txt'):
                sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
                sys.exit()
    featureFile = open(directory + '/learning/features.txt','w')
    masses = np.array([])
    intensities = np.array([])
    if options.annotate == 'yes':
        labelFile = open(directory + '/learning/labels.txt', 'w')
    else:
        labelFile = open(directory + '/learning/labels_null.txt', 'w')

if options.visualize == 'yes':
    if options.overwrite == 'yes':
        if not os.path.isdir(directory + '/visualization/'):
            os.mkdir(directory + '/visualization/')
    if options.overwrite == 'no':
        if not os.path.isdir(directory + '/visualization/'):
            os.mkdir(directory + '/visualization/')
        else:
            sys.stderr.write('ERROR: Directory already exists, exiting to avoid over-writing!!!!\n')
            sys.exit()
    visualFile = open(directory + '/visualization/spectacle-input.txt', 'w')
    visualCounter = 1
    visualTag = [100, '']

if options.bld > 0:
    if not os.path.isdir(directory + '/baseline/'):
        os.mkdir(directory + '/baseline/')
        denoisedFile = open(directory + '/baseline/baseline-denoised-spectrum-' +\
                                str(options.bld) + '.mgf2','w')
    elif options.overwrite == 'no':
        if os.path.exists(directory + '/baseline/baseline-denoised-spectrum-' +\
                            str(options.bld) + '.mgf2'):
            sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
            sys.exit()
        else:
            denoisedFile = open(directory + '/baseline/baseline-denoised-spectrum-' +\
                                    str(options.bld) + '.mgf2','w')
    else:
        denoisedFile = open(directory + '/baseline/baseline-denoised-spectrum-' +\
                                str(options.bld) + '.mgf2','w')
    blLines = []

for line in inputfile:
    line = line.rstrip()
    if options.annotate =='yes':
        annotatorfile.write(line)
    if options.visualize == 'yes':
            visualTag = [100, '']
#########################################    
    if line.startswith("BEGIN") == True:
        if options.bld > 0:
            denoisedFile.write(line + '\n')
        if options.test == "yes":
            testFile.write(line + '\n')
            cruxValues = []
        if options.genstats == "yes":
            numberOfPSMs = numberOfPSMs + 1
        if options.feature != 'no':
            PSMarray = []
        if options.annotate == 'yes':
            prefixArrays = []
            suffixArrays = []
            stringArrays = []  
#########################################
    elif line.startswith("END") == True:
        if options.bld > 0:
            blDenoise(blLines, denoisedFile)
            denoisedFile.write(line + '\n\n')
            blLines = []
	if options.annotate == 'yes':
            annotateEnding(sequences, annotatorfile, stringArrays)
        if options.test == "yes":
            testFile.write(line + "\n\n")
        if options.genstats == "yes":
            peptideCounter = np.append(peptideCounter, len(sequences))
            peakNumber = np.append(peakNumber, peakCounter)
            firstPeaks = np.append(firstPeaks, firstPeak)
            lastPeaks = np.append(lastPeaks, lastPeak)
        if options.feature != 'no':
            print 'peptide start:'
            extractFeatures(PSMarray, float(pepmass[0]))
            PSMarray = []
            print 'peptide end:'
        if options.visualize == 'yes':
            visualCounter = visualCounter + 1
            visualFile.write('\n')
        if options.pistats > 0:
            iontype1, iontype2 = appendIonTypes(iontag1, iontag2, iontype1, iontype2,
                                                ionMatrix1, ionMatrix2, piCharges[0]) 
            piCharges = []
#########################################
    elif line.startswith("SEQ") == True:
        if options.bld > 0:
            denoisedFile.write(line + '\n')
        sequences = line.rsplit('=')[1].rsplit(',')
        stuffIndex = 0
        if options.visualize == 'yes':
            visualSequence = line.rsplit('=')[1]
        if options.genstats == 'yes':
            peptideCounter = np.append(peptideCounter , len(sequences))
        for stuff in sequences:
            if options.annotate == 'yes':
                stringArrays.append([])
                for char in stuff:
                    stringArrays[stuffIndex].append('')
                reverseSequence = stuff[::-1]
                chosenPrefixArray = calculateMass(stuff ,
                                                         "prefix" , options.masstype )
                chosenSuffixArray = calculateMass(reverseSequence ,
                                                         "suffix" , options.masstype )
                prefixArrays.append(chosenPrefixArray)
                suffixArrays.append(chosenSuffixArray)
            if options.genstats == 'yes':
                peptideLength = np.append(peptideLength, len(stuff))            
            if options.test == "yes":
                testFile.write(stuff)
                testFile.write("\n")
            stuffIndex = stuffIndex + 1
#########################################
    elif line.startswith("CHARGE") == True:
        if options.bld > 0:
            denoisedFile.write(line + '\n')
        theCharges = line.rsplit('=')[1].rsplit(',')
        chargeArray = []
        for eachCharge in theCharges:
            eachCharge = eachCharge.lstrip('+')
            eachCharge = eachCharge.lstrip('-') 
            chargeArray.append(range(1,(int(eachCharge)+1)))
            if options.genstats == "yes":
                chargeNumber = np.append(chargeNumber , int(eachCharge))
            if options.pistats > 0:
                piCharges = np.append(piCharges , int(eachCharge))
        if options.test == "yes":
                cruxValues = getCruxPrediction(cruxValues, sequences, theCharges)
        if options.pistats > 0:
            if options.pistats > 0:
                # Does not work for multiple peptides per PSM
                ionMatrix1 = np.zeros((len(sequences[0]) - 1, piCharges[0]))
                ionMatrix2 = np.zeros((len(sequences[0]) - 1, piCharges[0]))
#########################################
    elif line.startswith("PEPMASS") == True:
        if options.bld > 0:
            denoisedFile.write(line + '\n')
        if options.genstats == 'yes':
            amIFirst = True
        if options.feature != 'no':
            pepmass = np.array(line.rsplit('=')[1].rsplit(','))
        if options.visualize == 'yes':
            visualFile.write('>spec' + str(visualCounter) + '\t' + line.rsplit('=')[1] +\
                                 '\t'+ visualSequence + '\n')
#########################################
    else:
        lines = line.rsplit(' ')
        if len(lines) > 1:
            if options.feature == 'yes':
                PSMarray.append(line)
                masses = np.append(masses, float(lines[0]))
                intensities = np.append(intensities, float(lines[1]))
            if options.visualize == 'yes':
                visualFile.write(lines[0] + '\t' + lines[1])
            if options. bld > 0:
                blLines.append([float(lines[0]), float(lines[1])])
            if options.test == "yes":
                cruxFinds = []
                annotatorFinds = []
                convertedValue = ""
                f2value = ""
                found = False
                for value in cruxValues:
                    if value[0] != '':
                        if abs(float(value[0]) - float(lines[0])) < (threshold +
                                                                     floatError):
                            cruxFinds.append(value[0])
                        if abs(float(value[0]) - float(lines[0])) < threshold:
                            cruxSym = symGen( found )
                            convertedValue = convertedValue + cruxSym +\
                                cruxConverter(value)
                            found = True
                cruxValue = line + convertedValue                    
            lines[1].rstrip()
            lines[0].rstrip()
            if options.genstats == 'yes':
                if amIFirst == True:
                    firstPeak = float(lines[0])
                    amIFirst = False
                peakCounter = peakCounter + 1
                lastPeak = float(lines[0])          
#########  THE SEQUENCE FOR-LOOP  #####################
            arrayIndex = 0
            amIHere = False
            wasFound = False
            for eachSequence in sequences:
                sequenceLength = len(eachSequence)
                for eachAcid in eachSequence:
                    if aminoAcids.find(eachAcid) < 0:
                        sequenceLength = sequenceLength - 1
                sym2 = peptideNumberGenerator(arrayIndex , sequences)
                for charges in chargeArray[arrayIndex]:
                    if options.annotate == 'yes':
                        annotateIon('b', '', '@', ('b-ch' + str(charges)) )
                        annotateIon('y', '', '@', ('y-ch' + str(charges))) 
                        if options.nl == "all":
                            annotateIon('b', '-h2o', '@',\
                                        ('b-ch' + str(charges) + '-h2o')) 
                            annotateIon('b', '-nh3', '@',\
                                        ('b-ch' + str(charges) + '-nh3')) 
                            annotateIon('b', '-h2o-h2o', '',\
                                        ('b-ch' + str(charges) + '-h2o-h2o')) 
                            annotateIon('b', '-h2o-nh3', '@',\
                                        ('b-ch' + str(charges) + '-h2o-nh3')) 
                            annotateIon('a', '', '',\
                                        ('a-ch' + str(charges) )) 
                            annotateIon('a', '-h2o', '',\
                                        ('a-ch' + str(charges) + '-h2o')) 
                            annotateIon('a', '-nh3', '',\
                                        ('a-ch' + str(charges) + '-nh3')) 
                            annotateIon('y', '-h2o', '@',\
                                        ('y-ch' + str(charges) + '-h2o')) 
                            annotateIon('y', '-nh3', '@',\
                                        ('y-ch' + str(charges) + '-nh3')) 
                            annotateIon('y', '-h2o-h2o', '',\
                                        ('y-ch' + str(charges) + '-h2o-h2o'))
                            annotateIon('y', '-h2o-nh3', '@',\
                                        ('y-ch' + str(charges) + '-h2o-nh3'))
                arrayIndex = arrayIndex + 1
            if wasFound == True:
                if options.feature == 'yes':
                    labelFile.write('+1\n')
                if options.sistats == 1:
                    signal = np.append(signal, float(lines[1]))
            else:
                if options.feature == 'yes':
                    labelFile.write('-1\n')
                if options.sistats == 1:
                    noise = np.append(noise, float(lines[1]))
            if options.test == "yes":
                testAnnotator()
            if options.visualize == 'yes':
                if visualTag[1] != '':
                    visualFile.write('\t' + visualTag[1] + '\t' + visualTag[2])
                visualFile.write('\n')
    if options.annotate == 'yes':
        annotatorfile.write("\n")
    
inputfile.close
if options.annotate == 'yes':
    annotatorfile.close
###### GENERAL-STATISTICS ######
if options.genstats == 'yes':
    if options.overwrite == 'no':
        if os.path.exists(directory + '/stats/general-statistics.txt'):
            sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
            sys.exit()
    genStatsFile = open(directory + '/stats/general-statistics.txt','w')
    genStatsFile.write("Number of PSMs in this file:  " + str(numberOfPSMs) + "\n\n")
    provideStatsFile(peptideCounter , 'Peptides per PSM:\n', genStatsFile)
    plotter(peptideCounter, 'Peptides per PSM', 'Peptide Number','Number', 'peptides_per_psm')
    provideStatsFile(peakNumber, 'Peaks per PSM:\n', genStatsFile)
    plotter(peakNumber, 'Number of Peaks per PSM', 'Peak Count', 'Number', 'number_of_peaks')
    provideStatsFile(chargeNumber, 'Charges:\n', genStatsFile)
    plotter(chargeNumber, 'Charge', 'Charge', 'Number', 'charge')
    provideStatsFile(peptideLength, 'Peptide Length:\n', genStatsFile)
    plotter(peptideLength, 'Peptide Length', 'Peptide Length', 'Number', 'peptide_length')
    provideStatsFile(firstPeaks, 'First Peak:\n', genStatsFile)
    plotter(firstPeaks, 'First Peak Position', 'M/Z', 'Number', 'first_peak')
    provideStatsFile(lastPeaks, 'Last Peak:\n', genStatsFile)
    plotter(lastPeaks, 'Last Peak Position', 'M/Z', 'Number', 'last_peak')
    genStatsFile.close()

if options.sistats != 0:
    if len(signal) != 0:
        statsArray = ['mean: ','median: ','standard deviation: ','minimum: ', 'maximum: ']
        siIndex = 0
        if options.sistats == 1:
            if options.overwrite == 'no':
                if os.path.exists(directory + '/stats/plots/ion_types' + '/signal_vs_noise.txt'):
                    sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
                    sys.exit()    
            siStatsFile = open(directory + '/stats/plots/ion_types' + '/signal_vs_noise.txt','w')
            provideStatsFile(noise, 'Noise:\n', siStatsFile)
            plotter(noise , 'NOISE' , '' , '' , '')
            provideStatsFile(signal, 'Signal:\n', siStatsFile)
            plotter(signal , 'SIGNAL vs NOISE' , 'Intensity' , 'Number' ,\
                        ('ion_types/signal_vs_noise' ))
            siStatsFile.close()

        if options.sistats > 1:
            #print 'Hello? '
            ID = generateStatsID(plotNumber)
            if options.overwrite == 'no':
                if os.path.exists(directory + '/stats/plots/ion_types/intensity_' + ID + '.txt'):
                    sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
                    sys.exit()
            siStatsFile = open(directory + '/stats/plots/ion_types/intensity_' + ID + '.txt','w')
            provideStatsFile(signal, ID + '\n', siStatsFile)
            plotter(signal , ID.upper() , 'Intensity' , 'Number' ,\
                        ('ion_types/intensity_' + ID) )
            siStatsFile.close()
    else:
        sys.stderr.write('ERROR: This ion was not found in the PSM\n')
        sys.exit()
if options.pistats > 0:
    if len(iontype1) != 0 and len(iontype2) !=0:
        if options.overwrite == 'no':
            if os.path.exists(directory + '/stats/plots/ion_types/correlations/' +\
                                       iontag1 + '_' + iontag2 + '.txt'):
                sys.stderr.write('ERROR: File already exists, exiting to avoid over-writing!!!!\n')
                sys.exit()
        piStatsFile = open(directory + '/stats/plots/ion_types/correlations/' +\
                               '/scatter_intensities_' + iontag1 + '_' + iontag2 + '.txt', 'w')
        correlationCoefficient = np.corrcoef(iontype1, iontype2)
        pyplot.scatter(iontype1,iontype2)
        pyplot.title(iontag1.upper() + ' versus ' + iontag2.upper())
        pyplot.xlabel('Intensity of ' + iontag1.upper()+ '\n' + '(' +\
                          str(correlationCoefficient[0][1]) + ')')
        pyplot.ylabel('Intensity of ' + iontag2.upper())
        pyplot.savefig( directory + '/stats/plots/ion_types/correlations' +\
                            '/scatter_intensities_' + iontag1 + '_' + iontag2)
        piStatsFile.write('CORRELATION COEFFICIENT: ' + str(correlationCoefficient[0][1]) + '\n')
        piStatsFile.write('Matching ION found: ' + str(counterX) + '\n')
        piStatsFile.write('Matching ION not found: ' + str(counterY) + '\n')
        piStatsFile.write('\nData:\n' + iontag1 + '\t' + iontag2 + '\n' )
        for x in range(0, len(iontype1)):
            piStatsFile.write(str(iontype1[x]) + '\t' + str(iontype2[x]) + '\n')
        piStatsFile.close()
###### FEATURE SELECTION ######
if options.feature != 'no':
    featureFile.close()
    labelFile.close()

if options.visualize == 'yes':
    subprocess.check_call(['perl', 'bin/spectacle.pl', directory +\
                               '/visualization/spectacle-input.txt','-html',\
                               directory + '/visualization/specrta.html', '-format', 'png',\
                               '-title'])
    subprocess.check_call(['mv spec*.png '+ directory + '/visualization/'], shell=True)
    

######## TEST RESULTS ########
if options.test == 'yes':
    if failed == False:
        testFile.write('\n\n\nThe Annotator PASSED!!!!  :)')
    else:
        testFile.write('\n\n\nThe Annotator FAILED!!!!  :(')

    testFile.close()
    





