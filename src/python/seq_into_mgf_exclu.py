import os
import csv 
import sys




def getIndices(filename, header):
    'Reads the header of the file and returns indices'
    scanIndex=-1
    chargeIndex=-1
    rankIndex=-1
    seqIndex=-1
    qvalueIndex=-1  
    header = header.split("\t")
    for i in range(len(header)):
        if header[i]=="scan": scanIndex=i
        elif header[i]=="charge": chargeIndex=i
        elif header[i]=="percolator q-value" && filename.find("percolator"):
            qvalueIndex=i
        elif header[i]=="q-ranker q-value" && filename.find("qranker"):
            qvalueIndex=i
        elif header[i]=="xcorr rank": rankIndex=i
        elif header[i]=="sequence": seqIndex=i

    return scanIndex,chargeIndex,rankIndex,seqIndex,qvalueIndex 



try:
    mgfFile = sys.argv[1]
    targetTxtFile = sys.argv[2]  #"../2009-09-29/mouse/P21-1-01/percolator.target.txt"
    outputdir = sys.argv[3]
    threshold = float(sys.argv[4])
    outputfile = outputdir+"/seq_output.mgf"
    targetBaseName = os.path.basename(targetTxtFile)
except:
    print "python seq_into_mgf.py <mgf_file> <target_txt_file> <output_dir>"
    sys.exit

targetTxtOpen =open(targetTxtFile, "r")
count = 0
    


indices = getIndices(targetBaseName, targetTxtOpen.readline())
scanIndex,chargeIndex,rankIndex,seqIndex,qvalueIndex = indices 


targetReader = csv.reader(targetTxtOpen, delimiter='\t')
matches = {}
# Iterate through the target file and store all sequences with thier 
# corresponding charge indexed through the scan number
for row in targetReader:
    rank = row[rankIndex]
    qvalue = row[qvalueIndex]
    if (rank == "1" and float(qvalue) < threshold):
        count+=1
        scan = row[scanIndex]
        charge = row[chargeIndex]
        seq = row[seqIndex]
        if (int(scan), charge) in matches:
            matches[(int(scan), charge)] += ","+seq
        else:
            matches[(int(scan), charge)] = seq
targetTxtOpen.close()        

    
mgfFileRead = open(mgfFile, "r")
outputWrite = open(outputfile, "w")
nextLine = mgfFileRead.readline()

# Iterate through mgf file file and check to see if the scan & charge number
# has a corresponding sequence and enter the peptide into the output file if
# it is
while(nextLine):
    if (nextLine.startswith("BEGIN IONS")):
        beforeString=nextLine
        nextLine = mgfFileRead.readline()
        spectraName = nextLine.split("=")[1]
        tempPair= spectraName.split("_")
        scanNumber = int(tempPair[1])
        charge = tempPair[2].strip('\n')
        if ((scanNumber, charge) in matches):
            resultPair = matches[(scanNumber, charge)]
            # write the title line and add the sequence to the next line
            outputWrite.write(beforeString+nextLine)
            nextLine = mgfFileRead.readline()
            nextLine = nextLine.strip("\n")+resultPair+"\n"
            outputWrite.write(nextLine)
            extra = resultPair.count(',')
            
            # if there are more than one sequence, we must write multiple
            # charge and pepmass so subsequent scripts can read it
            chargeLine = mgfFileRead.readline()
            pepMassLine = mgfFileRead.readline()
            charge = chargeLine.split('=')[1].strip('\n')
            pepMass = pepMassLine.split('=')[1].strip('\n')
            chargeString = "CHARGE="+charge
            pepMassString = "PEPMASS="+pepMass
            for i in range(0, extra):
                chargeString+=","+charge
                pepMassString+=","+pepMass
            chargeString +="\n"
            pepMassString += "\n"
            outputWrite.write(chargeString+pepMassString)
            nextLine = mgfFileRead.readline() 
            while (not nextLine.startswith("BEGIN IONS")):
                outputWrite.write(nextLine)
                nextLine = mgfFileRead.readline() 
    else:
        nextLine = mgfFileRead.readline()



    
