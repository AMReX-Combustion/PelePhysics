import numpy as np
import os

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False	

def readFile(folder,filename):
    fullFilename = folder + '/' +  filename
    f=open(fullFilename, 'r')
    Array = np.array(f.read().splitlines(),dtype=float)
    f.close()
    return Array

def readMultiColFile(filename,headerSize=0):
    f = open(filename,'r')
    lines=f.readlines()
    # Figure out number of column
    A=lines[headerSize].split()
    # Allocate array
    Array = np.zeros((len(lines)-headerSize,len(A)))
    counter = 0
    for line in lines[headerSize:]:
        A=line.split()
        for i in range(len(A)):
            Array[counter,i] = float(A[i])
        counter=counter+1
    f.close()
    return Array

def readKernelFile(folder,filename):
    fullFilename = folder + '/' + filename
    Array = np.loadtxt(fullFilename)
    return Array

def readPDFFile(folder,filename):
    fullFilename = folder + '/' + filename
    Array = np.loadtxt(fullFilename)
    return Array

def readOFfile(folder, filename,nHeader,nLength):
    Array = np.zeros(nLength)
    fullFilename = folder + '/' + filename
    f=open(fullFilename, 'r')
    lines = f.readlines()[nHeader:nLength+nHeader]
    f.close()
    it=0
    for line in lines:
        Array[it] = float(line)
        it = it + 1
    return Array

def readOFfile_decomp(folder,filename,nProc_decomp,nHeader,nLength,nCells_decomp):
    Array = np.zeros(nLength)
    it = 0
    for iproc in range(nProc_decomp):
        fullFilename = folder + '/processor' + str(iproc) + '/' + filename
        f=open(fullFilename, 'r')
        lines = f.readlines()[nHeader:int(nCells_decomp[iproc])+nHeader]
        f.close()
        for line in lines:
            Array[it] = float(line)
            it = it + 1
    return Array

def readOFVec(folder,filename,nHeader,nLength,indexRemove):
    Array1 = np.zeros(nLength)
    Array2 = np.zeros(nLength)
    Array3 = np.zeros(nLength)
    # Check that the field is not a uniform field
    fullFilename = folder + '/' + filename
    f=open(fullFilename, 'r')
    for i in range(20):
        line=f.readline()
    f.close()
    lineList=line.split()
    if len(lineList)==5 and lineList[1]=='uniform':
        # Uniform field
        val1 = float(lineList[2][1:])
        val2 = float(lineList[3])
        val3 = float(lineList[4][:-2])
        for it in range(nLength):
            Array1[it] = val1
            Array2[it] = val2
            Array3[it] = val3
    else:
        # Not a uniform field
        f=open(fullFilename, 'r')
        lines = f.readlines()[nHeader:nLength+nHeader]
        f.close()
        it=0
        for line in lines:
            lineRead = (line[1:-2]).split()
            Array1[it] = float(lineRead[0])
            Array2[it] = float(lineRead[1])
            Array3[it] = float(lineRead[2])
            it = it + 1

    Array1[indexRemove] = np.nan
    Array2[indexRemove] = np.nan
    Array3[indexRemove] = np.nan
    Array1 = Array1[np.logical_not(np.isnan(Array1))]
    Array2 = Array2[np.logical_not(np.isnan(Array2))]
    Array3 = Array3[np.logical_not(np.isnan(Array3))]

    return Array1,Array2,Array3

def readOFVec_decomp(folder,filename,nProc_decomp,nHeader,nLength,nCells_decomp,indexRemove):
    Array1 = np.zeros(nLength)
    Array2 = np.zeros(nLength)
    Array3 = np.zeros(nLength)
    it = 0
    for iproc in range(nProc_decomp):
        fullFilename = folder + '/processor' + str(iproc) + '/' + filename
        # Check that the field is not a uniform field
        f=open(fullFilename, 'r')
        for i in range(20):
            line=f.readline()
        f.close()
        lineList=line.split()
        if len(lineList)==5 and lineList[1]=='uniform':
            # Uniform field
            val1 = float(lineList[2][1:])
            val2 = float(lineList[3])
            val3 = float(lineList[4][:-2])
            for i in range(nLength):
                Array1[it] = val1
                Array2[it] = val2
                Array3[it] = val3
                it = it + 1
        else:
            # Not a uniform field
            f=open(fullFilename, 'r')
            lines = f.readlines()[nHeader:int(nCells_decomp[iproc])+nHeader]
            f.close()
            for line in lines:
                lineRead = (line[1:-2]).split()
                Array1[it] = float(lineRead[0])
                Array2[it] = float(lineRead[1])
                Array3[it] = float(lineRead[2])
                it = it + 1

    Array1[indexRemove] = np.nan
    Array2[indexRemove] = np.nan
    Array3[indexRemove] = np.nan
    Array1 = Array1[np.logical_not(np.isnan(Array1))]
    Array2 = Array2[np.logical_not(np.isnan(Array2))]
    Array3 = Array3[np.logical_not(np.isnan(Array3))]
    return Array1,Array2,Array3

def readOFScal(folder,filename,nHeader,nLength,indexRemove):
    Array = np.zeros(nLength)
    # Check that the field is not a uniform field
    fullFilename = folder + '/' + filename
    f=open(fullFilename, 'r')
    for i in range(20):
        line=f.readline()
    f.close()
    lineList=line.split()
    if len(lineList)==3 and lineList[1]=='uniform':
        # Uniform field
        val = float(lineList[2][:-1])
        for i in range(nLength):
            Array[i] = val
    else:
        # Not a uniform field
        f=open(fullFilename, 'r')
        lines = f.readlines()[nHeader:nLength+nHeader]
        f.close()
        it=0
        for line in lines:
            Array[it] = float(line)
            it = it + 1
    Array[indexRemove] = np.nan
    Array = Array[np.logical_not(np.isnan(Array))]
    return Array

def readOFScal_decomp(folder,filename,nProc_decomp,nHeader,nLength,nCells_decomp,indexRemove):
    Array = np.zeros(nLength)
    it = 0
    for iproc in range(nProc_decomp):
        fullFilename = folder + '/processor' + str(iproc) + '/' + filename
        # Check that the field is not a uniform field
        f=open(fullFilename, 'r')
        for i in range(20):
            line=f.readline()
        f.close()
        lineList=line.split()
        if len(lineList)==3 and lineList[1]=='uniform':
            # Uniform field
            val = float(lineList[2][:-1])
            for i in range(int(nCells_decomp[iproc])):
                Array[it] = val
                it = it + 1
        else:
            # Not a uniform field
            f=open(fullFilename, 'r')
            lines = f.readlines()[nHeader:int(nCells_decomp[iproc])+nHeader]
            f.close()
            for line in lines:
                Array[it] = float(line)
                it = it + 1

    Array[indexRemove] = np.nan
    Array = Array[np.logical_not(np.isnan(Array))]
    return Array


def readCylindCut(folder):
    folder = folder + '/cylindCut'
    listLog=os.listdir(folder)
    timeList=[]
    for name in listLog:
        if is_number(name):
            timeList.append(name)
    Output =  sorted(timeList, key = lambda x:float(x))
    timeList=list(Output)
    # Read the file in each directory
    data = []
    for time in timeList:
        if os.path.isfile(folder + '/' + time + '/surfaceFieldValue_0.00059988002.dat'):
            filename = folder + '/' + time + '/surfaceFieldValue_0.00059988002.dat'
        else:
            filename = folder + '/' + time + '/surfaceFieldValue.dat'
        data.append(readMultiColFile(filename,headerSize=4))
    # Concatenate the data properly
    dataFin = []
    for itime,time in enumerate(timeList):
        if itime<(len(timeList)-1):
            nextTime = float(timeList[itime+1])
            for idat in range(len(data[itime][:,0])):
                if data[itime][idat,0] > nextTime:
                    break
                else:
                    dataFin.append(data[itime][idat,:])
        else:
            for idat in range(len(data[itime][:,0])):
                dataFin.append(data[itime][idat,:])
    return np.array(dataFin)

    

def writeFile(filename,Array):
    f=open(filename, 'w+')
    if isinstance(Array,float):
        f.write( str(Array) + '\n'  )
    else:
        for Entry in Array:
            f.write( str(Entry) + '\n'  )
    f.close()

def writeMultiColFile(filename,header,Array):
    f=open(filename, 'w+')
    for i in range(len(Array[:,0])):
        for j in range(len(Array[0,:])):
            f.write( str(Array[i][j]) + ' '  )
        f.write(  '\n'  )
    f.close()
