#!/opt/local/bin/python

"""
Need to start Eclipse from command line (using executable inside app) to get
ADF env vars working right, when in unix-like system.
"""

import kf
import numpy as np
from random import randint
from os import name

"""
The two variables below are everything that will be fetched from CPs and Bader atoms
respectively.
To change the order things are printed in in the output, change the order here.
To exclude a piece of information, change it's "True" to "False".

Note that eigensystem stuff for CPs has it's own variables below CPInfoList.
"""

CPInfoList = [
              {'name':['x','y','z'],
               'kf_name':'CP coordinates',
               'collect':True},
              
              {'name':['rho'],
               'kf_name':'CP density at',
               'collect':True},
              
              {'name':['grad ' + i for i in ['x','y','z']],
               'kf_name':'CP density gradient at',
               'collect':True},
              
              {'name':['hess ' + i for i in ['xx','xy','xz','yy','yz','zz']],
               'kf_name':'CP density Hessian at',
               'collect':True}
              ]

CPEigVals = {'name':['eig val ' + str(i) for i in range(1,4)],
             'collect':True, 
             'data':[]}

CPEigVecs = {'name':['eig vec ' + str(i) + '.' + str(j) for i in range(1,4) for j in range(1,4)], 
             'collect':True,
             'data':[]}

AtomInfoList = [
                {'name':['charge'],
                 'kf_name':'Bader atomic charges',
                 'collect':True},
                
                {'name':['spin dens'],
                 'kf_name':'Bader atomic spin densities',
                 'collect':True},
                
                {'name':['dipole moment ' + i for i in ['x','y','z']],
                 'kf_name':'Bader atomic dipole moment',
                 'collect':True},
                
                {'name':['sph. tensor quadrupole moment ' + i for i in ['Q20','Q21c','Q21s','Q22c','Q22s']],
                 'kf_name':'Bader atomic quadrupole moment',
                 'collect':True},
                
                {'name':['laplacian'],
                 'kf_name':'Bader Laplacian',
                 'collect':True},
                
                {'name':['kinetic energy'],
                 'kf_name':'Bader Ts',
                 'collect':True},
                
                {'name':['virial kinetic energy'],
                 'kf_name':'Bader Tsl',
                 'collect':True},
                
                {'name':['local corr. kinetic energy'],
                 'kf_name':'Bader Tc',
                 'collect':True},
                
                {'name':['virial factor energy'],
                 'kf_name':'Bader EVF',
                 'collect':True},
                
                {'name':['ELF'],
                 'kf_name':'Bader ELF',
                 'collect':True}
                ]

IsWindows = name == 'nt'

# Takes the path to a Tape21 file WITH BADER CP INFO
# Returns a list of triplets with XYZ coords of each CP
CurrentDir = ''
SubDir = '/Output' + str(randint(1,999999))
PathCoords = []
T21BaseName = ''
ScriptDir = ''

def GetCPCoords(T21FileName):
    T21 = kf.kffile(T21FileName)
    FullList = T21.read("Properties", "CP coordinates")
    NumCPs = T21.read("Properties","CP number of")
    T21.close()
    
    CoordList = []
    
    for i in range(NumCPs):
        tmp = np.empty([3])
        k = 0
        for j in range(i,NumCPs * 3,NumCPs):
            tmp[k] = FullList[j]
            k += 1
        CoordList.append(tmp)
    
    return CoordList

HessInd = [
       [0,1,2],
       [1,3,4],
       [2,4,5]]
def GetFullHess(Hess):
    FullHess = np.empty([3,3])
    for i in range(3):
        for j in range(3):
            FullHess[i,j] = Hess[HessInd[i][j]]
    return FullHess

def GetCPEigVals(T21FileName):
    from numpy import linalg as la
    T21 = kf.kffile(T21FileName)
    Hess = T21.read("Properties", "CP density Hessian at")
    NumCPs = T21.read("Properties","CP number of")
    T21.close()
    
    EigValList = []
    
    for i in range(NumCPs):
        tmp = np.empty([6])
        k = 0
        for j in range(i, NumCPs * 6, NumCPs):
            tmp[k] = Hess[j]
            k += 1
        EVal, tmp1 = la.eigh(GetFullHess(tmp))
        EigValList.append(EVal)
        
    return EigValList

def GetCPInfo(T21FileName):
    from numpy import linalg as la
    global CPInfoList, CPEigVals, CPEigVecs
    T21 = kf.kffile(T21FileName)
    NumCPs = T21.read("Properties","CP number of")
    
    for i in [CPEigVals,CPEigVecs]:
        if i['collect']: i['data'] = []
    
    
    for InfoType in CPInfoList:
        if InfoType['collect']:
            tmp = T21.read("Properties",InfoType['kf_name'])
            InfoType['data'] = []
            for i in range(NumCPs):
                InfoType['data'].append([])
                for j in range(i, NumCPs * len(InfoType['name']), NumCPs):
                    InfoType['data'][i].append(tmp[j])
                if 'Hessian' in InfoType['kf_name']:
                    EVal, EVec = la.eigh(GetFullHess(InfoType['data'][i]))
                    CPEigVals['data'].append(EVal)
                    CPEigVecs['data'].append([EVec[y,x] for x in range(3) for y in range(3)])
    T21.close()
    return

def GetAtomInfo(T21FileName):
    global AtomInfoList
    T21 = kf.kffile(T21FileName)
    NumAtoms = len(T21.read("Properties","Bader atomic charges"))
    
    for InfoType in AtomInfoList:
        if InfoType['collect']:
            tmp = T21.read("Properties",InfoType['kf_name'])
            InfoType['data'] = []
            for i in range(NumAtoms):
                InfoType['data'].append([])
                for j in range(i, NumAtoms * len(InfoType['name']), NumAtoms):
                    InfoType['data'][i].append(tmp[j])
    T21.close()
    return
    
    
    
# Takes the path to a Tape21 file WITH BADER CP INFO
# Returns a list of CP type numbers
#The "type" numbers go as {1,2,3,4} = {n,c,b,r}
def GetCPTypes(T21FileName):
    TypeList = []
    T21 = kf.kffile(T21FileName)
    Types = T21.stringData("Properties","CP code number for (Rank,Signatu")
    Types = Types.split('\n')
    Types = Types[3:len(Types)-1]
    for i in Types:
        for j in i.split():
            TypeList.append(int(float(j)))
    
    T21.close()
    return TypeList

def GetCPsFromCoords(CPNum,NewCoords,NewTypes,OldCoords,OldTypes):
    NewCPNum = -1
    
    MinDist = 1e100
    MinInd = -1
    for i in range(len(NewCoords)):
        if NewTypes[i] == OldTypes[CPNum-1]:
            TmpDist = np.linalg.norm(np.subtract(NewCoords[i],OldCoords[CPNum-1]))
            if TmpDist < MinDist:
                MinDist = TmpDist
                MinInd = i
    if MinInd >= 0:
        NewCPNum = MinInd + 1
    
    return NewCPNum
    
    
CPTypeInd = ['N','C','B','R']
CPRankInd = [-3,3,-1,1]
def main():
    import ast
    import os
    import sys
    import subprocess
    global T21BaseName, ScriptDir, CurrentDir, CPInfoList, AtomInfoList
    
#     InputFilePath = "C:\\Users\\Haiiro\\Dropbox\\EclipseWin\\Python\\CPGobbler\\input.txt"
#     InputFilePath = "/Users/Haiiro/Safe/CP_Gobbler/input.txt"
    InputFilePath = "/Users/Haiiro/ADFdata/input.txt"
    
    
    
    if len(sys.argv) > 1:
        InputFilePath = sys.argv[1]
        

    InputFile = open(InputFilePath, 'r')
    InputFileContents = InputFile.readlines()
    InputFile.close()
    
    ScriptDir = os.path.dirname(os.path.realpath(__file__))
    
    InputDir = os.path.dirname(os.path.realpath(InputFilePath))
    
    FileNum = 0
    
    CurrentDir = os.path.dirname(InputFilePath)
    OutBRCCSVName = CurrentDir + SubDir + "/CPs.csv"
    OutNCSVName = CurrentDir + SubDir + "/BaderAtoms.csv"
    
    if IsWindows:
        subprocess.Popen('mkdir "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()
    else:
        subprocess.Popen('mkdir -p "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()
    
    
    BRCHeaders = "Filename,InCP#,CP#,CP Rank,CP Type,"
    NHeaders = "Filename,InCP#,CP#,x,y,z,"
    
    FirstCP = True
    FirstBaderAtom = True
    
    for aLine in InputFileContents:
        FileNum += 1
        if len(aLine) > 0 and aLine[0] == '#':
            continue
        print "Running file " + str(FileNum) + " of " + str(len(InputFileContents))
        
        aLine = aLine.split('\t')
        if len(aLine) > 1:
            InputT21FileName = aLine[0]
            TargetCPs = ast.literal_eval(aLine[1])
            AllCPs = (TargetCPs[0] < 0)
            if AllCPs: print "Getting all CP information"
            
            T21BaseName = os.path.basename(InputT21FileName)
                
            CurrentDir = os.path.dirname(InputT21FileName)
                
            if not os.path.exists(CurrentDir):
                if IsWindows:
                    InputT21FileName = InputDir + '\\' + T21BaseName
                else:
                    InputT21FileName = InputDir + '/' + T21BaseName
            
            if '?' in InputT21FileName or '*' in InputT21FileName:
                if IsWindows:
                    T21FileNames = subprocess.Popen("dir /b " + InputT21FileName, shell=True, stdout=subprocess.PIPE).stdout.read()
                else:
                    T21FileNames = subprocess.Popen('ls ' + InputT21FileName, shell=True, stdout=subprocess.PIPE).stdout.read()
                T21FileNames = T21FileNames.split('\n')
            else:
                T21FileNames = [InputT21FileName]
                
            OrigCPCoords = []
            OrigCPTypes = []
            
            for T21FileName in T21FileNames:
                if len(T21FileName) > 2:
                    if IsWindows:
                        T21FileName = (CurrentDir + '\\' + T21FileName).replace('\r','')
                    if True:#os.path.exists(T21FileName):
                        T21BaseName = os.path.basename(T21FileName)
                        
                        if '.' in T21BaseName:
                            T21BaseName = T21BaseName.rpartition('.')[0]
                            
                            
                        T21 = kf.kffile(T21FileName)
                        ExitStatus = T21.stringData("General", "termination status")
                        if not 'NORMAL TERMINATION' in ExitStatus:
                            print "T21 file reports abnormal termination: %s. skipping..." % ExitStatus
                            continue
                        NumCPs = T21.read("Properties", "CP number of")
                        if NumCPs is None or NumCPs <= 0:
                            print "No CPs found, skipping..."
                            continue
                        T21.close()
                            
                        GetCPInfo(T21FileName)
                        GetAtomInfo(T21FileName)
                    
                        CPCoords = GetCPCoords(T21FileName)
                        CPTypes = GetCPTypes(T21FileName)
#                         CPEigVals = GetCPEigVals(T21FileName)
                        T21 = kf.kffile(T21FileName)
#                         CPDens = T21.read("Properties", "CP density at")
#                         BaderCharges = T21.read("Properties", "Bader atomic charges")
                        if AllCPs:
                            NumCPs = T21.read("Properties", "CP number of")
                            TargetCPs = [i+1 for i in range(NumCPs)]
                            
                        
                        AtomOrderIndex = T21.read("Geometry","atom order index")
                        NumAtoms = len(AtomOrderIndex) / 2
                        T21.close()
                        
                        if len(OrigCPCoords) == 0: OrigCPCoords = CPCoords
                        if len(OrigCPTypes) == 0: OrigCPTypes = CPTypes
                        
                        for CPNum in TargetCPs:
                            CurCPNum = GetCPsFromCoords(CPNum, CPCoords, CPTypes, OrigCPCoords, OrigCPTypes) if not AllCPs else CPNum
                            CurCPType = CPTypes[CurCPNum-1]
                            
                            print "Getting CP " + str(CurCPNum) + " from " + T21BaseName
                            
                            if FirstCP:
                                FirstCP = False
                                if IsWindows:
                                    subprocess.Popen('mkdir "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()
                                else:
                                    subprocess.Popen('mkdir -p "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()
                                CPCSV = open(OutBRCCSVName, 'a')
                                CPCSV.write(BRCHeaders)
                                for i in CPInfoList:
                                    if i['collect']:
                                        for j in i['name']:
                                            CPCSV.write(j + ',')
                                        if 'Hessian' in i['kf_name']:
                                            if CPEigVals['collect']:
                                                for j in CPEigVals['name']:
                                                    CPCSV.write(j + ',')
                                            if CPEigVecs['collect']:
                                                for j in CPEigVecs['name']:
                                                    CPCSV.write(j + ',')
                                CPCSV.write('\n')
                            for i in [T21FileName,CPNum,CurCPNum,CPRankInd[CPTypes[CurCPNum-1]-1],CPTypeInd[CPTypes[CurCPNum-1]-1]]:
                                CPCSV.write(str(i) + ",")
                            for i in CPInfoList:
                                if i['collect']:
                                    for j in i['data'][CurCPNum-1]:
                                        CPCSV.write(str(j) + ",")
                                    if 'Hessian' in i['kf_name']:
                                            if CPEigVals['collect']:
                                                for j in CPEigVals['data'][CurCPNum-1]:
                                                    CPCSV.write(str(j) + ',')
                                            if CPEigVecs['collect']:
                                                for j in CPEigVecs['data'][CurCPNum-1]:
                                                    CPCSV.write(str(j) + ',')
#                                 for i in CPCoords[CurCPNum-1]:
#                                     CPCSV.write(str(i) + ",")
#                                 CPCSV.write(str(CPDens[CurCPNum-1]) + ",")
#                                 for i in CPEigVals[CurCPNum-1]:
#                                     CPCSV.write(str(i) + ",")
                            CPCSV.write("\n")
                                
                            if CurCPType == 1:
                                if FirstBaderAtom:
                                    FirstBaderAtom = False
                                    if IsWindows:
                                        subprocess.Popen('mkdir "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()
                                    else:
                                        subprocess.Popen('mkdir -p "' + CurrentDir + SubDir + '"', shell=True, stdout=subprocess.PIPE).stdout.read()                                    
                                    BACSV = open(OutNCSVName, 'a')
                                    BACSV.write(NHeaders)
                                    for i in AtomInfoList:
                                        if i['collect']:
                                            for j in i['name']:
                                                BACSV.write(j + ',')
                                    BACSV.write("\n")
#                                 Fetch the input atom number of the current nuclear CP.
#                                 I guess this isn't necessary, so it's commented out.
#                                 AtomNum = AtomOrderIndex[CurCPNum-1 + NumAtoms]
                                AtomNum = CurCPNum
                                for i in [T21FileName,CPNum,AtomNum]:
                                    BACSV.write(str(i) + ",")
                                for i in CPCoords[CurCPNum-1]:
                                    BACSV.write(str(i) + ",")
                                for i in AtomInfoList:
                                    if i['collect']:
                                        for j in i['data'][AtomNum-1]:
                                            BACSV.write(str(j) + ",")
#                                 BACSV.write(str(CPDens[CurCPNum-1]) + ",")
#                                 BACSV.write(str(BaderCharges[CurCPNum-1]) + ",")
                                BACSV.write("\n")
        
        
    if not FirstCP:
        CPCSV.close()
    if not FirstBaderAtom:
        BACSV.close()
    
    print "Finished"
    return


main()