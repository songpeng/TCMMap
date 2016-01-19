#!/usr/bin/python

# Analysis the results of AL-Mode-1
# Songpeng Zu
# 2015-01-08

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import rdMolDescriptors

class ALanalog:
    def __int__(self):
        self.name = ""
        self.smilesstr = ""
        self.IC50 = 0.0
        self.logP = 0.0
        self.TPSA = 0.0
    def getSMILEs(self,strng):
        self.smilesstr = strng
    def getMolName(self,namestr):
        self.name = namestr
    def getIC50(self,value):
        self.IC50 = value
    def getlogP(self):
        if not not self.smilesstr:
            self.logP = Descriptors.MolLogP(Chem.MolFromSmiles(self.smilesstr))
        else:
            self.logP = 0.0
            print "Current ALanalog's SMILES is EMPTY!\n"
    def getTPSA(self):
        if not not self.smilesstr:
            self.TPSA = Descriptors.TPSA(Chem.MolFromSmiles(self.smilesstr))
        else:
            self.TPSA = 0.0
            print "Current ALanalog's SMILES is EMPTY!\n"
    def printAnalogStr(self):
        if not not self.smilesstr:
            tmpstr = self.name + "\t" + self.smilesstr + "\t" + str(self.IC50) + "\t" +\
                    str(self.logP) + "\t" + str(self.TPSA) + "\n"
            print tmpstr
        else:
            tmpstr = ""
            print "Current Analog's SMILES is EMPTY!\n"
        return tmpstr


def GetAnalogFromFile(fname = "/home/zusongpeng/lab/albiflorin/Analog_1/AL-mode-1-IC50muM"):
    tlist = []
    with open(fname) as fid:
        for line in fid:
            line = line.strip()
            (tname,tsmiles,tIC50) = line.split()
            tmp = ALanalog()
            tmp.getSMILEs(tsmiles)
            tmp.getMolName(tname)
            tmp.getIC50(float(tIC50))
            tmp.getlogP()
            tmp.getTPSA()
            tlist.append(tmp)
    return tlist

def DrawMaplogP(Analog, imgfpre):
    if not not Analog.smilesstr:
        mol = Chem.MolFromSmiles(Analog.smilesstr)
        contribs = rdMolDescriptors._CalcCrippenContribs(mol)
        fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], colorMap='jet', contourLines=10)
        fname = imgfpre + Analog.name + '.png'
        fig.savefig(fname,bbox_inches = "tight")
    else:
        raise "EmptySMILES!"

def DrawMapTPSA(Analog, imgfpre):
    if not not Analog.smilesstr:
        mol = Chem.MolFromSmiles(Analog.smilesstr)
        contribs = rdMolDescriptors._CalcTPSAContribs(mol)
        fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x in contribs],colorMap='jet', contourLines=10)
        fname = imgfpre + Analog.name + '.png'
        fig.savefig(fname, bbox_inches = "tight")
    else:
        raise "EmptySMILES!"

def DrawGasteigerPartialCharge(Analog,imgfpre):
    if not not Analog.smilesstr:
        mol = Chem.MolFromSmiles(Analog.smilesstr)
        AllChem.ComputeGasteigerCharges(mol)
        contribs = [float(mol.GetAtomWithIdx(i).GetProp('_GasteigerCharge')) for i in range(mol.GetNumAtoms())]
        fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,contribs,colorMap='jet', contourLines=10)
       # fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,contribs,contourLines=10)
        fname = imgfpre + Analog.name + '.png'
        fig.savefig(fname, bbox_inches = "tight")
    else:
        raise "EmptySMILES!"

def main(resultfname):
    tlist = GetAnalogFromFile()
    imgfpre1 = "/home/zusongpeng/lab/albiflorin/Analog_1/logP_"
    [DrawMaplogP(x,imgfpre1) for x in tlist]
    
    imgfpre2 = "/home/zusongpeng/lab/albiflorin/Analog_1/TPSA_"
    [DrawMapTPSA(x,imgfpre2) for x in tlist]
    
    imgfpre3 = "/home/zusongpeng/lab/albiflorin/Analog_1/GasteigerPartialCharges_"
    [DrawGasteigerPartialCharge(x,imgfpre3) for x in tlist]
    
    with open(resultfname,"w") as fid:
        fid.write("AnalogName\tSMILTES\tIC50muM\tlogP\tTPSA\n")
        [fid.write(x.printAnalogStr()) for x in tlist]

if __name__ == "__main__":
    resultfname = "/home/zusongpeng/lab/albiflorin/Analog_1/Result_Analog_1.txt"
    main(resultfname)





