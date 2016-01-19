#!/usr/bin/python

# Find drugs (approved and experimental from DrugBank) that has the substructure: benzoic acid.
# Songpeng Zu
# 2015-01-05 10:40:52 CST

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
 

def rdictfile(fname):
    dictn = {}
    with open(fname) as fid:
        for line in fid:
            line = line.strip()
            (key,val) = line.split("\t")
            dictn[key] = val
    return dictn

def SMARTmatchSMILE(patt,smiledict):
    dictn = {}
    SMARTpatt = Chem.MolFromSmarts(patt)
    for key, val in smiledict.iteritems():
        try:
            if  Chem.MolFromSmiles(val).HasSubstructMatch(SMARTpatt):
                dictn[key] = val
        except:
            print "The drug cannot be kekulized.\t"
        finally:
            print "Current Drug is: ", key+"\n"
    return dictn

def SMILEmatchSMILE(patt='c1ccccc1C(=O)O',smiledict={}):
    dictn = {}
    SMARTpatt = Chem.MolFromSmiles(patt)
    for key, val in smiledict.iteritems():
        try:
            if  Chem.MolFromSmiles(val).HasSubstructMatch(SMARTpatt):
                dictn[key] = val
        except:
            print "The drug cannot be kekulized.\t"
        finally:
            print "Current Drug is: ", key+"\n"
    return dictn

def DrawMultiMatch(smiledict,savefile='benacidmatch.png',molsPerRow=5,subImgSize=(200,200)):
    m = [Chem.MolFromSmiles(value) for value in smiledict.itervalues()]
    tmp = [t1.SetProp("_Name",t2) for t1, t2 in zip(m,smiledict.iterkeys())]
    tmp = [AllChem.Compute2DCoords(mi) for mi in m]
    
    img = Draw.MolsToGridImage(m,molsPerRow,subImgSize,legends=[x.GetProp("_Name") for x in m])
    img.save(savefile)








