#!/usr/bin/env python

## Date. 210804
## FilterSpeciese(Human)
## ECFP fingerprint for SoM atom and neighbor atom

## Date. 211101
## Add SDF data

from __future__ import print_function
import copy

# Import modules required 
import sys, time, os, ast, pickle, copy
from rdkit import Chem
from rdkit.Chem import AllChem
from glob import glob
from edge_fingerprint import generate_fp_oddt, generate_fp, atom_featurize, get_properties
import numpy as np

class trans_data():
    def __init__(self):
        # self.filename = filename
        # self.metaboliteN = '_'
        self.sdfns = []; self.pmols = []; self.soms = []; self.ys = []
        self.fps1 = []; self.fps2 = []; self.fps3 = []; self.fps4 = []; self.fps5 = []
        self.fps6 = []; self.fps7 = []; self.fps8 = []; self.fps9 = []; self.fps10 = []
        
    def loadMolFile(self, filename):
        if os.path.splitext(filename)[1] == '.sdf':
            idx2item = {}; idx2mol = {}
            supplier = Chem.SDMolSupplier(str(filename)
                                        ,sanitize=False
                                        ,strictParsing=False)
            
            metabolic_steps = 0
            for mol in supplier:
                mol_name = mol.GetProp('_Name')
                idx2mol[int(mol_name)] = Chem.RemoveHs(mol)
                metabolic_steps += 1
            ids = list(idx2mol.keys())
            ids.sort()

            # Get connectivity info. from parent mol's properties
            connec_pairs = supplier[0].GetProp('Connectivity')
            connec_pairs = ast.literal_eval(connec_pairs) # Trans type: str to list
            
            print(filename)
            for connec_pair in connec_pairs:
                p = connec_pair[0] # Parent
                m_pair = connec_pair[1:] # Metabolite
                pmol = idx2mol[p]
                pmol = Chem.rdmolops.AddHs(pmol)

                # pmol = Chem.RemoveHs(pmol) #2011.11.10
                patoms = [patom for patom in pmol.GetAtoms()]

                for m in m_pair:
                    mmol = idx2mol[m]
                    mmol_prop = mmol.GetPropsAsDict()
                    if mmol_prop['Species'] == 'Human': 
                        # print("@@@@@@@@@@@ SOM: ",mmol_prop['SOM'])
                        somlist = ast.literal_eval(mmol_prop['SOM'])
                        
                        # select maximum numboer of som atoms
                        if len(somlist) > 2:
                            somlist = somlist[:2]
                        
                        for atom in pmol.GetAtoms():
                            atomIdx = atom.GetIdx()
                            if atomIdx in somlist:
                                y = 'SoM'
                            else:
                                y = 'NoSoM'

                            yield [filename, pmol, atomIdx, y, \
                                [(atomIdx, \
                                atom_featurize(atom.GetSymbol(), \
                                atom.GetIsAromatic(), \
                                atom.GetHybridization()),
                                atom.GetDegree())], \
                                generate_fp(pmol, atomIdx, 1), \
                                generate_fp(pmol, atomIdx, 2), \
                                generate_fp(pmol, atomIdx, 3), \
                                generate_fp(pmol, atomIdx, 4), \
                                generate_fp(pmol, atomIdx, 5), \
                                generate_fp(pmol, atomIdx, 6), \
                                generate_fp(pmol, atomIdx, 7), \
                                generate_fp(pmol, atomIdx, 8), \
                                generate_fp(pmol, atomIdx, 9)]

    def sdf_to_list(self, genrator):
        for g in genrator:
            if len(g) == 14:
                sdfn = g[0]; pm = g[1]; som = g[2]
                y = g[3]; fp1 = g[4]; fp2 = g[5]; fp3 = g[6]; fp4 = g[7]; fp5 = g[8]
                fp6 = g[9]; fp7 = g[10]; fp8 = g[11]; fp9 = g[12]; fp10 = g[13]
                
                self.sdfns.append(sdfn)
                self.pmols.append(pm)
                self.soms.append(som)
                self.ys.append(y)
                self.fps1.append(fp1)
                self.fps2.append(fp2)
                self.fps3.append(fp3)
                self.fps4.append(fp4)
                self.fps5.append(fp5)
                self.fps6.append(fp6)
                self.fps7.append(fp7)
                self.fps8.append(fp8)
                self.fps9.append(fp9)
                self.fps10.append(fp10)
                
    def sdf_to_dump(self, dump_fn):
        f = open(dump_fn,'wb')
        pickle.dump(self.sdfns, f)
        pickle.dump(self.pmols, f)
        pickle.dump(self.soms, f)
        pickle.dump(self.ys, f)
        pickle.dump(self.fps1, f)
        pickle.dump(self.fps2, f)
        pickle.dump(self.fps3, f)
        pickle.dump(self.fps4, f)
        pickle.dump(self.fps5, f)
        pickle.dump(self.fps6, f)
        pickle.dump(self.fps7, f)
        pickle.dump(self.fps8, f)
        pickle.dump(self.fps9, f)
        pickle.dump(self.fps10, f)
        f.close()
        print(">> End")

class select_molecules():
    def __init__(self) -> None:
        self.selection = []
        self.select_reaction = {'R1NHC(=O)R2 to R1NH2 + R2COOH':[] # More than 300 : All of these reaction are included in Phase I.
                            ,'R1CR2 to R1C(OH)R2':[]
                            ,'R1NCR2 to R1NH + R2=O':[]
                            ,'R1OCR2 to R1OH + R2=O':[]
                            ,'RC to RC-OH':[]
                            ,'RC-OH to RC=O':[]
                            ,'Arene to Arene-OH':[]}

    def get_sdf_lst(self,dir_path):
        return glob(dir_path)
        
def load_to_dump(dump_fn):
    f = open(dump_fn,'rb')
    sdfns = pickle.load(f)
    pmols = pickle.load(f)
    soms = pickle.load(f)
    ys = pickle.load(f)
    fps1 = pickle.load(f)
    fps2 = pickle.load(f)
    fps3 = pickle.load(f)
    fps4 = pickle.load(f)
    fps5 = pickle.load(f)
    fps6 = pickle.load(f)
    fps7 = pickle.load(f)
    fps8 = pickle.load(f)
    fps9 = pickle.load(f)
    fps10 = pickle.load(f)
    f.close()
    return sdfns, pmols, soms, ys,\
            fps1, fps2, fps3, fps4, fps5, \
            fps6, fps7, fps8, fps9, fps10
        
def do1():
    l = trans_data()
    # fn = "D:\\script\\4Metabolism\\Data\\SDF\\*.sdf" # fn = "D:\script\\4Metabolism\Editor\*.sdf"
    fn = "../../../Data/SDF/*.sdf"

    # dump_file = "D:\\script\\4Metabolism\\Data\\SDF_to_data_211101.dump" # sdf filename, rdkit_mol type, Site of Metabolism
    dump_file = "../../../Data/SDF_to_data_211101.dump"

    for f in glob(fn):
        if f.split("\\")[-1][0] in ['K','S','Y']:
            l.sdf_to_list(l.loadMolFile(f))
    print(" > > > > ", len(set(l.sdfns))) # Print Original SDF file name & path
    print(len(l.sdfns))
    l.sdf_to_dump(dump_file)

def do2():
    # sdfns, pmols, soms, ys, fps1, fps2, fps3, fps4, fps5, fps6, fps7, fps8, fps9, fps10 = load_to_dump("D:\\script\\4Metabolism\\Data\\SDF_to_data_211101.dump")
    sdfns, pmols, soms, ys, fps1, fps2, fps3, fps4, fps5, fps6, fps7, fps8, fps9, fps10 = load_to_dump("../../../data/SDF_to_data_211101.dump")
    
    print(" @ In do2 function >> ")
    print(sdfns[0])
    print(pmols[0])
    print(soms[0])
    print(ys[0])
    print(fps1[0])
    print(fps2[0])
    print(fps3[0])
    print(fps4[0])
    print(fps5[0])
    print(fps6[0])
    print(fps7[0])
    print(fps8[0])
    print(fps9[0])
    print(fps10[0])


def do3():
    fns = []; mols = []
    soms = []; ys = []
    fps1 = []; fps2 = []; fps3 = []; fps4 = []; fps5 = []
    fps6 = []; fps7 = []; fps8 = []; fps9 = []; fps10 = []

    # for fn in glob("D:\\script\\4Metabolism\\Data\\CYP_DB_Zaretzki_XenoSite\\*.sdf"):
    for fn in glob("../../../data/CYP_DB_Zaretzki_XenoSite/*.sdf"):
        
        print(fn)
        for m in Chem.SDMolSupplier(fn):
            mprop = m.GetPropsAsDict()
            # m = Chem.RemoveHs(m) #2011.11.10
            
            somdict = {}
            for key in mprop:
                if 'SOM' in key:
                    somdict[key] = mprop[key]
            
            for skey in somdict:
                somlist = [int(som)-1 for som in str(somdict[skey]).split(' ')]
                
                # select maximum numboer of som atoms
                if len(somlist) > 2:
                    somlist = somlist[:2]

                for atom in m.GetAtoms():
                    atomIdx = atom.GetIdx()
                    if atomIdx in somlist:
                        y = 'SoM'
                    else:
                        y = 'NoSoM'
            
                    fns.append(fn)
                    mols.append(m)
                    soms.append(atomIdx)
                    ys.append(y)
                    fps1.append([(atomIdx, \
                                atom_featurize(atom.GetSymbol(), \
                                atom.GetIsAromatic(), \
                                atom.GetHybridization()),
                                atom.GetDegree())])
                    fps2.append(generate_fp(m, atomIdx, 1))
                    fps3.append(generate_fp(m, atomIdx, 2))
                    fps4.append(generate_fp(m, atomIdx, 3))
                    fps5.append(generate_fp(m, atomIdx, 4))
                    fps6.append(generate_fp(m, atomIdx, 5))
                    fps7.append(generate_fp(m, atomIdx, 6))
                    fps8.append(generate_fp(m, atomIdx, 7))
                    fps9.append(generate_fp(m, atomIdx, 8))
                    fps10.append(generate_fp(m, atomIdx, 9))

    # dump_fn = "D:\\script\\4Metabolism\\Data\\Zretzki_to_data_211101.dump" # sdf filename, rdkit_mol type, Site of Metabolism
    dump_fn = "../../../data/Zretzki_to_data_211101.dump"

    f = open(dump_fn,'wb')
    pickle.dump(fns, f)
    pickle.dump(mols, f)
    pickle.dump(soms, f)
    pickle.dump(ys, f)
    pickle.dump(fps1, f)
    pickle.dump(fps2, f)
    pickle.dump(fps3, f)
    pickle.dump(fps4, f)
    pickle.dump(fps5, f)
    pickle.dump(fps6, f)
    pickle.dump(fps7, f)
    pickle.dump(fps8, f)
    pickle.dump(fps9, f)
    pickle.dump(fps10, f)
    f.close()
    print(">> End")

def do4():
    # sdfns, pmols, soms, ys, fps1, fps2, fps3, fps4, fps5, fps6, fps7, fps8, fps9, fps10 = load_to_dump("D:\\script\\4Metabolism\\Data\\Zretzki_to_data_211101.dump")
    sdfns, pmols, soms, ys, fps1, fps2, fps3, fps4, fps5, fps6, fps7, fps8, fps9, fps10 = load_to_dump("../../../data/Zretzki_to_data_211101.dump")
    
    print(" @ In do4 function >> ")
    print(sdfns[0])
    print(pmols[0])
    print(soms[0])
    print(ys[0])
    print(fps1[0])
    print(fps2[0])
    print(fps3[0])
    print(fps4[0])
    print(fps5[0])
    print(fps6[0])
    print(fps7[0])
    print(fps8[0])
    print(fps9[0])
    print(fps10[0])
    
def do5(): # Generate ATF map & y_map template
    # sdf_sdfns, sdf_pmols, sdf_soms, sdf_ys, sdf_fps1, sdf_fps2, sdf_fps3, sdf_fps4, sdf_fps5, sdf_fps6, sdf_fps7, sdf_fps8, sdf_fps9, sdf_fps10 = load_to_dump("../data/SDF_to_data_211101.dump")
    sdf_sdfns, sdf_pmols, sdf_soms, sdf_ys, sdf_fps1, sdf_fps2, sdf_fps3, sdf_fps4, sdf_fps5, sdf_fps6, sdf_fps7, sdf_fps8, sdf_fps9, sdf_fps10 = load_to_dump("../../../data/SDF_to_data_211101.dump")
    # zret_sdfns, zret_pmols, zret_soms, zret_ys, zret_fps1, zret_fps2, zret_fps3, zret_fps4, zret_fps5, zret_fps6, zret_fps7, zret_fps8, zret_fps9, zret_fps10 = load_to_dump("../data/Zretzki_to_data_211101.dump")
    zret_sdfns, zret_pmols, zret_soms, zret_ys, zret_fps1, zret_fps2, zret_fps3, zret_fps4, zret_fps5, zret_fps6, zret_fps7, zret_fps8, zret_fps9, zret_fps10 = load_to_dump("../../../data/Zretzki_to_data_211101.dump")
    
    sdf_fps_all_1 = [sdf_fp for sdf_fp in sdf_fps1]; sdf_fps_all_2 = [sdf_fp for sdf_fp in sdf_fps2]
    sdf_fps_all_3 = [sdf_fp for sdf_fp in sdf_fps3]; sdf_fps_all_4 = [sdf_fp for sdf_fp in sdf_fps4]
    sdf_fps_all_5 = [sdf_fp for sdf_fp in sdf_fps5]; sdf_fps_all_6 = [sdf_fp for sdf_fp in sdf_fps6]
    sdf_fps_all_7 = [sdf_fp for sdf_fp in sdf_fps7]; sdf_fps_all_8 = [sdf_fp for sdf_fp in sdf_fps8]
    sdf_fps_all_9 = [sdf_fp for sdf_fp in sdf_fps9]; sdf_fps_all_10 = [sdf_fp for sdf_fp in sdf_fps10]

    zret_fps_all_1 = [zret_fp for zret_fp in zret_fps1]; zret_fps_all_2 = [zret_fp for zret_fp in zret_fps2]
    zret_fps_all_3 = [zret_fp for zret_fp in zret_fps3]; zret_fps_all_4 = [zret_fp for zret_fp in zret_fps4]
    zret_fps_all_5 = [zret_fp for zret_fp in zret_fps5]; zret_fps_all_6 = [zret_fp for zret_fp in zret_fps6]
    zret_fps_all_7 = [zret_fp for zret_fp in zret_fps7]; zret_fps_all_8 = [zret_fp for zret_fp in zret_fps8]
    zret_fps_all_9 = [zret_fp for zret_fp in zret_fps9]; zret_fps_all_10 = [zret_fp for zret_fp in zret_fps10]
    
    
    sdf_fps_all_1.extend(zret_fps_all_1); sdf_fps_all_2.extend(zret_fps_all_2)
    sdf_fps_all_3.extend(zret_fps_all_3); sdf_fps_all_4.extend(zret_fps_all_4)
    sdf_fps_all_5.extend(zret_fps_all_5); sdf_fps_all_6.extend(zret_fps_all_6)
    sdf_fps_all_7.extend(zret_fps_all_7); sdf_fps_all_8.extend(zret_fps_all_8)
    sdf_fps_all_9.extend(zret_fps_all_9); sdf_fps_all_10.extend(zret_fps_all_10)

    # print(">> ", len(sdf_fps_all))
    
    atfmap1 = []
    for sdf_fp in sdf_fps_all_1:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap1:
                atfmap1.append(atf)
                
    atfmap2 = []
    for sdf_fp in sdf_fps_all_2:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap2:
                atfmap2.append(atf)

    atfmap3 = []
    for sdf_fp in sdf_fps_all_3:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap3:
                atfmap3.append(atf)

    atfmap4 = []
    for sdf_fp in sdf_fps_all_4:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap4:
                atfmap4.append(atf)

    atfmap5 = []
    for sdf_fp in sdf_fps_all_5:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap5:
                atfmap5.append(atf)

    atfmap6 = []
    for sdf_fp in sdf_fps_all_6:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap6:
                atfmap6.append(atf)
                
    atfmap7 = []
    for sdf_fp in sdf_fps_all_7:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap7:
                atfmap7.append(atf)

    atfmap8 = []
    for sdf_fp in sdf_fps_all_8:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap8:
                atfmap8.append(atf)

    atfmap9 = []
    for sdf_fp in sdf_fps_all_9:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap9:
                atfmap9.append(atf)

    atfmap10 = []
    for sdf_fp in sdf_fps_all_10:
        for fp in sdf_fp:
            atf = fp[1]+"_"+str(fp[2])
            if fp[1][0] != 'H' and not atf in atfmap10:
                atfmap10.append(atf)

    print(len(atfmap1))
    print(atfmap1)
    print(len(atfmap10))
    print(atfmap10)

    sdf_ys_all = [sdf_y for sdf_y in sdf_ys]
    zret_ys_all = [zret_y for zret_y in zret_ys]
    
    sdf_ys_all.extend(zret_ys_all)

    print(">> ", len(sdf_ys_all))
    
    y_map = []
    for sdf_y in sdf_ys_all:
        # if sdf_y[0] != "H" and sdf_y not in y_map:
        if sdf_y not in y_map:
            y_map.append(sdf_y)
    
    print(len(y_map))
    print(y_map)

class atfmap():
    def __init__(self):
        self.map_dicts = {'C_SP3_4':0, 'O_SP2_1':0, 'N_SP3_3':0, 'O_SP2_2':0, \
                        'N_SP2_3':0, 'C_SP2_3':0, 'S_SP3_2':0, 'CAr_SP2_3':0, \
                        'C_SP_2':0, 'NAr_SP2_3':0, 'NAr_SP2_2':0, 'F_SP3_1':0, \
                        'Cl_SP3_1':0, 'N_SP2_2':0, 'O_SP3_2':0, 'SAr_SP2_2':0, \
                        'S_SP3_4':0, 'N_SP3_4':0, 'N_SP_1':0, 'OAr_SP2_2':0, \
                        'Br_SP3_1':0, 'S_SP2_1':0, 'S_SP3_3':0, 'I_SP3_1':0, \
                        'P_SP3_4':0, 'O_SP3_1':0, 'N_SP2_1':0, 'S_SP3D_4':0, \
                        'C_SP3_1':0, 'C_SP2_2':0, 'C_SP3_2':0, 'CAr_SP2_2':0, \
                        'C_SP3_3':0, 'N_SP3_2':0, 'N_SP3_1':0, 'B_SP2_3':0, \
                        'C_SP2_1':0, 'C_SP_1':0} # Get form do5 function
        
class y_map():
    def __init__(self):
        self.ymap = {'SoM':0, 'NoSoM':0}
                    
def do6():
    # sdf_sdfns, sdf_pmols, sdf_soms, sdf_ys, sdf_fps1, sdf_fps2, sdf_fps3, sdf_fps4, sdf_fps5, sdf_fps6, sdf_fps7, sdf_fps8, sdf_fps9, sdf_fps10 = load_to_dump("D:\\script\\4Metabolism\\Data\\SDF_to_data_211101.dump")
    sdf_sdfns, sdf_pmols, sdf_soms, sdf_ys, sdf_fps1, sdf_fps2, sdf_fps3, sdf_fps4, sdf_fps5, sdf_fps6, sdf_fps7, sdf_fps8, sdf_fps9, sdf_fps10 = load_to_dump("../../../data/SDF_to_data_211101.dump")
    # zret_sdfns, zret_pmols, zret_soms, zret_ys, zret_fps1, zret_fps2, zret_fps3, zret_fps4, zret_fps5, zret_fps6, zret_fps7, zret_fps8, zret_fps9, zret_fps10 = load_to_dump("D:\\script\\4Metabolism\\Data\\Zretzki_to_data_211101.dump")
    zret_sdfns, zret_pmols, zret_soms, zret_ys, zret_fps1, zret_fps2, zret_fps3, zret_fps4, zret_fps5, zret_fps6, zret_fps7, zret_fps8, zret_fps9, zret_fps10 = load_to_dump("../../../data/Zretzki_to_data_211101.dump")

    m = atfmap(); ydict = y_map()
    all_sdfns = []; all_pmols = []; all_soms = []; all_ys = []
    all_fps1 = []; all_fps2 = []; all_fps3 = []; all_fps4 = []; all_fps5 = []
    all_fps6 = []; all_fps7 = []; all_fps8 = []; all_fps9 = []; all_fps10 = []

    for sdfns, pmols, som, y, fps1, fps2, fps3, fps4, fps5 \
        , fps6, fps7, fps8, fps9, fps10 in zip(sdf_sdfns, sdf_pmols, sdf_soms, sdf_ys \
                                                , sdf_fps1, sdf_fps2, sdf_fps3, sdf_fps4, sdf_fps5 \
                                                , sdf_fps6, sdf_fps7, sdf_fps8, sdf_fps9, sdf_fps10):

        ym = copy.deepcopy(ydict.ymap)
        ym[y] += 1
        y_list = list(ym.values())
        
        mmap1 = copy.deepcopy(m.map_dicts)
        mmap2 = copy.deepcopy(m.map_dicts)
        mmap3 = copy.deepcopy(m.map_dicts)
        mmap4 = copy.deepcopy(m.map_dicts)
        mmap5 = copy.deepcopy(m.map_dicts)
        mmap6 = copy.deepcopy(m.map_dicts)
        mmap7 = copy.deepcopy(m.map_dicts)
        mmap8 = copy.deepcopy(m.map_dicts)
        mmap9 = copy.deepcopy(m.map_dicts)
        mmap10 = copy.deepcopy(m.map_dicts)
        
        if len(fps1):
            for fps in fps1:
                atf1 = fps[1]+"_"+str(fps[2])
                if atf1 in mmap1:
                    mmap1[atf1] += 1 
            
        if len(fps2):
            for fps in fps2:
                atf2 = fps[1]+"_"+str(fps[2])
                if atf2 in mmap2:
                    mmap2[atf2] += 1 
            
        if len(fps3):
            for fps in fps3:
                atf3 = fps[1]+"_"+str(fps[2])
                if atf3 in mmap3:
                    mmap3[atf3] += 1 
            
        if len(fps4):
            for fps in fps4:
                atf4 = fps[1]+"_"+str(fps[2])
                if atf4 in mmap4:
                    mmap4[atf4] += 1 
            
        if len(fps5):
            for fps in fps5:
                atf5 = fps[1]+"_"+str(fps[2])
                if atf5 in mmap5:
                    mmap5[atf5] += 1 
            
        if len(fps6):
            for fps in fps6:
                atf6 = fps[1]+"_"+str(fps[2])
                if atf6 in mmap6:
                    mmap6[atf6] += 1 
            
        if len(fps7):
            for fps in fps7:
                atf7 = fps[1]+"_"+str(fps[2])
                if atf7 in mmap7:
                    mmap7[atf7] += 1 
            
        if len(fps8):
            for fps in fps8:
                atf8 = fps[1]+"_"+str(fps[2])
                if atf8 in mmap8:
                    mmap8[atf8] += 1 
            
        if len(fps9):
            for fps in fps9:
                atf9 = fps[1]+"_"+str(fps[2])
                if atf9 in mmap9:
                    mmap9[atf9] += 1 
            
        if len(fps10):
            for fps in fps10:
                atf10 = fps[1]+"_"+str(fps[2])
                if atf10 in mmap10:
                    mmap10[atf10] += 1 
            
        # if sum(mmap.values()):
        all_sdfns.append(sdfns)
        all_pmols.append(pmols)
        all_soms.append(som)
        all_ys.append(y_list)
        all_fps1.append(mmap1); all_fps2.append(mmap2)
        all_fps3.append(mmap3); all_fps4.append(mmap4)
        all_fps5.append(mmap5); all_fps6.append(mmap6)
        all_fps7.append(mmap7); all_fps8.append(mmap8)
        all_fps9.append(mmap9); all_fps10.append(mmap10)

    for sdfns, pmols, som, y, fps1, fps2, fps3, fps4, fps5 \
        , fps6, fps7, fps8, fps9, fps10 in zip(zret_sdfns, zret_pmols, zret_soms, zret_ys \
                                                , zret_fps1, zret_fps2, zret_fps3, zret_fps4, zret_fps5 \
                                                , zret_fps6, zret_fps7, zret_fps8, zret_fps9, zret_fps10):
        
        ym = copy.deepcopy(ydict.ymap)
        ym[y] += 1
        y_list = list(ym.values())

        mmap1 = copy.deepcopy(m.map_dicts)
        mmap2 = copy.deepcopy(m.map_dicts)
        mmap3 = copy.deepcopy(m.map_dicts)
        mmap4 = copy.deepcopy(m.map_dicts)
        mmap5 = copy.deepcopy(m.map_dicts)
        mmap6 = copy.deepcopy(m.map_dicts)
        mmap7 = copy.deepcopy(m.map_dicts)
        mmap8 = copy.deepcopy(m.map_dicts)
        mmap9 = copy.deepcopy(m.map_dicts)
        mmap10 = copy.deepcopy(m.map_dicts)
        
        if len(fps1):
            for fps in fps1:
                    atf1 = fps[1]+"_"+str(fps[2])
                    if atf1 in mmap1:
                        mmap1[atf1] += 1 
        
        if len(fps2):
            for fps in fps2:
                atf2 = fps[1]+"_"+str(fps[2])
                if atf2 in mmap2:
                    mmap2[atf2] += 1 
            
        if len(fps3):
            for fps in fps3:
                atf3 = fps[1]+"_"+str(fps[2])
                if atf3 in mmap3:
                    mmap3[atf3] += 1 
            
        if len(fps4):
            for fps in fps4:
                atf4 = fps[1]+"_"+str(fps[2])
                if atf4 in mmap4:
                    mmap4[atf4] += 1 
            
        if len(fps5):
            for fps in fps5:
                atf5 = fps[1]+"_"+str(fps[2])
                if atf5 in mmap5:
                    mmap5[atf5] += 1 
            
        if len(fps6):
            for fps in fps6:
                atf6 = fps[1]+"_"+str(fps[2])
                if atf6 in mmap6:
                    mmap6[atf6] += 1 
            
        if len(fps7):
            for fps in fps7:
                atf7 = fps[1]+"_"+str(fps[2])
                if atf7 in mmap7:
                    mmap7[atf7] += 1 
            
        if len(fps8):
            for fps in fps8:
                atf8 = fps[1]+"_"+str(fps[2])
                if atf8 in mmap8:
                    mmap8[atf8] += 1 
            
        if len(fps9):
            for fps in fps9:
                atf9 = fps[1]+"_"+str(fps[2])
                if atf9 in mmap9:
                    mmap9[atf9] += 1 
            
        if len(fps10):
            for fps in fps10:
                atf10 = fps[1]+"_"+str(fps[2])
                if atf10 in mmap10:
                    mmap10[atf10] += 1 

        # if sum(mmap.values()):
        all_sdfns.append(sdfns)
        all_pmols.append(pmols)
        all_soms.append(som)
        all_ys.append(y_list)
        all_fps1.append(mmap1); all_fps2.append(mmap2)
        all_fps3.append(mmap3); all_fps4.append(mmap4)
        all_fps5.append(mmap5); all_fps6.append(mmap6)
        all_fps7.append(mmap7); all_fps8.append(mmap8)
        all_fps9.append(mmap9); all_fps10.append(mmap10)

    # dump_fn = "D:\\script\\4Metabolism\\Data\\all_data_211101.dump"
    dump_fn = "../../../data/all_data_211101.dump"
    f = open(dump_fn,'wb')
    pickle.dump(all_sdfns, f)
    pickle.dump(all_pmols, f)
    pickle.dump(all_soms, f)
    pickle.dump(all_ys, f)
    pickle.dump(all_fps1, f); pickle.dump(all_fps2, f)
    pickle.dump(all_fps3, f); pickle.dump(all_fps4, f)
    pickle.dump(all_fps5, f); pickle.dump(all_fps6, f)
    pickle.dump(all_fps7, f); pickle.dump(all_fps8, f)
    pickle.dump(all_fps9, f); pickle.dump(all_fps10, f)
    f.close()
    print(">> End")
    print(len(all_sdfns))
    print(len(all_pmols))
    print(len(all_soms))
    print(len(all_ys))
    print(len(all_fps1))
    print(len(all_fps2))
    print(len(all_fps3))
    print(len(all_fps4))
    print(len(all_fps5))
    print(len(all_fps6))
    print(len(all_fps7))
    print(len(all_fps8))
    print(len(all_fps9))
    print(len(all_fps10))
    print(all_fps1[0])

def do7():
    # read dump file
    # dump_fn = "D:\\script\\4Metabolism\\Data\\all_data_211101.dump"
    dump_fn = "../../../data/all_data_211101.dump"
    
    f = open(dump_fn,'rb')
    all_sdfns = pickle.load(f)
    all_pmols = pickle.load(f)
    all_soms = pickle.load(f)
    all_ys = pickle.load(f)
    all_fps1 = pickle.load(f); all_fps2 = pickle.load(f)
    all_fps3 = pickle.load(f); all_fps4 = pickle.load(f)
    all_fps5 = pickle.load(f); all_fps6 = pickle.load(f)
    all_fps7 = pickle.load(f); all_fps8 = pickle.load(f)
    all_fps9 = pickle.load(f); all_fps10 = pickle.load(f)
    f.close()
    print("@@ Load Dump file")
    print("Data Index: ",len(all_sdfns))
    
    print(len(all_fps1[0]))
    
def do8():
    # Generate training data
    # dump_fn = "D:\\script\\4Metabolism\\Data\\all_data_211101.dump"
    dump_fn = "../../../data/all_data_211101.dump"
    
    f = open(dump_fn,'rb')
    all_sdfns = np.array(pickle.load(f))
    all_pmols = np.array(pickle.load(f))
    all_soms = np.array(pickle.load(f))
    all_ys = np.array(pickle.load(f))
    all_fps1 = pickle.load(f); all_fps2 = pickle.load(f)
    all_fps3 = pickle.load(f); all_fps4 = pickle.load(f)
    all_fps5 = pickle.load(f); all_fps6 = pickle.load(f)
    all_fps7 = pickle.load(f); all_fps8 = pickle.load(f)
    all_fps9 = pickle.load(f); all_fps10 = pickle.load(f)
    f.close()

    # List to NumPy array
    all_fps1 = np.array(all_fps1); all_fps2 = np.array(all_fps2)
    all_fps3 = np.array(all_fps3); all_fps4 = np.array(all_fps4)
    all_fps5 = np.array(all_fps5); all_fps6 = np.array(all_fps6)
    all_fps7 = np.array(all_fps7); all_fps8 = np.array(all_fps8)
    all_fps9 = np.array(all_fps9); all_fps10 = np.array(all_fps10)
    
    all_fps = [np.array([]) for _ in range(len(all_fps1))]
    all_bit_fps = [np.array([]) for _ in range(len(all_fps1))]    
    data_index = len(all_fps)

    for i in range(data_index):
        all_fps[i] = np.array(all_fps[i])
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps1[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps2[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps3[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps4[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps5[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps6[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps7[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps8[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps9[i].values())), axis=None)
        all_fps[i] = np.concatenate((all_fps[i], list(all_fps10[i].values())), axis=None)
        
        bit_fps = [np.array([]) for _ in range(data_index)]

        klst = []
        for j in all_fps[i]:
            j = int(j)
            for k in range(32):
                if k < j:
                    klst.append(1)
                else:
                    klst.append(0)
        bit_fps[j] = np.concatenate((bit_fps[j], klst), axis=None)
        bit_fps = np.array(all_bit_fps)
        all_bit_fps.append(bit_fps)
    all_bit_fps = np.array(all_bit_fps)
    print(">>>>>>>>>>>>>>>>>>>>>>>>>> ",all_bit_fps.shape)

    # new_dump_fn = "D:\\script\\4Metabolism\\Data\\all_bit_fp_211101.dump"
    new_dump_fn = "../../../Data/all_bit_fp_211126.dump"
    
    f = open(new_dump_fn,'wb')
    pickle.dump(all_bit_fps, f)
    pickle.dump(all_ys, f)
    f.close()

if __name__=='__main__':
    # do1() # Generate Our SDF data
    # do2() # read dump file
    # do3() # Generate Zaretzki's data
    # do4() # read dump file
    # # # ===============================
    # do5() # Generate ATF map template
    # do6() # Generate curated data
    # do7() # read dump file

    do8() # Generate training data
    # ===============================
    