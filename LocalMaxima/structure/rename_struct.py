#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 21:47:22 2019

@author: noel
"""
import Bio.PDB as struct
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
import sys
import copy

class rename(object):
    def __init__(self, d, f):
        self.dir1 = d
        self.file = f
        pdb_parser = PDBParser(QUIET = True)
        self.strct = pdb_parser.get_structure('main', self.dir1+self.file)
        self.new_s = struct.StructureBuilder.Structure("Renamed S")
        
    def model_error_collision(self, entity, entity_id):
        print("Error: "+entity+" "+str(entity_id)+" cannot be renamed to"+
              "a "+entity.lower()+" that already exist. Name collision.")
        print("       Exit without completion.")
        sys.exit(1)
        
    def chain_error_collision(self, entity, chain_id, model_id):
        print("Error: "+entity+" "+chain_id+" in model "+model_id+" cannot be "+
              "renamed to a "+entity.lower()+" that already exist. Name collision.")
        print("       Exit without completion.")
        sys.exit(1)
        
    def residue_error_collision(self, entity, residue_id, chain_id, model_id):
        print("Error: "+entity+" "+residue_id+" in model "+model_id+" and chain "+chain_id+
              " cannot be renamed to a "+entity.lower()+" that already exist. Name collision.")
        print("Exit without completion.")
        sys.exit(1)
        
    def atom_error_collision(self, entity, atom_id, residue_id, chain_id, model_id):
        print("Error: "+entity+" "+atom_id+" in model "+model_id+", chain "+chain_id+" and residue"+
              residue_id+" cannot be renamed to a "+entity.lower()+" that already exist. Name collision.")
        print("Exit without completion.")
        sys.exit(1)
    
    def error_missmatch_ast_nonast(self, entity):
        print("Error: "+entity+" selection to be ranemed, and new "+entity.lower()+
              "s' names missmatch.")
        print("       "+entity+" to be rename is '*', and new name is not. A one ")
        print("       to one correspondance is needed and "+entity.lower()+" must")
        print("       be explictly declared, no '*'.")
        print("       Exit without completion.")
        sys.exit(1)
    
    def error_missmatch_nonast_ast(self, entity):
        print("Error: New "+entity.lower()+" name and "+entity.lower()+"s'")
        print("       selection to be renamed missmatch. New "+entity.lower())
        print("       is '*', and "+entity.lower()+" to be renamed is not. A one")
        print("       to one correspondance is needed and "+entity.lower()+"s must")
        print("       be explict declared, no '*'.")
        print("       Exit without completion.")
        sys.exit(1)
    
    def error_equal_entity_selection(self, entity):
        print("Error: "+entity+"s to be renamed and new "+entity.lower()+"s' names")
        print("       listed must be equal in size. A one to one correspondance is")
        print("       needed and "+entity.lower()+"s must be explicitily declared, no '*'")
        print("       Exit without completion.")
        sys.exit(1)
                       
    def gen_closer_lower_atom(self, i, j, res1_list, entity_to_add):
        closer_lower = {}
        max_res = max(res1_list)
        for kk in entity_to_add[i.get_id()][j.get_id()]:
            if kk > max_res:
                if max_res in closer_lower:
                    pass
                else:
                    closer_lower[max_res] = []
                closer_lower[max_res].append(kk)
            else:
                for ll in range(len(res1_list)-1):
                    if kk > res1_list[ll] and kk < res1_list[ll+1]:
                        if res1_list[ll] in closer_lower:
                            pass
                        else:
                            closer_lower[res1_list[ll]] = []
                        closer_lower[res1_list[ll]].append(kk)
        return closer_lower
    
    def gen_closer_lower_resi(self, i, j, res1_list, residues_to_add):
        closer_lower = {}
        max_res = max(res1_list)
        for kk in residues_to_add[i.get_id()][j.get_id()]:
            if kk.get_id()[1] > max_res:
                if max_res in closer_lower:
                    pass
                else:
                    closer_lower[max_res] = []
                closer_lower[max_res].append(kk.get_id()[1])
            else:
                for ll in range(len(res1_list)-1):
                    if kk.get_id()[1] > res1_list[ll] and kk.get_id()[1] < res1_list[ll+1]:
                        if res1_list[ll] in closer_lower:
                            pass
                        else:
                            closer_lower[res1_list[ll]] = []
                        closer_lower[res1_list[ll]].append(kk.get_id()[1])
        return closer_lower
    
    def gen_closer_lower_chain(self, i, chn1_list, chains_to_add):
        closer_lower = {}
        max_res = max(chn1_list)
        for kk in chains_to_add[i.get_id()]:
            if kk.get_id() > max_res:
                if max_res in closer_lower:
                    pass
                else:
                    closer_lower[max_res] = []
                closer_lower[max_res].append(kk.get_id())
            else:
                for ll in range(len(chn1_list)-1):
                    if kk.get_id()[1] > chn1_list[ll] and kk.get_id() < chn1_list[ll+1]:
                        if chn1_list[ll] in closer_lower:
                            pass
                        else:
                            closer_lower[chn1_list[ll]] = []
                        closer_lower[chn1_list[ll]].append(kk.get_id())
        return closer_lower
    # TODO f03_0001 has old algorithm.                                  
    def f01_0001(self, se1, se2):
        if se1.sel_dic['atoms'] == se2.sel_dic['atoms']:
            pass
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                new_m = struct.Model.Model(i.get_id())
                for j in i.get_chains():
                    new_c = struct.Chain.Chain(j.get_id())
                    for k in j.get_residues():
                        for l in k.get_atoms():
                            if l.get_id() in se2.sel_dic['atoms']:
                                self.atom_error_collision("Atom", str(l.get_id()), 
                                                                  str(k.get_id()),
                                                                  str(j.get_id()), 
                                                                  str(i.get_id()))
                        for l in k.get_atoms():
                            if l.get_id() in se1.sel_dic['atoms']:
                                indx = se1.sel_dic['atoms'].index(l.get_id())
                                atm2 = se2.sel_dic['atoms'][indx]
                                l.id = atm2
                                l.name = atm2
                                l.fullname = atm2
                        new_r = copy.deepcopy(k)
                        new_r.detach_parent()
                        new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
    # TODO f02_0010 has old algorithm. 
    def f02_0010(self, se1, se2):
        if se1.sel_dic['residues'] == se2.sel_dic['residues']:
            pass
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                new_m = struct.Model.Model(i.get_id())
                for j in i.get_chains():
                    new_c = struct.Chain.Chain(j.get_id())
                    # These loops check that the renaming destinations are available
                    res1_list = []
                    for k in j.get_residues():
                        res1_list.append(k.get_id()[1])
                        if k.get_id()[1] in se2.sel_dic['residues']:
                            self.residue_error_collision("Residue", str(k.get_id()),
                                                                    str(j.get_id()), 
                                                                    str(i.get_id()))
                    # if destinations are all available, the code is still running, then
                    # remove residues from se1 and store them in a list for adding later
                    residues_to_add = {}
                    residues_to_add[i.get_id()] = {}
                    residues_to_add[i.get_id()][j.get_id()] = []
                    residues_to_rem = {}
                    residues_to_rem[i.get_id()] = {}
                    residues_to_rem[i.get_id()][j.get_id()] = []       
                    for k in j.get_residues():
                        if k.get_id()[1] in se1.sel_dic['residues']:
                            mov_r = copy.deepcopy(k)
                            # so you want to move from this residue.
                            # note not where it comes from, but where it goes.
                            indx = se1.sel_dic['residues'].index(k.get_id()[1])
                            res2 = se2.sel_dic['residues'][indx]
                            mov_r.id = (k.get_id()[0],res2,k.get_id()[2])
                            residues_to_add[i.get_id()][j.get_id()].append(mov_r)
                            residues_to_rem[i.get_id()][j.get_id()].append(k.get_id())
                    # In order to know where to add residues, the closer resdies index
                    # that is less than the residues to be added is store in a dictionary
                    closer_lower = self.gen_closer_lower_resi(i, j, res1_list, residues_to_add)
                    for k in j.get_residues():
                        new_r = copy.deepcopy(k)
                        new_r.detach_parent()
                        if k.get_id() not in residues_to_rem[i.get_id()][j.get_id()]:
                            new_c.add(new_r)
                        if k.get_id()[1] in closer_lower.keys():
                            # Ranks and gives indexes of closer lower to guarantee order of output
                            # TODO the following lines work, but sorting might not be necessary.
                            for aa in [ai[0] for ai in sorted(enumerate(closer_lower[k.get_id()[1]]), key=lambda x:x[1])]:
                                new_c.add(residues_to_add[i.get_id()][j.get_id()][aa])
                    new_m.add(new_c)
                self.new_s.add(new_m)
    
    def f03_0011(self, se1, se2):
        if len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                for j in i.get_chains():
                    for k in j.get_residues():
                        if k.get_id()[1] in se2.sel_dic['residues']:
                            for l in k.get_atoms():
                                if l.get_id() in se2.sel_dic['atoms']:
                                    #print("Error", str(l.get_id()), str(k.get_id()),str(j.get_id()), str(i.get_id()))                            
                                    self.atom_error_collision("Atom", str(l.get_id()), 
                                                                      str(k.get_id()),
                                                                      str(j.get_id()), 
                                                                      str(i.get_id()))
            resid_resname_segid_1 = {}
            resid_resname_segid_2 = {}
            strc_dic = {}
            for i in self.strct.get_models():
                resid_resname_segid_1[i.get_id()] = {}
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid_2[i.get_id()] = {}
                for j in i.get_chains():
                    resid_resname_segid_1[i.get_id()][j.get_id()] = {}
                    if j.get_id() in strc_dic[i.get_id()]:
                        pass
                    else:
                        strc_dic[i.get_id()][j.get_id()] = {}
                        resid_resname_segid_2[i.get_id()][j.get_id()] = {}
                    for k in j.get_residues():
                        resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                        resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                        if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                        if k.get_id()[1] in se1.sel_dic['residues']:
                            indx = se1.sel_dic['residues'].index(k.get_id()[1])
                            res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])
                            if res2 in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][res2] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][res2] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][res2].append((res2,k.resname,k.segid))
                            for l in k.get_atoms():
                                if l.get_id() in se1.sel_dic['atoms']:
                                    mov_a = copy.deepcopy(l)
                                    mov_a.detach_parent()
                                    indx = se1.sel_dic['atoms'].index(l.get_id())
                                    atm2 = se2.sel_dic['atoms'][indx]
                                    mov_a.id = atm2
                                    mov_a.name = atm2
                                    mov_a.fullname = atm2
                                    strc_dic[i.get_id()][j.get_id()][res2].append(mov_a)
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                        else:
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    if j in resid_resname_segid_1[i].keys() and\
                       j in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            if k in resid_resname_segid_1[i][j].keys() and\
                               k in resid_resname_segid_2[i][j].keys():
                                res_info_1 = resid_resname_segid_1[i][j][k]
                                res_info_2 = resid_resname_segid_2[i][j][k]
                                if res_info_1[0][1] == res_info_2[0][1]:
                                    new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                                else:
                                    new_r = struct.Residue.Residue(res_info_2[0][0],res_info_1[0][1],res_info_2[0][2])
                                for l in strc_dic[i][j][res_info_2[0][0]]:
                                    new_r.add(l)
                                new_c.add(new_r)
                            elif k not in resid_resname_segid_1[i][j].keys() and\
                                 k in resid_resname_segid_2[i][j].keys():
                                res_info_2 = resid_resname_segid_2[i][j][k]
                                new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                                for l in strc_dic[i][j][res_info_2[0][0]]:
                                    new_r.add(l)
                                new_c.add(new_r)
                            elif k in resid_resname_segid_1[i][j].keys() and\
                                 k not in resid_resname_segid_2[i][j].keys():
                                     print("WARNING: Anticipated situation not encountered before.")
                                     print("         Needs to be coded.")
                    elif j not in resid_resname_segid_1[i].keys() and\
                         j in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            res_info_2 = resid_resname_segid_2[i][j][k]
                            new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                            for l in strc_dic[i][j][res_info_2[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    elif j in resid_resname_segid_1[i].keys() and\
                         j not in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            res_info_1 = resid_resname_segid_1[i][j][k]
                            new_r = struct.Residue.Residue(res_info_1[0][0],res_info_1[0][1],res_info_1[0][2])
                            for l in strc_dic[i][j][res_info_1[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
    # TODO f04_0100 has old algorithm. 
    def f04_0100(self, se1, se2):
        if se1.sel_dic['chains'] == se2.sel_dic['chains']:
            pass
        else:
            for i in self.strct.get_models():
                new_m = struct.Model.Model(i.get_id())
                chn1_list = []
                for j in i.get_chains():
                    chn1_list.append(j.get_id())
                    if j.get_id() in se2.sel_dic['chains']:
                        self.chain_error_collision("Chain", str(j.get_id()), 
                                                            str(i.get_id()))
                    # if destinations are all available, the code is still running
                chains_to_add = {}
                chains_to_add[i.get_id()] = []
                chains_to_rem = {}
                chains_to_rem[i.get_id()] = []
                for j in i.get_chains():
                    if j.get_id() in se1.sel_dic['chains']:
                        new_c = copy.deepcopy(j)
                        new_c.detach_parent()
                        indx = se1.sel_dic['chains'].index(j.get_id())
                        chn2 = se2.sel_dic['chains'][indx]
                        new_c.id = chn2
                        chains_to_add[i.get_id()].append(new_c)
                        chains_to_rem[i.get_id()].append(j.get_id())
                closer_lower = self.gen_closer_lower_chain(i, chn1_list, chains_to_add)
                for j in i.get_chains():                                  
                    new_c = copy.deepcopy(j)
                    new_c.detach_parent()
                    if j.get_id() not in chains_to_rem[i.get_id()]:
                        new_m.add(new_c)
                    if j.get_id() in closer_lower.keys():
                        for aa in [ai[0] for ai in sorted(enumerate(closer_lower[j.get_id()]), key=lambda x:x[1])]:
                            new_m.add(chains_to_add[i.get_id()][aa])
                self.new_s.add(new_m)

    def f05_0101(self, se1, se2):
        for i in self.strct.get_models():
            for j in i.get_chains():
                if j.get_id() in se2.sel_dic['chains']:
                    for k in j.get_residues():
                        for l in k.get_atoms():
                            if l.get_id() in se2.sel_dic['atoms']:
                                self.atom_error_collision("Atom", str(l.get_id()), 
                                                                  str(k.get_id()),
                                                                  str(j.get_id()), 
                                                                  str(i.get_id()))        
        resid_resname_segid_1 = {}
        resid_resname_segid_2 = {}
        strc_dic = {}
        for i in self.strct.get_models():
            resid_resname_segid_1[i.get_id()] = {}
            if i.get_id() in strc_dic:
                pass
            else:
                strc_dic[i.get_id()] = {}
                resid_resname_segid_2[i.get_id()] = {}
            for j in i.get_chains():
                resid_resname_segid_1[i.get_id()][j.get_id()] = {}
                if j.get_id() in strc_dic[i.get_id()]:
                    pass
                else:
                    strc_dic[i.get_id()][j.get_id()] = {}
                    resid_resname_segid_2[i.get_id()][j.get_id()] = {}
                if j.get_id() in se1.sel_dic['chains']:
                    indx = se1.sel_dic['chains'].index(j.get_id())
                    chn2 = se2.sel_dic['chains'][indx]
                    if chn2 in strc_dic[i.get_id()]:
                        pass
                    else:
                        strc_dic[i.get_id()][chn2] = {}
                        resid_resname_segid_2[i.get_id()][chn2] = {}
                    for k in j.get_residues():
                        strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                        resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                        resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                        resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                        resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                        if k.get_id() in strc_dic[i.get_id()][chn2]: 
                            pass
                        else:
                            strc_dic[i.get_id()][chn2][k.get_id()] = []
                            resid_resname_segid_2[i.get_id()][chn2][k.get_id()] = []
                            resid_resname_segid_2[i.get_id()][chn2][k.get_id()].append((k.get_id(),k.resname,k.segid))
                        for l in k.get_atoms():
                            if l.get_id() in se1.sel_dic['atoms']:
                                mov_a = copy.deepcopy(l)
                                mov_a.detach_parent()
                                indx = se1.sel_dic['atoms'].index(l.get_id())
                                atm2 = se2.sel_dic['atoms'][indx]
                                mov_a.id = atm2
                                mov_a.name = atm2
                                mov_a.fullname = atm2
                                strc_dic[i.get_id()][chn2][k.get_id()].append(mov_a)
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    for k in j.get_residues():
                        resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                        resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))                   
                        if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                        for l in k.get_atoms():
                            strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
        for ii in sorted(strc_dic.keys()):
            for jj in sorted(strc_dic[ii].keys()):
                for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                    if len(strc_dic[ii][jj][kk]) == 0:
                        del(strc_dic[ii][jj][kk])
                if len(strc_dic[ii][jj]) == 0:
                    del(strc_dic[ii][jj])
            if len(strc_dic[ii]) == 0:
                del(strc_dic[ii])
        new_m = struct.Model.Model(i.get_id())
        for i in sorted(strc_dic.keys()):
            new_m = struct.Model.Model(i)
            for j in sorted(strc_dic[i].keys()):
                if j in resid_resname_segid_1[i].keys() and\
                   j in resid_resname_segid_2[i].keys():
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        if k in resid_resname_segid_1[i][j].keys() and\
                           k in resid_resname_segid_2[i][j].keys():
                            res_info_1 = resid_resname_segid_1[i][j][k]
                            res_info_2 = resid_resname_segid_2[i][j][k]
                            if res_info_1[0][1] == res_info_2[0][1]:
                                new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                            else:
                                new_r = struct.Residue.Residue(res_info_2[0][0],res_info_1[0][1],res_info_2[0][2])
                            for l in strc_dic[i][j][res_info_2[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                        elif k not in resid_resname_segid_1[i][j].keys() and\
                             k in resid_resname_segid_2[i][j].keys():
                            res_info_2 = resid_resname_segid_2[i][j][k]
                            new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                            for l in strc_dic[i][j][res_info_2[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                        elif k in resid_resname_segid_1[i][j].keys() and\
                             k not in resid_resname_segid_2[i][j].keys():
                                 print("WARNING: Anticipated situation not encountered before.", k)
                                 print("         Needs to be coded.")
                elif j not in resid_resname_segid_1[i].keys() and\
                     j in resid_resname_segid_2[i].keys():
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info_2 = resid_resname_segid_2[i][j][k]
                        new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                        for l in strc_dic[i][j][res_info_2[0][0]]:
                            new_r.add(l)
                        new_c.add(new_r)
                elif j in resid_resname_segid_1[i].keys() and\
                     j not in resid_resname_segid_2[i].keys():
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info_1 = resid_resname_segid_1[i][j][k]
                        new_r = struct.Residue.Residue(res_info_1[0][0],res_info_1[0][1],res_info_1[0][2])
                        for l in strc_dic[i][j][res_info_1[0][0]]:
                            new_r.add(l)
                        new_c.add(new_r)
                new_m.add(new_c)
            self.new_s.add(new_m)
    
    def f06_0110(self, se1, se2):
        if len(se1.sel_dic['chains']) != len(se2.sel_dic['chains']):
            self.error_equal_entity_selection("Chain")
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        else:
            if se1.sel_dic['chains'] == se2.sel_dic['chains']:
                for i in self.strct.get_models():
                    for j in i.get_chains():
                        if j.get_id() in se1.sel_dic['chains']:
                            for k in j.get_residues():
                                if k.get_id()[1] in se2.sel_dic['residues']:
                                    self.atom_error_collision("Atom", str(k.get_id()),
                                                                      str(j.get_id()), 
                                                                      str(i.get_id()))
            resid_resname_segid_1 = {}
            resid_resname_segid_2 = {}
            strc_dic = {}
            for i in self.strct.get_models():
                resid_resname_segid_1[i.get_id()] = {}
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid_2[i.get_id()] = {}
                for j in i.get_chains():
                    resid_resname_segid_1[i.get_id()][j.get_id()] = {}
                    if j.get_id() in strc_dic[i.get_id()]:
                        pass
                    else:
                        strc_dic[i.get_id()][j.get_id()] = {}
                        resid_resname_segid_2[i.get_id()][j.get_id()] = {}
                    if j.get_id() in se1.sel_dic['chains']:
                        indx = se1.sel_dic['chains'].index(j.get_id())
                        chn2 = se2.sel_dic['chains'][indx]
                        if chn2 in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][chn2] = {}
                            resid_resname_segid_2[i.get_id()][chn2] = {}
                        for k in j.get_residues():
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id()[1] in se1.sel_dic['residues']:
                                indx = se1.sel_dic['residues'].index(k.get_id()[1])
                                res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])
                                if res2 in strc_dic[i.get_id()][chn2]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][chn2][res2] = []
                                    resid_resname_segid_2[i.get_id()][chn2][res2] = []
                                    resid_resname_segid_2[i.get_id()][chn2][res2].append((res2,k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][chn2][res2].append(l)
                            else:
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                    else:
                        for k in j.get_residues():
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))                   
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            # TODO the following loops are for assigning the right residue when
            # an atom is moved from one residue to another. In 0110, no single
            # atoms are remove, only all at once, so it is not necessary.
            # Check and simplify. It works, but seems unnecessary, but what if you add
            # a chain into another chain with repeated residues and no conflicting
            # atom numbers? Then, it is not uneccesary. Test!
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    if j in resid_resname_segid_1[i].keys() and\
                       j in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            if k in resid_resname_segid_1[i][j].keys() and\
                               k in resid_resname_segid_2[i][j].keys():
                                res_info_1 = resid_resname_segid_1[i][j][k]
                                res_info_2 = resid_resname_segid_2[i][j][k]
                                if res_info_1[0][1] == res_info_2[0][1]:
                                    new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                                else:
                                    new_r = struct.Residue.Residue(res_info_2[0][0],res_info_1[0][1],res_info_2[0][2])
                                for l in strc_dic[i][j][res_info_2[0][0]]:
                                    new_r.add(l)
                                new_c.add(new_r)
                            elif k not in resid_resname_segid_1[i][j].keys() and\
                                 k in resid_resname_segid_2[i][j].keys():
                                res_info_2 = resid_resname_segid_2[i][j][k]
                                new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                                for l in strc_dic[i][j][res_info_2[0][0]]:
                                    new_r.add(l)
                                new_c.add(new_r)
                            elif k in resid_resname_segid_1[i][j].keys() and\
                                 k not in resid_resname_segid_2[i][j].keys():
                                     print("WARNING: Anticipated situation not encountered before.")
                                     print("         Needs to be coded.")
                    elif j not in resid_resname_segid_1[i].keys() and\
                         j in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            res_info_2 = resid_resname_segid_2[i][j][k]
                            new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                            for l in strc_dic[i][j][res_info_2[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    elif j in resid_resname_segid_1[i].keys() and\
                         j not in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            res_info_1 = resid_resname_segid_1[i][j][k]
                            new_r = struct.Residue.Residue(res_info_1[0][0],res_info_1[0][1],res_info_1[0][2])
                            for l in strc_dic[i][j][res_info_1[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
                
    def f07_0111(self, se1, se2):
        if len(se1.sel_dic['chains']) != len(se2.sel_dic['chains']):
            self.error_equal_entity_selection("Chain")
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                for j in i.get_chains():
                    if j.get_id() in se2.sel_dic['chains']:
                        for k in j.get_residues():
                            if k.get_id()[1] in se2.sel_dic['residues']:
                                for l in k.get_atoms():
                                    if l.get_id() in se2.sel_dic['atoms']:
                                        self.atom_error_collision("Atom", str(l.get_id()), 
                                                                          str(k.get_id()),
                                                                          str(j.get_id()), 
                                                                          str(i.get_id()))
            resid_resname_segid_1 = {}
            resid_resname_segid_2 = {}
            strc_dic = {}
            for i in self.strct.get_models():
                resid_resname_segid_1[i.get_id()] = {}
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid_2[i.get_id()] = {}
                for j in i.get_chains():
                    resid_resname_segid_1[i.get_id()][j.get_id()] = {}
                    if j.get_id() in strc_dic[i.get_id()]:
                        pass
                    else:
                        strc_dic[i.get_id()][j.get_id()] = {}
                        resid_resname_segid_2[i.get_id()][j.get_id()] = {}
                    if j.get_id() in se1.sel_dic['chains']:
                        indx = se1.sel_dic['chains'].index(j.get_id())
                        chn2 = se2.sel_dic['chains'][indx]
                        if chn2 in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][chn2] = {}
                            resid_resname_segid_2[i.get_id()][chn2] = {}
                        for k in j.get_residues():
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id()[1] in se1.sel_dic['residues']:
                                indx = se1.sel_dic['residues'].index(k.get_id()[1])
                                res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])
                                if res2 in strc_dic[i.get_id()][chn2]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][chn2][res2] = []
                                    resid_resname_segid_2[i.get_id()][chn2][res2] = []
                                    # NOTE only add res2 to resid_resname_segid_2 if not already in there.
                                    #      This avoids duplicates that can make it harder to prevent
                                    #      that an atom added to an existing residue, with atoms already in it,
                                    #      renames a residue that already exist.
                                    resid_resname_segid_2[i.get_id()][chn2][res2].append((res2,k.resname,k.segid))
                                for l in k.get_atoms():
                                    if l.get_id() in se1.sel_dic['atoms']:
                                        mov_a = copy.deepcopy(l)
                                        mov_a.detach_parent()
                                        indx = se1.sel_dic['atoms'].index(l.get_id())
                                        atm2 = se2.sel_dic['atoms'][indx]
                                        mov_a.id = atm2
                                        mov_a.name = atm2
                                        mov_a.fullname = atm2
                                        strc_dic[i.get_id()][chn2][res2].append(mov_a)
                                    else:
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                            else:
                                # NOTE this is redundant
                                #if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                #    pass
                                #else:
                                #    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                #    resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                                #resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                    else:
                        # Note this is redundant
                        #if j.get_id() in strc_dic[i.get_id()]:
                        #    pass
                        #else:
                        #    strc_dic[i.get_id()][j.get_id()] = {}
                        #    resid_resname_segid_2[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid_1[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))                   
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid_2[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    if j in resid_resname_segid_1[i].keys() and\
                       j in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            if k in resid_resname_segid_1[i][j].keys() and\
                               k in resid_resname_segid_2[i][j].keys():
                                res_info_1 = resid_resname_segid_1[i][j][k]
                                res_info_2 = resid_resname_segid_2[i][j][k]
                                if res_info_1[0][1] == res_info_2[0][1]:
                                    new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                                else:
                                    new_r = struct.Residue.Residue(res_info_2[0][0],res_info_1[0][1],res_info_2[0][2])
                                for l in strc_dic[i][j][res_info_2[0][0]]:
                                    new_r.add(l)
                                new_c.add(new_r)
                            elif k not in resid_resname_segid_1[i][j].keys() and\
                                 k in resid_resname_segid_2[i][j].keys():
                                res_info_2 = resid_resname_segid_2[i][j][k]
                                new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                                for l in strc_dic[i][j][res_info_2[0][0]]:
                                    new_r.add(l)
                                new_c.add(new_r)
                            elif k in resid_resname_segid_1[i][j].keys() and\
                                 k not in resid_resname_segid_2[i][j].keys():
                                     print("WARNING: Anticipated situation not encountered before.")
                                     print("         Needs to be coded.")
                    elif j not in resid_resname_segid_1[i].keys() and\
                         j in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            res_info_2 = resid_resname_segid_2[i][j][k]
                            new_r = struct.Residue.Residue(res_info_2[0][0],res_info_2[0][1],res_info_2[0][2])
                            for l in strc_dic[i][j][res_info_2[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    elif j in resid_resname_segid_1[i].keys() and\
                         j not in resid_resname_segid_2[i].keys():
                        new_c = struct.Chain.Chain(j)
                        for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                            res_info_1 = resid_resname_segid_1[i][j][k]
                            new_r = struct.Residue.Residue(res_info_1[0][0],res_info_1[0][1],res_info_1[0][2])
                            for l in strc_dic[i][j][res_info_1[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
    # TODO f08_0111 to f15_1111 has old algorithm. 
    def f08_1000(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        for k in j.get_residues():
                            for l in k.get_atoms():
                                self.model_error_collision("Model", str(i.get_id()))
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in strc_dic[mdl2]:
                            pass
                        else:
                            strc_dic[mdl2][j.get_id()] = {}
                            resid_resname_segid[mdl2][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id() in strc_dic[mdl2][j.get_id()]:
                                pass
                            else:
                                strc_dic[mdl2][j.get_id()][k.get_id()] = []
                                resid_resname_segid[mdl2][j.get_id()][k.get_id()] = []
                            resid_resname_segid[mdl2][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[mdl2][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)

    def f09_1001(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        for k in j.get_residues():
                            for l in k.get_atoms():
                                if l.get_id() in se2.sel_dic['atoms']:
                                    self.atom_error_collision("Atom", str(l.get_id()), 
                                                                      str(k.get_id()),
                                                                      str(j.get_id()), 
                                                                      str(i.get_id()))
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in strc_dic[mdl2]:
                            pass
                        else:
                            strc_dic[mdl2][j.get_id()] = {}
                            resid_resname_segid[mdl2][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id() in strc_dic[mdl2][j.get_id()]:
                                pass
                            else:
                                strc_dic[mdl2][j.get_id()][k.get_id()] = []
                                resid_resname_segid[mdl2][j.get_id()][k.get_id()] = []
                            resid_resname_segid[mdl2][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                if l.get_id() in se1.sel_dic['atoms']:
                                    mov_a = copy.deepcopy(l)
                                    mov_a.detach_parent()
                                    indx = se1.sel_dic['atoms'].index(l.get_id())
                                    atm2 = se2.sel_dic['atoms'][indx]
                                    mov_a.id = atm2
                                    mov_a.name = atm2
                                    mov_a.fullname = atm2
                                    strc_dic[mdl2][j.get_id()][k.get_id()].append(mov_a)
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
    
    def f10_1010(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        for k in j.get_residues():
                            if k.get_id()[1] in se2.sel_dic['residues']:
                                for l in k.get_atoms():
                                    self.atom_error_collision("Atom", str(l.get_id()), 
                                                                      str(k.get_id()),
                                                                      str(j.get_id()), 
                                                                      str(i.get_id()))
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in strc_dic[mdl2]:
                            pass
                        else:
                            strc_dic[mdl2][j.get_id()] = {}
                            resid_resname_segid[mdl2][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id()[1] in se1.sel_dic['residues']:
                                indx = se1.sel_dic['residues'].index(k.get_id()[1])
                                res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])                    
                                if res2 in strc_dic[mdl2][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[mdl2][j.get_id()][res2] = []
                                    resid_resname_segid[mdl2][j.get_id()][res2] = []
                                resid_resname_segid[mdl2][j.get_id()][res2].append((res2,k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[mdl2][j.get_id()][res2].append(l)
                            else:
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
                    
    def f11_1011(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                    if i.get_id() in se2.sel_dic['models']:
                        for j in i.get_chains():
                            for k in j.get_residues():
                                if k.get_id()[1] in se2.sel_dic['residues']:
                                    for l in k.get_atoms():
                                        if l.get_id() in se2.sel_dic['atoms']:
                                            self.atom_error_collision("Atom", str(l.get_id()), 
                                                                              str(k.get_id()),
                                                                              str(j.get_id()), 
                                                                              str(i.get_id()))
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in strc_dic[mdl2]:
                            pass
                        else:
                            strc_dic[mdl2][j.get_id()] = {}
                            resid_resname_segid[mdl2][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            if k.get_id()[1] in se1.sel_dic['residues']:
                                indx = se1.sel_dic['residues'].index(k.get_id()[1])
                                res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])                    
                                if res2 in strc_dic[mdl2][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[mdl2][j.get_id()][res2] = []
                                    resid_resname_segid[mdl2][j.get_id()][res2] = []
                                resid_resname_segid[mdl2][j.get_id()][res2].append((res2,k.resname,k.segid))
                                for l in k.get_atoms():
                                    if l.get_id() in se1.sel_dic['atoms']:
                                        mov_a = copy.deepcopy(l)
                                        mov_a.detach_parent()
                                        indx = se1.sel_dic['atoms'].index(l.get_id())
                                        atm2 = se2.sel_dic['atoms'][indx]
                                        mov_a.id = atm2
                                        mov_a.name = atm2
                                        mov_a.fullname = atm2
                                        strc_dic[mdl2][j.get_id()][res2].append(mov_a)
                                    else:
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                            else:
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)

    def f12_1100(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['chains']) != len(se2.sel_dic['chains']):
            self.error_equal_entity_selection("Chain")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        if j.get_id() in se2.sel_dic['chains']:
                            for k in j.get_residues():
                                for l in k.get_atoms():
                                    self.atom_error_collision("Atom", str(l.get_id()), 
                                                                      str(k.get_id()),
                                                                      str(j.get_id()), 
                                                                      str(i.get_id()))
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in se1.sel_dic['chains']:
                            indx = se1.sel_dic['chains'].index(j.get_id())
                            chn2 = se2.sel_dic['chains'][indx]
                            if chn2 in strc_dic[mdl2]:
                                pass
                            else:
                                strc_dic[mdl2][chn2] = {}
                                resid_resname_segid[mdl2][chn2] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                if k.get_id() in strc_dic[mdl2][chn2]:
                                    pass
                                else:
                                    strc_dic[mdl2][chn2][k.get_id()] = []
                                    resid_resname_segid[mdl2][chn2][k.get_id()] = []
                                resid_resname_segid[mdl2][chn2][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[mdl2][chn2][k.get_id()].append(l)
                        else:
                            if j.get_id() in strc_dic[i.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()] = {}
                                resid_resname_segid[i.get_id()][j.get_id()] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
    
    def f13_1101(self, se1, se2):
        if se1.sel_dic['models'] == se2.sel_dic['models']:
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['chains']) != len(se2.sel_dic['chains']):
            self.error_equal_entity_selection("Chain")
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        if j.get_id() in se2.sel_dic['chains']:
                            for k in j.get_residues():
                                for l in k.get_atoms():
                                    if l.get_id()[1] in se2.sel_dic['atoms']:
                                        self.atom_error_collision("Atom", str(l.get_id()), 
                                                                          str(k.get_id()),
                                                                          str(j.get_id()), 
                                                                          str(i.get_id()))
            # Ok so far
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in se1.sel_dic['chains']:
                            indx = se1.sel_dic['chains'].index(j.get_id())
                            chn2 = se2.sel_dic['chains'][indx]
                            if chn2 in strc_dic[mdl2]:
                                pass
                            else:
                                strc_dic[mdl2][chn2] = {}
                                resid_resname_segid[mdl2][chn2] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                if k.get_id() in strc_dic[mdl2][chn2]:
                                    pass
                                else:
                                    strc_dic[mdl2][chn2][k.get_id()] = []
                                    resid_resname_segid[mdl2][chn2][k.get_id()] = []
                                resid_resname_segid[mdl2][chn2][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    if l.get_id() in se1.sel_dic['atoms']:
                                        mov_a = copy.deepcopy(l)
                                        mov_a.detach_parent()
                                        indx = se1.sel_dic['atoms'].index(l.get_id())
                                        atm2 = se2.sel_dic['atoms'][indx]
                                        mov_a.id = atm2
                                        mov_a.name = atm2
                                        mov_a.fullname = atm2
                                        strc_dic[mdl2][chn2][k.get_id()].append(mov_a)
                                    else:
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                        else:
                            if j.get_id() in strc_dic[i.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()] = {}
                                resid_resname_segid[i.get_id()][j.get_id()] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            # Remove empty entities
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)

    def f14_1110(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['chains']) != len(se2.sel_dic['chains']):
            self.error_equal_entity_selection("Chain")
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        if j.get_id() in se2.sel_dic['chains']:
                            for k in j.get_residues():
                                if k.get_id()[1] in se2.sel_dic['residues']:
                                    for l in k.get_atoms():
                                        self.atom_error_collision("Atom", str(l.get_id()), 
                                                                          str(k.get_id()),
                                                                          str(j.get_id()), 
                                                                          str(i.get_id()))
            # Ok so far
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in se1.sel_dic['chains']:
                            indx = se1.sel_dic['chains'].index(j.get_id())
                            chn2 = se2.sel_dic['chains'][indx]
                            if chn2 in strc_dic[mdl2]:
                                pass
                            else:
                                strc_dic[mdl2][chn2] = {}
                                resid_resname_segid[mdl2][chn2] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                if k.get_id()[1] in se1.sel_dic['residues']:
                                    indx = se1.sel_dic['residues'].index(k.get_id()[1])
                                    res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])
                                    if res2 in strc_dic[mdl2][chn2]:
                                        pass
                                    else:
                                        strc_dic[mdl2][chn2][res2] = []
                                        resid_resname_segid[mdl2][chn2][res2] = []
                                    resid_resname_segid[mdl2][chn2][res2].append((res2,k.resname,k.segid))
                                    for l in k.get_atoms():
                                        strc_dic[mdl2][chn2][res2].append(l)
                                else:
                                    if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                        pass
                                    else:
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                        resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                    for l in k.get_atoms():
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                        else:
                            if j.get_id() in strc_dic[i.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()] = {}
                                resid_resname_segid[i.get_id()][j.get_id()] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            # Remove empty entities
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)

    def f15_1111(self, se1, se2):
        if len(se1.sel_dic['models']) != len(se2.sel_dic['models']):
            self.error_equal_entity_selection("Model")
        elif len(se1.sel_dic['chains']) != len(se2.sel_dic['chains']):
            self.error_equal_entity_selection("Chain")
        elif len(se1.sel_dic['residues']) != len(se2.sel_dic['residues']):
            self.error_equal_entity_selection("Residue")
        elif len(se1.sel_dic['atoms']) != len(se2.sel_dic['atoms']):
            self.error_equal_entity_selection("Atom")
        else:
            for i in self.strct.get_models():
                if i.get_id() in se2.sel_dic['models']:
                    for j in i.get_chains():
                        if j.get_id() in se2.sel_dic['chains']:
                            for k in j.get_residues():
                                if k.get_id()[1] in se2.sel_dic['residues']:
                                    for l in k.get_atoms():
                                        if l.get_id() in se2.sel_dic['atoms']:
                                            self.atom_error_collision("Atom", str(l.get_id()), 
                                                                              str(k.get_id()),
                                                                              str(j.get_id()), 
                                                                              str(i.get_id()))
            resid_resname_segid = {}
            strc_dic = {}
            for i in self.strct.get_models():
                if i.get_id() in strc_dic:
                    pass
                else:
                    strc_dic[i.get_id()] = {}
                    resid_resname_segid[i.get_id()] = {}
                if i.get_id() in se1.sel_dic['models']:
                    indx = se1.sel_dic['models'].index(i.get_id())
                    mdl2 = se2.sel_dic['models'][indx]
                    if mdl2 in strc_dic:
                        pass
                    else:
                        strc_dic[mdl2] = {}
                        resid_resname_segid[mdl2] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        if j.get_id() in se1.sel_dic['chains']:
                            indx = se1.sel_dic['chains'].index(j.get_id())
                            chn2 = se2.sel_dic['chains'][indx]
                            if chn2 in strc_dic[mdl2]:
                                pass
                            else:
                                strc_dic[mdl2][chn2] = {}
                                resid_resname_segid[mdl2][chn2] = {}
                            for k in j.get_residues():
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                if k.get_id()[1] in se1.sel_dic['residues']:
                                    indx = se1.sel_dic['residues'].index(k.get_id()[1])
                                    res2 = (k.get_id()[0],se2.sel_dic['residues'][indx],k.get_id()[2])
                                    if res2 in strc_dic[mdl2][chn2]:
                                        pass
                                    else:
                                        strc_dic[mdl2][chn2][res2] = []
                                        resid_resname_segid[mdl2][chn2][res2] = []
                                    resid_resname_segid[mdl2][chn2][res2].append((res2,k.resname,k.segid))
                                    for l in k.get_atoms():
                                        if l.get_id() in se1.sel_dic['atoms']:
                                            mov_a = copy.deepcopy(l)
                                            mov_a.detach_parent()
                                            indx = se1.sel_dic['atoms'].index(l.get_id())
                                            atm2 = se2.sel_dic['atoms'][indx]
                                            mov_a.id = atm2
                                            mov_a.name = atm2
                                            mov_a.fullname = atm2
                                            strc_dic[mdl2][chn2][res2].append(mov_a)
                                        else:
                                            strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                                else:
                                    if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                        pass
                                    else:
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                        resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                    for l in k.get_atoms():
                                        strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                        else:
                            if j.get_id() in strc_dic[i.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()] = {}
                                resid_resname_segid[i.get_id()][j.get_id()] = {}
                            for k in j.get_residues():                       
                                if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                    pass
                                else:
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                    resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                                for l in k.get_atoms():
                                    strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
                else:
                    if i.get_id() in strc_dic:
                        pass
                    else:
                        strc_dic[i.get_id()] = {}
                        resid_resname_segid[i.get_id()] = {}
                    for j in i.get_chains():
                        if j.get_id() in strc_dic[i.get_id()]:
                            pass
                        else:
                            strc_dic[i.get_id()][j.get_id()] = {}
                            resid_resname_segid[i.get_id()][j.get_id()] = {}
                        for k in j.get_residues():
                            if k.get_id() in strc_dic[i.get_id()][j.get_id()]:
                                pass
                            else:
                                strc_dic[i.get_id()][j.get_id()][k.get_id()] = []
                                resid_resname_segid[i.get_id()][j.get_id()][k.get_id()] = []
                            resid_resname_segid[i.get_id()][j.get_id()][k.get_id()].append((k.get_id(),k.resname,k.segid))
                            for l in k.get_atoms():
                                strc_dic[i.get_id()][j.get_id()][k.get_id()].append(l)
            # Remove empty entities
            for ii in sorted(strc_dic.keys()):
                for jj in sorted(strc_dic[ii].keys()):
                    for kk in sorted([ss for ss in strc_dic[ii][jj].keys()]):
                        if len(strc_dic[ii][jj][kk]) == 0:
                            del(strc_dic[ii][jj][kk])
                    if len(strc_dic[ii][jj]) == 0:
                        del(strc_dic[ii][jj])
                if len(strc_dic[ii]) == 0:
                    del(strc_dic[ii])
            new_m = struct.Model.Model(i.get_id())
            for i in sorted(strc_dic.keys()):
                new_m = struct.Model.Model(i)
                for j in sorted(strc_dic[i].keys()):
                    new_c = struct.Chain.Chain(j)
                    for k in sorted([ss for ss in strc_dic[i][j].keys()]):
                        res_info = resid_resname_segid[i][j][k]
                        if len(res_info) == 0:
                            pass
                        else:
                            new_r = struct.Residue.Residue(res_info[0][0],res_info[0][1],res_info[0][2])
                            for l in strc_dic[i][j][res_info[0][0]]:
                                new_r.add(l)
                            new_c.add(new_r)
                    new_m.add(new_c)
                self.new_s.add(new_m)
                    
    def check_renaming_expressions(self, se1, se2):
        # Check that new model names do not exist
        all_models_ast = False
        all_chains_ast = False
        all_residu_ast = False
        all_atoms_ast = False
        if se1.sel_dic['models'][0] == '*' and se2.sel_dic['models'][0] == '*':
            all_models_ast = True
        elif se1.sel_dic['models'][0] == '*' and se2.sel_dic['models'][0] != '*':
            self.error_missmatch_ast_nonast("Model")
        elif se1.sel_dic['models'][0] != '*' and se2.sel_dic['models'][0] == '*':
            self.error_missmatch_nonast_ast("Model")
        else:
            # check_models(strct, se1, se2)
            pass
         
        # Check that new chain names do not exist
        if se1.sel_dic['chains'][0] == '*' and se2.sel_dic['chains'][0] == '*':
            all_chains_ast = True
        elif se1.sel_dic['chains'][0] == '*' and se2.sel_dic['chains'][0] != '*':
            self.error_missmatch_ast_nonast("Chain")
        elif se1.sel_dic['chains'][0] != '*' and se2.sel_dic['chains'][0] == '*':
            self.error_missmatch_nonast_ast("Chain")
        else:
            # check_chains(strct, se1, se2) 
            pass
                
        # Check that new residue names do not exist
        if se1.sel_dic['residues'][0] == '*' and se2.sel_dic['residues'][0] == '*':
            all_residu_ast = True
        elif se1.sel_dic['residues'][0] == '*' and se2.sel_dic['residues'][0] != '*':
            self.error_missmatch_ast_nonast("Residue")
        elif se1.sel_dic['residues'][0] != '*' and se2.sel_dic['residues'][0] == '*':
            self.error_missmatch_nonast_ast("Residue")
        else:
            pass
            #check_residues(strct, se1, se2)
        
        # Check that new atom names do not exist
        if se1.sel_dic['atoms'][0] == '*' and se2.sel_dic['atoms'][0] == '*':
            all_atoms_ast = True
        elif se1.sel_dic['atoms'][0] == '*' and se2.sel_dic['atoms'][0] != '*':
            self.error_missmatch_ast_nonast("Atom")
        elif se1.sel_dic['atoms'][0] != '*' and se2.sel_dic['atoms'][0] == '*':
            self.error_missmatch_nonast_ast("Atom")
        else:
            pass
        if all_models_ast and all_chains_ast and all_residu_ast and all_atoms_ast:
            print("Error: There must be at least an entity that is not an asterisc '*'.")
            print("       A structural expression such as 'm[*]c[*]r[*]a[*]' can't ")
            print("       be used for both entities to be rename and new entity names.")
            print("       It would be like reneaming all entities to their current name.")
            print("       This is valid but pointless, and it's flagged as an error FYI.")
            sys.exit(1)
        elif all_models_ast and all_chains_ast and all_residu_ast and not all_atoms_ast: # 0001
            self.f01_0001(se1, se2)
        elif all_models_ast and all_chains_ast and not all_residu_ast and all_atoms_ast: # 0010
            self.f02_0010(se1, se2)
        elif all_models_ast and all_chains_ast and not all_residu_ast and not all_atoms_ast: # 0011
            self.f03_0011(se1, se2)
        elif all_models_ast and not all_chains_ast and all_residu_ast and all_atoms_ast: # 0100
            self.f04_0100(se1, se2)
        elif all_models_ast and not all_chains_ast and all_residu_ast and not all_atoms_ast: # 0101
            self.f05_0101(se1, se2)
        elif all_models_ast and not all_chains_ast and not all_residu_ast and all_atoms_ast: # 0110
            self.f06_0110(se1, se2)
        elif all_models_ast and not all_chains_ast and not all_residu_ast and not all_atoms_ast: # 0111
            self.f07_0111(se1, se2)
        elif not all_models_ast and all_chains_ast and all_residu_ast and all_atoms_ast: # 1000
            self.f08_1000(se1, se2)
        elif not all_models_ast and all_chains_ast and all_residu_ast and not all_atoms_ast: #  1001
            self.f09_1001(se1, se2)
        elif not all_models_ast and all_chains_ast and not all_residu_ast and all_atoms_ast: # 1010
            self.f10_1010(se1, se2)
        elif not all_models_ast and all_chains_ast and not all_residu_ast and not all_atoms_ast: # 1011
            self.f11_1011(se1, se2)
        elif not all_models_ast and not all_chains_ast and all_residu_ast and all_atoms_ast: # 1100
            self.f12_1100(se1, se2)
        elif not all_models_ast and not all_chains_ast and all_residu_ast and not all_atoms_ast: # 1101
            self.f13_1101(se1, se2)
        elif not all_models_ast and not all_chains_ast and not all_residu_ast and all_atoms_ast: # 1110
            self.f14_1110(se1, se2)
        elif not all_models_ast and not all_chains_ast and not all_residu_ast and not all_atoms_ast: # 1111
            self.f15_1111(se1, se2)
            
    def output_strct(self, path):
        io = PDBIO()
        io.set_structure(self.strct)
        io.save(path)
        
    def output_new_s(self, path):
        io = PDBIO()
        io.set_structure(self.new_s)
        io.save(path)
