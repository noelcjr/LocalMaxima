# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 23:26:00 2016
@author: noel
"""
import sys
import os
import Bio.PDB as struct
import LocalMaxima.utilities as utl
#import string
import LocalMaxima.structure.MolecularRigidManipulations as MRM
from LocalMaxima.structure.XYZ_Formats import crd
from jinja2 import Environment, PackageLoader
from LocalMaxima.ThirdPartySoftware import RunInputs
import LocalMaxima.utilities.ColRows as CR
import LocalMaxima.sequence.sequence as SQ

class GenInputs(object):
    def __init__(self):
        self.run = RunInputs.RunCharmm()
        self.env = Environment(loader=PackageLoader('LocalMaxima.ThirdPartySoftware.charmm',
                                                    'templates'))
    def print_template(self,outFile,template):
        for i in template:
            outFile.write(i)
    def gen_parameters(self):
        template = self.env.get_template('parameters.template')
        return template.render()
    # TODO Paths in templates must be setup during installation to point in the
    #      right direction. 
    # param_absolute_path="/home//nnoel/Code/EntropyMaxima/em/params/charmm27.ff/")
    
    def stream_nina_radii(self):
        template = self.env.get_template('stream_nina_radii.template')
        return template.render()
    # param_absolute_path="/home//nnoel/Code/EntropyMaxima/em/params/")
    
    def fix_sequence(self,chain_info):
        template = self.env.get_template('fix_sequence.template')
        return template.render(
               a_sequence=chain_info)
    
    def disulfide_patch(self,disu_bond):
        template = self.env.get_template('disulfides.template')
        return template.render(
               chain_num_chain_num=disu_bond)
    
    def wmain(self,value):
        template = self.env.get_template('wmain.template')
        return template.render(val=value)
    
    def write_psf(self,file_name):
        template = self.env.get_template('write_psf.template')
        return template.render(out_filename=file_name)
    
    def write_psf_xplor(self,file_name):
        template = self.env.get_template('write_psf_xplor.template')
        return template.render(out_filename=file_name)
    
    def read_psf(self,file_name):
        template = self.env.get_template('read_psf.template')
        return template.render(in_filename=file_name)
    
    def read_psf2(self,file_name):
        template = self.env.get_template('read_psf2.template')
        return template.render(in_filename=file_name)
    
    def write_pdb(self,file_name):
        template = self.env.get_template('write_pdb.template')
        return template.render(out_filename=file_name)
    
    def write_ic(self,file_name):
        template = self.env.get_template('write_ic.template')
        return template.render(out_filename=file_name)
    
    def write_crd(self,file_name):
        template = self.env.get_template("write_crd.template")
        return template.render(out_filename=file_name)
    
    def write_output(self,file_name):
        template = self.env.get_template("output_file.template")
        return template.render(out_filename=file_name)
    
    def read_crd(self,file_name):
        template = self.env.get_template("read_crd.template")
        return template.render(in_filename=file_name)
    
    def read_crd2(self,file_name):
        template = self.env.get_template("read_crd2.template")
        return template.render(in_filename=file_name)
    
    def build_heavy_atoms(self):
        template = self.env.get_template("build_heavy.template")
        return template.render()
    
    def translate_xyz(self, chains, zdir):
        template = self.env.get_template("translate_xyz.template")
        return template.render(chain=chains, z=zdir)
    
    def build_hydrogens(self):
        template = self.env.get_template("build_hydrogens.template")
        return template.render()
    
    def mmgbsa_ca(self):
        template = self.env.get_template("mmgbsa_ca.template")
        return template.render()
    
    def mmgbsa(self):
        template = self.env.get_template("mmgbsa.template")
        return template.render()

    def mmgbsa2(self, chains, zdir):
        template = self.env.get_template("mmgbsa2.template")
        return template.render(chain=chains, z=zdir)
    
    def mmgbsa_ca_z(self):
        template = self.env.get_template("mmgbsa_ca_z.template")
        return template.render()
     
    def stop(self):
        template = self.env.get_template("stop.template")
        return template.render()
    
    def bomlev(self, bomlev):
        template = self.env.get_template("bomlev.template")
        return template.render(lev=bomlev)
    
    def minimization1(self,outfile,input_list):
        template = self.env.get_template("minimization1.template")
        return template.render(
            OUTFILE=outfile,
            A=input_list[0][0],
            ROTA="resid "+input_list[0][1]+" : "+input_list[0][2],
            LINA="resid "+input_list[0][3]+" : "+input_list[0][4],
            CENA="resid "+input_list[0][5]+" : "+input_list[0][6],
            B = input_list[1][0],
            ROTB="resid "+input_list[1][1]+" : "+input_list[1][2],
            LINB="resid "+input_list[1][3]+" : "+input_list[1][4],
            CENB="resid "+input_list[1][5]+" : "+input_list[1][6])
    
    def gen_protein(self,Chn,Name):
        template = self.env.get_template("gen_sequence.template")
        return template.render(chn=Chn,name=Name)
    
    def gen_protein_ter(self,Chn,Name,First='first none',Last='last none'):
        template = self.env.get_template("gen_sequence_ter.template")
        return template.render(chn=Chn,name=Name,first=First,last=Last)
    
    def gen_setup_one(self):
        template = self.env.get_template('setup_one.inp.template')
        a_pair = {
            'id': 'A',
            'seq': "A.SEQ",
            'first': 'none',
            'last': 'none',
            'inp': "A_FIXRES.INP"
        }
        b_pair = {
            'id': 'B',
            'seq': "B.SEQ",
            'first': 'none',
            'last': 'none',
            'inp': "B_FIXRES.INP"
        }
    
        return template.render(
            param_absolute_path='/home/noel/Code/EntropyMaxima/em/params/charmm27.ff/',
            a_sequence=a_pair,
            b_sequence=b_pair,
            in_filename='INFILE',
            out_filename='OUTFILE')
    
    def get_folding_mmgbsa(self, options, dir_path = ''):
        """ Gets folding MMGBSA without Component ANalysis (i.e. --CA)
            It runs a generated CHARMM script.
        """
        xyz_name = os.path.basename(options.xyz).split('.')
        psf_name = os.path.basename(options.psf).split('.')[0]
        outFile = open(dir_path+"mmgbsa.inp", 'w')
        
    def get_binding_mmgbsa(self, options, dir_path = ''):
        """ Gets binding MMGBSA without Component ANalysis (i.e. --CA)
            It runs a generated CHARMM script.
        """
        xyz_name = os.path.basename(options.xyz).split('.')
        psf_name = os.path.basename(options.psf).split('.')[0]
        outFile = open(dir_path+"mmgbsa.inp", 'w')
        if options.mov == '':
            print("ERROR: Chain to move is an empty string. Exiting Now.")
            sys.exit(1)
        str1 = "segid "
        first_chain = True
        chains = options.mov.strip().split(',')
        for i in chains:
            if first_chain:
                str1 += i
                first_chain = False
            else:
                str1 += " .or. segid "+i
        self.print_template(outFile,self.gen_parameters().split('/n'))
        self.bomlev('-5')
        self.print_template(outFile,self.read_psf(psf_name.lower()))
        self.print_template(outFile,self.write_output(xyz_name[0].lower()))
        if xyz_name[1] == 'crd':
            self.print_template(outFile,self.read_crd(xyz_name[0].lower()))
        elif xyz_name[1] == 'pdb':
            pass
        else:
            print("ERROR: only pdb ir crd xoordinate files allowed.")
            sys.exit(1)
        self.print_template(outFile,self.stream_nina_radii())
        self.print_template(outFile,self.mmgbsa())
        self.print_template(outFile,self.mmgbsa2(str1,options.trans))
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm("mmgbsa.inp","mmgbsa.out")
        os.system("rm mmgbsa.inp mmgbsa.out")

    def get_gbsa_folding(self, options, dir_path = ''):
        """
        Gets the GB and SA valiues for atoms. SA values only for Heavy, non
        hydrogen, atoms. Files are output to files.

        Parameters
        ----------
        options : TYPE
            DESCRIPTION.
        dir_path : TYPE, optional
            DESCRIPTION. The default is ''.

        Returns
        -------
        None.

        """
        # NOTE: Z0 is folded gbsa.
        #       Z500 is unfolded gbsa.
        xyz_name = os.path.basename(options.xyz).split('.')
        psf_name = os.path.basename(options.psf).split('.')[0]
        outFile = open(dir_path+"gbsa_z0.inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        self.bomlev('-1')
        self.print_template(outFile,self.read_psf(psf_name.lower()))
        if xyz_name[1] == 'crd':
            self.print_template(outFile,self.read_crd(xyz_name[0].lower()))
        elif xyz_name == 'pdb':
            pass
        else:
            print("ERROR: only pdb or crd xoordinate files allowed.")
            sys.exit(1)
        self.print_template(outFile,self.stream_nina_radii())
        self.print_template(outFile,self.mmgbsa_ca())
        #self.print_template(outFile,self.write_psf_xplor('temp_abcxyz'))
        self.print_template(outFile,self.write_pdb('temp'))
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm("gbsa_z0.inp","gbsa_z0.out")
        os.system("mv fort.90 sasa_z0.txt")
        os.system("mv fort.91 gb_z0.txt")
        
        outFile = open(dir_path+"gbsa_z500.inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        self.bomlev('-1')
        self.print_template(outFile,self.read_psf(psf_name.lower()))
        if xyz_name[1] == 'crd':
            self.print_template(outFile,self.read_crd(xyz_name[0].lower()))
        elif xyz_name == 'pdb':
            pass
        else:
            print("ERROR: only pdb or crd xoordinate files allowed.")
            sys.exit(1)
        self.print_template(outFile,self.stream_nina_radii())
        self.print_template(outFile,self.mmgbsa_ca())
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm("gbsa_z500.inp","gbsa_z500.out")
        os.system("mv fort.90 sasa_z500.txt")
        os.system("mv fort.91 gb_z500.txt")
    
    def get_gbsa_binding(self, options, dir_path = ''):
        """
        Gets the GB and SA valiues for atoms. SA values only for Heavy, non
        hydrogen, atoms. Files are output to files.

        Parameters
        ----------
        options : TYPE
            DESCRIPTION.
        dir_path : TYPE, optional
            DESCRIPTION. The default is ''.

        Returns
        -------
        None.

        """
        xyz_name = os.path.basename(options.xyz).split('.')
        psf_name = os.path.basename(options.psf).split('.')[0]
        outFile = open(dir_path+"gbsa_z0.inp", 'w')
        if options.mov == '':
            print("ERROR: Chain to move is an empty string. Exiting Now.")
            sys.exit(1)
        str1 = "segid "
        first_chain = True
        chains = options.mov.strip().split(',')
        for i in chains:
            if first_chain:
                str1 += i
                first_chain = False
            else:
                str1 += " .or. segid "+i
        self.print_template(outFile,self.gen_parameters().split('/n'))
        self.bomlev('-1')
        self.print_template(outFile,self.read_psf(psf_name.lower()))
        if xyz_name[1] == 'crd':
            self.print_template(outFile,self.read_crd(xyz_name[0].lower()))
        elif xyz_name == 'pdb':
            pass
        else:
            print("ERROR: only pdb ir crd xoordinate files allowed.")
            sys.exit(1)
        self.print_template(outFile,self.stream_nina_radii())
        self.print_template(outFile,self.mmgbsa_ca())
        self.print_template(outFile,self.write_psf_xplor('temp'))
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm("gbsa_z0.inp","gbsa_z0.out")
        os.system("mv fort.90 sasa_z0.txt")
        os.system("mv fort.91 gb_z0.txt")
        
        outFile = open(dir_path+"gbsa_z500.inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        self.bomlev('-1')
        self.print_template(outFile,self.read_psf(psf_name.lower()))
        if xyz_name[1] == 'crd':
            self.print_template(outFile,self.read_crd(xyz_name[0].lower()))
        elif xyz_name == 'pdb':
            pass
        else:
            print("ERROR: only pdb ir crd xoordinate files allowed.")
            sys.exit(1)
        self.print_template(outFile,self.translate_xyz(str1,options.trans))
        self.print_template(outFile,self.stream_nina_radii())
        self.print_template(outFile,self.mmgbsa_ca())
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm("gbsa_z500.inp","gbsa_z500.out")
        os.system("mv fort.90 sasa_z500.txt")
        os.system("mv fort.91 gb_z500.txt")
    
    def prepare_pdb_for_charmm(self, pdbin, prot_ends, dir_path = ''):
        if not os.path.exists(pdbin):
            print("Error: Invalid path in --structure option.")
            sys.exit(1)
        if (os.path.basename(pdbin).split('.')[1] != 'pdb'):
            print("""ERROR: A PDB file with lower case 'pdb' suffix is required.
                  Program will exit without results.""")
            sys.exit(1)
        input_name = os.path.basename(pdbin).split('.')[0]
        if dir_path == '':
            output_name = input_name.lower()
        else:
            if dir_path[-1] == '/':
                pass
            else:
                dir_path += '/'
            output_name = dir_path+input_name.lower()
        reduceopt="-HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000"
        os.system("reduce {} {} > {} &".format(reduceopt,pdbin,output_name+"r.pdb"))
        print("Ran Reduce\n")
        if not os.path.exists(output_name+"r.pdb"):
            print("Error: OS System call to run reduce didn't output a PDB file.")
            sys.exit(1)
        print("Passed\n")
        pdb_parser = struct.PDBParser(QUIET=True)
        self.sp = pdb_parser.get_structure(input_name, pdbin)
        count_models = 0
        #for i in self.structure.get_models():
        for i in self.sp.get_models():
            count_models += 1
        if count_models != 1:
            print("ERROR: This function works on a PDB file with only one model.")
            print("       Extract PDBs from CIF if needed. Simulations can only be done")
            print("       from a single structure and not multiple structure models.")
            sys.exit(1)
        self.count_ions = 0
        ''' Some amino acids can come with alternative names.'''
        self.num_residues = {}
        self.check_unusual_names()
        '''
         This fixes when there are multiple chains and the residue number gets
         reset at the beginning of each chain. With this fix, residue numbers
         will be renumbered
        '''
        coor = crd("")
        coor.write_crd(self.sp,output_name)
        # TODO we need way to check prot_ends match number of chains in 
        #      structure if ends are already present
        self.chains_info = {}
        # Gen SEQ files for each chain
        self.gen_seq(dir_path, prot_ends)
        # Gen FIXERS to keep offset amino acid numbers
        self.gen_fixers(dir_path)
        print(dir_path+"setup_"+input_name+".inp")
        outFile = open(dir_path+"setup_"+input_name+".inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        for i in self.sp.get_chains():
            self.print_template(outFile,self.fix_sequence(self.chains_info[i.id]))
        self.print_template(outFile,self.wmain('charge'))
        self.print_template(outFile,self.write_psf(input_name.lower()+'r'))
        self.print_template(outFile,self.write_ic(input_name.lower()+'r'))
        self.print_template(outFile,self.read_crd(input_name.lower()))
        self.print_template(outFile,self.build_heavy_atoms())
        self.print_template(outFile,self.build_hydrogens())
        self.print_template(outFile,self.wmain('charge'))
        self.print_template(outFile,self.write_crd(output_name+'r'))
        self.print_template(outFile,self.write_pdb(output_name+'r'))
        self.print_template(outFile,self.write_psf_xplor(output_name+'r'))
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm(dir_path+"setup_"+input_name+".inp",dir_path+"setup_"+input_name+".out")
        fixPPDB = CR.ColumnsRows(output_name+"r.pdb")
        fixPPDB.fix_pdb_from_CHARMM()
    
    def gen_unfolded_xtruct_file_wrapper(self, options):
        self.gen_unfolded_xtruct_file(options.file, options.name, options.terminals)
        
    def gen_unfolded_xtruct_file(self, options_file, options_name,options_terminals=''):
        self.sequences = {}
        if (os.path.basename(options_file).split('.')[1].lower() == 'pdb'):
            self.file_not_found_ERROR(options_file, "--file")
            pdb_parser = struct.PDBParser(QUIET=True)
            self.sp = pdb_parser.get_structure("For Seq", options_file)
            for i in self.sp.get_models():
                if i.id == 0:
                    for j in i.get_chains():
                        self.sequences[j.get_id()] = []
                        for k in j.get_residues():
                            if k.get_resname() == 'HIS':
                                self.sequences[j.get_id()].append('HSE')
                            else:
                                self.sequences[j.get_id()].append(k.get_resname())
            for i in self.sequences.keys():
                print("Generating unfolded structure for "+i+": ", self.sequences[i])
                print("in "+options_name+i.lower()+"_u.pdb")
            count_models = 0
            for i in self.sp.get_models():
                count_models += 1
            if count_models > 1:
                print("WARNING: Multiple models in the structure. Only one unfolded")
                print("         structure per chain will be generated.")
            ''' Some amino acids can come with alternative names.'''
            self.num_residues = {}
            self.check_unusual_names()
            name = options_name.lower()
            first_model = True
            ter_dict = self.process_terminals(options_terminals)
            for i in self.sp.get_models():
                if first_model:
                    for j in i.get_chains():
                        j_lbl = str(j.get_id()).lower()+"_u"
                        l = len(self.sequences[j.get_id()])
                        self.gen_unfolded_charmm(name+j_lbl,l,j.get_id(),ter_dict)
                    first_model = False
        else:
            self.file_not_found_ERROR(options_file, "--file")
            a = SQ.sequence()
            a.add_Fasta_sequence(options_file)
            if len(a.chain_sequences) == 0:
                print("ERROR: File did not contain any sequences.")
                sys.exit(1)
            for i in a.chain_sequences.keys():
                chain_id = i.split('|')[0].split(':')[1]
                self.sequences[chain_id] = []
                for j in a.chain_sequences[i]:
                    for k in j:
                        if j == 'H':
                            self.sequences[chain_id].append('HSE')
                        else:
                            self.sequences[chain_id].append(utl.AA_Info.residueDict1_1[j])
            ter_dict = self.process_terminals(options_terminals)
            for i in self.sequences.keys():
               name = options_name.lower()+i.lower()+"_u"
               print("Generating unfolded structure for "+i+": ", self.sequences[i])
               print("in "+name)
               l = len(self.sequences[i])
               self.gen_unfolded_charmm(name,l,i,ter_dict)
        os.system("rm *.inp *.out")

    def gen_unfolded_xtruct_seq_wrapper(self, options):
        if not options.chn_id:
            print("ERROR: You need a value for --chn_id")
            print("       when --seq option is used.")
            sys.exit(1)
        self.gen_unfolded_xtruct_seq(options.seq, options.name, 
                                     options.chn_id, options.terminals)
        
    def gen_unfolded_xtruct_seq(self, options_seq, options_name, 
                                options_chn_id, options_terminals):
        self.sequences = {}
        # A chain id is assigned: A
        self.sequences[options_chn_id] = []
        if options_seq.find(',') == -1:
            # Either a three leter single AA or sequence in short notation
            seq_len = len(options_seq)
            if seq_len > 3:
                for i in options_seq:
                    if i.upper() not in utl.AA_Info.aa:
                        print("ERROR: "+str(i)+" Amino acid not found in parameters")
                        sys.exit(1)
                    else:                            
                        if utl.AA_Info.residueDict1_1[i] == 'HIS':
                            self.sequences[options_chn_id].append('HSE')
                        else:
                            self.sequences[options_chn_id].append(utl.AA_Info.residueDict1_1[i])
            elif seq_len == 3:
                if options_seq.upper() in utl.AA_Info.residueDict1_2:
                    if options_seq.upper() == 'HIS':
                        self.sequences[options_chn_id].append('HSE')
                    else:
                        self.sequences[options_chn_id].append(options_seq.upper())
                else:
                    print("ERROR: "+str(options_seq.upper())+" Amino acid not\
                          found in parameters")
                    sys.exit(1)
            else:
                if len(options_seq) <= 2:
                    print("ERROR: Need to enter a three letter identifier")
                    print("       for sequences of three amino acids or")
                    print("       in single letter format.")
                    sys.exit(1)       
        else:            
            # check amino acids are valid
            for i in options_seq.split(','):
                if i.upper() not in utl.AA_Info.residueDict1_2:
                    print("ERROR: Amino acid not found in parameters")
                    sys.exit(1)
                else:
                    if i.upper() == 'HIS':
                        self.sequences[options_chn_id].append('HSE')
                    else:
                        self.sequences[options_chn_id].append(i.upper())
        name = options_name.lower()+options_chn_id.lower()+"_u"
        l = len(self.sequences[options_chn_id])
        ter_dict = self.process_terminals(options_terminals)
        if len(ter_dict) > 1:
            print("""ERROR: Only one unfolded chain of the protein can be
                 generated at the same time. Terminal information for only
                 one residue is allowed. More than one present.""")
            print("Number of terminal entries:"+str(len(ter_dict)))
            print("Exit without results.")
            sys.exit(1)
        self.gen_unfolded_charmm(name,l,options_chn_id,ter_dict)
        
    def gen_unfolded_charmm(self,name,l,chn_id,ter_dict):
        # TODO Incorporate ter_dict into this function
        outFile = open(name+".inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        outFile.write('\n')
        outFile.write("read sequence card\n")
        outFile.write("* gen seq for chain "+chn_id+"\n")
        outFile.write("*\n")
        outFile.write(str(l)+"\n")
        count = 0
        step = 10
        while count < l:
            for k in self.sequences[chn_id][count:(count+step)]:
                outFile.write(k+' ')
            count += step
            outFile.write("\n")
        outFile.write("\n")
        if len(ter_dict) == 0:
            self.print_template(outFile,self.gen_protein(chn_id,name))
        else:
            if chn_id in ter_dict:
                if ter_dict[chn_id]['first'] == 'none' and ter_dict[chn_id]['last'] == 'none':
                    self.print_template(outFile,self.gen_protein(chn_id,name))
                elif ter_dict[chn_id]['first'] != 'none' and ter_dict[chn_id]['last'] != 'none':
                    NTER = 'first '+ter_dict[chn_id]['first']
                    CTER = 'last '+ter_dict[chn_id]['last']
                    self.print_template(outFile,self.gen_protein_ter(chn_id,name,NTER,CTER))
                else:
                    print("ERROR: Unfolded Protein can only be built with no terminals, or")
                    print("       with both terminals. Use 'prepr pdb' to add terms.")
                    sys.exit(1)
            else:
                print("ERROR: No terminal information for chain Id "+chn_id)
                print("       Program terminated without results.")
                sys.exit(1)
        outFile.close()
        self.run.run_charmm(name+".inp",name+".out")
        print(name)
        fixPPDB = CR.ColumnsRows(name+".pdb")
        fixPPDB.fix_pdb_from_CHARMM()
        os.system("rm *.inp *.out")
    
    def prepare_pdb_for_charmm(self, pdbin, prot_ends, dir_path = ''):
        if not os.path.exists(pdbin):
            print("Error: Invalid path in --structure option.")
            sys.exit(1)
        if (os.path.basename(pdbin).split('.')[1] != 'pdb'):
            print("""ERROR: A PDB file with lower case 'pdb' suffix is required.
                  Program will exit without results.""")
            sys.exit(1)
        input_name = os.path.basename(pdbin).split('.')[0]
        if dir_path == '':
            output_name = input_name.lower()
        else:
            if dir_path[-1] == '/':
                pass
            else:
                dir_path += '/'
            output_name = dir_path+input_name.lower()
        reduceopt="-HIS -FLIP -OH -ROTEXOH -BUILD -OCC0.0 -H2OOCC0.0 -H2OB1000"
        os.system("reduce {} {} > {} &".format(reduceopt,pdbin,output_name+"r.pdb"))
        print("Ran Reduce\n")
        if not os.path.exists(output_name+"r.pdb"):
            print("Error: OS System call to run reduce didn't output a PDB file.")
            sys.exit(1)
        print("Passed\n")
        pdb_parser = struct.PDBParser(QUIET=True)
        self.sp = pdb_parser.get_structure(input_name, pdbin)
        count_models = 0
        #for i in self.structure.get_models():
        for i in self.sp.get_models():
            count_models += 1
        if count_models != 1:
            print("ERROR: This function works on a PDB file with only one model.")
            print("       Extract PDBs from CIF if needed. Simulations can only be done")
            print("       from a single structure and not multiple structure models.")
            sys.exit(1)
        self.count_ions = 0
        ''' Some amino acids can come with alternative names.'''
        self.num_residues = {}
        self.check_unusual_names()
        '''
         This fixes when there are multiple chains and the residue number gets
         reset at the beginning of each chain. With this fix, residue numbers
         will be renumbered
        '''
        coor = crd("")
        coor.write_crd(self.sp,output_name)
        # TODO we need way to check prot_ends match number of chains in 
        #      structure if ends are already present and detect automatically
        self.chains_info = {}
        # Gen SEQ files for each chain
        self.gen_seq(dir_path, prot_ends)
        # Gen FIXERS to keep offset amino acid numbers
        self.gen_fixers(dir_path)
        print(dir_path+"setup_"+input_name+".inp")
        outFile = open(dir_path+"setup_"+input_name+".inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        for i in self.sp.get_chains():
            self.print_template(outFile,self.fix_sequence(self.chains_info[i.id]))
        self.print_template(outFile,self.wmain('charge'))
        self.print_template(outFile,self.write_psf(input_name.lower()+'r'))
        self.print_template(outFile,self.write_ic(input_name.lower()+'r'))
        self.print_template(outFile,self.read_crd(input_name.lower()))
        self.print_template(outFile,self.build_heavy_atoms())
        self.print_template(outFile,self.build_hydrogens())
        self.print_template(outFile,self.wmain('charge'))
        self.print_template(outFile,self.write_crd(output_name+'r'))
        self.print_template(outFile,self.write_pdb(output_name+'r'))
        self.print_template(outFile,self.write_psf_xplor(output_name+'r'))
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm(dir_path+"setup_"+input_name+".inp",dir_path+"setup_"+input_name+".out")
        fixPPDB = CR.ColumnsRows(output_name+"r.pdb")
        fixPPDB.fix_pdb_from_CHARMM()
        
    def simple_minimization(self):
        pass
    
    def minimization_1(self,input_file,input_info):
        input_name = os.path.basename(input_file).split('.')[0]
        dir_name = os.path.dirname(input_file)
        output_name = input_name.lower()
        outFile = open("min1_"+input_name+".inp", 'w')
        self.print_template(outFile,self.gen_parameters().split('/n'))
        self.print_template(outFile,self.stream_nina_radii())
        if len(dir_name) == 0:
            self.print_template(outFile,self.read_psf(input_name))
            self.print_template(outFile,self.read_crd2(input_name))
        else:
            self.print_template(outFile,self.read_psf(dir_name+"/"+input_name))
            self.print_template(outFile,self.read_crd2(dir_name+"/"+input_name))
        input_list = input_info.split(':')
        input_list = [i.replace('-', ',') for i in input_list]
        input_list = [i.split(',') for i in input_list]
        self.print_template(outFile,self.minimization1(output_name,input_list))
        self.print_template(outFile,self.stop())
        outFile.close()
        self.run.run_charmm("min1_"+input_name+".inp","min1_"+input_name+".out")
        
    def check_unusual_names(self):
        HC = MRM.Histidin_Correction()
        for i in self.sp.get_models():
            for j in i.get_chains():
                count = 0
                for k in j.get_residues():
                    count += 1
                    if k.resname == 'HIS':
                        newname = HC.rename_HIS(k)
                        if not newname == 'HIS':
                            for l in k.get_atoms():
                                k.resname = newname
                    elif k.resname == ' NA':
                        k.resname = 'SOD'
                        for l in k.get_atoms():
                            l.id = 'SOD'
                            l.fullname = ' SOD'
                    elif k.resname == '  A':
                        k.resname = 'ADE'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == ' DA':
                        k.resname = 'ADE'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == '  G':
                        k.resname = 'GUA'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == ' DG':
                        k.resname = 'GUA'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == '  C':
                        k.resname = 'CYT'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == ' DC':
                        k.resname = 'CYT'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == '  T':
                        k.resname = 'THY'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == ' DT':
                        k.resname = 'THY'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                    elif k.resname == '  U':
                        k.resname = 'URA'
                        for l in k.get_atoms():
                            if l.id == 'OP1':
                                l.id = 'O1P'
                                l.fullname = ' O1P'
                            elif l.id == 'OP2':
                                l.id = 'O2P'
                                l.fullname = ' O2P'
                self.num_residues[j.id] = count
                
    def gen_seq(self, dir_path, prot_ends):
        temp = self.process_terminals(prot_ends)
        print(temp)
        # TODO: Can an unfolded protein be generated by a .SEQ file? 
        for i in self.sp.get_chains():
            print("Chain:"+i.get_id())
            file_name = i.id + '.SEQ'
            self.chains_info[i.id] = {}
            self.chains_info[i.id]['id'] = i.id
            self.chains_info[i.id]['seq'] = file_name
            self.chains_info[i.id]['first'] = temp[i.id]['first']
            self.chains_info[i.id]['last'] = temp[i.id]['last']
            outFile = open(dir_path+file_name, 'w')
            outFile.write('* Sequence of Chain ' + i.id + '\n')
            outFile.write('* Name: Molecule\n')
            outFile.write('* Title: Generated by prepr protein\n')
            outFile.write('*\n')
            x = 5 - len(str(self.num_residues[i.id]))
            outFile.write(' ' * x + str(self.num_residues[i.id]) + '\n')
            for j in i.get_residues():
                outFile.write(j.resname.strip(' ') + '\n')
            outFile.close()
            
    def process_terminals(self, prot_ends):
        temp = {}
        if prot_ends != None:
            chain_terminals = prot_ends.split(':')
            for i in chain_terminals:
                temp2 = i.split(',')
                temp[temp2[0]] = {}
                temp[temp2[0]]['first'] = temp2[1]
                temp[temp2[0]]['last'] = temp2[2]
        return temp
        
    def gen_fixers(self, dir_path):
        offset = 0
        residue_count = 1
        for i in self.sp.get_chains():
            file_name = i.id + '_FIXRES.INP'
            self.chains_info[i.id]['inp'] = file_name
            outFile = open(dir_path+file_name, 'w')
            outFile.write("! CHARMM Script to renumber RESID's after reading from SEQ file\n")
            outFile.write('! Generated by pdb_to_crd.py\n')
            outFile.write('\n')
            outFile.write('! First pass .... assigning temporary RESID\n')
            count = 1
            for j in self.sp.get_residues():
                if j.get_parent().id == i.id:
                    if count == 1:
                        offset = j.get_id()[1] - 1
                    outFile.write('! ' + str(offset + count) + ' , ' + str(count) + '\n')
                    outFile.write('rename resid x' + str(count) + ' sele segid ' + i.id + ' .and. resid ' + str(
                        count) + ' end\n')
                    count += 1
                    residue_count += 1
            outFile.write('! First pass complete.\n')
            outFile.write('\n')
            outFile.write('! Second pass .... re-assigning correct RESID\n')
            count = 1
            for j in self.sp.get_residues():
                if j.get_parent().id == i.id:
                    outFile.write(
                        'rename resid ' + str(offset + count) + ' sele segid ' + i.id + ' .and. resid x' + str(
                            count) + ' end\n')
                    count += 1
            outFile.write('! Second pass complete.\n')
            outFile.write('\n')
            outFile.write('! End of generated script\n')
            outFile.close()

    def file_not_found_ERROR(self, option, command):
        if not os.path.exists(option):
            print("Error: Invalid path in "+command+" option.")
            sys.exit(1)
