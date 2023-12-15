import os

class RunCharmm(object):
    def __init__(self):
        pass
    
    def run_charmm(self, command_file, out_file,charmm_dist="charmm"):
        string_format  = charmm_dist+" < {} > {}"
        command = string_format.format(command_file, out_file)
        os.system(command)
