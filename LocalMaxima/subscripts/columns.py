import argparse
import LocalMaxima.utilities.ColRows as CR

def register_parser(subparsers):
    parser = subparsers.add_parser('columns', usage=usage(), description=description())
    add_arguments(parser)

def add_arguments(parser):
    parser.add_argument("--structure", metavar="FILE", 
                        help="Input structure file in PDB or CIF format.", required=True)
    parser.add_argument('--frm', type=int, action="store", default=73, 
                        help="This option takes an integer value for initial column. Default is 73.")
    parser.add_argument('--to', type=int, action="store", default=22, 
                        help="This option takes an integer value for destination column. Default is 22.")
    parser.set_defaults(func=run)

def run(options):
    structure_file = CR.ColumnsRows(options.structure)
    structure_file.fix_pdb_from_CHARMM(options.to, options.frm)

def description():
        return '''It fixes a PDB file that was output by CHARMM. The chain identifier is placed by CHARMM in a column
        that is different from what Biopython and Entropy Maxima can handle. The chain identifier is swapped by default
        from column 73 to column 22 in the PDB file, but other columns swaps are allowed for possible formating
        variations. Columns are indexed from 1 to N.
        '''

def usage():
    return '\nmodifr columns --input 2GIA.pdb \n' \
           'modifr columns --input 2GIA.pdb --frm 72 --to 21\
            \n or type strctr fixpdb --help'

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description=description())
    add_arguments(arg_parser)
    args = arg_parser.parse_args()
    args.func(args)
