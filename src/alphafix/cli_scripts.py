from argparse import ArgumentParser
from alphafix._version import __version__
from alphafix.tools import alpha_check, alpha_fix


def alpha_check_cli():
    parser = ArgumentParser(
        description="Check how well a PDB file matches its UniProt sequences.")
    parser.add_argument("-i", "--inpdb",
                        help="Input PDB file.", required=True)
    parser.add_argument("-u", "--uniprot_ids", nargs='*', required=False,
                        help="List of UniProt IDs for the input structure.")
    parser.add_argument("-l", "--log", help="Log file for output.")
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()
    if not args.inpdb:
        parser.print_help()
        return
    try:
        log = alpha_check(args.inpdb, args.uniprot_ids)
    except Exception as e:
        print(f"{e}")
        return
    if args.log:
        with open(args.log, 'w') as log_file:
            log_file.write(log)
    else:
        print(log)


def alpha_fix_cli():
    parser = ArgumentParser(
        description="Fix missing residues in a PDB file using Alphafold.")
    parser.add_argument("-i", "--inpdb",
                        help="Input PDB file.", required=True)
    parser.add_argument("-o", "--outpdb",
                        help="Fixed PDB file.", required=True)
    parser.add_argument("-c", "--chains", nargs='*', required=False,
                        help="List of chain IDs to include in the output.")
    parser.add_argument("-u", "--uniprot_ids", nargs='*', required=False,
                        help="List of UniProt IDs for the input structure.")
    parser.add_argument("-l", "--log", help="Log file.")
    parser.add_argument("-n", "--no_trim", action="store_true",
                        help="Don't trim the fixed PDB file"
                        " to match the input.")
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()
    if not args.inpdb or not args.outpdb:
        parser.print_help()
        return
    try:
        out_pdb, log = alpha_fix(args.inpdb, args.uniprot_ids,
                                 chains=args.chains, trim=not args.no_trim)
    except Exception as e:
        print(f"{e}")
        return
    out_pdb.save(args.outpdb)
    if args.log:
        with open(args.log, 'w') as log_file:
            log_file.write(log)
    else:
        print(log)
