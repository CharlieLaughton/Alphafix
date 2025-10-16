# Alphafix

Alphafix uses data from the Alphafold database of protein models to remediate incomplete protein structures - e.g. crystal structures missing atoms in flexible loops or sidechains.

```bash
# Remediate chain A from PDB entry 5FDR:
% alpha_fix -i 5fdr.pdb -c A -o 5fdr_A_fixed.pdb
*** ALPHAFIX V. 0.1 ***
Comparison of protein chain A with uniprot entry Q07820:
  Missing residues: M1-E170 (Ignored)
  Mutation:         E171->D
  Missing atoms:    D195: CG, OD1, OD2
  Missing atoms:    K197: CG, CD, CE, NZ
  Missing atoms:    R201: CG, CD, NE, NH1, NH2, CZ
  Missing atoms:    E325: CG, CD, OE1, OE2
  Missing residues: G326-R350 (Ignored)


Optimisation 1:
Restrained minimization performed using OpenMM
  force constant: 1.0 kcal/mol/Å²
  maximum cycles: 200
  initial energy: 570.52 kcal/mol
  final energy: -265.45 kcal/mol
  number of restrained atoms: 1233

Optimisation 2:
Restrained minimization performed using OpenMM
  force constant: 1.0 kcal/mol/Å²
  maximum cycles: 200
  initial energy: 129.55 kcal/mol
  final energy: -273.78 kcal/mol
  number of restrained atoms: 1255

%
```

## How it works

1. UniProt entries matching each chain in the protein are identified.
2. The corresponding Alphafold database structures are downloaded.
3. The alphafold structures are least-squares fitted to the input structure on the basis of a sequence alignment, and conserved and missing atoms identified.
4. Restrained energy minimisation of the alphafold structures drives conserved atoms towards their input file counterparts.
5. Conserved atoms in the optimised alphafold structure have their coordinates replaced with those of their input file counderparts.
6. A second round of restrained optimisation is performed.

## Requirements

Alphafix requires that Ambertools and blastp (with the Swissprot database files) are installed. The easiest way to do this is via conda. Note it is neccessary to roll the version of numpy back to 2.2 to avoid problems with ambertools::

```bash
% conda install conda-forge::ambertools bioconda::blast -y
% conda install numpy==2.2
% export BLASTDB="/path/to/blast/database/folder"
% update_blastdb.pl --decompress swissprot
```

## Installation

Alphafix can be installed via pip:

```bash
pip install git+https://github.com/CharlieLaughton/Alphafix.git
```
Insalling Alphafix makes two command-line tools available:

 - `alpha_fix`: remediates protein structures in PDB files
  - `alpha_check`: reports on remediations required to fix each protein chain in a PDB file.

## Usage

```bash
% alpha_fix -h
usage: alpha_fix [-h] -i INPDB -o OUTPDB [-c [CHAINS ...]]
                 [-u [UNIPROT_IDS ...]] [-l LOG] [-n] [--version]

Fix missing residues in a PDB file using Alphafold.

options:
  -h, --help            show this help message and exit
  -i, --inpdb INPDB     Input PDB file.
  -o, --outpdb OUTPDB   Fixed PDB file.
  -c, --chains [CHAINS ...]
                        List of chain IDs to include in the output.
  -u, --uniprot_ids [UNIPROT_IDS ...]
                        List of UniProt IDs for the input structure.
  -l, --log LOG         Log file for Alphafold output.
  -n, --no_trim         Don't trim the fixed PDB file to match the input.
  --version             show program's version number and exit



% alpha_check -h
usage: alpha_check [-h] -i INPDB [-u [UNIPROT_IDS ...]] [-l LOG] [--version]

Check how well a PDB file matches its UniProt sequences.

options:
  -h, --help            show this help message and exit
  -i, --inpdb INPDB     Input PDB file.
  -u, --uniprot_ids [UNIPROT_IDS ...]
                        List of UniProt IDs for the input structure.
  -l, --log LOG         Log file for output.
  --version             show program's version number and exit

%
```

## Author

Charlie Laughton charles.laughton@nottingham.ac.uk
