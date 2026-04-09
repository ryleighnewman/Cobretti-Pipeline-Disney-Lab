#!/usr/bin/python3
"""
  __    __     ______     ______     ______     __         ______     ______
 /\ "-./  \   /\  __ \   /\  ___\   /\  ___\   /\ \       /\  __ \   /\  == \
 \ \ \-./\ \  \ \ \/\ \  \ \___  \  \ \___  \  \ \ \____  \ \  __ \  \ \  __<
  \ \_\ \ \_\  \ \_____\  \/\_____\  \/\_____\  \ \_____\  \ \_\ \_\  \ \_____\
   \/_/  \/_/   \/_____/   \/_____/   \/_____/   \/_____/   \/_/\/_/   \/_____/

                                   COBRETTI
               "Does the job nobody wants" -Cobra (1986) trailer
Automates redundant RNA bioinformatic tasks, saving valuable time and resources

                  Contact: Snake Peterson - jakemp@iastate.edu
                  or the Moss Lab: https://github.com/moss-lab


Required modules for Pronto/HPC users:
    py-biopython
    python

Use:
    python cobretti.py -stage [#] -email [email@iastate.edu]
    Example stages: 1A, 2C, 3A (see manual for all options)

Manual:
    https://github.com/moss-lab/Cobretti/blob/main/cobretti_manual.docx
"""

import argparse
from Bio import Entrez
from Bio import SeqIO
import logging
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

current_directory = Path.cwd()

# Program default locations
COBRETTI = Path('/blue/mdisney/ryleighnewman/cobretti_project/Cobretti/cobretti.py')

SCANFOLD_1 = Path('/apps/cobretti/1.0/scripts/ScanFold/ScanFold.py')
CMBUILDER = Path('/apps/cobretti/1.0/scripts/labtools/cm-builder')
RSCAPE = Path('/apps/cobretti/1.0/programs/rscape/bin/R-scape')
PERL = Path('/apps/cobretti/1.0/programs/viennarna/lib/perl5/')
RNAFRAMEWORK = Path('/apps/cobretti/1.0/programs/RNAFramework/lib/')
KNOTTY = Path('/apps/cobretti/1.0/programs/knotty/Knotty/knotty')
HFOLD = Path('/apps/cobretti/1.0/programs/iterative-hfold/Iterative-HFold/Iterative-HFold')
INFERNAL = Path('/apps/cobretti/1.0/programs/infernal')



SCANFOLD_2 = Path('/apps/cobretti/1.0/programs/ScanFold2.0/ScanFoldBothForInforna.py')
SCANFOLD_2_ENV = Path('/apps/cobretti/1.0/envs/SF2')

R2R = Path('/apps/cobretti/1.0/programs/r2r/bin/r2r')
#PERL = Path('/work/LAS/wmoss-lab/programs/lib/perl5/')
#RNAFRAMEWORK = Path('/work/LAS/wmoss-lab/programs/RNAFramework/lib/')
#KNOTTY = Path('/work/LAS/wmoss-lab/programs/knotty/knotty')
#HFOLD = Path('/work/LAS/wmoss-lab/programs/hfold/HFold_iterative')
SIMRNA = Path('/work/LAS/wmoss-lab/programs/SimRNA')
QRNAS = Path('/work/LAS/wmoss-lab/programs/qrnas/QRNA')
QRNAS_FF = Path('/work/LAS/wmoss-lab/programs/qrnas/forcefield')
ARES = Path('/work/LAS/wmoss-lab/programs/ares')
ARES_ENV = Path('/work/LAS/wmoss-lab/programs/envs/ares')
FPOCKET = Path('/work/LAS/wmoss-lab/programs/fpocket2/bin/fpocket')


def main():
    parser = argparse.ArgumentParser()

    # Required inputs
    parser.add_argument('-stage', type=str, default='0',
                        help='input stage number to run')
    parser.add_argument('-email', type=str, default='',
                        help='input email address')

    # Optional inputs
    parser.add_argument('-i', type=str, default='',
                        help='input list of fasta accession numbers (stage 1A)')
    parser.add_argument('-seq', type=Path, default=current_directory,
                        help='input alternate location for fasta sequences')
    parser.add_argument('-db', type=Path, default=current_directory,
                        help='input alternate location for BLAST sequence files')
    parser.add_argument('-dbn', type=Path, default=current_directory,
                        help='input alternate location for fasta sequences')
    parser.add_argument('-pk', type=Path, default=current_directory,
                        help='input alternate location for pseudoknot motifs')
    parser.add_argument('-set_dbn_extend', type=bool, default=True,
                        help='extend pseudoknot motifs True/False (stage 1B)')

    # Program location inputs
    parser.add_argument('-cobretti', type=Path, default=COBRETTI,
                        help='input location of cobretti.py')
    parser.add_argument('-sf', type=Path, default=SCANFOLD_1,
                        help='input location of ScanFold.py')
    parser.add_argument('-sf2', type=Path, default=SCANFOLD_2,
                        help='input location of ScanFold2.0.py')
    parser.add_argument('--sf2env', type=Path, default=SCANFOLD_2_ENV,
                        help='input location of ScanFold2.0 conda environment')
    parser.add_argument('-cmb', type=Path, default=CMBUILDER,
                        help='input location of cm-builder')
    parser.add_argument('-rs', type=Path, default=RSCAPE,
                        help='input location of R-Scape')
    parser.add_argument('-r2r', type=Path, default=R2R,
                        help='input location of R2R')
    parser.add_argument('-perl', type=Path, default=PERL,
                        help='input location of Perl directory (perl5) (for R-Scape)')
    parser.add_argument('-rf', type=Path, default=RNAFRAMEWORK,
                        help='input location of RNAFramework directory (RNAFramework/lib/) (for R-Scape)')
    parser.add_argument('-ky', type=Path, default=KNOTTY,
                        help='input location of Knotty (knotty)')
    parser.add_argument('-hf', type=Path, default=HFOLD,
                        help='input location of Iterative HFold (HFold_iterative)')
    parser.add_argument('-sim', type=Path, default=SIMRNA,
                        help='input location of SimRNA directory')
    parser.add_argument('-qrnas', type=Path, default=QRNAS,
                        help='input location of QRNAS')
    parser.add_argument('-qrnasff', type=Path, default=QRNAS_FF,
                        help='input location of QRNAS force field')
    parser.add_argument('-ares', type=Path, default=ARES,
                        help='input location of ARES')
    parser.add_argument('-aresenv', type=Path, default=ARES_ENV,
                        help='input location of ARES conda environment')
    parser.add_argument('-fpocket', type=Path, default=FPOCKET,
                        help='input location of fpocket')

    # Add args
    args = parser.parse_args()
    fasta_list = args.i
    stage = args.stage
    email = args.email
    sequence_directory = args.seq
    database_directory = args.db
    dbn_directory = args.dbn
    pk_directory = args.pk
    set_dbn_extend = args.set_dbn_extend
    cobretti_program = args.cobretti
    scanfold_program = args.sf
    scanfold2_program = args.sf2
    scanfold2_environment = args.sf2env
    cmbuilder_program = args.cmb
    rscape_program = args.rs
    r2r_program = args.r2r
    perl_program = args.perl
    rnaframework_directory = args.rf
    knotty_program = args.ky
    hfold_program = args.hf
    simrna_directory = args.sim
    qrnas_program = args.qrnas
    qrnas_ff = args.qrnasff
    ares_program = args.ares
    ares_environment = args.aresenv
    fpocket_program = args.fpocket

    # Exit if required inputs missing
    if email == '':
        logging.error('No email provided, use "-email" to define string. Exiting...')
        sys.exit()
    if stage == '0':
        logging.error('No stage specified, use "-stage" to define string (e.g., "1A"). Exiting...')
        sys.exit()
    # If stage/email provided, execute stage functions
    logging.info(f'Cobretti stage {stage} run in progress...')
    # Checksum to see if any stage(s) ran
    is_stage = False

    if stage in ('1A', '1AA', '1AB', '1A1', '1A1A'):
        # Prepare files for ScanFold run
        is_stage = True
        logging.info('Preparing sequence files...')
        if sequence_directory == current_directory:
            # Make a sequence directory if there is none
            Path.mkdir(Path.joinpath(current_directory, 'sequences'), exist_ok=True)
            sequence_directory = Path.joinpath(current_directory, 'sequences')
        if fasta_list != '':
            # Read sequences in textfile, pull sequences from Entrez database, place in sequence directory
            fasta_build(fasta_list, sequence_directory, email)
        else:
            # Move sequences from current directory to sequence directory
            for filepath in current_directory.glob('*.fa*'):
                try:
                    with open(filepath, 'r') as readfile:
                        # Standardize naming convention, uses internal name from ">" line NOT the filename
                        fasta_lines = readfile.readlines()
                        fasta_name = fasta_lines[0].rstrip().replace('.', '-').replace('_', '-').strip('>').split(' ')[
                                         0] + '.fasta'
                    shutil.move(filepath, Path.joinpath(sequence_directory, fasta_name))
                except:
                    logging.error(f'Unable to move file: {filepath.name}')

    if stage in ('1A', '1AA'):
        # Run ScanFold 2.0
        is_stage = True
        logging.info('Running ScanFold 2.0...')
        scanfold2_prep(sequence_directory, scanfold2_program, scanfold2_environment, email)
        scanfold2_run(current_directory)

    if stage in ('1A1','1A1A'):
        # Legacy code to run ScanFold 1.0
        is_stage = True
        logging.info('Running ScanFold 1.0...')
        scanfold_prep(sequence_directory, scanfold_program, email)
        scanfold_run(current_directory)

    if stage in ('1A', '1AB', '1A1'):
        # Run NCBI BLAST using Pronto databases, save results to database directory
        is_stage = True
        logging.info('Running BLAST...')
        if database_directory == current_directory:
            # Make a database directory if there is none
            Path.mkdir(Path.joinpath(current_directory, 'databases'), exist_ok=True)
            database_directory = Path.joinpath(current_directory, 'databases')
        blast_prep(sequence_directory, database_directory, email)
        blast_run(current_directory)

    if stage in ('1B', '1BA', '1BB', '1BC', '1BC1'):
        # Build/set locations of files
        is_stage = True
        logging.info('Preparing files for PK/cm-builder...')
        if sequence_directory == current_directory:
            sequence_directory = folder_check(sequence_directory, 'sequences')
        if database_directory == current_directory:
            database_directory = folder_check(database_directory, 'databases')
        if dbn_directory == current_directory:
            # Make a dot-bracket notation directory if there is none
            Path.mkdir(Path.joinpath(current_directory, 'motifs'), exist_ok=True)
            dbn_directory = Path.joinpath(current_directory, 'motifs')
        if pk_directory == current_directory:
            # Make a pseudoknot directory if there is none
            Path.mkdir(Path.joinpath(current_directory, 'pk_motifs'), exist_ok=True)
            pk_directory = Path.joinpath(current_directory, 'pk_motifs')

    if stage in ('1B', '1BA'):
        # Move BLAST files and edit databases to make them more readable
        is_stage = True
        logging.info('Optimizing BLAST databases...')
        blast_cleanup(current_directory, database_directory, dbn_directory)

    if stage in ('1B', '1BB'):
        # Refold motifs using Knotty and HFold
        is_stage = True
        if set_dbn_extend:
            # Extend all motifs for pseudoknot refolding
            logging.info('Extending motifs...')
            extend_motifs(sequence_directory, dbn_directory)
        logging.info('Checking for pseudoknots...')
        pk_fold(knotty_program, hfold_program)
        logging.info('Cleaning up pseudoknot results...')
        pk_cleanup()
        logging.info('Nesting pseudoknots...')
        pk_breakdown(pk_directory)

    if stage in ('1B', '1BC', '1BC1'):
        # Prepare cm-builder files, stage 1BC1 will NOT run them
        is_stage = True
        logging.info('Preparing cm-builder scripts...')
        cmbuilder_prep(sequence_directory, database_directory, pk_directory, cmbuilder_program, email, rscape_program,
                       r2r_program,
                       cobretti_program, perl_program, rnaframework_directory)

    if stage in ('1B', '1BC'):
        # Run cm-builder shell scripts
        is_stage = True
        logging.info('Running cm-builder scripts...')
        cmbuilder_run(current_directory)

    if stage in ('1C'):
        # Clean up cm-builder runs and run R-Scape
        is_stage = True
        logging.info('Organizing files and running R-Scape...')
        cmbuilder_cleanup(current_directory)

    if stage in ('1CA'):
        # Substage called by R-Scape script to organize remaining files
        is_stage = True
        logging.info('R-Scape complete, organizing results...')
        cmbuilder_cleanup2(current_directory)

    if stage in ('2A', '2AA'):
        # Prepare SimRNA scripts
        is_stage = True
        logging.info('Preparing SimRNA scripts...')
        simrna_prep(current_directory, simrna_directory, email)

    if stage in ('2A', '2AB'):
        # Run SimRNA scripts
        is_stage = True
        logging.info('Running SimRNA scripts...')
        simrna_run(current_directory)

    if stage in ('2B', '2BA'):
        # Organize SimRNA files
        is_stage = True
        logging.info('Organizing SimRNA files...')
        simrna_cleanup(current_directory)

    if stage in ('2B', '2BB'):
        # Prepare QRNAS scripts
        is_stage = True
        logging.info('Preparing QRNAS scripts...')
        qrnas_prep(current_directory, qrnas_program, qrnas_ff, email)

    if stage in ('2B', '2BC'):
        # Run QRNAS scripts
        is_stage = True
        logging.info('Running QRNAS scripts...')
        qrnas_run(current_directory)

    if stage in ('2C', '2CA'):
        # Organize QRNAS files
        is_stage = True
        logging.info('Organizing QRNAS files...')
        qrnas_cleanup(current_directory)

    if stage in ('2C', '2CB'):
        # Prepare ARES scripts
        is_stage = True
        logging.info('Preparing ARES scripts...')
        ares_prep(current_directory, ares_program, ares_environment, email)

    if stage in ('2C', '2CC'):
        # Run ARES scripts
        is_stage = True
        logging.info('Running ARES scripts...')
        ares_run()

    if stage in ('2D'):
        # Final stage cleanup and fpocket run
        # Separated stage to avoid ARES reading duplicate PDB files in fpocket subdirectories
        is_stage = True
        logging.info('Running fpocket...')
        fpocket_run(current_directory, fpocket_program)
        logging.info('Organizing fpocket results...')
        fpocket_cleanup(current_directory)
        fpocket_read(current_directory)

    if stage in ('3A', '3AA'):
        # Prepare DOCK 6 scripts
        is_stage = True
        logging.info('Preparing DOCK 6 scripts...')
        dock6_prep()

    if stage in ('3A', '3AB'):
        # Run DOCK 6 scripts
        is_stage = True
        logging.info('Running DOCK 6 scripts...')
        dock6_run(current_directory)

    if stage in ('3B'):
        # Run AnnapuRNA
        is_stage = True
        logging.info('Running AnnapuRNA...')
        annapurna_run(current_directory)
    if not is_stage:
        logging.error('Incorrect stage specified, use "-stage" to define string (e.g., "1A"). Exiting...')
        sys.exit()


def folder_check(working_directory, subdirectory):
    # Check if subdirectory exists; if so, return as path; if not, return directory as path
    if Path.exists(Path.joinpath(working_directory, subdirectory)):
        return Path.joinpath(working_directory, subdirectory)
    else:
        return working_directory


def shell_build_start(filename, job, email, time=3, nodes=1, mem=50, tasks=1, notify='ALL'):
    # Build shell scripts with Pronto/HPC settings and modules
    with open(filename, 'w', newline='\n') as writefile:
        writefile.writelines('#!/bin/bash -l\n')
        writefile.writelines('#SBATCH --partition=hpg-default\n')  # Partition to submit to
        writefile.writelines(f'#SBATCH --time={time}-00:00:00\n')  # Time limit for this job in days
        writefile.writelines(f'#SBATCH --nodes={nodes}\n')  # Nodes to be used for this job during runtime
        writefile.writelines(f'#SBATCH --account=scripps-dept\n') #scripps dept resources
        writefile.writelines(f'#SBATCH --qos=scripps-dept-b\n') #scripps dept resources
        if mem != 0:
            writefile.writelines(f'#SBATCH --mem={mem}G\n')  # Optional memory allocation
        writefile.writelines(f'#SBATCH --ntasks-per-node={tasks}\n')
        writefile.writelines(f'#SBATCH --job-name={job}\n')  # Name of this job in work queue
        writefile.writelines(f'#SBATCH --mail-user={email}\n')  # Email to send notifications to
        writefile.writelines(f'#SBATCH --mail-type={notify}\n')  # Email notification type (BEGIN, END, FAIL, ALL)
        writefile.writelines('\n')
        # Add Pronto modules
        if job.startswith('scanfold2'):
            writefile.writelines('#SBATCH --export=NONE\n\n')
            writefile.writelines('module purge\n')
            writefile.writelines('module load miniconda3\n')
        elif job.startswith('scanfold1'):
            #writefile.writelines('module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core\n')
            writefile.writelines('module load cobretti\n')
            #writefile.writelines('module load python/3.6.5-fwk5uaj\n')
        elif job.startswith('blast'):
            writefile.writelines('module load ncbi_blast\n')
            writefile.writelines('echo ${NCBI_BLAST_DB_PATH}: /data/reference/blast/db\n')
            writefile.writelines('module load cobretti\n')
            #writefile.writelines('export BLASTDB=${NCBI_BLAST_DB_PATH}\n')
        elif job.startswith('cmbuilder'):
            writefile.writelines('module load cobretti\n')
            #writefile.writelines('module load infernal\n')
        elif job.startswith('rscape'):
            writefile.writelines('module load gcc ghostscript gnuplot\n')
            #writefile.writelines('module load ghostscript\n')
            #writefile.writelines('module load gnuplot\n')
        elif job.startswith('qrnas'):
            writefile.writelines('module load gcc\n')
        elif job.startswith('ares'):
            writefile.writelines('#SBATCH --export=NONE\n\n')
            writefile.writelines('module purge\n')
            writefile.writelines('module load miniconda3\n')
        else:
            writefile.writelines('module load py-biopython\n')
            writefile.writelines('module load python\n')
        writefile.writelines('\n')


def fasta_build(filename, sequence_directory, email):
    # read a text file of accession numbers, download associated fasta files
    Path.mkdir(sequence_directory, exist_ok=True)
    with open(filename, 'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
            logging.debug(f'Name: {line}')
            with open(Path.joinpath(sequence_directory, line + '.fa'), 'w', newline='\n') as writefile:
                writefile.writelines(f'>{line}\n')
                entrez_sequence = fasta_find(line, email).replace('T', 'U')
                writefile.writelines(f'{entrez_sequence}\n')
                logging.debug(f'Seq: {entrez_sequence}')


def fasta_find(accession_number, email):  # Download fasta sequences via Entrez
    Entrez.email = email  # Required for Entrez
    logging.debug(f'Running Entrez for {accession_number}')
    entrez_search = Entrez.efetch(db='nucleotide', id=accession_number, rettype='fasta')
    entrez_record = SeqIO.read(entrez_search, 'fasta')
    return entrez_record.seq


def scanfold2_prep(sequence_directory, scanfold2_directory, scanfold2_environment, email):
    # Make shell scripts for ScanFold 2.0 runs
    for filepath in sequence_directory.glob('*.fa*'):
        with open(filepath, 'r') as readfile:
            lines = readfile.readlines()
            # Redundant header replace in case a stage was skipped
            name = lines[0].rstrip().replace('.', '-').replace('_', '-').strip('>').split(' ')[0]
            logging.debug(f'ScanFold 2.0 generic filename: {name}')
        shell_build_start(f'scanfold2_{name}.sh', f'scanfold2_{name}', email, mem=10, notify='END,FAIL')
        with open(f'scanfold2_{name}.sh', 'a', newline='\n') as writefile:
            writefile.writelines(f'conda activate {scanfold2_environment}\n')
            writefile.writelines('wait;\n')
            writefile.writelines(f'python {scanfold2_directory} {filepath} --folder_name {name} --global_refold -f -2 -w 120 &\n')
            writefile.writelines('wait;\n')


def scanfold2_run(working_directory):
    # Run ScanFold 2.0 shell scripts
    for filepath in working_directory.glob('scanfold2_*.sh'):
        logging.debug(f'ScanFold 2 shell path: {filepath}')
        os.system(f'sbatch {filepath.name}')


def scanfold_prep(sequence_directory, scanfold_directory, email):
    # Makes shell script(s) for ScanFold run(s)
    for filepath in sequence_directory.glob('*.fa*'):
        with open(filepath, 'r') as readfile:
            lines = readfile.readlines()
            # Redundant header replace in case a stage was skipped
            name = lines[0].rstrip().replace('.', '-').replace('_', '-').strip('>').split(' ')[0]
            logging.debug(f'ScanFold generic filename: {name}')
        shell_build_start(f'scanfold1_{name}.sh', f'scanfold1_{name}', email, mem=10, notify='END,FAIL')
        with open(f'scanfold1_{name}.sh', 'a', newline='\n') as writefile:
            writefile.writelines(f'python {scanfold_directory} {filepath} --name {name} --global_refold -f -2 -w 120 &\n')
            writefile.writelines('wait;\n')


def scanfold_run(working_directory):
    # Run ScanFold shell scripts
    for filepath in working_directory.glob('scanfold1_*.sh'):
        logging.debug(f'ScanFold shell path: {filepath}')
        os.system(f'sbatch {filepath.name}')


def blast_prep(sequence_directory, database_directory, email):
    # BLAST search using local databases, output to database directory
    for filepath in sequence_directory.glob('*.fa*'):
        for record in SeqIO.parse(str(filepath), 'fasta'):
            name = record.name.replace('.', '-').replace('_', '-')
            logging.debug(f'BLAST sequence name: {name}')
            blast_tempfile = f'{name}_db_raw.txt'
            tempfile_location = Path.joinpath(database_directory, blast_tempfile)
            shell_build_start(f'blast_{name}.sh', f'blast_{name}', email, tasks=8, notify='END,FAIL')
            with open(f'blast_{name}.sh', 'a', newline='\n') as writefile:
                writefile.writelines(f'blastn -task blastn -db nt -query {filepath} -out {tempfile_location} -outfmt '
                                     f'"6 sacc sseq" -max_target_seqs 2500 -num_threads 8 -max_hsps 1;\n')
                writefile.writelines('wait;')


def blast_run(working_directory):
    # Run BLAST shell scripts   
    for filepath in working_directory.glob('blast_*.sh'):
        logging.debug(f'BLAST shell path: {filepath}')
        os.system(f'sbatch {filepath.name}')


def blast_cleanup(working_directory, database_directory, dbn_directory):
    # Organize BLAST output files and move sequences to a single line
    folders = ['sh', 'out']
    for folder in folders:
        Path.mkdir(Path.joinpath(working_directory, folder), exist_ok=True)
    sh_directory = Path.joinpath(working_directory, 'sh')
    out_directory = Path.joinpath(working_directory, 'out')
    # Move ScanFold shell scripts
    for filepath in working_directory.glob('scanfold*_*.sh'):
        try:
            shutil.move(filepath, Path.joinpath(sh_directory, filepath.name))
        except:
            logging.error(f'Unable to move file: {filepath.name}')
    # Move BLAST shell scripts
    for filepath in working_directory.glob('blast*.sh'):
        try:
            shutil.move(filepath, Path.joinpath(sh_directory, filepath.name))
        except:
            logging.error(f'Unable to move file: {filepath.name}')
    # Move SLURM output files
    for filepath in working_directory.glob('slurm*.out'):
        try:
            shutil.move(filepath, Path.joinpath(out_directory, filepath.name))
        except:
            logging.error(f'Unable to move file: {filepath.name}')
    # Edit database files to proper FASTA format
    for filepath in database_directory.glob('*_db_raw.txt'):
        stem = filepath.name.split('_db_raw.txt')[0]
        database_filename = stem + '_db.fa'
        with open(Path.joinpath(database_directory, database_filename), 'w', newline='\n') as writefile:
            with open(filepath, 'r') as readfile:
                lines = readfile.readlines()
                for line in lines:
                    accession = line.rstrip().split('\t')[0]
                    sequence = line.rstrip().split('\t')[1]
                    sequence = sequence.replace('T', 'U').replace('-', '')
                    writefile.writelines(f'>{accession}\n')
                    writefile.writelines(f'{sequence}\n')
    # Make copies of all motifs and consolidate
    for filepath in working_directory.rglob('*_motif_*.dbn'):
        try:
            shutil.copy(filepath, Path.joinpath(dbn_directory, filepath.name))
        except:
            logging.error(f'Unable to copy file: {filepath.name}')


def extend_motifs(sequence_directory, dbn_directory, dbn_writefile='extended.dbn'):
    # Find all DBN files iin a directory and extend them using sequence files
    open(dbn_writefile, 'w', newline='\n')
    for filepath in dbn_directory.glob('*_motif_*.dbn'):
        try:
            logging.debug(f'Extending {filepath.name}')
            dbn_extend(sequence_directory, filepath, dbn_writefile)
        except:
            logging.error(f'Failed to extend file: {filepath.name}')


def dbn_extend(sequence_directory, dbn_filepath, dbn_writefile):
    # Extends ScanFold DBN motifs on 5' and 3' end
    extension = 30
    with open(dbn_filepath, 'r') as readfile:
        dbn_lines = readfile.readlines()
        tag = dbn_filepath.name.split('_')[0]
        dbn_header = str(dbn_lines[0]).rstrip()
        dbn_sequence = str(dbn_lines[1]).rstrip()
        dbn_structure = str(dbn_lines[2]).rstrip()
    for filepath in sequence_directory.glob(tag + '*.fa*'):
        logging.debug(f'Matching sequence: {filepath.name}')
        for record in SeqIO.parse(str(filepath), 'fasta'):
            fasta_sequence = str(record.seq).replace('T', 'U')
            sequence_start = fasta_sequence.find(dbn_sequence)
            sequence_length = len(dbn_sequence)
            sequence_end = sequence_start + sequence_length
            logging.debug(f'Sequence start: {sequence_start}')
            logging.debug(f'Sequence length: {sequence_length}')
            logging.debug(f'Sequence end: {sequence_end}')
            # If motif is near the end of the sequence, don't extend past the end
            if sequence_start < extension:
                dbn_start = 0
                five_prime_filler = sequence_start
            else:
                dbn_start = sequence_start - extension
                five_prime_filler = extension
            if sequence_end + extension > len(fasta_sequence):
                dbn_end = len(fasta_sequence)
                three_prime_filler = dbn_end - sequence_end
            else:
                dbn_end = sequence_end + extension
                three_prime_filler = extension
            # Append results, adding dots to unstructured extensions to maintain proper alignment
            with open(dbn_writefile, 'a+', newline='\n') as writefile:
                writefile.writelines(f'>{dbn_header}\n')
                writefile.writelines(f'{fasta_sequence[dbn_start:dbn_end]}\n')
                writefile.writelines(f'{"." * five_prime_filler}{dbn_structure}{"." * three_prime_filler}\n')


def pk_fold(knotty_program, hfold_program, dbn_readfile='extended.dbn', pk_writefile='tmppk.txt'):
    # Fold DBN motifs with Knotty, HFold and output results to a textfile
    with open(pk_writefile, 'w', newline='\n') as writefile, open(dbn_readfile, 'r') as readfile:
        i = 0
        nucleotides = {'A', 'C', 'G', 'U'}
        left_brackets = {'(', '[', '{', '<'}
        right_brackets = {')', ']', '}', '>'}
        ambiguous_nucleotide_codes = {
            "N": ('A', 'C', 'G', 'U'),
            "B": ('C', 'G', 'U'),
            "D": ('A', 'G', 'U'),
            "H": ('A', 'C', 'U'),
            "V": ('A', 'C', 'G'),
            "K": ('G', 'U'),
            "M": ('A', 'C'),
            "R": ('A', 'G'),
            "S": ('C', 'G'),
            "W": ('A', 'U'),
            "Y": ('C', 'U')}
        lines = [line for line in readfile.readlines() if line.strip()]
        for line_index, line in enumerate(lines):
            if (line_index + 3) > len(lines):
                # Stop when there aren't enough lines left to contain a DBN
                break
            logging.debug(f'First character in line {line_index}: {line[0]}')
            if line[0] == '>':
                dbn_header = str(lines[i]).rstrip()
                dbn_sequence = str(lines[i + 1]).rstrip()
                dbn_structure = str(lines[i + 2]).rstrip()
                left_pos = 0
                right_pos = len(dbn_sequence)
                if dbn_sequence[0] == '>' or dbn_structure[0] == '>':
                    # Skip if the sequence or structure contain a header
                    logging.error(f'Misaligned DBN on line {line_index}: {line}')
                    continue
                if not dbn_sequence[0].isalnum():
                    # Skip if the sequence contains special characters
                    logging.error(f'Sequence contains non-canonical nucleotides on line {line_index}: {line}')
                    continue
                if not (dbn_structure[0] == '.' or dbn_structure[0] in left_brackets or dbn_structure[
                        0] in right_brackets):
                    # Skip if the structure is non-standard
                    logging.error(f'Non-canonical pairing found on line {line_index}: {line}')
                    continue
                for k, j in enumerate(dbn_sequence):
                    # Shorten the 5' end until an unambiguous nucleotide or base pair is encountered
                    if j in nucleotides or dbn_structure[k] in left_brackets:
                        left_pos = k
                        break
                for k, j in enumerate(dbn_sequence[::-1]):
                    # Shorten the 3' end until an unambiguous nucleotide or base pair is encountered
                    if j in nucleotides or (dbn_structure[::-1])[k] in right_brackets:
                        right_pos -= k
                        break
                # Adjust sequence/structure length
                dbn_sequence = dbn_sequence[left_pos:right_pos]
                dbn_structure = dbn_structure[left_pos:right_pos]
                # Replace ambiguous nucleotides with a random canonical nucleotide
                # Knotty and HFold will error out unless ambiguous nucleotides are replaced
                new_sequence = ""
                for k in dbn_sequence:
                    if k in ambiguous_nucleotide_codes:
                        new_sequence += random.choice(ambiguous_nucleotide_codes[k])
                    else:
                        new_sequence += k
                dbn_sequence = new_sequence
                # Write ScanFold motif
                logging.debug(f'ScanFold starting sequence: {dbn_sequence}')
                writefile.writelines(f'{dbn_header} ScanFold\n')
                writefile.writelines(f'{dbn_sequence}\n')
                writefile.writelines(f'{dbn_structure}\n')
                # Write Knotty motif
                writefile.writelines(f'{dbn_header} Knotty\n')
                # Knotty writes its own sequence line
                knotty_input = f'{knotty_program} {dbn_sequence}'
                logging.debug(f'Knotty input: {knotty_input}')
                writefile.flush()
                # Flushing the buffer helps keep DBNs properly aligned
                knotty = subprocess.run(knotty_input, stdout=writefile, shell=True)
                # Write HFold motif
                writefile.writelines(f'{dbn_header} HFold\n')
                # HFold writes its own sequence line
                # HFold only structure (blank)
                blank_structure = str('_' * len(dbn_sequence))
                # HFold constraints require '_' instead of '.' to fold properly
                hfold_blank_input = f'{hfold_program} --s "{dbn_sequence}" --r "{blank_structure}"'
                logging.debug(f'HFold input: {hfold_blank_input}')
                writefile.flush()
                hfold_only = subprocess.run(hfold_blank_input, stdout=writefile, shell=True)
                # Write HFold with ScanFold constraints
                writefile.writelines(f'{dbn_header} HFold ScanFold\n')
                constraints = dbn_structure.replace('.', '_')
                hfold_input = f'{hfold_program} --s "{dbn_sequence}" --r "{constraints}"'
                logging.debug(f'HFold w/ ScanFold input: {hfold_input}')
                writefile.flush()
                hfold = subprocess.run(hfold_input, stdout=writefile, shell=True)
                i += 1
            else:
                i += 1


def pk_cleanup(pk_readfile='tmppk.txt', pk_writefile='pk_clean.txt'):
    # Clean up pseudoknot file for easier visual comparison and breakdown
    i = 0
    with open(pk_writefile, 'w', newline='\n') as writefile, open(pk_readfile, 'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
            if line.rstrip().endswith(' HFold ScanFold'):
                writefile.writelines(lines[i].rstrip() + '\n')
                try:
                    sequence = lines[i + 5].rstrip().split(' ')
                    writefile.writelines(sequence[1] + '\n')
                except:
                    writefile.writelines('HFold ScanFold Sequence ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                try:
                    structure = lines[i + 6].rstrip().split(' ')
                    writefile.writelines(structure[1] + '\n')
                except:
                    writefile.writelines('HFold ScanFold Structure ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                writefile.writelines('\n')
                i += 1
            elif line.rstrip().endswith(' Knotty'):
                writefile.writelines(lines[i].rstrip() + '\n')
                try:
                    sequence = lines[i + 1].rstrip().split(' ')
                    writefile.writelines(sequence[1] + '\n')
                except:
                    writefile.writelines('Knotty Sequence ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                try:
                    structure = lines[i + 2].rstrip().split(' ')
                    writefile.writelines(structure[1] + '\n')
                except:
                    writefile.writelines('Knotty Structure ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                i += 1
            elif line.rstrip().endswith(' HFold'):
                writefile.writelines(lines[i].rstrip() + '\n')
                try:
                    sequence = lines[i + 5].rstrip().split(' ')
                    writefile.writelines(sequence[1] + '\n')
                except:
                    writefile.writelines('HFold Sequence ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                try:
                    structure = lines[i + 6].rstrip().split(' ')
                    writefile.writelines(structure[1] + '\n')
                except:
                    writefile.writelines('HFold Structure ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                i += 1
            elif line.rstrip().endswith(' ScanFold'):
                writefile.writelines(lines[i].rstrip() + '\n')
                try:
                    writefile.writelines(lines[i + 1].rstrip() + '\n')
                except:
                    writefile.writelines('ScanFold Sequence ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                try:
                    writefile.writelines(lines[i + 2].rstrip() + '\n')
                except:
                    writefile.writelines('ScanFold Structure ERROR\n')
                    logging.error(f'PK cleanup error: {lines[i]}')
                i += 1
            else:
                i += 1


def pk_breakdown(pk_directory, pk_readfile='pk_clean.txt'):
    # Takes clean DBN file (pk_cleanup) and converts into separate DBN files
    # Breaks pseudoknots into nested (non-pseudoknotted) motifs
    i = 0
    with open(pk_readfile, 'r') as readfile:
        lines = readfile.readlines()
        for line in lines:
            if line[0] == '>':
                dbn_header = lines[i].rstrip()
                # Preserve naming conventions
                if lines[i + 1].endswith('ERROR') or lines[i + 2].endswith('ERROR'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = dbn_header
                    pk_filename = pk_filename[0] + '_error.dbn'
                elif dbn_header.endswith(' HFold ScanFold'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>' + pk_filename[0] + ' HFold ScanFold'
                    pk_filename = pk_filename[0] + '_hfold_scanfold.dbn'
                elif dbn_header.endswith(' Knotty'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>' + pk_filename[0] + ' Knotty'
                    pk_filename = pk_filename[0] + '_knotty.dbn'
                elif dbn_header.endswith(' HFold'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>' + pk_filename[0] + ' HFold'
                    pk_filename = pk_filename[0] + '_hfold.dbn'
                elif dbn_header.endswith(' ScanFold'):
                    pk_filename = dbn_header.strip('>').split('_coordinates')
                    pk_header = '>' + pk_filename[0] + ' ScanFold'
                    pk_filename = pk_filename[0] + '_scanfold.dbn'
                pk_sequence = lines[i + 1].rstrip()
                pk_structure = lines[i + 2].rstrip()
                # Output error files, no need to check for pseudoknots
                if pk_sequence.endswith('ERROR') or pk_sequence.endswith('source') or pk_structure.endswith(
                        'ERROR') or pk_structure.endswith('HFold'):
                    with open(Path.joinpath(pk_directory, pk_filename), 'w', newline='\n') as writefile:
                        writefile.writelines(pk_header + '\n')
                        writefile.writelines(pk_sequence + '\n')
                        writefile.writelines(pk_structure + '\n')
                    i += 1
                # Check for pseudoknots, if none then output results
                elif not pk_istrue(pk_structure):
                    # No pseudoknot is present, so default all non-pseudoknot structures to open/close parentheses
                    pk_structure = pk_structure.replace(
                        '[', '(').replace(']', ')').replace('{', '(').replace('}', ')').replace('<', '(').replace('>',
                                                                                                                  ')')
                    with open(Path.joinpath(pk_directory, pk_filename), 'w', newline='\n') as writefile:
                        # Remove unpaired sequence from motifs prior to writing
                        left_pos = 0
                        right_pos = len(pk_sequence)
                        left_brackets = ('(', '[', '{', '<')
                        right_brackets = (')', ']', '}', '>')
                        while True:
                            if left_pos == right_pos:
                                left_pos = 0
                                break
                            elif pk_structure[left_pos] in left_brackets:
                                break
                            else:
                                left_pos += 1
                        while True:
                            if left_pos == right_pos:
                                right_pos = len(pk_sequence)
                                break
                            elif pk_structure[right_pos - 1] in right_brackets:
                                break
                            else:
                                right_pos -= 1
                        writefile.writelines(pk_header + '\n')
                        writefile.writelines(pk_sequence[left_pos:right_pos] + '\n')
                        writefile.writelines(pk_structure[left_pos:right_pos] + '\n')
                    i += 1
                else:
                    # Pseudoknot is present, will need to be broken down
                    pk_splitter(pk_directory, pk_filename, pk_header, pk_sequence, pk_structure)
                    i += 1
            else:
                i += 1


def pk_istrue(structure):
    # Will check for presence of pseudoknots, returning either True or False
    # Able to search 4 levels deep (),[],{},<> which is assumed to be enough depth for <200 nt structures
    # First eliminate any sequences with only one bracket type (can't be a pseudoknot structure)
    bracket1 = structure.count('(')
    bracket2 = structure.count('[')
    bracket3 = structure.count('{')
    bracket4 = structure.count('<')
    if bracket1 >= 0 and bracket2 == 0 and bracket3 == 0 and bracket4 == 0:
        return False
    elif bracket1 == 0 and bracket2 > 0 and bracket3 == 0 and bracket4 == 0:
        return False
    elif bracket1 == 0 and bracket2 == 0 and bracket3 > 0 and bracket4 == 0:
        return False
    elif bracket1 == 0 and bracket2 == 0 and bracket3 == 0 and bracket4 > 0:
        return False
    # Anything left will be checked for non-nested structures
    else:
        ijcoords = {}
        icoords = []
        jcoords = []
        klcoords = {}
        kcoords = []
        lcoords = []
        mncoords = {}
        mcoords = []
        ncoords = []
        opcoords = {}
        ocoords = []
        pcoords = []
        column = 0
        for character in structure:
            # Build lists of coordinates (i, j, etc.)
            if character == '(':
                icoords.append(column)
                column += 1
            elif character == ')':
                jcoords.append(column)
                column += 1
            elif character == '[':
                kcoords.append(column)
                column += 1
            elif character == ']':
                lcoords.append(column)
                column += 1
            elif character == '{':
                mcoords.append(column)
                column += 1
            elif character == '}':
                ncoords.append(column)
                column += 1
            elif character == '<':
                ocoords.append(column)
                column += 1
            elif character == '>':
                pcoords.append(column)
                column += 1
            else:
                column += 1
        for jcoord in jcoords:
            # Finds matching i,j coordinates by taking the left-most (lowest) j value
            # and pairing it with the right-most (highest) i value that is less than j
            # then removes that value from the list
            for icoord in icoords[::-1]:
                if int(icoord) < int(jcoord):
                    ijcoords.update({icoord: jcoord})
                    icoords.remove(icoord)
                    break
        for lcoord in lcoords:
            for kcoord in kcoords[::-1]:
                if int(kcoord) < int(lcoord):
                    klcoords.update({kcoord: lcoord})
                    kcoords.remove(kcoord)
                    break
        for ncoord in ncoords:
            for mcoord in mcoords[::-1]:
                if int(mcoord) < int(ncoord):
                    mncoords.update({mcoord: ncoord})
                    mcoords.remove(mcoord)
                    break
        for pcoord in pcoords:
            for ocoord in ocoords[::-1]:
                if int(ocoord) < int(pcoord):
                    opcoords.update({ocoord: pcoord})
                    ocoords.remove(ocoord)
                    break
        for i, j in ijcoords.items():
            # Checks for presence of pseudoknots in all possible base pair bracket configurations
            for k, l in klcoords.items():
                if i < k < j < l:
                    return True
                elif k < i < l < j:
                    return True
                else:
                    continue
            for m, n in mncoords.items():
                if i < m < j < n:
                    return True
                elif m < i < n < j:
                    return True
                else:
                    continue
            for o, p in opcoords.items():
                if i < o < j < p:
                    return True
                elif o < i < p < j:
                    return True
                else:
                    continue
        for k, l in klcoords.items():
            for m, n in mncoords.items():
                if k < m < l < n:
                    return True
                elif m < k < n < l:
                    return True
                else:
                    continue
            for o, p in opcoords.items():
                if k < o < l < p:
                    return True
                elif o < k < p < l:
                    return True
                else:
                    continue
        for m, n in mncoords.items():
            for o, p in opcoords.items():
                if m < o < n < p:
                    return True
                elif o < m < p < n:
                    return True
                else:
                    continue
        # If the script reaches this line, no pseudoknots were found (structure is nested)
        return False


def pk_splitter(dbn_directory, dbn_filename, dbn_header, dbn_sequence, dbn_structure):
    # Takes DBN data and breaks it down into multiple nested structures
    bracket1 = dbn_structure.count('(')
    bracket2 = dbn_structure.count('[')
    bracket3 = dbn_structure.count('{')
    bracket4 = dbn_structure.count('<')
    count = 0
    dbn_filename = dbn_filename.rsplit('.dbn')
    dbn_filename = dbn_filename[0]
    if bracket1 > 0:
        # Add to count so filenames don't override each other
        count += 1
        dbn_filename1 = str(dbn_filename) + '_' + str(count) + '.dbn'
        # Defaults all extracted structures to open/close parentheses
        dbn_structure1 = dbn_structure.replace('[', '.').replace(']', '.').replace('{', '.').replace('}', '.').replace(
            '<', '.').replace('>', '.')
        with open(Path.joinpath(dbn_directory, dbn_filename1), 'w', newline='\n') as writefile:
            writefile.writelines(dbn_header + ' ' + str(count) + '\n')
            # Remove unpaired sequence from motifs prior to writing
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = '('
            right_brackets = ')'
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure1[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure1[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence1 = dbn_sequence[left_pos:right_pos]
            dbn_structure1 = dbn_structure1[left_pos:right_pos]
            writefile.writelines(dbn_sequence1 + '\n')
            writefile.writelines(dbn_structure1)
    if bracket2 > 0:
        count += 1
        dbn_filename2 = str(dbn_filename) + '_' + str(count) + '.dbn'
        dbn_structure2 = dbn_structure.replace('(', '.').replace(')', '.').replace('{', '.').replace('}', '.').replace(
            '<', '.').replace('>', '.')
        dbn_structure2 = dbn_structure2.replace('[', '(').replace(']', ')')
        with open(Path.joinpath(dbn_directory, dbn_filename2), 'w', newline='\n') as writefile:
            writefile.writelines(dbn_header + ' ' + str(count) + '\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = '('
            right_brackets = ')'
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure2[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure2[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence2 = dbn_sequence[left_pos:right_pos]
            dbn_structure2 = dbn_structure2[left_pos:right_pos]
            writefile.writelines(dbn_sequence2 + '\n')
            writefile.writelines(dbn_structure2)
    if bracket3 > 0:
        count += 1
        dbn_filename3 = str(dbn_filename) + '_' + str(count) + '.dbn'
        dbn_structure3 = dbn_structure.replace('[', '.').replace(']', '.').replace('(', '.').replace(')', '.').replace(
            '<', '.').replace('>', '.')
        dbn_structure3 = dbn_structure3.replace('{', '(').replace('}', ')')
        with open(Path.joinpath(dbn_directory, dbn_filename3), 'w', newline='\n') as writefile:
            writefile.writelines(dbn_header + ' ' + str(count) + '\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = '('
            right_brackets = ')'
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure3[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure3[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence3 = dbn_sequence[left_pos:right_pos]
            dbn_structure3 = dbn_structure3[left_pos:right_pos]
            writefile.writelines(dbn_sequence3 + '\n')
            writefile.writelines(dbn_structure3)
    if bracket4 > 0:
        count += 1
        dbn_filename4 = str(dbn_filename) + '_' + str(count) + '.dbn'
        dbn_structure4 = dbn_structure.replace('[', '.').replace(']', '.').replace('{', '.').replace('}', '.').replace(
            '(', '.').replace(')', '.')
        dbn_structure4 = dbn_structure4.replace('<', '(').replace('>', ')')
        with open(Path.joinpath(dbn_directory, dbn_filename4), 'w', newline='\n') as writefile:
            writefile.writelines(dbn_header + ' ' + str(count) + '\n')
            left_pos = 0
            right_pos = len(dbn_sequence)
            left_brackets = '('
            right_brackets = ')'
            while True:
                if left_pos == right_pos:
                    left_pos = 0
                    break
                elif dbn_structure4[left_pos] in left_brackets:
                    break
                else:
                    left_pos += 1
            while True:
                if left_pos == right_pos:
                    right_pos = len(dbn_sequence)
                    break
                elif dbn_structure4[right_pos - 1] in right_brackets:
                    break
                else:
                    right_pos -= 1
            dbn_sequence4 = dbn_sequence[left_pos:right_pos]
            dbn_structure4 = dbn_structure4[left_pos:right_pos]
            writefile.writelines(dbn_sequence4 + '\n')
            writefile.writelines(dbn_structure4)


def cmbuilder_prep(sequence_directory, database_directory, dbn_directory, cmbuilder_program, email, rscape_program,
                   r2r_program, cobretti_program, perl_program, rnaframework_directory, cm_writefile='cmbuilder'):
    # Reads all DBN files and outputs a shell script for cm-builder, based on Van's code
    count = 1
    current_size = 0
    max_size = 10  # Number of runs per shell script
    max_db_size = 8000000000  # Files over 8GB are run separately to aid with HPC out of memory errors
    for dbn_filepath in dbn_directory.glob('*_error.dbn'):
        # Exit if any motifs with known errors are present
        logging.error(f'Error with DBN {dbn_filepath}, check error and then run stage 1BC')
        sys.exit()
    for dbn_filepath in dbn_directory.glob('*.dbn'):
        logging.debug(f'DBN file: {dbn_filepath.name}')
        gene = dbn_filepath.name.split('_')[0]
        logging.debug(f'Motif count: {current_size}')
        if current_size == 0:
            # Build initial shell scripts
            shell_build_start(f'{cm_writefile}_{count}.sh', f'{cm_writefile}_{count}', email, mem=100, tasks=20,
                              notify='FAIL')
            with open(f'{cm_writefile}_{count}.sh', 'a', newline='\n') as writefile:
                # Add RNAFramework libraries
                writefile.writelines(f'export PERL5LIB={perl_program}:{rnaframework_directory}:${{PERL5LIB}}\n')

        with open(f'{cm_writefile}_{count}.sh', 'a', newline='\n') as writefile:
            # Find matching sequence, structure, database files
            for seq_filepath in sequence_directory.glob(gene + '*.fa*'):
                logging.debug(f'SEQ file: {seq_filepath.name} MATCH')
                for db_filepath in database_directory.glob(gene + '*_db.fa*'):
                    logging.debug(f'DB file: {db_filepath.name} MATCH')
                    db_size = int(os.path.getsize(db_filepath))
                    logging.debug(f'DB size: {db_size}\tMax size: {max_db_size}')
                    # Check size and separate out large database files
                    if db_size > max_db_size:
                        # Check oversize number until it doesn't exist and create that file
                        oversize_db_count = 1
                        while os.path.isfile(f'{cm_writefile}_{gene}_{oversize_db_count}.sh'):
                            oversize_db_count += 1
                        shell_build_start(f'{cm_writefile}_{gene}_{oversize_db_count}.sh',
                                          f'{cm_writefile}_{gene}_{oversize_db_count}', email,
                                          mem=100, tasks=20, notify='FAIL')
                        with open(f'{cm_writefile}_{gene}_{oversize_db_count}.sh', 'a',
                                  newline='\n') as oversize_writefile:
                            oversize_writefile.writelines(
                                f'perl {cmbuilder_program} -s {seq_filepath} -m {dbn_filepath} -d {db_filepath} -c 4 '
                                f'-M ./tmp{count} -k -t 1 &\n')
                            oversize_writefile.writelines('wait;\n')
                    else:
                        writefile.writelines(
                            f'perl {cmbuilder_program} -s {seq_filepath} -m {dbn_filepath} -d {db_filepath} -c 4 -M '
                            f'./tmp{count} -k -t 1 &\n')
                        current_size += 1
        if current_size >= max_size:
            # Close out shell scripts when full, increase shell count and reset size count
            with open(f'{cm_writefile}_{count}.sh', 'a', newline='\n') as writefile:
                writefile.writelines('wait;\n')
            current_size = 0
            count += 1
    with open(f'{cm_writefile}_{count}.sh', 'a', newline='\n') as writefile:
        # Close out the last shell script
        writefile.writelines('wait;\n')
    shell_build_start('rscape.sh', 'rscape', email, time=1, notify='END,FAIL')
    # Prepare R-Scape shell script
    with open('rscape.sh', 'a', newline='\n') as writefile:
        rscape_runs = 10
        writefile.writelines(f'for f in *.stockholm; do {rscape_program} -s --ntree {rscape_runs} $f; done\n')
        #writefile.writelines(f'for f in *.stockholm; do {rscape_program} -s $f; done\n')
        # Back up Stockholm files, then delete lines that error out R2R
        writefile.writelines('sed -i.bak "/#=GF R2R*/d" *.sto\n')
        writefile.writelines(
            f'for g in *.sto; do {r2r_program} --disable-usage-warning $g $(basename $g sto)pdf; done\n')
        # Write all R2R outputs to a single PDF
        writefile.writelines('gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=All.Rscape.pdf -dBATCH *.R2R.pdf\n')
        writefile.writelines('wait;\n')
        # Since code is linear, automatically start the next substage
        writefile.writelines(f'module load cobretti\n') ##########################################################
        writefile.writelines(f'python {cobretti_program} -stage 1CA -email {email} &\n')
        writefile.writelines('wait;\n')


def cmbuilder_run(working_directory):
    # Run cm-builder shell scripts
    for filepath in working_directory.glob('cmbuilder_*.sh'):
        logging.debug(f'cm-builder shell path: {filepath}')
        os.system(f'sbatch {filepath.name}')


def cmbuilder_cleanup(working_directory):
    # Clean up working directory, run R-Scape
    folders = ['cm', 'dbn', 'out', 'sh']
    for folder in folders:
        # Create the directory if it doesn't exist
        Path.mkdir(Path.joinpath(working_directory, folder), exist_ok=True)
    for filepath in working_directory.glob('*'):
        for folder in folders:
            if filepath.name != 'rscape.sh':
                # Ignore the R-Scape shell script as it will be needed later
                try:
                    if os.path.isdir(filepath):
                        # Delete temporary directories
                        if filepath.name.startswith('tmp'):
                            shutil.rmtree(filepath)
                    else:
                        # Move files based on extension
                        if filepath.name.endswith(folder):
                            shutil.move(filepath, Path.joinpath(working_directory, folder, filepath.name))
                except:
                    logging.error(f'Unable to move file: {filepath.name}')
    os.system('sbatch rscape.sh')
    # Automatically starts cmbuilder_cleanup2() upon completion of shell script


def cmbuilder_cleanup2(working_directory, cm_writefile='covariance.txt'):
    # Move remaining files into folders, write covariance power and output to file
    covariance_read(cm_writefile, working_directory)
    folders = ['bak', 'cm', 'cov', 'dbn', 'out', 'pdf', 'power', 'ps', 'sh', 'sto', 'stockholm', 'surv', 'svg']
    for folder in folders:
        Path.mkdir(Path.joinpath(working_directory, folder), exist_ok=True)
    for filepath in working_directory.glob('*.*'):
        for folder in folders:
            # Ignore directories and R-Scape PDF output
            if not filepath.name.endswith('Rscape.pdf'):
                if filepath.name.endswith(folder):
                    try:
                        shutil.move(filepath, Path.joinpath(working_directory, folder, filepath.name))
                    except:
                        logging.error(f'Unable to move file: {filepath.name}')


def covariance_read(cm_writefile, working_directory):
    # Read power files and output all power results to single file - Adapted from Collin's script
    with open(cm_writefile, 'w', newline='\n') as writefile:
        writefile.writelines(
            f'Name\tBPs\tavg_substitutions\texpected BPs\tSTD_DEV\tObserved BPs\t# of BP power 0-0.1\t# of BP power '
            f'0.1-0.25\t# of BP power >=0.25\t0-0.1 BP info\t0.1-0.25 BP info\t>=0.25 BP info\n')
        for filepath in working_directory.glob('*.power'):
            with open(filepath, 'r') as readfile:
                name = filepath.name.replace("_1.power", "")
                low_cv = []
                mid_cv = []
                high_cv = []
                lines = readfile.readlines()
                for line in lines:
                    if len(line) <= 1:
                        pass
                    elif line.startswith('# BPAIRS observed to covary'):
                        obv_bp = line.rstrip().split(' ')[5]
                    elif line.startswith('# BPAIRS expected to covary'):
                        ex_bp = line.rstrip().split(' ')[5]
                        std_dev = line.rstrip().split(' ')[5]
                    elif line.startswith('# avg substitutions per BP'):
                        subs = line.rstrip().split(' ')[6]
                    elif line.startswith('# BPAIRS'):
                        bps = line.rstrip().split(' ')[2]
                    elif line.startswith('     *'):
                        power = line.rsplit()
                        if float(power[-1]) >= 0.25:
                            high_cv.append((int(power[1]), int(power[2]), float(power[-1])))
                        elif float(power[-1]) >= 0.1:
                            mid_cv.append((int(power[1]), int(power[2]), float(power[-1])))
                        else:
                            low_cv.append((int(power[1]), int(power[2]), float(power[-1])))
                    else:
                        pass
                low_num = len(low_cv)
                mid_num = len(mid_cv)
                high_num = len(high_cv)
                writefile.writelines(f'{name}\t{bps}\t{subs}\t{ex_bp}\t{std_dev}\t{obv_bp}\t{low_num}\t'
                                     f'{mid_num}\t{high_num}\t{low_cv}\t{mid_cv}\t{high_cv}\n')


def simrna_prep(working_directory, simrna_directory, email, config_file='None', instances=8, replicas=10,
                simrna_writefile='simrna'):
    # Writes shell scripts for optimized SimRNA runs
    simrna_program = Path.joinpath(simrna_directory, 'SimRNA')
    simrna_data = Path.joinpath(simrna_directory, 'data')
    simrna_clustering = Path.joinpath(simrna_directory, 'clustering')
    simrna_trafl2pdb = Path.joinpath(simrna_directory, 'SimRNA_trafl2pdbs')
    if config_file == 'None':
        simrna_config = Path.joinpath(simrna_directory, 'config.dat')
    else:
        simrna_config = config_file
    data_link = Path.joinpath(working_directory, 'data')
    if not Path.exists(data_link):
        os.symlink(simrna_data, data_link)
    file_count = 1
    for filepath in working_directory.glob('*.dbn'):
        dbn_headers = Path(filepath).stem
        # Split sequence and structure into separate inputs so SimRNA can read them
        sequence_tempfile = dbn_headers + '.sequence'
        structure_tempfile = dbn_headers + '.structure'
        rmsd = [5.0, 7.0, 10.0]
        with open(filepath, 'r') as readfile:
            with open(sequence_tempfile, 'w', newline='\n') as writefile:
                line = readfile.readlines()
                writefile.write(line[1])
                sequence_length = len(line[1])
                rmsd.append(float(sequence_length / 10.0))
                logging.debug(f'RMSD: {rmsd}')
            with open(structure_tempfile, 'w', newline='\n') as writefile:
                # Break down non-nested structures
                pk_count = line[2].count('[') + line[2].count('{') + line[2].count('<')
                if pk_count == 0:
                    # If there are no pseudoknots, output structure
                    writefile.write(line[2])
                else:
                    # If there are pseudoknots, break structures down by bracket
                    newline = line[2].replace('[', '.').replace(']', '.').replace('{', '.').replace('}',
                                                                                                    '.').replace(
                        '<', '.').replace('>', '.').rstrip()
                    writefile.write(newline + '\n')
                    if line[2].count('[') != 0:
                        newline = line[2].replace('(', '.').replace(')', '.').replace('{', '.').replace('}',
                                                                                                        '.').replace(
                            '<', '.').replace('>', '.').replace('[', '(').replace(']', ')').rstrip()
                        writefile.write(newline + '\n')
                    if line[2].count('{') != 0:
                        newline = line[2].replace('(', '.').replace(')', '.').replace('[', '.').replace(']',
                                                                                                        '.').replace(
                            '<', '.').replace('>', '.').replace('{', '(').replace('}', ')').rstrip()
                        writefile.write(newline + '\n')
                    if line[2].count('<') != 0:
                        newline = line[2].replace('(', '.').replace(')', '.').replace('{', '.').replace('}',
                                                                                                        '.').replace(
                            '[', '.').replace(']', '.').replace('<', '(').replace('>', ')').rstrip()
                        writefile.write(newline + '\n')
        shell_build_start(f'{simrna_writefile}_{file_count}.sh', f'{simrna_writefile}_{file_count}',
                          email, tasks=50, notify='END,FAIL')
        with open(f'{simrna_writefile}_{file_count}.sh', 'a', newline='\n') as writefile:
            # Prepare shell scripts, multiple instances changes the naming convention
            if instances == 1:
                writefile.writelines(
                    f'{simrna_program} -c {simrna_config} -E {replicas} -s {sequence_tempfile} -S {structure_tempfile}'
                    f' -o {dbn_headers} -R {random.randint(0, 9999999999)} >& {dbn_headers}.log &\n')
            elif instances < 1:
                logging.error('Number of instances is less than one')
            elif instances > 100:
                logging.error('Number of instances is too big')
            else:
                for instance_count in range(instances):
                    writefile.writelines(
                        f'{simrna_program} -c {simrna_config} -E {replicas} -s '
                        f'{sequence_tempfile} -S {structure_tempfile} -o {dbn_headers}')
                    if instance_count < 10:
                        writefile.writelines(
                            f'_0{instance_count} -R {random.randint(0, 9999999999)}'
                            f' >& {dbn_headers}_0{instance_count}.log &\n')
                    else:
                        writefile.writelines(
                            f'_{instance_count} -R {random.randint(0, 9999999999)}'
                            f' >& {dbn_headers}_{instance_count}.log &\n')
            writefile.writelines('wait;\n')
            if instances != 1:
                # Concatenate all instances into a single trajectory file
                writefile.writelines(f'cat {dbn_headers}_??')
                if replicas != 1:
                    writefile.writelines('_??')
                else:
                    # SimRNA dds '119' to single replicas
                    writefile.writelines('_119')
                writefile.writelines(f'.trafl > {dbn_headers}_all.trafl &\n')
                writefile.writelines('wait;\n')
                # Cluster lowest 1% energy frames
                writefile.writelines(f'{simrna_clustering} {dbn_headers}_all.trafl 0.01 ')
                for value in rmsd:
                    # Add RMSD values for clustering
                    writefile.writelines(f'{value} ')
                writefile.writelines(f' >& {dbn_headers}_clust.log &\n')
                writefile.writelines('wait;\n')
                for value in rmsd:
                    for rank in range(1, 6):
                        # Output all-atom PDB files for top 5 clusters for each RMSD
                        writefile.writelines(
                            f'{simrna_trafl2pdb} {dbn_headers}_01_01-000001.pdb '
                            f'{dbn_headers}_all_thrs{value}0A_clust0{rank}.trafl 1 AA &\n')
                writefile.writelines('wait;\n')
            else:
                if replicas != 1:
                    writefile.writelines(f'cat {dbn_headers}_??.trafl > {dbn_headers}_all.trafl &\n')
                    writefile.writelines('wait;\n')
                writefile.writelines(f'{simrna_clustering} {dbn_headers}')
                if replicas != 1:
                    writefile.writelines('_all')
                else:
                    writefile.writelines('_119')
                writefile.writelines('.trafl 0.01 ')
                for value in rmsd:
                    writefile.writelines(f'{value} ')
                writefile.writelines(f' >& {dbn_headers}_clust.log &\n')
                writefile.writelines('wait;\n')
                for value in rmsd:
                    for rank in range(1, 6):
                        writefile.writelines(
                            f'{simrna_trafl2pdb} {dbn_headers}_119-000001.pdb '
                            f'{dbn_headers}_119_thrs{value}0A_clust0{rank}.trafl 1 AA &\n')
                writefile.writelines('wait;\n')
        file_count += 1


def simrna_run(working_directory):
    for filepath in working_directory.glob('simrna*.sh'):
        os.system(f'sbatch {filepath.name}')


def simrna_cleanup(working_directory):
    folders = ['out', 'sh', 'log', 'trafl', 'motifs', 'bonds', 'ss_detected', 'models', 'cluster_models']
    for folder in folders:
        # Create directories as necessary
        Path.mkdir(Path.joinpath(working_directory, folder), exist_ok=True)
    for filepath in working_directory.glob('*.*'):
        # Ignore directories, move files that don't match directory names
        if filepath.name.endswith('_clust.log') or filepath.name.endswith('_AA.pdb') or (
                filepath.name.endswith('.ss_detected') and filepath.name.__contains__('clust')):
            try:
                shutil.move(filepath, Path.joinpath(working_directory, 'cluster_models', filepath.name))
            except:
                logging.error(f'Unable to move file: {filepath.name}')
        elif (filepath.name.endswith('.structure') or filepath.name.endswith('.sequence') or filepath.name.endswith(
                '.dbn')):
            try:
                shutil.move(filepath, Path.joinpath(working_directory, 'motifs', filepath.name))
            except:
                logging.error(f'Unable to move file: {filepath.name}')
        elif filepath.name.endswith('.pdb'):
            try:
                shutil.move(filepath, Path.joinpath(working_directory, 'models', filepath.name))
            except:
                logging.error(f'Unable to move file: {filepath.name}')
        for folder in folders:
            # Move remaining files
            if filepath.name.endswith(folder):
                try:
                    shutil.move(filepath, Path.joinpath(working_directory, folder, filepath.name))
                except:
                    logging.error(f'Unable to move file: {filepath.name}')


def qrnas_prep(working_directory, qrnas_program, qrnas_ff_directory, email):
    # Prepare QRNAS shell scripts
    count = 1
    current_size = 0
    max_size = 10
    cluster_directory = folder_check(working_directory, 'cluster_models')
    qrnas_directory = Path.joinpath(working_directory, 'qrnas_models')
    Path.mkdir(qrnas_directory, exist_ok=True)
    for filepath in cluster_directory.glob('*_AA.pdb'):
        output_file = filepath.name.replace('.pdb', '_QRNAS.pdb')
        if current_size == 0:
            # Build initial shell scripts
            shell_build_start(f'qrnas_{count}.sh', f'qrnas_{count}', email, time=7, tasks=10, notify='END,FAIL')
            with open(f'qrnas_{count}.sh', 'a', newline='\n') as writefile:
                # Make and use temporary local directories to improve runtime
                writefile.writelines('mkdir $TMPDIR/qrnas_${SLURM_JOB_ID}\n')
                writefile.writelines('cd $TMPDIR/qrnas_${SLURM_JOB_ID}\n')
                writefile.writelines(f'export QRNAS_FF_DIR={qrnas_ff_directory}\n')
        with open(f'qrnas_{count}.sh', 'a', newline='\n') as writefile:
            writefile.writelines(f'{qrnas_program} -i {filepath} -o {output_file} &\n')
            current_size += 1
        if current_size >= max_size:
            # Close out shell script when full, reset count
            with open(f'qrnas_{count}.sh', 'a', newline='\n') as writefile:
                writefile.writelines('wait;\n')
                # Copy files from temporary directory
                writefile.writelines('cp -r $TMPDIR/qrnas_${SLURM_JOB_ID} %s\n' % qrnas_directory)
                count += 1
                current_size = 0
    with open('qrnas_' + str(count) + '.sh', 'a', newline='\n') as writefile:
        # Close out final shell script, copy files from temporary directory
        writefile.writelines('wait;\n')
        writefile.writelines('cp -r $TMPDIR/qrnas_${SLURM_JOB_ID}/. %s\n' % qrnas_directory)


def qrnas_run(working_directory):
    # Run QRNAS scripts in working directory
    for filepath in working_directory.glob('qrnas*.sh'):
        os.system(f'sbatch {filepath.name}')


def qrnas_cleanup(working_directory):
    # Move QRNAS files from working directory
    folders = ['out', 'sh', 'log', 'qrnas_models']
    # Create directories as necessary
    for folder in folders:
        Path.mkdir(Path.joinpath(working_directory, folder), exist_ok=True)
    for filepath in working_directory.glob('*.*'):
        if filepath.name.endswith('_AA_QRNAS.pdb'):
            # Move QRNAS models
            try:
                shutil.move(filepath, Path.joinpath(working_directory, 'qrnas_models', filepath.name))
            except:
                logging.error(f'Unable to move file: {filepath.name}')
        for folder in folders:
            # Move everything else
            if filepath.name.endswith(folder):
                try:
                    shutil.move(filepath, Path.joinpath(working_directory, folder, filepath.name))
                except:
                    logging.error(f'Unable to move file: {filepath.name}')


def ares_prep(working_directory, ares_directory, ares_environment, email):
    # ARES runs on ALL PDBs within the working directory including subdirectories
    # If running ARES alone, isolate PDB files in a subdirectory to avoid running unoptimized models
    qrnas_directory = folder_check(working_directory, 'qrnas_models')
    shell_build_start('ares.sh', 'ares', email, nodes=8, mem=10, notify='END,FAIL')
    with open('ares.sh', 'a', newline='\n') as writefile:
        writefile.writelines(f'conda activate {ares_environment}\nwait;\n')
        writefile.writelines(f'cd {ares_directory}\n')
        writefile.writelines(
            f'python -m ares.predict {qrnas_directory} data/epoch=0-step=874.ckpt {working_directory}/ares.csv -f pdb '
            f'--nolabels\n')


def ares_run():
    os.system('sbatch ares.sh')


def fpocket_run(working_directory, fpocket_program):
    # Runs fpocket, no shell scripts needed
    qrnas_directory = folder_check(working_directory, 'qrnas_models')
    for filepath in qrnas_directory.glob('*_QRNAS.pdb'):
        os.system(f'{fpocket_program} -f {filepath}')


def fpocket_cleanup(working_directory):
    # Move pocket results into their own directory
    Path.mkdir(Path.joinpath(working_directory, 'pockets'), exist_ok=True)
    pocket_directory = Path.joinpath(working_directory, 'pockets')
    qrnas_directory = folder_check(working_directory, 'qrnas_models')
    for filename in os.listdir(qrnas_directory):  # Move pocket directories to pocket folder
        if filename.endswith('AA_QRNAS_out'):
            shutil.move(Path.joinpath(qrnas_directory, filename), Path.joinpath(pocket_directory, filename))


def fpocket_read(working_directory):
    # Read the fpocket results and collate data
    with open('pockets.txt', 'w', newline='\n') as writefile:
        writefile.writelines('Model Fold RMSD Cluster Pocket Score Druggability\n')
        pocket_directory = folder_check(working_directory, 'pockets')
        for filepath in pocket_directory.rglob('*_info.txt'):
            with open(filepath, 'r') as readfile:
                header = filepath.name.split('_')
                model = header[0]
                fold = ''
                rmsd = ''
                cluster = ''
                pocket_number = 0
                pocket_score = 0.0
                # Decrypt filename to extract RMSD, cluster info
                for i in range(1, 20):
                    # Find the 'all', may have shifted
                    try:
                        if header[i] == 'all':
                            # Everything before the 'all'
                            for j in range(1, i):
                                fold += header[j]
                        if header[i].startswith('thrs'):
                            # RMSD, cluster values
                            rmsd = header[i].split('thrs')[1]
                            cluster = header[i + 1].split('-')[0]
                    except:
                        logging.error(f'info.txt file format does not match: {filepath.name}')
                lines = readfile.readlines()
                for line in lines:
                    # Read the info.txt file, extract pocket data
                    if line.startswith('Pocket'):
                        pocket_number = int(line.split()[1])
                    if line.startswith('	Score'):
                        pocket_score = float(line.split()[2])
                    if line.startswith('	Druggability'):
                        drug_score = float(line.split()[3])
                        # Write the results once the druggability score is found
                        writefile.writelines(
                            f'{model} {fold} {rmsd} {cluster} {pocket_number} {pocket_score} {drug_score}\n')
                        # Empty the scores between writes
                        pocket_number = 0
                        pocket_score = 0.0
                        drug_score = 0.0


def dock6_prep():
    # TODO read SimRNA pdb and pocket pdb (might need dummy molecule?)
    # TODO Convert PDB to mol2 format, add hydrogens and charges
    # TODO make surface file (Google "Chimera Programmer's Guide"
    # TODO Find and center on pocket location
    # TODO ZINC15 database(s), copy database of choice to working directory? Or just reference copy
    # TODO Write shell script
    # TODO Write INSPH file
    logging.error('DOCK 6 is a WIP')


def dock6_run(working_directory):
    # Run DOCK 6
    for filepath in working_directory.glob('dock*.sh'):
        os.system(f'sbatch {filepath.name}')


def annapurna_run(working_directory):
    # TODO AnnapuRNA directory, read current directory, find dock results
    # TODO run AnnapuRNA on results
    # TODO module load rnaframe/2.7.2
    # TODO annapurna_env (might not need for shell)
    # TODO run_annapurna_env annapurna.py -r 1AJU.pdb -l ARG.sdf -m kNN_modern -o output --groupby
    logging.error('AnnapuRNA is a WIP')


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stdout,
                        format='%(asctime)s %(levelname)s %(process)d: %(message)s',
                        level=logging.INFO)
    main()
    logging.info('Cobretti run completed successfully!')
