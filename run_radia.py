#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_radia.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import division, print_function
from subprocess import call
from datetime import datetime as dt
from collections import defaultdict

import sys
import os
import prepare
import pi_errors

def test_database(params_db_value, db_name, radia_pkg_path, vcf_location,
                  data_subfolder, db_mapper):
    '''
    This test is used to search for a database (db_name) in the following
    locations in order of preference:
    1. In the parent directory of the radia executable (radia_path)
        radia is packaged as radia
                               |- data
                               |- scripts
                               |   |-radia.py
                               |   |-etc.
                               |
                               |- etc.
       The existence of the data_subfolder will be queried for within the data
       folder.
    2. VCF_LOCATION - provided by the user in the arguments will be queried for
       a folder called data_subfolder

    If the data_subfolder is found in neither location, the script returns a
    signal suggesting the data needs to be downloaded

    Module arguments
    db_name - The name of the database
    params_db_value - the value for the database db_name from params
    radia_path - The full path to the radia package (2 levels above radia.py)
    vcf_location - The value of params.vcf_location
    data_subfolder - The subfolder in the radia package where db_name is found
    db_mapper - an object of type collections.defaultdict() which will hold the
                value (or 'DOWNLOAD' if not present) of the db_name's path
    '''
    #  The default value if this
    if params_db_value == 'PACKAGED':
        #  Is it packaged with radia?
        if os.path.exists('/'.join([radia_pkg_path, data_subfolder])):
            db_mapper[db_name] = '/'.join([radia_pkg_path, data_subfolder])
        #  Is it in VCF_LOCATION?
        elif os.path.exists('/'.join([vcf_location, data_subfolder])):
            db_mapper[db_name] = '/'.join([vcf_location, data_subfolder])
        # If not either of the previous, download it
        else:
            db_mapper[db_name] = 'DOWNLOAD'
    else:
        if os.path.exists(params_db_value):
            db_mapper[db_name] = os.path.abspath(params_db_value)
        else:
            db_mapper[db_name] = 'DOWNLOAD'
    return None


def download_databases(db_map, logfile):
    '''
    Downloads the data folder from the radia package on git using sparse
    checkout (http://stackoverflow.com/questions/600079/is-there-any-way-to \
    -clone-a-git-repositorys-sub-directory-only) and sets the values of unmapped
    databases to the downloaded ones.
    '''
    #  Download the database
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Downloading radia data folder from github.', file=logfile)
    call(['git', 'init', 'temp_radia_dbs'])
    os.chdir('temp_radia_dbs')
    call(['git', 'remote', 'add', '-f', 'origin',
          'https://github.com/aradenbaugh/radia.git'])
    call(['git', 'config', 'core.sparseCheckout', 'true'])
    with open('.git/info/sparse-checkout', 'a') as sparse_config_file:
        print('data/', file=sparse_config_file)
    call(['git', 'pull', 'origin', 'master'])
    os.chdir('../')
    temp_db_dir = os.path.realpath('temp_radia_dbs/data')
    #  set the individual datbases to the downloaded ones
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Populating database map.', file=logfile)
    for db_name, val in db_map.items():
        if val == 'DOWNLOAD':
            if db_name == 'blacklist':
                db_map[db_name] = '/'.join([temp_db_dir,
                                            'hg19/blacklists/1000Genomes/' +
                                            'phase1/'])
            elif db_name == 'retrogenes':
                db_map[db_name] = '/'.join([temp_db_dir, 'hg19/retroGenes/'])
            elif db_name == 'pseudogenes':
                db_map[db_name] = '/'.join([temp_db_dir, 'hg19/peudoGenes/'])
            elif db_name == 'broad_targets':
                db_map[db_name] = '/'.join([temp_db_dir, 'hg19/broadTargets/'])
            elif db_name == 'rna_blacklist':
                db_map[db_name] = '/'.join([temp_db_dir,
                                            'rnaGeneBlacklist.tab'])
            elif db_name == 'rna_family_blacklist':
                db_map[db_name] = '/'.join([temp_db_dir,
                                            'rnaGeneFamilyBlacklist.tab'])
            else:
                if logfile.name != '<stderr>':
                    logfile.close()
                raise SystemExit('Shouldn\'t ever get here')
        else:
            pass
    return None


def process_parameters(params):
    '''
    This module conducts the error handling for all parmeters passed to the
    program.
    '''
    #  Does the provided radia binary provided exist?
    params.radia_executable = pi_errors.test_param_value(
        params.radia_executable, 'radia', '--radia', params.logfile)
    #  Setup filterRadia.py
    params.filter_radia_executable = '/'.join([os.path.split(
        params.radia_executable)[0], 'filterRadia.py'])
    params.filter_radia_executable = pi_errors.test_param_value(
        params.filter_radia_executable, 'filterradia', '--radia',
        params.logfile)
    #  Test input files
    params.tum_d_file = pi_errors.test_param_value(params.tum_d_file,
                                                   'Tumor DNA',
                                                   '--tum_dna_file',
                                                   params.logfile)
    params.norm_d_file = pi_errors.test_param_value(params.norm_d_file,
                                                    'Normal DNA',
                                                    '--norm_dna_file',
                                                    params.logfile)
    if params.tum_r_file is not None:
        params.tum_r_file = pi_errors.test_param_value(params.tum_r_file,
                                                       'Tumor RNA',
                                                       '--tum_rna_file',
                                                       params.logfile)

    #  If you don't have a reference, you need twoBitToFasta
    if params.index_location is None:
        params.tbtf_executable = pi_errors.test_param_value(
            params.tbtf_executable, 'twoBitToFa', '--twoBitToFa',
            params.logfile)
    #  Are dnsnp or cosmic vcf required?
    if params.dbsnp_file == 'DOWNLOAD' or params.cosmic_file == 'DOWNLOAD' or \
            params.genome_fasta == 'DOWNLOAD':
        # Ensure the vcf storage location has been provided
        if params.vcf_location is None:
            raise pi_errors.ParameterError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': --vcf_location cannot be empty if either --cosmic, ' + \
                '--dbsnp, or --genome_fasta are empty.', params.logfile)
        else:
            params.vcf_location = os.path.abspath(params.vcf_location)
        # Download dbsnp file if required
        if params.dbsnp_file == 'DOWNLOAD':
            if os.path.exists('/'.join([params.vcf_location, '00-All.vcf'])):
                params.dbsnp_file = '/'.join([params.vcf_location,
                                              '00-All.vcf'])
            else:
                params.dbsnp_file = prepare.download_vcf('dbsnp', params)
        # Download cosmic file if required
        if params.cosmic_file == 'DOWNLOAD':
            if os.path.exists('/'.join([params.vcf_location,
                                        'Cosmic_sorted.vcf'])):
                params.cosmic_file = '/'.join([params.vcf_location,
                                               'Cosmic_sorted.vcf'])
            else:
                params.cosmic_file = prepare.download_vcf('cosmic', params)
        # Download genome fasta if required
        if params.genome_fasta == 'DOWNLOAD' or not \
                os.path.exists(params.genome_fasta):
            if os.path.exists(''.join([params.vcf_location, '/',
                                       params.genome_version, '.fa'])):
                params.genome_fasta = ''.join([params.vcf_location, '/',
                                               params.genome_version, '.fa'])
            else:
                params.genome_fasta = prepare.get_genome(
                    params.genome_version, params.vcf_location,
                    params.twoBitToFa_executable, sys.stderr)
        else:
            params.genome_fasta = os.path.abspath(params.genome_fasta)
    #  Set up the value for rna_fasta
    if params.rna_fasta == 'GENOME_FASTA':
        params.rna_fasta = params.genome_fasta
    else:
        params.rna_fasta = pi_errors.test_param_value(params.rna_fasta,
                                                      'RNA Fasta',
                                                      '--rna_fasta',
                                                      params.logfile)
    #  Ensure the other databases are set up correctly
    #  The package path is 2 levels above the
    radia_pkg_path = os.path.split(os.path.split(params.radia_executable)[0])[0]
    database_map = defaultdict()
    test_database(params.blacklist, 'blacklist', radia_pkg_path,
                  params.vcf_location, 'data/hg19/blacklists/1000Genomes/' + \
                  'phase1/', database_map)
    test_database(params.retrogenes, 'retrogenes', radia_pkg_path,
                  params.vcf_location, 'data/hg19/retroGenes/', database_map)
    test_database(params.pseudogenes, 'pseudogenes', radia_pkg_path,
                  params.vcf_location, 'data/hg19/peudoGenes/', database_map)
    test_database(params.broad_targets, 'broad_targets', radia_pkg_path,
                  params.vcf_location, 'data/hg19/broadTargets/',
                  database_map)
    test_database(params.rna_blacklist, 'rna_blacklist', radia_pkg_path,
                  params.vcf_location, 'data/rnaGeneBlacklist.tab',
                  database_map)
    test_database(params.rna_family_blacklist, 'rna_family_blacklist',
                  radia_pkg_path, params.vcf_location,
                  'data/rnaGeneFamilyBlacklist.tab', database_map)
    #  If any of the above were returned as 'DOWNLOAD' then download the radia
    #  data folder to a temp directory from git and set the values for the
    #  invalid ones.
    if len([db for db, val in database_map.items() if val == 'DOWNLOAD']) > 0:
        download_databases(database_map, params.logfile)
    # if the -C all option was specified, expand params.chromosome
    if params.chromosome == 'all':
        params.chromosome = [''.join(['chr', str(i)]) for i in \
                             range(1, 23)+['X', 'Y']]
    return database_map


def main():
    """
    This wrapper script will run the tool radia within the  radia docker
    container for the precision immuno project. The wrapper requires
    1. radia
    2. snpeff (if the --use_snpeff flag is used)
    3. twoBitToFa from the kent tools library (For extracting the reference
            genome in case indexing is required)
    4. lftp for downloading the cosmic vcf

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.

    If you want to use a genome build other than hg19 then download cosmic
    and dbsnp vcfs manually and pass them to this program. This program
    currently only works with hg19 due to how cosmic and NCBI's ownload pages
    work. Other options may be made available in the future.

    The vcfs for the various databases are assumed to be in the parent folders
    of the radia executable (../data/*).  If the data isn't found is the parent
    directories, the program will search VCF_LOCATION before throwing a warning
    and continuing without the said database.  Only dbsnp and cosmic are
    downloaded if not present.
    """
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'radia', 'radia_calls')
    # params ERROR handling and processing
    database_map = process_parameters(params)

    # Move to working directory before doing I/O intensive alignment
    os.chdir(params.working_dir)

    # Call the program
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Starting radia run.', file=params.logfile)
    for chrom in params.chromosome:
        radia_call = [params.radia_executable] #  Base radia call
        radia_call.extend([params.out_prefix, chrom])
        radia_call.extend(['-n', params.norm_d_file])
        radia_call.extend(['-t', params.tum_d_file])
        if params.tum_r_file is not None:
            radia_call.extend(['-r', params.tum_r_file])
        radia_call.append(''.join(['--rnaTumorFasta=', params.rna_fasta]))
        radia_call.extend(['-f', params.genome_fasta])
        radia_call.extend(['-o', ''.join([params.out_prefix, '_', chrom, '.vcf']
                                        )])
        radia_call.extend(['-i', params.genome_version])
        radia_call.extend(['-m', params.genome_fasta])
        radia_call.extend(['-d', params.data_source])
        radia_call.extend(['-q', params.seq_platform])
        radia_call.extend(['--disease', 'params.disease'])
        radia_call.extend(['-l', '\"INFO\"'])
        radia_call.extend(['-g', ''.join([params.out_prefix, '_', chrom, '.log']
                                        )])
        return_value = call(radia_call)
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': radia failed.', params.logfile)
    # Call radia filtering
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Radia completed. Running FilterRadia now.', file=params.logfile)
    for chrom in params.chromosome:
        filter_radia_call = [params.radia_executable] #  Base filter radia call
        filter_radia_call.extend([params.out_prefix, chrom])
        filter_radia_call.append(''.join([params.out_prefix, '_', chrom,
                                          '_filtered.vcf']))
        filter_radia_call.append(params.working_dir)
        filter_radia_call.append(os.path.split(params.radia_executable)[0])
        filter_radia_call.extend(['-b', database_map['blacklist']])
        filter_radia_call.extend(['-d'])
        filter_radia_call.extend(['-r', database_map['retrogenes']])
        filter_radia_call.extend(['-p', database_map['pseudogenes']])
        filter_radia_call.extend(['-c'])
        filter_radia_call.extend(['-t', database_map['broad_targets']])
        if params.use_snpeff:
            filter_radia_call.extend(['-s', params.snpeff_jar])
            filter_radia_call.extend(['-e', params.genome_version])
            if not params.no_canonical:
                filter_radia_call.append(['--canonical'])
        else:
            filter_radia_call.append('--noSnpEff')
        filter_radia_call.extend(['--rnaGeneBlckFile',
                                  database_map['rna_blacklist']])
        filter_radia_call.extend(['--rnaGeneFamilyBlckFile',
                                  database_map['rna_family_blacklist']])
        filter_radia_call.extend(['-f', params.genome_fasta])
        filter_radia_call.extend(['-l', '\"INFO\"'])
        filter_radia_call.extend(['-g', ''.join([params.out_prefix, '_', chrom,
                                                 '_filter.log'])])
        return_value = call(radia_call)
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': FilterRadia failed.', params.logfile)
        call(['rm', ''.join([params.out_prefix, '_', chrom, '.log']),
              ''.join([params.out_prefix, '_', chrom, '.vcf'])])


    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Radia run completed. Finishing up...', file=params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
