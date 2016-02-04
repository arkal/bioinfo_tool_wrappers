#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_snpeff.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import division, print_function
from subprocess import call
from datetime import datetime as dt

import sys
import os
import prepare
import pi_errors

def process_parameters(params):
    '''
    This module conducts the error handling for all parmeters passed to the
    program.
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') +
          ': Processing input parameters.', file=params.logfile)
    # Does the input vcf file exist?
    if not os.path.exists(''.join([params.file_path, '/', params.file_prefix,
                                   '.vcf'])):
        raise pi_errors.InputFileError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Please provide a valid input file using --file_prefix',
            params.logfile)
    #  The memory option for java should be of the form Xmx10G or Xmx10M
    if not (params.java_Xmx.endswith('G') or params.java_Xmx.endswith('M')):
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Please use a suitable value for --Xmx.', params.logfile)
    params.java_executable = pi_errors.test_param_value(params.java_executable,
                                                        'java',
                                                        '--java',
                                                        params.logfile)
    #  Does the provided snpeff binary provided exist?
    params.snpeff_jar = pi_errors.test_param_value(params.snpeff_jar,
                                                   'snpeff',
                                                   '--snpeff_jar',
                                                   params.logfile)
    params.use_snpeff_db = False
    #  Does the user want a snpEff packaged database?
    if params.config_file == 'PACKAGED':
        params.use_snpeff_db = True
        #  Has the snpeff reference to be used been provided?
        if params.reference_name == None:
            raise pi_errors.ParameterError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': --snp_reference is required if --config=PACKAGED.',
                params.logfile)
    # If a custom databse is desired, does it need to be created?
    if params.index_location is None:
        #  If the user has provided the location to the parent directory of data
        #  directory, make DATA_DIRECTORY point to data. If they have provided
        #  the link to data, make INDEX_DESTINATION point to the parent and
        #  DATA_DIRECTORY point to data.
        if os.path.split(params.index_destination.rstrip('/'))[1] != 'data':
            params.data_directory = '/'.join([params.index_destination, 'data'])
        else:
            params.data_directory = params.index_destination
            params.index_destination = \
                params.index_destination.rstrip('/').rstrip('/data')
        #  Create the data directory if needed
        if not os.path.exists(params.data_directory):
            prepare.py_mkdir(params.data_directory)
        #  If we're using a custom databse, thre is nothing more to do
        if params.use_snpeff_db:
            return None
        #  Initialise the reference name
        params.reference_name = ''.join([params.genome_version, '_custom'])
        # make a variable to gold GENOME_VERSION_custom
        genome_folder = '/'.join([params.data_directory, params.reference_name])
        prepare.py_mkdir(genome_folder)
        #  If the genome fasta isn't provided or is provided a wrong value,
        #  download it
        if params.genome_fasta == 'DOWNLOAD' or not \
                os.path.exists(params.genome_fasta):
            #  Does the provided tbtf binary point to a valid file?
            params.tbtf_executable = pi_errors.test_param_value(
                params.tbtf_executable, 'twoBitToFa', '--twoBitToFa',
                params.logfile)
            params.genome_fasta = prepare.get_genome(
                params.genome_version, genome_folder,
                params.tbtf_executable, params.logfile)
            #  Rename genome fasta
            call(['mv', params.genome_fasta, '/'.join([genome_folder,
                                                       'sequences.fa'])])
        else:
            params.genome_fasta = os.path.abspath(params.genome_fasta)
            #  Link sequencesfa to genome fasta
            call(['ln', '-s', '-T', params.genome_fasta,
                  '/'.join([genome_folder, 'sequences.fa'])])
        #  Download the gencode GTF file
        params.gtf_file = prepare.get_gtf(params.genome_version,
                                          genome_folder,
                                          params.logfile)
        #  Rename gtf file
        call(['mv', params.gtf_file, '/'.join([genome_folder, 'genes.gtf'])])
    # If it has been provided, set up the config file
    else:
        #  If the user has provided the location to the parent directory of data
        #  directory, make DATA_DIRECTORY point to data. If they have provided
        #  the link to data, make INDEX_LOCATION point to the parent and
        #  DATA_DIRECTORY point to data.
        if os.path.split(params.index_location.rstrip('/'))[1] != 'data':
            params.data_directory = '/'.join([params.index_location, 'data'])
        else:
            params.data_directory = params.index_location
            params.index_location = \
                params.index_location.rstrip('/').rstrip('/data')
        #  If we're using a custom databse, thre is nothing more to do
        if params.use_snpeff_db:
            return None
        #  If the config file hasn't been provided, is it in INDEX_LOCATION
        #  AND does GENOME_VERSION_custom exist (i.e. was it created by
        #  this script?)
        if params.config_file is None:
            params.config_file = pi_errors.test_param_value(
                '/'.join([params.index_location, 'snpEff.config']),
                'snpEff.config', '--config', params.logfile)
            #  Dummy variable to ensure the GENOME_VERSION_custom exists
            _ = pi_errors.test_param_value(
                ''.join([params.data_directory, '/', params.genome_version,
                         '_custom']), '_'.join([params.genome_version, 'custom'
                                               ]), '--snpeff_reference and' +
                '--config', params.logfile)
            params.reference_name = '_'.join([params.genome_version, 'custom'])
        #  If a config file has been provided, does it point to a legit file and
        #  has the reference name also been provided?
        else:
            params.config_file = pi_errors.test_param_value(
                params.config_file, 'snpEff config file', '--config',
                params.logfile)
            if params.reference_name is None:
                raise pi_errors.ParameterError(
                    dt.now().strftime('%I:%M %p %b %d, %Y') + \
                    ': --snpeff_reference is required if --config points to' + \
                    ' a custom file.', params.logfile)
    return None



def create_custom_config_file(params):
    '''
    This module will create a custom snpEff.config file pointing to an existing
    or potential target database location. The data used to create the config
    file include INDEX_DESTINATION, GENOME_VERSION, GENOME_FASTA, GTF_FILE,
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') +
          ': creating custom config file', file=params.logfile)
    #  Make the config file
    params.config_file = '/'.join([params.index_destination, 'snpEff.config'])
    with open(params.config_file, 'w') as my_conf_file, \
            open('/'.join([os.path.split(params.snpeff_jar)[0], 'snpEff.config']
                         ), 'r') as pkg_conf_file:
        ###### index_dest or data_dir?
        print('datadir = ', params.index_destination, sep='', file=my_conf_file)
        for line in pkg_conf_file:
            if line.startswith('#'):
                continue
            if line.startswith('database_repository') or \
                    line.startswith('versions_url') or \
                    line.startswith('lof.ignoreProteinCodingAfter') or \
                    line.startswith('lof.ignoreProteinCodingBefore') or \
                    line.startswith('lof.deleteProteinCodingBases') or \
                    line.startswith('codon.Standard'):
                print(line, file=my_conf_file)
            if line.startswith('codon.Vertebrate_Mitochondrial'):
                print(line, file=my_conf_file)
                break
        print(params.genome_version, '_custom.genome : Human', sep='',
              file=my_conf_file)
        print(params.genome_version, '_custom.M.codonTable : Vertebrate_'
              'Mitochondrial', sep='', file=my_conf_file)
    return None


def build_snpeff_database(params):
    '''
    This module will build a custom snpeff database for the provided
    config file. If creating froma  custom config file, the database will be
    created as

    INDEX_DESTINATION
        |-snpEff.config
        |-data
            |-GENOME_VERISON_custom
                |-genes.gtf
                |-sequences.fa
                |-snpEffectPredictor.bin
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') +
          'Building snpEff databse', file=params.logfile)
    #  Build the database
    build_db_call = [params.java_executable, ''.join(['-Xmx', params.java_Xmx]),
                     '-jar'] #  Base java call
    build_db_call.append(params.snpeff_jar)  # snpEff
    build_db_call.append('build')
    if params.use_snpeff_db:
        ###### index_dest or data_dir?
        build_db_call.extend(['-dataDir', params.index_destination])
    else:
        build_db_call.extend(['-c', params.config_file])  # custom-specific
        build_db_call.append('-gtf22')
    build_db_call.append(params.reference_name)
    print(build_db_call, file=params.logfile)
    return_value = call(build_db_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': snpeff build failed.', params.logfile)
    return None

def run_snpeff(params):
    '''
    This module takes the processed parameters and runs the snpeff tool
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Starting' +
          ' snpEff run.', file=params.logfile)
    snpeff_call = [params.java_executable, ''.join(['-Xmx', params.java_Xmx]),
                   '-jar']
    snpeff_call.append(params.snpeff_jar)  # snpEff
    snpeff_call.append('eff')
    #  If we are using snpeff database, we need to override the datadir in the
    #  packaged config file
    if params.use_snpeff_db:
        #  If location has been provided, use it
        if params.index_destination is None:
            snpeff_call.extend(['-dataDir', params.index_location])
        #  Else use the -download flag to downlaod it
        else:
            snpeff_call.append('-download')
            snpeff_call.extend(['-dataDir', params.index_destination])
    #  If we are using a custom database, we need to provide the config file
    else:
        snpeff_call.extend(['-c', params.config_file])
    if not params.use_intergenic:
        snpeff_call.append('-no-intergenic')
    if not params.use_downstream:
        snpeff_call.append('-no-downstream')
    if not params.use_upstream:
        snpeff_call.append('-no-upstream')
    if not params.no_canonical:
        snpeff_call.append('-canon')
    snpeff_call.append('-noStats')
    snpeff_call.append('-snp')
    snpeff_call.append(params.reference_name)
    snpeff_call.append(''.join([params.file_path, '/', params.file_prefix,
                                '.vcf']))
    print(snpeff_call, file=params.logfile)
    with open(''.join([params.out_prefix, '_snpeffed.vcf']), 'w') as \
            snpeff_file:
        return_value = call(snpeff_call, stdout=snpeff_file,
                            stderr=params.logfile)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': snpEff failed.', params.logfile)
    return None

def main():
    """
    This wrapper script will run the tool snpeff for the precision immuno
    project. The wrapper requires
    1. snpeff

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.

    A standard snpeff database looks like
    <database folder>
        |-data
           |
           |-<reference name>
                 |-genes.gtf
                 |-sequences.fa
                 |-snpEffectPredictor.bin

    When using --index_location, provide <database folder> as mentioned above.
    DO NOT ADD THE THE data folder at the end of the file path.  If the database
    was made using this script, the snpEff.config file will be in
    <database folder> and will be used to run snpeff using --snpeff_reference
     GENOME_VERISON_custom.  If the database was created manually, and
    snpEff.config is not in <database folder>, it must be explicitly provided
    with --config and --snpeff_reference flags.

    NOTE: If you want to use a snpEff packaged database, use --config=PACKAGED
    and --snpeff_reference=<snpeff reference name>
    If used along with --require_index (instead of --index_location), the
    script download the database
    """
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'snpeff', 'snpeffed')

    # params ERROR handling
    process_parameters(params)

    # Create databses if required
    if params.index_location is None:
        #  If we are not using a snpeff packaged database, we need to create a
        #  config file
        if not params.use_snpeff_db:
            create_custom_config_file(params)
            build_snpeff_database(params)
    # Move to working directory before doing I/O intensive work
    os.chdir(params.working_dir)
    # Run Snpeff
    run_snpeff(params)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Snpeff ' +
          'run finished. Finishing up...', file=params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()

if __name__ == '__main__':
    sys.exit(main())
