#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_rsem.py

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
    #  Does the rsem path point to a valid file?
    params.rsem_path = pi_errors.test_param_value(params.rsem_path,
                                                  'rsem binaries',
                                                  '--rsem_path',
                                                  params.logfile)
    #  Does the file prefix point to a bam or a fastq file?
    if not (os.path.exists(''.join([params.file_path, '/', params.file_prefix,
                                    '.bam'])) or
            os.path.exists(''.join([params.file_path, '/', params.file_prefix,
                                    '_1.fastq']))):
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' + \
            'Please check input files. Neither bam nor fastq found at:\n' + \
            params.file_path, params.logfile)
    #  If the input is not a bam file, bowtie will be required
    if not os.path.exists(''.join([params.file_path, '/', params.file_prefix,
                                   '.bam'])):
        fasta_input = True
        bowtie_status = [1, 1]  #  Assume both exist [bt1, bt2]
        #  Try bowtie first then bowtie 2.  If both exist, bowtie_path and
        #  bowtie_version will take the values for bowtie2.  If only one exists,
        #  then they will take the values from the one that exists.  Since
        #  test_param_value closes the logfile, the except cases in both try
        #  instances reopen params.logfile.
        try:
            params.bowtie_executable = pi_errors.test_param_value(
                params.bowtie_executable, 'bowtie', '--bowtie',
                params.logfile)
        except pi_errors.ParameterError:
            bowtie_status[0] = 0  #  bt1 doesn't exist
            if params.logfile.name != '<stderr>':
                params.logfile = open(params.logfile.name, params.logfile.mode,
                                      0)
        else:
            bowtie_path, bowtie_version = \
                os.path.split(params.bowtie_executable)
        try:
            params.bowtie2_executable = pi_errors.test_param_value(
                params.bowtie2_executable, 'bowtie2', '--bowtie2',
                params.logfile)
        except pi_errors.ParameterError:
            bowtie_status[1] = 0  #  bt2 doesn't exist
            if params.logfile.name != '<stderr>':
                params.logfile = open(params.logfile.name, params.logfile.mode,
                                      0)
        else:
            bowtie_path, bowtie_version = \
                os.path.split(params.bowtie2_executable)
        if bowtie_status == [0, 0]:  # This means that both aren't present
            raise pi_errors.ParameterError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': bowtie binaries not found on $PATH. Please specify' + \
                ' one explicitly with --bowtie or --bowtie2.', params.logfile)
    else:
        fasta_input = False
        bowtie_path, bowtie_version = None, None # Set these to None

    # Check for indexes. If the user has specified that indexes need to
    # be created then do so.
    if params.index_location is None:
        #  If you need indexes, you need twoBitToFa
        params.tbtf_executable = pi_errors.test_param_value(
            params.tbtf_executable, 'twoBitToFa', '--twoBitToFa',
            params.logfile)
        index_path = rsem_index(''.join([params.rsem_path,
                                         '/rsem-prepare-reference']),
                                         fasta_input,
                                         (bowtie_path, bowtie_version),
                                         params)
    else:
        index_path = os.path.abspath(params.index_location)

    try:
        #  One of the files crucial to the index is GENOME_FASTA_PREFIX.chrlist
        index_prefix = [x.split('.')[0] for x in os.listdir(index_path) if \
                        x.endswith('chrlist')][0]
    except IndexError:
        if params.logfile.name != '<stderr>':
            params.logfile.close()
        raise SystemExit('ERROR ' + dt.now().strftime('%I:%M %p %b %d, %Y') +
                         ': Indexes not found at specified location')
    return bowtie_path, bowtie_version, index_path, index_prefix


def rsem_index(rsem_index_executable, fasta_input, bowtie_info, params):
    '''
    This module will create the rsem indexes at params.index_destination using
    RSEM_INDEX_EXECUTABLE. If FASTA_INPUT = True, it will use the bowtie version
    to make bowtie indexes as well.
    bowtie_info is a tuple of (bowtie_path, bowtie_version)
    params contains
    index_destination - Folder to store the indexes
    n - number of cores to use
    genome_fasta - path to genomic fasta file. Can also specify DOWNLOAD.
    genome_version - hg19/hg38
    logfile - Open file handle to a log file
    RETURN VALUES
    index_path - Path to directory where nidexes were stored
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Creating rsem references...', file=params.logfile)
    index_path = os.path.abspath(params.index_destination)
    #  If the directory doesn't exist, create it
    if not os.path.exists(index_path):
        prepare.py_mkdir(index_path)
    if params.genome_fasta == 'DOWNLOAD':
        params.genome_fasta = prepare.get_genome(params.genome_version,
                                                 index_path,
                                                 params.tbtf_executable,
                                                 params.logfile)
    else:
        params.genome_fasta = pi_errors.test_param_value(params.genome_fasta,
                                                         'Genomic Fasta',
                                                         '--genome_fasta',
                                                         params.logfile)
    #  If the gtf file is required, download it
    gencode_file = prepare.get_gtf(params.genome_version, index_path,
                                   params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ':    ' +
          'Running rsem-prepare-reference on fasta reference.',
          file=params.logfile)
    rsem_prepref_call = [rsem_index_executable] # base call
    rsem_prepref_call.extend(['--gtf', gencode_file]) # gtf file
    if fasta_input:
        rsem_prepref_call.extend([''.join(['--', bowtie_version]),
                                  ''.join(['--', bowtie_version, '-path']),
                                  bowtie_path])
    else:
        rsem_prepref_call.append('--no-bowtie')
    rsem_prepref_call.append(params.genome_fasta)
    rsem_prepref_call.extend([''.join([index_path, '/',
                                       params.genome_version])])
    print(rsem_prepref_call, file=params.logfile)
    return_value = call(rsem_prepref_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Indexing Failed', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Indexing completed.', file=params.logfile)
    return index_path


def rsem_calculate_expression(rsem_calexp_executable, bowtie_info,
                              index_info, params):
    '''
    This module will process the bam of fastq files pointed to by 
    PARAMS.FILE_PREFIX using RSEM_CALEXP_EXECUTABLE.
    bowtie_info is a tuple of (bowtie_path, bowtie_version)
    index_info is a tuple of (index_path, index_prefix)
    params contains:
    logfile - Open file handle to a log file
    n - number of cores to use
    file_path, file_prefix - path and prefix to input bam/fastq files
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Calculating gene expression using rsem...', file=params.logfile)
    index_path, index_prefix = index_info
    bowtie_path, bowtie_version  = bowtie_info
    rsem_calcexp_call = [rsem_calexp_executable] # base call
    rsem_calcexp_call.extend(['--paired-end'])
    rsem_calcexp_call.extend(['-p', str(params.n)])
    if not fasta_input:
        rsem_calcexp_call.extend(['--bam'])
        rsem_calcexp_call.extend([''.join([params.file_path, '/',
                                           params.file_prefix, '.bam'])])
        rsem_calcexp_call.extend(['--no-bam-output'])
    else:
        rsem_calcexp_call.extend(['--output-genome-bam'])
        rsem_calcexp_call.extend([''.join(['--', bowtie_version]),
                                  ''.join(['--', bowtie_version, '-path']),
                                  bowtie_path])
        rsem_calcexp_call.extend([''.join([params.file_path, '/',
                                           params.file_prefix, '_1.fastq']),
                                  ''.join([params.file_path, '/',
                                           params.file_prefix, '_2.fastq'])])
    rsem_calcexp_call.extend(['/'.join([index_path, index_prefix])])
    rsem_calcexp_call.extend([params.file_prefix])
    print(rsem_calcexp_call, file=params.logfile)
    return_value = call(rsem_calcexp_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError('ERROR ' + \
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': RSEM failed', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Gene ' +
          'expression calculated. Finishing up...', file=params.logfile)


def main():
    '''
    This wrapper script will run the gene expression analysis for RNA-seq data.
    The input material can be ither a bam file containing RNA-seq reads mapped
    to the genome (Using STAR or bowtie), or the RNA-Seq fastq files themselves.
    The wrapper can produce the rsem references if required. The wrapper
    requires the following
    1. rsem
    2. bowtie/bowtie2 - Only if input is fastq
    3. twoBitToFa from the kent tools library (For extracting the reference
                genome in case indexing is required)

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    '''
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'rsem', 'RSEM_quant')

    # params ERROR handling
    bowtie_path, bowtie_version, index_path, index_prefix = \
                                                      process_parameters(params)
    # Move to working directory before doing I/O intensive alignment
    os.chdir(params.working_dir)

    # Process with rsem
    rsem_calculate_expression(''.join([params.rsem_path,
                                       '/rsem-calculate-expression']),
                              (bowtie_path, bowtie_version),
                              (index_path, index_prefix), params)

    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
