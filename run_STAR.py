#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_STAR.py

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
    # Do the STAR executables point to valid files?
    params.STAR_executable = pi_errors.test_param_value(
        params.STAR_executable, 'STAR', '--STAR', params.logfile)
    params.STAR_executable = pi_errors.test_param_value(
        params.STARlong_executable, 'STARlong', '--STARlong', params.logfile)
    #  If Indexing is required, does twoBitToFa point to a valid file?
    if params.index_location is None:
        params.tbtf_executable = pi_errors.test_param_value(
            params.tbtf_executable, 'twoBitToFa', '--twoBitToFa',
            params.logfile)
    #  Obtain read length from the second non-empty line of the file (seq line).
    with open(''.join([params.file_path, '/', params.file_prefix,
                       '_1.fastq']), 'r') as fastq_1_file:
        line_num = 0
        for line in fastq_1_file:
            line = line.strip()
            if len(line) == 0:
                continue
            line_num += 1
            if line_num == 2:
                read_length = len(line)
                break
    #  Set up the executable based on the read length
    if int(read_length) > 300:
        star_executable = params.STARlong_executable
    else:
        star_executable = params.STAR_executable
    # Check for indexes. If the user has specified that indexes need to
    # be created then do so.
    if params.index_location is None:
        index_path = star_indexing(star_executable, read_length, params)        
    else:
        params.index_location = os.path.abspath(params.index_location)
        if not os.path.exists(''.join([params.index_location, '/SA'])):
            raise pi_errors.InputFileError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Index file not found', params.logfile)
        else:
            index_path = params.index_location
    return star_executable, index_path

def star_indexing(star_executable, read_length, params):
    '''
    This module indexes a genome using STAR_EXECUTABLE using READ_LENGTH to set
    edge size.
    params contains
    index_destination - The location where the index should be stored
    logfile - Open file handle to a log file
    genome_version - hg19/hg38
    n - number of cores to use
    tbtf_executable - path to twoBitToFa
    RETURN VALUES
    index_path - path ot the directory where indexes were stored
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Indexing fasta...', file=params.logfile)
    params.index_destination = os.path.abspath(params.index_destination)
    if not os.path.exists(params.index_destination):
        prepare.py_mkdir(params.index_destination)
    edge_size = max(50, int(round(read_length / 50, 0) * 50)) # minimum edge
                                                              # size = 50
    index_path = ''.join([params.index_destination, '/STAR_',
                          str(edge_size), '_references'])
    if not os.path.exists(index_path): # make reference based on edge size
        prepare.py_mkdir(index_path)
    genome_fasta = prepare.get_genome(params.genome_version, index_path,
                                      params.tbtf_executable,
                                      params.logfile)
    gencode_file = prepare.get_gtf(params.genome_version, index_path,
                                   params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ':    ' +
          'Running STAR index on fasta reference.', file=params.logfile)
    starindex_call = [star_executable] # Base call
    starindex_call.extend(['--runThreadN', str(params.n)]) # Threads
    starindex_call.extend(['--runMode', 'genomeGenerate']) # Indexing module
    starindex_call.extend(['--genomeDir', index_path]) # index directory
    starindex_call.extend(['--genomeFastaFiles', genome_fasta]) # Genomic fa
    starindex_call.extend(['--sjdbGTFfile', gencode_file]) # gencode annots
    starindex_call.extend(['--sjdbOverhang', str(read_length)]) # edge size
    print(starindex_call, file=params.logfile)
    return_value = call(starindex_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Indexing Failed', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Indexing completed.', file=params.logfile)
    return index_path

def star_alignment(star_executable, index_path, params):
    '''
    This module will align the reads to the indexes stored in INDEX_PATH using
    STAR_EXECUTABLE.
    params contains
    n - Number of cores to use
    out_prefix - Prefix to use for output files
    file_prefix - Prefix to input fastq files
    logfile - Open file handle to a log file
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Aligning' +
          ' reads using STAR...', file=params.logfile)
    staralign_call = [star_executable]  # Base
    staralign_call.extend(['--runThreadN', str(params.n)]) # Threads
    staralign_call.extend(['--genomeDir', index_path]) # index directory
    staralign_call.extend(['--outFileNamePrefix', params.out_prefix])
    staralign_call.extend(['--readFilesIn',
                           ''.join([params.file_path, '/', params.file_prefix,
                                    '_1.fastq']),
                           ''.join([params.file_path, '/', params.file_prefix,
                                    '_2.fastq'])])
    staralign_call.extend(['--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD'])
    staralign_call.extend(['--outSAMtype', 'BAM', 'SortedByCoordinate'])
    staralign_call.extend(['--quantMode', 'TranscriptomeSAM'])
    staralign_call.extend(['--outSAMunmapped', 'Within'])
    return_value = call(staralign_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Alignment failed', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Alignment completed. Finishing up...', file=params.logfile)


def main():
    '''
    This wrapper script will run the entire alignment pipeline for RNA from
    alignment of fastqs (to both the genome, and the transcriptome), to sorting
    and indexing. The wrapper can even download and produce STAR references if
    required. The wrapper requires
    1. STAR (For aligning reads)
    2. twoBitToFa from the kent tools library (For extracting the reference
                genome in case indexing is required)

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    '''
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'star', 'STAR_alignment')

    #  Params ERROR handling
    #  Do the STAT executables point to valid files?
    star_executable, index_path = process_parameters(params)

    # Move to working directory before doing I/O intensive alignment
    os.chdir(params.working_dir)

    # Align reads to sam file
    star_alignment(star_executable, index_path, params)

    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
