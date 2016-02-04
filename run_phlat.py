#!/usr/bin/env python2.7
'''
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_phlat.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
'''
from __future__ import division, print_function
from subprocess import call
from datetime import datetime as dt

import sys
import os
import prepare
import pi_errors


def main():
    '''
    This wrapper script will run the tool PHLAT within the phlat docker
    container for the precision immuno project. The wrapper requires:
    1. PHLAT.py
    2. bowtie2
    3. gdown.pl (For donwloading the PHLAT Index - available from
       https://raw.githubusercontent.com/Nanolx/patchimage/master/tools/gdown.pl)

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    '''
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'phlat', 'MHC_typing')

    # params ERROR handling
    if not params.phlat_executable.endswith('PHLAT.py'):
        params.phlat_executable = '/'.join([params.phlat_executable,
                                            'PHLAT.py'])
    params.phlat_executable = pi_errors.test_param_value(
        params.phlat_executable, 'PHLAT', '--phlat', params.logfile)
    params.bowtie2_executable = pi_errors.test_param_value(
        params.bowtie2_executable, 'bowtie2', '--bowtie2', params.logfile)
    phlat_dir = os.path.split(os.path.split(params.phlat_executable)[0])[0]
    params.gdownpl_executable = pi_errors.test_param_value(
        params.gdownpl_executable, 'gdown.pl', '--gdownpl', params.logfile)

    if params.index_location is None:
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') +': ' +
              'Downloading Indexes...', file=params.logfile)
        params.index_destination = os.path.abspath(params.index_destination)
        if not os.path.exists(params.index_destination):
            prepare.py_mkdir(params.index_destination)
        getindex_call = [params.gdownpl_executable, 'https://drive.google.com' +
                         '/uc?export=download&confirm=yAjx&id=0Bz-w5tutuZIYY3' +
                         'h5YlMzTjhnbGM', ''.join([params.index_destination,
                                                   '/index4phlat.tar.gz'])]
        print(getindex_call, file=params.logfile)
        return_value = call(getindex_call)
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Could not download indexes. Try manually downloading.',
                params.logfile)
        extract_call = ['tar', '-C', params.index_destination, '-zxvf',
                        '/'.join([params.index_destination,
                                  'index4phlat.tar.gz'])]
        return_value = call(extract_call)
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Index4phlat could not be extracted.', params.logfile)
        else:
            call(['rm', '/'.join([params.index_destination,
                                  'index4phlat.tar.gz'])])
        index_path = '/'.join([params.index_destination, 'index4phlat'])
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Indexes Downloaded.', file=params.logfile)
    else:
        params.index_location = os.path.abspath(params.index_location)
        if not os.path.exists(''.join([params.index_location,
                                       '/ucsc.artHLA.1.bt2'])):
            raise pi_errors.InputFileError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Index file not found.', params.logfile)
        else:
            index_path = params.index_location

    # Move to working directory before doing I/O intensive alignment
    os.chdir(params.working_dir)

    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') +':    ' +
          'Begining MHC Haplotyping', file=params.logfile)
    system_call = ['/usr/bin/env', 'python2.7', '-O', params.phlat_executable]
    system_call.extend(['-1', ''.join([params.file_path, '/',
                                       params.file_prefix, '_1.fastq'])]) # Fq1
    system_call.extend(['-2', ''.join([params.file_path, "/",
                                       params.file_prefix, '_2.fastq'])]) # Fq2
    system_call.extend(['-index', index_path]) # Index files
    system_call.extend(['-b2url', params.bowtie2_executable]) # Bowtie2
    system_call.extend(['-tag', ''.join([params.out_prefix])]) # DNA/RNA
    system_call.extend(['-e', phlat_dir]) # Phlat directory home
    system_call.extend(['-o', params.outdir]) # Output directory
    system_call.extend(['-p', str(params.n)]) # Number of threads

    # Call the program
    return_value = call(system_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': MHC Haplotyping failed.', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Alignment completed. Finishing up...', file=params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()

if __name__ == '__main__':
    sys.exit(main())
