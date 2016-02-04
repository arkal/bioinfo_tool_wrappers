#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_mutect.py

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


def main():
    """
    This wrapper script will run the tool mutect within the  mutect docker
    container for the precision immuno project. The wrapper requires
    1. mutect
    2. java (For running mutect)
    3. twoBitToFa from the kent tools library (For extracting the reference
            genome in case indexing is required)
    4. lftp for downloading the cosmic vcf

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    """
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'mutect', 'mutect_calls')
    # params ERROR handling
    if not (params.java_Xmx.endswith('G') or params.java_Xmx.endswith('M')):
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Please use a suitable value for --Xmx.', params.logfile)
    params.java_executable = pi_errors.test_param_value(params.java_executable,
                                                        'java',
                                                        '--java',
                                                        params.logfile)
    params.mutect_jar = pi_errors.test_param_value(params.mutect_jar,
                                                   'Mutect jar',
                                                   '--mutect_jar',
                                                   params.logfile)
    #  If Indexing is required, does twoBitToFa point to a valid file?
    if params.index_location is None:
        params.tbtf_executable = pi_errors.test_param_value(
            params.tbtf_executable, 'twoBitToFa', '--twoBitToFa',
            params.logfile)
    #  Do the dnsnp and cosmic vcfs exist?
    if params.dbsnp_file == 'DOWNLOAD' or params.cosmic_file == 'DOWNLOAD':
        #  First ensure the vcf storage location has been provided
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
    if params.genome_fasta == 'DOWNLOAD':
        if params.vcf_location is None:
            #  If params.vcf_location is None, set it to the output directory
            params.vcf_location = params.outdir
        #  Does the fasta exist in the vcf_location directory?
        if os.path.exists(''.join([params.vcf_location, '/',
                                   params.genome_version, '.fa'])):
            params.genome_fasta = ''.join([params.vcf_location, '/',
                                           params.genome_version, '.fa'])
        else:
            params.genome_fasta = prepare.get_genome(params.genome_version,
                                                     params.vcf_location,
                                                     params.tbtf_executable,
                                                     params.logfile)
    else:
        params.genome_fasta = pi_errors.test_param_value(params.genome_fasta,
                                                         'Genomic Fasta',
                                                         '--genome_fasta',
                                                         params.logfile)

    # Move to working directory before doing I/O intensive work
    os.chdir(params.working_dir)

    # Call the program
    mutect_call = [params.java_executable, ''.join(['-Xmx', params.java_Xmx]),
                   '-jar'] #  Base java call
    mutect_call.append(params.mutect_jar)
    mutect_call.extend(['-T', 'MuTect'])
    mutect_call.extend(['-R', params.genome_fasta])
    mutect_call.extend(['--cosmic', params.cosmic_file])
    mutect_call.extend(['--dbsnp', params.dbsnp_file])
    mutect_call.extend(['--input_file:normal', params.norm_d_file])
    mutect_call.extend(['--input_file:tumor', params.tum_d_file])
    mutect_call.extend(['--out', ''.join([params.out_prefix, '.out'])])
    return_value = call(mutect_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': MuTect failed.', params.logfile)

    with open(''.join([params.out_prefix, '.out']), 'r') as mutect_file, \
            open(''.join([params.out_prefix, 'non_rejected.out']), 'w') as \
            nr_file:
        for line in mutect_file:
            line = line.strip()
            if line.startswith('#'):
                print(line, file=nr_file)
                continue
            if line.startswith('contig'):
                print('#', line, sep='', file=nr_file)
                continue
            line = line.split('\t')
            if line[50] == 'REJECT':
                continue
            else:
                print(line, sep='\t', file=nr_file)

    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Mutect run completed. Finishing up...', file=params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
