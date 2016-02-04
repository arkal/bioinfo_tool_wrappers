#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_cutadapt.py

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
    This wrapper script will run the tool cutadapt within the  cutadapt docker
    container for the precision immuno project. The wrapper requires
    1. cutadapt
    2. GNU sed (Tested on version 4.2.1)

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    """
    #  Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'cutadapt', 'adapter_fixed')

    # params ERROR handling
    params.cutadapt_executable = pi_errors.test_param_value(
        params.cutadapt_executable, 'cutadapt', '--cutadapt', params.logfile)
    if not(set(params.fwd_3pr_adapter).issubset(set("ACTGN")) and \
           set(params.rev_3pr_adapter).issubset(set("ACTGN"))):
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Adapter sequences can only contain A, C, T, G, and N.',
            params.logfile)

    #  Move to working directory before doing I/O intensive work
    os.chdir(params.working_dir)

    #  Remvove adapter contamination
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + \
          ': Trimming adapters using cutadapt.', file=params.logfile)
    cutadapt_call = [params.cutadapt_executable] # base call
    cutadapt_call.extend(['-a', params.fwd_3pr_adapter])  # Fwd read 3' adapter
    cutadapt_call.extend(['-A', params.rev_3pr_adapter])  # Rev read 3' adapter
    cutadapt_call.extend(['-m', '35'])  # Minimum size of read
    cutadapt_call.extend(['-o', ''.join([params.file_prefix,
                                         '_cutadapt_1.fastq'])])
    cutadapt_call.extend(['-p', ''.join([params.file_prefix,
                                         '_cutadapt_2.fastq'])])
    cutadapt_call.append(''.join([params.file_path, '/', params.file_prefix,
                                  '_1.fastq']))
    cutadapt_call.append(''.join([params.file_path, '/', params.file_prefix,
                                  '_2.fastq']))
    print(' '.join(cutadapt_call), file=params.logfile)
    return_value = call(cutadapt_call)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': cutadapt failed', params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
