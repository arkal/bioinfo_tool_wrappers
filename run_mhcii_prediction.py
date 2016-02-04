#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_MHC_prediction.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import division, print_function
from subprocess import call
from datetime import datetime as dt
from collections import Counter

import re
import sys
import os
import prepare
import pi_errors

def process_parameters(params):
    '''
    This module conducts the error handling for all parmeters passed to the
    program.
    '''
    #  Are predictor executables set up correctly?
    if os.path.split(params.mhc_executable)[1] == 'mhc_II_binding.py':
        params.mhc_executable = pi_errors.test_param_value(
            params.mhc_executable, 'mhc_II_binding.py', '--mhc_predictor',
            params.logfile)
        params.netmhciipan_executable = pi_errors.test_param_value(
            params.netmhciipan_executable, 'netMHCIIpan', '--netMHCIIpan',
            params.logfile)
    else:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': --mhc_predictor has to be predict_binding.py', params.logfile)
    #  List of acceptable prediction methods.  Used to ensure the correct
    #  prediction method has been provided.
    prediction_methods = set(['comblib', 'consensus3', 'IEDB_recommended',
                              'NetMHCIIpan', 'nn_align', 'smm_align',
                              'tepitope', ])
    if params.pred_meth not in prediction_methods:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': --prediction_method has to be one of ' + \
            ''.join(prediction_methods), params.logfile)
    #  For MHCII, peplen can only be 15.
    params.peplen = [15]
    #  Ensure the --file_prefix has been set up right
    if params.file_prefix.endswith('.faa'):
        pepfilename = pi_errors.test_param_value(
            '/'.join([params.file_path, params.file_prefix]), 'Input file',
            '--file_prefix', params.logfile)
    else:
        #  More that 1 peptide provided and file_prefix points to a single file
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': --file_prefix must point to a valid .faa file.', params.logfile)
    #  Test the input alleles for validity
    #  First, if an input .allele file is provided, parse it to obtain the list
    #  of alleles
    if len(params.alleles) == 1 and params.alleles[0].endswith(".alleles"):
        #  This means the user has provided a .alleles file
        pi_errors.test_param_value(params.alleles[0], params.alleles[0],
                                   '--alleles', params.logfile)
        with open(params.alleles[0], 'r') as allele_file:
            params.alleles = []
            for line in allele_file:
                params.alleles.append(line.strip())
    #  Once the .allele file has been parsed, params.alleles now either contains
    #  the parsed alleles, or the list of alleles provided by the user.  Now we
    #  need to ensure that all alleles have been formatted in the required
    #  format Eg. HLA-DRB1*01:04, HLA-DQA1*03:01/DQB1*03:02
    #  IMPORTANT: HLA-DP isn't implemented yet -- none of the predictors use it
    MHCII = re.compile(
        r'HLA-D((RB\d\*\d{2}:\d{2})|'
        r'([PQ]A\d\*\d{2}:\d{2}/D[PQ]B\d\*\d{2}:\d{2}))')

    for allele in params.alleles:
        #  For each allele in params.alleles, cehck if it matches the regex for
        #  a generic MHCI molecule
        if not MHCII.match(allele):
            return pi_errors.ParameterError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Alleles for mhcii must be in the form HLA-DRB1*XX:XX or' + \
                ' HLA-DQA1*XX:XX/DQB1*XX:XX where X is in [0, 9]',
                params.logfile)
    return pepfilename

def main():
    """
    This wrapper script will run the IEDB tools within the  cutadapt docker
    container for the precision immuno project. The wrapper requires
    1. IEDB tools for MHCI prediction - http://tools.iedb.org/mhci
    2. netMHCIIpan - In case IEDB tools fails
    3. python

    This script requires an input file (--file_prefix) that contains
                        (2 * PEPLEN - 1)-mer
    fasta records for analysis. FILE_PREFIX must be a .faa file.

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    """
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'mhc', 'mhci_predictions')

    #  Params ERROR handling
    #  peplen_filenames is a dictionary with peptide length as key and the full
    #  path to the filename associated with the peplen as the value.
    pepilename = process_parameters(params)

    # Move to working directory before doing I/O intensive work
    os.chdir(params.working_dir)

    #  set up the different allele regexes
    strip_allele_regex = re.compile(r'[\*:/-]') # For strip allele
    dpqa_allele_regex_1 = re.compile(r'[\*:]')  # For DPA and DQA if netMHCIIpan 
    dpqa_allele_regex_2 = re.compile(r'/')      # is used
    for allele in params.alleles:
        for peptide_length in params.peplen:
            #  Setup the output file
            #  Strip allele converts HLA-DRB1*15:01 to HLA_DRB1_15_01 and
            #  HLA-DQA1*01:02/DQB1*03:02 to HLA_DQA1_01_02_DQB1_03_02
            strip_allele = re.sub(strip_allele_regex, '_', allele)
            #  Setup the call
            mhc_ii_call = ['python', params.mhc_executable] #  base call
            mhc_ii_call.append(params.pred_meth) #  prediction method
            mhc_ii_call.append(allele) #  Allele
            mhc_ii_call.append(pepfilename)
            mhc_outfile_name = ''.join([params.out_prefix, '_', allele, '.tsv'])
            with open(mhc_outfile_name, 'w') as mhc_outfile:
                return_value = call(mhc_ii_call, stdout=mhc_outfile,
                                    stderr=params.logfile)
            if return_value != 0:
                print('WARNING: IEDBtools failed.  Attempting netMHCIIpan',
                      file=params.logfile)
                #  netmHCIIpan needs a different formatting for allele
                #  HLA-DQA1*01:02/DQB1*03:02 should be HLA-DQA10102-DQB10302
                #  HLA-DRB1*15:01 should be DRB1_1501. DP and DQ are similar
                if allele.startswith('HLA-DQ') or allele.startswith('HLA-DP'):
                    allele = re.sub(dpqa_allele_regex_1, '', allele)
                    allele = re.sub(dpqa_allele_regex_2, '-', allele)
                else:
                    allele = strip_allele[4:] # Easier than starting from allele
                netmhc_ii_call = [params.netmhciipan_executable]
                netmhc_ii_call.extend(['-a', allele]) #  Allele
                netmhc_ii_call.extend(['-xls', 1])
                netmhc_ii_call.extend(['-xlsfile', mhc_outfile_name])
                netmhc_ii_call.extend(['-f', pepfilename])
                return_value = call(mhc_ii_call, stderr=params.logfile)
                if return_value != 0:
                    raise pi_errors.MyRuntimeError(
                        dt.now().strftime('%I:%M %p %b %d, %Y') + \
                        ': MHCII prediction failed.', params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
