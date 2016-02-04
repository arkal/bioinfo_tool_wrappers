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

import sys
import os
import prepare
import pi_errors

def process_parameters(params):
    '''
    This module conducts the error handling for all parmeters passed to the
    program.
    '''
    #  Is mhc_executable set up correctly?
    if os.path.split(params.mhc_executable)[1] == 'predict_binding.py':
        params.mhc_executable = pi_errors.test_param_value(
            params.mhc_executable, 'predict_binding.py', '--mhc_predictor',
            params.logfile)
    else:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': --mhc_predictor has to be predict_binding.py', params.logfile)
    #  List of acceptable prediction methods.  Used to ensure the correct
    #  prediction method has been provided.
    prediction_methods = set(["ann", "comblib_sidney2008", "consensus", 
                              "IEDB_recommended", "netmhcpan", "smm",
                              "smmpmbec"])
    if params.pred_meth not in prediction_methods:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': --prediction_method has to be one of ' + \
            ''.join(prediction_methods), params.logfile)
    #  Test the value of peplen. For MHCI it can only be 9, 10 or 11.
    if params.peplen == "mhc dependent":
        params.peplen = [9, 10]
    elif set([int(x) for x in params.peplen]).difference(set([9, 10, 11])):
        #  This is true iff params.peplen contains an item not in [9, 10, 11].
        #  The values are stored in params as strings and need to be converted
        #  to ints hence the list comprehension.
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': --peplen must be provided space separated values in [9, 10, 11]',
            params.logfile)
    else:
        #  This means the user has passed acceptable values to --peplen
        params.peplen = [int(x) for x in params.peplen]
    #  Once we have the peptide lengths, we need to ensure the --file_prefix has
    #  been set up right
    peplen_filenames = Counter()  #  The 
    if len(params.peplen) == 1 and params.file_prefix.endswith('.faa'):
        #  There is only 1 peplen and the full filename has been provided
        peplen_filenames[params.peplen[0]] = pi_errors.test_param_value(
            '/'.join(params.file_path, params.file_prefix]), 'Input file',
            '--file_prefix', params.logfile)
    elif not params.file_prefix.endswith('.faa')):
        #  1 or more peptide lengths provided and file_prefix is a prefix and
        #  not a path to a .faa file
        for peplen in params.peplen:
            peplen_filenames[params.peplen[0]] = pi_errors.test_param_value(
                ''.join([params.file_path, '/', params.file_prefix, '_',
                         str(params.peplen), '_mer_snpeffed.faa']),
                'Input file', '--file_prefix', params.logfile)
    else:
        #  More that 1 peptide provided and file_prefix points to a single file
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': If more than 1 peplen is provided, --file_prefix must be a ' + \
            'prefix to the individual .faa files and not the path to a ' + \
            'single .faa file itself.', params.logfile)
    #  Test the input alleles for validity
    #  First, if an input .allele file is provided, parse it to obtain the list
    #  of alleles
    if length(params.alleles) == 1 and params.alleles[0].endswith(".alleles"):
        #  This means the user has provided a .alleles file
        pi_errors.test_param_value(params.alleles[0], params.alleles[0],
            '--alleles', params.logfile)
        with open(params.alleles[0], 'r') as allele_file:
            params.alleles=[]
            for line in allele_file:
                params.alleles.append(line.strip())
    #  Once the .allele file has been parsed, params.alleles now either contains
    #  the parsed alleles, or the list of alleles provided by the user.  Now we
    #  need to ensure that all alleles have been formatted in the required
    #  format Eg. HLA-A*01:04
    MHCI = re.compile('HLA-[ABC]\*\d{2}:\d{2}')
    for allele in params.alleles:
        #  For each allele in params.alleles, cehck if it matches the regex for
        #  a generic MHCI molecule
        if not MHCI.match(allele):
            return pi_errors.ParameterError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Alleles for mhci must be in the form HLA-[A/B/C]*XX:XX ' + \
                ' where X is in [0, 9]', params.logfile)
    return peplen_filenames

def main():
    """
    This wrapper script will run the IEDB tools within the  cutadapt docker
    container for the precision immuno project. The wrapper requires
    1. IEDB tools for MHCI prediction - http://tools.iedb.org/mhci
    2. python

    This script requires an input file (--file_prefix) that contains
                        (2 * PEPLEN - 1)-mer
    fasta records for analysis.  If you are looking into multiple peplens,
     --file_prefix should be expandable to:
                FILE_PREFIX_PEPLEN_mer_snpeffed.faa
    for each value of PEPLEN.  If you specify only 1 peplen, if the value of
    FILE_PREFIX ends with .faa then the value is used as is.  If not, the prefix
    needs to expand as above.
    
    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    """
    # Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'mhc', 'mhci_predictions')

    #  Params ERROR handling
    #  peplen_filenames is a dictionary with peptide length as key and the full
    #  path to the filename associated with the peplen as the value.
    peplen_filenames = process_parameters(params)

    # Move to working directory before doing I/O intensive work
    os.chdir(params.working_dir)

    #  set up the strip allele regex
    strip_allele_regex = re.compile('[\*:-]')
    for allele in params.alleles:
        for peptide_length in params.peplen:
            #  Setup the output file
            strip_allele = re.sub(strip_allele_regex, '_', allele)
            #  Setup the call
            mhc_i_call = ['python', params.mhc_executable] #  base call
            mhc_i_call.append(params.pred_meth) #  prediction method
            mhc_i_call.append(allele) #  Allele
            mhc_i_call.append(peptide_length)
            mhc_i_call.append(peplen_filenames[peptide_length])
            mhc_outfile_name = ''.join([params.out_prefix, '_', allele, '.tsv'])
            with open(mhc_outfile_name, 'w') as mhc_outfile:
                return_value = call(mhc_i_call, stdout=mhc_outfile,
                                    stderr=params.logfile)
            if return_value != 0:
                raise pi_errors.MyRuntimeError(
                    dt.now().strftime('%I:%M %p %b %d, %Y') + \
                    ': MHCI prediction failed.', params.logfile)
    
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
