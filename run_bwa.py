#!/usr/bin/env python2.7
"""
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/run_bwa.py

Program info can be found in the docstring of the main function.
Details can also be obtained by running the script with -h .
"""
from __future__ import division, print_function
from subprocess import call, check_output
from datetime import datetime as dt

import sys
import os
import prepare
import pi_errors


def main():
    '''
    This wrapper script will run the entire alignment pipeline for genomic DNA
    (WGS or WXS) from alignment of fastqs, to sorting, indexing, and Read Group
    incorporation. The wrapper can even download and produce bwa references if
    required. The wrapper requires
    1. bwa (For aligning reads)
    2. java (For picard)
    3. picard tools (For read groups)
    4. samtools (For sam/bam manipulation)
    5. twoBitToFa from the kent tools library (For extracting the reference
            genome in case indexing is required)

    Unless specified, the program will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.
    '''
    #  Parse the arguments using prepare.parse_args()
    params = prepare.parse_args(main.__doc__, 'bwa', 'bwa_alignment')

    #  Params ERROR handling
    #  The memory option for java should be of the form Xmx10G or Xmx10M
    if not (params.java_Xmx.endswith('G') or params.java_Xmx.endswith('M')):
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Please use a suitable value for --Xmx.', params.logfile)
    params.bwa_executable = pi_errors.test_param_value(params.bwa_executable,
                                                       'bwa',
                                                       '--bwa',
                                                       params.logfile)
    params.samtools_executable = pi_errors.test_param_value(
        params.samtools_executable, 'samtools', '--samtools', params.logfile)
    params.java_executable = pi_errors.test_param_value(params.java_executable,
                                                        'java',
                                                        '--java',
                                                        params.logfile)
    #  If Indexing is required, does twoBitToFa point to a valid file?
    if params.index_location is None:
        params.tbtf_executable = pi_errors.test_param_value(
            params.tbtf_executable, 'twoBitToFa', '--twoBitToFa',
            params.logfile)
    if not params.picard_jar.endswith('jar'):
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Please specify a valid jar file for picard!', params.logfile)
    else:
        params.picard_jar = pi_errors.test_param_value(params.picard_jar,
                                                       'picard',
                                                       '--picard_jar',
                                                       params.logfile)

    if params.RGID is None:
        params.RGID = params.file_prefix

    #read_group = ''.join(['\'@RG\\tID:', params.RGID, '\\tPL:ILLUMINA\\tSM:',
    #                      params.sample_type, '\''])
    # Check for indexes. If the user has specified that indexes need to
    # be created then do so.
    if params.index_location is None:
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Indexing fasta...', file=params.logfile)
        if not os.path.exists(params.index_destination):
            prepare.py_mkdir(params.index_destination)
        index_path = params.index_destination
        genome_fasta = prepare.get_genome(params.genome_version, index_path,
                                          params.twoBitToFa_executable,
                                          params.logfile)
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ':    ' +
              'Running BWA index on fasta reference.', file=params.logfile)
        return_value = call([params.bwa_executable, 'index', genome_fasta])
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': bwa index failed.', params.logfile)
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ':    ' +
              'Running samtools faidx.', file=params.logfile)
        return_value = call([params.samtools_executable, 'faidx', genome_fasta])
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': samtools faidx failed', params.logfile)
        index_prefix = genome_fasta
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Indexing completed.', file=params.logfile)
    else:
        if params.index_location.endswith('.fa'):
            assert os.path.exists(params.index_location), 'Index file not found'
            index_prefix = params.index_location
        else:
            fastas = [x for x in os.listdir(params.index_location) if
                      x.endswith(".fa")]
            if len(fastas) == 1:
                index_prefix = "".join([params.index_location, '/', fastas[0]])
            elif len(fastas) == 0:
                raise pi_errors.InputFileError(
                    dt.now().strftime('%I:%M %p %b %d, %Y') + \
                    ': No valid fasta found in provided index folder',
                    params.logfile)
            else:
                raise pi_errors.InputFileError(
                    dt.now().strftime('%I:%M %p %b %d, %Y') + \
                    ':Multiple fastas found in provided index folder. Try ' + \
                    'running with --index_location /path/to/file/filename.fa',
                    params.logfile)

    # Move to working directory before doing I/O intensive alignment
    os.chdir(params.working_dir)

    # Align reads to sam file
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Aligning' +
          ' reads to reference.', file=params.logfile)
    bwa_call = [params.bwa_executable, 'mem'] # base call
    bwa_call.extend(['-t', str(params.n)])  # Number of threads
    #bwa_call.extend(['-R', read_group])  # Read group
    bwa_call.append(index_prefix)  # bwa index
    bwa_call.append(''.join([params.file_path, '/', params.file_prefix,
                             '_1.fastq']))
    bwa_call.append(''.join([params.file_path, '/', params.file_prefix,
                             '_2.fastq']))
    print(' '.join(bwa_call), file=params.logfile)
    with open(''.join([params.file_prefix, '.sam']), 'w') as samfile, \
            open(''.join([params.file_prefix, '_bwa_log.txt']), 'w') as logfile:
        return_value = call(bwa_call, stdout=samfile, stderr=logfile)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': bwa mem failed.', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Alignment completed. Converting to bam', file=params.logfile)
    # Convert the sam to a bam file
    with open(''.join([params.file_prefix, '.bam']), 'w') as bamfile:
        call([params.samtools_executable, 'view', '-bS',
              ''.join([params.file_prefix, '.sam'])], stdout=bamfile)
    call(['rm', ''.join([params.file_prefix, '.sam'])])
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': bam file' +
          ' created. Preparing file for inserting RG into header.',
          file=params.logfile)
    # Fix PG line
    sam_header = check_output([params.samtools_executable, 'view', '-H',
                               ''.join([params.file_prefix, '.bam'])])
    sam_header = sam_header.strip().split('\n')  # Strip whitespace and separate
    pg_line = sam_header[-1].split('\t')  # Grab @PG line + split by tab
    # Then remove the CL field form the PG line
    sam_header[-1] = '\t'.join([x for x in pg_line if not x.startswith('CL')])
    with open(''.join([params.file_prefix, '_sam.header']), 'w') as hdr_file:
        print('\n'.join(sam_header), file=hdr_file)
    with open(''.join([params.file_prefix, '_fixPG.bam']), 'w') as \
              fixpg_bamfile:
        return_value = call([params.samtools_executable, 'reheader',
                             ''.join([params.file_prefix, '_sam.header']),
                             ''.join([params.file_prefix, '.bam'])],
                            stdout=fixpg_bamfile)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': samtools reheader failed', params.logfile)
    call(['rm', ''.join([params.file_prefix, '.bam']),
          ''.join([params.file_prefix, '_sam.header'])])
    # Sort and Index the _fixPG.bam file
    return_value = call([params.samtools_executable, 'sort',
                         ''.join([params.file_prefix, '_fixPG.bam']),
                         ''.join([params.file_prefix, '_fixPG_sorted'])])
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': samtools sort failed.', params.logfile)
    return_value = call([params.samtools_executable, 'index',
                         ''.join([params.file_prefix, '_fixPG_sorted.bam'])])
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': samtools index failed.', params.logfile)
    call(['rm', ''.join([params.file_prefix, '_fixPG.bam'])])
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Inserting @RG tag into header.', file=params.logfile)
    # Reheader the indexed _fixPG_sorted.bam to prepare for mutect
    picard_call = [params.java_executable, ''.join(['-Xmx', params.java_Xmx]),
                   '-jar'] #  Base java call
    picard_call.append(params.picard_jar)  # picard
    picard_call.append('AddOrReplaceReadGroups')  # module
    picard_call.append('CREATE_INDEX=true')
    picard_call.append(''.join(['I=', params.file_prefix, '_fixPG_sorted.bam']))
    picard_call.append(''.join(['O=', params.file_prefix,
                                '_fixPG_sorted_reheader.bam']))
    picard_call.append('SO=coordinate')
    picard_call.append('ID=1')
    picard_call.append(''.join(['LB=', params.file_prefix]))
    picard_call.append('PL=ILLUMINA')
    picard_call.append('PU=12345')
    picard_call.append(''.join(['SM=', params.sample_type]))
    with open(''.join([params.file_prefix, '_picard_log.txt']), 'w') as logfile:
        return_value = call(picard_call, stdout=logfile)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': picard AddOrReplaceReadGroups failed.', params.logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': @RG ' +
          'inserted. Indexing bam', file=params.logfile)
    # Index _fixPG_sorted_reheader.bam file
    return_value = call([params.samtools_executable, 'index',
                         ''.join([params.file_prefix,
                                  '_fixPG_sorted_reheader.bam'])])
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': samtools index failed.', params.logfile)
    # Remove intermediate files
    call(['rm', ''.join([params.file_prefix, '_fixPG_sorted.bam']),
          ''.join([params.file_prefix, '_fixPG_sorted.bam.bai'])])

    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Alignment completed. Finishing up...', params.logfile)
    # Move files from temp directory to outdir
    prepare.move_output(params)
    print('RESULT ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Process ' +
          'completed', file=params.logfile)
    params.logfile.close()


if __name__ == "__main__":
    sys.exit(main())
