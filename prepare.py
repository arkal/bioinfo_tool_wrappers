#!/usr/bin/env python2.7
'''
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/prepare.py

This file contains the modules used by multiple programs in the docker_scripts
folder. The modules contained withing include:
py_which(): A python implementation of the bash 'which' command,
get_genome(): A module to download the genome to a specified directory.
get_gtf(): A module to download the gtf file for the specified genome to a
           specified directory.
'''
from __future__ import division, print_function
from subprocess import call, check_call, check_output, CalledProcessError
from datetime import datetime as dt

import sys
import argparse
import os
import gzip
import pi_errors

class MyUniversalHelpFormatter(argparse.HelpFormatter):
    '''
    This formatter formats both the description and argument defaults formatting
    of the argparse help string.
    '''
    def _fill_text(self, text, width, indent):
        '''
        This module was taken from the ArgumentDefaultsHelpFormatter class
        within argparse.  It deals with the formatting of arguments in that it
        appends the default value to teh description of each argument.
        '''
        return ''.join(indent + line for line in text.splitlines(True))
    def _get_help_string(self, action):
        '''
        This module was taken from the RawDescriptionHelpFormatter class
        within argparse.  It deals with the formatting of the description string
        and allows properly formatted descriptions ot be printed without line
        wrapping.
        '''
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help

# Define the constant for current gencode release version for hg38
GENCODE_HG38_CURR = 23

def allowed_chrom_list():
    '''
    This returns the list of all options the -C flag in radia can accept.
    '''
    chromlist = [str(i) for i in range(1, 23) +['X', 'Y']]
    chromlist += [''.join(['chr', str(i)]) for i in range(1, 23)+['X', 'Y']]
    chromlist += ['all']
    return chromlist


def parse_args(program_doc, program, out_postfix):
    '''
    This function is the all-encompassing argument parser for all the functions
    in the docker scripts project.

    Unless specified, the parser will look for default executables on $PATH.
    The program DOES NOT look for jar files and they are required to be
    passed during execution.

    The program currently only indexes using hg19. Other options will be
    made available in the future.
    '''
    # Universal required argument parser
    universal_required = argparse.ArgumentParser(add_help=False)
    univ_req = universal_required.add_argument_group(
        ''.join(['Mandatory arguments required by most scripts in the package'])
        )
    univ_req.add_argument('--file_prefix', '-f', dest='prefix',
                          help='Path to  input file with file_prefix. Eg. ' +
                          '~/abc/def/ghi if the input files are ghi_1.fastq ' +
                          'and ghi_2.fastq, or ghi.bam, or ghi.vcf. For ' +
                          'snp callers this is the prefix name used in the ' +
                          'ouput file', type=str,
                          required=True)
    # Universal optional argument parser
    universal_optional = argparse.ArgumentParser(add_help=False)
    univ_opt = universal_optional.add_argument_group(
        ''.join(['Optional arguments shared by all scripts in the package']))
    univ_opt.add_argument('--out_prefix', dest='out_prefix', help='Prefix ' +
                          'to use on output files.', type=str, required=False,
                          default='FILE_PREFIX')
    univ_opt.add_argument('--outdir', '-o', dest='outdir', help='Desired ' +
                          'directory to store output in.', type=str,
                          required=False, default='INPUT_DIR')
    univ_opt.add_argument('--working_directory', '-W', dest='working_dir',
                          help='Optional directory to work in. This is ' +
                          'useful in a cluster environment to prevent excess ' +
                          'I/O on a single file system (Use /scratch for ' +
                          'working and the results will be moved to the ' +
                          'appropriate directory on completion.', type=str,
                          required=False, default='OUTDIR')
    univ_opt.add_argument('--genome_version', '-G', dest='genome_version',
                          help='What genome version should be used for the ' +
                          'analysis?',
                          choices=['hg19', 'hg38', 'GRCh37', 'GRCh38'],
                          type=str, required=False, default='hg19')
    univ_opt.add_argument('--genome_fasta', '-F', dest='genome_fasta',
                          help='Path to genomic fasta file. If this is ' +
                          'not provided, the value of GENOME_VERSION ' +
                          'will be used to download the file. The file will ' +
                          'be stored in order of preference in the index/vcf' +
                          ' location followed by the output directory.',
                          type=str, required=False, default='DOWNLOAD')
    univ_opt.add_argument('--num_threads', '-n', dest='n', help='Number ' +
                          'of threads to use.', type=int, required=False,
                          default=4)
    univ_opt.add_argument('--Logfile', '-L', dest='logfile', help='File ' +
                          'where script logs whould be written.', type=str,
                          required=False, default='STDERR')
    univ_opt.add_argument('--use_docker', dest='use_docker', help='Use ' +
                          'docker images?', type=bool,
                          required=False, default='True')
    univ_opt.add_argument('--dockerhub', dest='dockerhub', help='Hub where ' +
                          'docker images may be obtained.', type=str,
                          required=False, default='aarjunrao')
    # Indexing options parser
    index_general = argparse.ArgumentParser(add_help=False)
    req_index_grp = index_general.add_argument_group('Indexing related options')
    req_index_grp.add_argument('--index_location', dest='index_location',
                               help='If already indexed/referenced, where are' +
                               'the indexes/references stored?', type=str,
                               required=False, default=None)
    req_index_grp.add_argument('--require_index', dest='index_destination',
                               help='If no indexes/references exist, specify ' +
                               'a directory to store the indexes/references.',
                               type=str, required=False, default=None)
    # twoBitToFa options parser
    tbtf_option = argparse.ArgumentParser(add_help=False)
    tbtf_group = tbtf_option.add_argument_group('twoBitToFa options')
    tbtf_group.add_argument('--twoBitToFa', dest='tbtf_executable',
                            help='Path to twoBitToFa.', type=str,
                            required=False, default=py_which('twoBitToFa'))
    # Java options parser
    java_options = argparse.ArgumentParser(add_help=False)
    java_group = java_options.add_argument_group('Java related options')
    java_group.add_argument('--java', '-J', dest='java_executable',
                            help='Path to java.', type=str,
                            required=False, default=py_which('java'))
    java_group.add_argument('--Xmx', '-X', dest='java_Xmx', help='Virtual' +
                            ' memory to allocate to java for picard reheader.' +
                            ' Eg 2G, 100M, 100G, etc.', type=str,
                            required=False, default='20G')
    # Parser for bowtie
    bowtie_parser = argparse.ArgumentParser(add_help=False)
    bt_group = bowtie_parser.add_argument_group('Bowtie')
    bt_group.add_argument('--bowtie', '-b', dest='bowtie_executable',
                          help='Path to bowtie.', type=str, required=False,
                          default=py_which('bowtie'))
    # Parser for bowtie2
    bowtie2_parser = argparse.ArgumentParser(add_help=False)
    bt2_group = bowtie2_parser.add_argument_group('Bowtie2')
    bt2_group.add_argument('--bowtie2', '-B', dest='bowtie2_executable',
                           help='Path to bowtie2.', type=str, required=False,
                           default=py_which('bowtie2'))
    # Parser for SNP callers
    snpcallers_options = argparse.ArgumentParser(add_help=False)
    snpcallers_group = snpcallers_options.add_argument_group(
        'Options related to SNP callers')
    snpcallers_group.add_argument('--tum_dna_file', '-T', dest='tum_d_file',
                                  help='Path to tumor bam WGS/WXS bam file.',
                                  type=str, required=True)
    snpcallers_group.add_argument('--norm_dna_file', '-N', dest='norm_d_file',
                                  help='Path to normal WGS/WXS bam file.',
                                  type=str, required=True)
    snpcallers_group.add_argument('--vcf_location', dest='vcf_location',
                                  help='Link to the folder where the cosmic ' +
                                  'and dbsnp vcfs are stored. If any of the ' +
                                  'required vcfs are not provided, this is ' +
                                  'where they will be queried for before ' +
                                  'proceeding to download them. This '
                                  'argument is mandatory if any required vcf ' +
                                  'is not provided.', type=str,
                                  required=False, default=None)
    snpcallers_group.add_argument('--cosmic', '-c', dest='cosmic_file',
                                  help='Path to cosmic vcf file. Required ' +
                                  'if Cosmic_sorted.vcf not present in ' +
                                  'VCF_LOCATION.', type=str, required=False,
                                  default='DOWNLOAD')
    snpcallers_group.add_argument('--cosmic_version', dest='cosmic_version',
                                  help='Version of cosmic to use if cosmic' +
                                  'need to be downloaded.', type=int,
                                  required=False, default='73')
    snpcallers_group.add_argument('--Sanger_username', '-U', dest='username',
                                  help='Username at the sanger website used ' +
                                  'to download the cosmic vcf. This is '+
                                  'required only if the --cosmic flag isn\'t ' +
                                  'provided AND CosmicCodingMuts.vcf isn\'t ' +
                                  'found in VCF_LOCATION. A sanger account ' +
                                  'can be made at https://cancer.sanger.ac.uk' +
                                  '/cosmic/register', type=str, required=False,
                                  default=None)
    snpcallers_group.add_argument('--Sanger_password', '-P', dest='password',
                                  help='Username at the sanger website used ' +
                                  'to download the cosmic vcf. USERNAME AND ' +
                                  'PASSWORD IS NOT STORED BY THIS PROGRAM.',
                                  type=str, required=False, default=None)
    snpcallers_group.add_argument('--dbsnp', '-s', dest='dbsnp_file',
                                  help='Path to dbsnp vcf file. Required ' +
                                  'if 00-All.vcf not present in VCF_LOCATION.',
                                  type=str, required=False, default='DOWNLOAD')
    snpcallers_group.add_argument('--dbsnp_genome_patch',
                                  dest='dbsnp_genome_patch', help='Patch of ' +
                                  'the genome for which the dbsnp vcf should ' +
                                  'be downloaded. Defaults to 13 for hg19 and' +
                                  ' 2 for hg38. Verify that the patch is ' +
                                  'current yourself at ' +
                                  'ftp://ftp.ncbi.nih.gov/snp/organisms/',
                                  type=str, required=False, default='LATEST')
    snpcallers_group.add_argument('--dbsnp_version', dest='dbsnp_version',
                                  help='Version of dbsnp to use if dbsnp' +
                                  'need to be downloaded.', type=int,
                                  required=False, default=144)
    # MHC caller dependent parsers
    mhc_options = argparse.ArgumentParser(add_help=False)
    mhc_group = mhc_options.add_argument_group(
            'Options related to MHC predictors')
    mhc_group.add_argument('--mhc_predictor', dest='mhc_executable', type=str,
                           help='path to predict_binding.py for MHCI or ' +
                           'mhc_II_binding.py for MHCII.', required=True)
    mhc_group.add_argument('--prediction_method', dest='pred_meth',
                           type=str, help='IEDB prediction method to use.',
                           required=False, default='IEDB_recommended')
    mhc_group.add_argument('--alleles', dest='alleles', type=str,
                           help='allele(s) to predict over. This can be ' +
                           'passed as a list of space-separated values, OR ' +
                           'as a ".alleles" file which is a text file ' +
                           'containing one allele per line.', required=True,
                           nargs='+')
    # Define parsers for individual programs
    if program == 'cutadapt':
        cutadapt_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional])
        cutadapt_group = cutadapt_parser.add_argument_group(
            ''.join(['Options for run_cutadapt.py']))
        cutadapt_group.add_argument('-c', '--cutadapt',
                                    dest='cutadapt_executable', help='Path ' +
                                    'to cutadapt executable.', type=str,
                                    required=False,
                                    default=py_which('cutadapt'))
        cutadapt_group.add_argument('-a', '--fwd_3pr_adapter',
                                    dest='fwd_3pr_adapter', help='Adapter' +
                                    ' sequence to trim from the 3\' end of' +
                                    ' the forward read in the pair', type=str,
                                    required=False, default='AGATCGGAAGAG')
        cutadapt_group.add_argument('-A', '--rev_3pr_adapter',
                                    dest='rev_3pr_adapter', help='Adapter' +
                                    ' sequence to trim from the 3\' end of' +
                                    ' the reverse read in the pair', type=str,
                                    required=False, default='AGATCGGAAGAG')
        params = cutadapt_parser.parse_args()
    # Parser specific to run_bwa.py
    elif program == 'bwa':
        bwa_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    index_general, tbtf_option, java_options])
        bwa_group = bwa_parser.add_argument_group('Options for run_bwa.py')
        bwa_group.add_argument('--sample_type', '-s', dest='sample_type',
                               help='Is the sample Tumor or Normal DNA?',
                               type=str, required=True)
        bwa_group.add_argument('--RGID', '-R', dest='RGID', help='Read ' +
                               'Group information to write in the BAM file.',
                               type=str, required=False, default=None)
        bwa_group.add_argument('--bwa', '-B', dest='bwa_executable',
                               help='Path to bwa executable.', type=str,
                               required=False, default=py_which('bwa'))
        bwa_group.add_argument('--samtools', '-S',
                               dest='samtools_executable', help='Path to' +
                               ' samtools.', type=str, required=False,
                               default=py_which('samtools'))
        bwa_group.add_argument('--picard_jar', '-P', dest='picard_jar',
                               help='Path to the picard jar file.',
                               type=str, required=True)
        params = bwa_parser.parse_args()
    # Parser specific to run_rsem.py
    elif program == 'rsem':
        rsem_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required,
                                    universal_optional, index_general,
                                    tbtf_option, bowtie_parser, bowtie2_parser])
        rsem_group = rsem_parser.add_argument_group('Options for run_rsem.py')
        rsem_group.add_argument('--rsem_path', '-R', dest='rsem_path',
                                help='Path to the rsem binaries', type=str,
                                required=False,
                                default=py_where('rsem-prepare-reference'))
        params = rsem_parser.parse_args()
    elif program == 'phlat':
        phlat_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    index_general, bowtie2_parser])
        phlat_group = phlat_parser.add_argument_group(
            'Options for run_phlat.py')
        phlat_group.add_argument('--phlat', '-p', dest='phlat_executable',
                                 type=str, help='Path to PHLAT.py.',
                                 required=False, default=py_which('PHLAT.py'))
        phlat_group.add_argument('--gdownpl', '-g', dest='gdownpl_executable',
                                 type=str, help='Path to the gdown.pl script',
                                 required=False, default=py_which('gdown.pl'))
        params = phlat_parser.parse_args()
    elif program == 'mutect':
        mutect_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_optional, index_general,
                                    java_options, snpcallers_options])
        mutect_group = mutect_parser.add_argument_group(
            'Options for run_mutect.py')
        mutect_group.add_argument('--mutect', '-M', dest='mutect_jar',
                                  help='Path to mutect jar file,', type=str,
                                  required=True)
        params = mutect_parser.parse_args()
    elif program == 'radia' or program == 'snpeff':
        # snpEff options
        snpeff_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    index_general, tbtf_option, java_options])
        snpeff_group = snpeff_parser.add_argument_group(
            'Options for run_snpeff.py')
        snpeff_group.add_argument('--snpeff', dest='snpeff_jar', help='Path ' +
                                  'to snpeff jar file,', type=str,
                                  required=True)
        snpeff_group.add_argument('--intergenic', dest='use_intergenic',
                                  help='Look at intergenic regions?',
                                  action='store_true', required=False,
                                  default=False)
        snpeff_group.add_argument('--downstream', dest='use_downstream',
                                  help='Look at downstream regions?',
                                  action='store_true', required=False,
                                  default=False)
        snpeff_group.add_argument('--upstream', dest='use_upstream',
                                  help='Look at upstream regions?',
                                  action='store_true', required=False,
                                  default=False)
        snpeff_group.add_argument('--no_canonical', dest='no_canonical',
                                  help='Look at all transcripts for a gene ' +
                                  'instead of just canonical?',
                                  action='store_true', required=False,
                                  default=False)
        snpeff_group.add_argument('--config', dest='config_file',
                                  help='Config file to point snpeff to the' +
                                  'database. Remember, file paths change when' +
                                  ' you run a tool in docker.', type=str,
                                  required=False, default=None)
        snpeff_group.add_argument('--snpeff_reference', dest='reference_name',
                                  help='SPECIFIC snpEff provided database to ' +
                                  'use.', type=str,
                                  required=False, default=None)
        if program == 'radia':
            radia_parser = argparse.ArgumentParser(
                description=program_doc,
                formatter_class=MyUniversalHelpFormatter,
                add_help=True, parents=[snpcallers_options, snpeff_parser],
                conflict_handler='resolve')
            radia_group = radia_parser.add_argument_group(
                'Options for run_radia.py')
            radia_group.add_argument('--tum_rna_file', '-R', dest='tum_r_file',
                                     help='Path to tumor RNA-Seq bam file.',
                                     type=str, required=False, default=None)
            radia_group.add_argument('--rna_fasta', dest='rna_fasta',
                                     help='Path to rna fasta file to pass as' +
                                     ' a reference to radia.', type=str,
                                     required=False, default='GENOME_FASTA')
            radia_group.add_argument('--chromosome', '-C', dest='chromosome',
                                     help='Chromosome to process.',
                                     choices=allowed_chrom_list(), type=str,
                                     required=False, default='all')
            radia_group.add_argument('--data_source', dest='data_source',
                                     help='The source of the data - used in ' +
                                     'the sample VCF meta tag.', required=False,
                                     default=''.join([os.getenv('USER'), '@',
                                                      check_output(['echo '+ \
                                                                    '$HOSTNAME'
                                                                   ],
                                                                   shell=True
                                                                  ).strip()]))
            radia_group.add_argument('--radia', dest='radia_executable',
                                     type=str, help='Path to radia.',
                                     required=False, default=py_which('radia.py'
                                                                     ))
            radia_group.add_argument('--sequencingPlatform',
                                     dest='seq_platform', help='The ' +
                                     'sequencing platform - used in the ' +
                                     'sample VCF meta tag', type=str,
                                     required=False, default='Illumina')
            radia_group.add_argument('--disease', dest='disease', help='What' +
                                     ' is the disease under study?', type=str,
                                     required=False, default='CANCER')
            radia_group.add_argument('--blacklist', '-b', dest='blacklist',
                                     help='Path to 1000 genomes blacklist ' +
                                     'file. Defaults to the one packaged in ' +
                                     'radia.', type=str, required=False,
                                     default='PACKAGED')
            radia_group.add_argument('--retro', '-r', dest='retrogenes',
                                     help='Path to retrogenes file. Defaults' +
                                     ' to the one packaged in radia.', type=str,
                                     required=False, default='PACKAGED')
            radia_group.add_argument('--pseudo', '-p', dest='pseudogenes',
                                     help='Path to pseudogenes file. ' +
                                     'Defaults to the one packaged in radia.',
                                     type=str, required=False,
                                     default='PACKAGED',)
            radia_group.add_argument('--broad', '-B', dest='broad_targets',
                                     help='Path to broad targets file. ' +
                                     'Defaults to the one packaged in radia.',
                                     type=str, required=False,
                                     default='PACKAGED')
            radia_group.add_argument('--rnablacklist', dest='rna_blacklist',
                                     help='Path to rna blacklist file. ' +
                                     'Defaults to the one packaged in radia.',
                                     type=str, required=False,
                                     default='PACKAGED')
            radia_group.add_argument('--rnafamilyblacklist',
                                     dest='rna_family_blacklist', help='Path' +
                                     ' to rna family blacklist file. Defaults' +
                                     ' to the one packaged in radia.', type=str,
                                     required=False, default='PACKAGED')
            radia_group.add_argument('--usesnpEff',
                                     dest='use_snpeff', help='Should the ' +
                                     'output VCF be snpEff\'ed?',
                                     action='store_true', required=False,
                                     default=False)
            params = radia_parser.parse_args()
        else:
            params = snpeff_parser.parse_args()
    elif program == 'star':
        star_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    index_general, tbtf_option])
        star_group = star_parser.add_argument_group('Options for run_star.py')
        star_group.add_argument('--STAR', '-S', dest='STAR_executable',
                                type=str, help='Path to STAR.',
                                required=False, default=py_which('STAR'))
        star_group.add_argument('--STARlong', '-X', dest='STARlong_executable',
                                type=str, help='Path to STARlong.',
                                required=False, default=py_which('STARlong'))
        params = star_parser.parse_args()
    elif program == 'kallisto':
        kallisto_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    index_general, tbtf_option])
        kallisto_group = kallisto_parser.add_argument_group(
            'Options for run_kallisto.py')
        # Fill in program body here.
        params = kallisto_parser.parse_args()
    elif program == 'mhci':
        mhci_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    mhc_options])
        mhci_group = mhci_parser.add_argument_group(
            'Options for run_mhci_prediction.py')
        mhci_group.add_argument('--peplen', dest='peplen', type=str,
                               help='Size(s) of peptides to predict for.',
                               required=False, default='mhc dependent',
                               nargs='+')
        params = mhci_parser.parse_args()
    elif program == 'mhcii':
        mhcii_parser = argparse.ArgumentParser(
            description=program_doc,
            formatter_class=MyUniversalHelpFormatter,
            add_help=True, parents=[universal_required, universal_optional,
                                    mhc_options])
        mhcii_group = mhcii_parser.add_argument_group(
            'Options for run_mhcii_prediction.py')
        mhcii_group.add_argument('--netMHCIIpan', dest='netmhciipan_executable',
                               type=str, help='Path to netMHCIIpan. ' + \
                               'Required if IEDBtools fails', required=False,
                               default=py_which('netMHCIIpan'))
        params = mhcii_parser.parse_args()
    else:
        print('This module has to be passed one of:', file=sys.stderr)
        print('\tcutadapt', file=sys.stderr)
        print('\tbwa', file=sys.stderr)
        print('\trsem', file=sys.stderr)
        print('\tphlat', file=sys.stderr)
        print('\tstar', file=sys.stderr)
        print('\tmutect', file=sys.stderr)
        print('\tradia', file=sys.stderr)
        print('\tsnpeff', file=sys.stderr)
        print('\tkallisto', file=sys.stderr)
        print('\tmhc', file=sys.stderr)
        raise pi_errors.ParameterError('', sys.stderr)

    #  General argument handling
    #  First set up the Log file
    if params.logfile == 'STDERR':
        params.logfile = sys.stderr
    else:
        try:
            params.logfile = open(params.logfile, 'w', 0) #  Unbuffered
        except IOError:
            raise pi_errors.InputFileError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                'Could not open log file. Please check the path and try again.',
                sys.stderr)
    #  Ensure either an index is provided or the user has specified where to
    #  store the produced one.
    if program != 'cutadapt' and params.index_location is None and \
            params.index_destination is None:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Please specify either the index location, or a directory ' +
            'wherein the indexes may be stored or located.', params.logfile)
    #  Split input file into filename and file path
    params.file_path, params.file_prefix = \
        os.path.split(os.path.abspath(params.prefix))
    params.file_prefix = params.file_prefix.rstrip("_")
    #  If the user has provided the full file name
    if params.file_prefix.endswith('.vcf') or \
            params.file_prefix.endswith('.bam') or \
            params.file_prefix.endswith('.fastq') or \
            params.file_prefix.endswith('.fa'):
        params.file_prefix, file_ext = os.path.splitext(params.file_prefix)
        if file_ext == '.fastq':
            params.file_prefix = params.file_prefix[:-2] # discard the _1 or _2
    #  Setup the output directory
    if params.outdir == 'INPUT_DIR':
        params.outdir = params.file_path
    params.outdir = ''.join([os.path.abspath(params.outdir), '/',
                             params.file_prefix, '_', out_postfix])
    if not os.path.exists(params.outdir):
        py_mkdir(params.outdir)
    #  Setup the working directory
    if params.working_dir == 'OUTDIR':
        params.working_dir = params.outdir
    else:
        params.working_dir = os.path.abspath(params.working_dir)
    params.working_dir = ''.join([params.working_dir, '/', params.file_prefix,
                                  '_temp'])
    if not os.path.exists(params.working_dir):
        py_mkdir(params.working_dir)
    #  Modify genome versions to UC notation
    if params.genome_version == 'GRCh37':
        params.genome_version = 'hg19'
    elif params.genome_version == 'GRCh38':
        params.genome_version = 'hg38'
    #  If out_prefix hasn't been provided, set it up. Mutect is slightly
    #  different form the rest.
    if params.out_prefix == 'FILE_PREFIX':
        if program == 'mutect':
            params.out_prefix = '_'.join([os.path.split(params.tum_d_file)[1],
                                          os.path.split(params.norm_d_file)[1],
                                          'mutect', 'results'])
        else:
            params.out_prefix = params.file_prefix
    return params


def py_which(prog):
    '''
    This function accepts a program name as a string and returns the default
    executable found on the user's PATH.
    '''
    prog_exists = True
    with open(os.devnull, 'w') as dev_null:
        try:
            check_call(['which', prog], stderr=dev_null, stdout=dev_null)
        except CalledProcessError:
            prog_exists = False
            prog_exec = '404 - Not found'
    if prog_exists:
        prog_exec = check_output(['which', prog]).strip()
    return prog_exec


def get_genome(genome_version, index_path, tbtf_executable, logfile):
    '''
    This function accepts a genome version to download, and a directory for
    storing the downloaded 2bit file from ucsc genome browser. The program also
    accepts the path to the twoBitToFa executable to convert the downloaded 2bit
    file to fasta, and a logfile to write progress and error messages. The
    function returns the full path to the downloaded fasta.
    '''
    if genome_version not in ['hg19', 'hg38']:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': get_genome was passed an incompatible genome version.', logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ':    Downl' +
          'oading fasta reference in 2bit format.', file=logfile)
    with open(''.join([index_path, '/', genome_version, '.2bit']), 'w') as \
          twobit_f, open(os.devnull, 'w') as dev_null:
        return_value = call(['curl', '-L', '-s',
                             ''.join(['http://hgdownload.cse.ucsc.edu/golden',
                                      'Path/', genome_version, '/bigZips/',
                                      genome_version, '.2bit'])],
                            stdout=twobit_f, stderr=dev_null)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Could not download 2bit reference file.', logfile)
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ':    2bit ' +
          'downloaded... Converting to fasta.', file=logfile)
    with open(''.join([index_path, '/', 'chr_names.list']), 'w') as chr_list:
        for chrom in range(1, 23) + ['X', 'Y', 'M']:
            print(''.join(['chr', str(chrom)]), file=chr_list)
    return_value = call([tbtf_executable,
                         ''.join(['-seqList=', index_path, '/',
                                  'chr_names.list']),
                         ''.join([index_path, '/', genome_version, '.2bit']),
                         ''.join([index_path, '/', genome_version, '.fa'])])
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': 2bit -> Fa conversion failed.', logfile)
    call(['rm', ''.join([index_path, '/', genome_version, '.2bit']),
          ''.join([index_path, '/', 'chr_names.list'])])
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Genomic ' +
          'fasta downloaded.', file=logfile)
    return ''.join([index_path, '/', genome_version, '.fa'])

def get_gtf(genome_version, index_path, logfile):
    '''
    This function will download the appropriate gtf file for a provided genome
    version, and will store it in the specified location. It returns the full
    path to the downloaded gtf file.

    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Downl' +
          'oading gtf file.', file=logfile)
    if genome_version == 'hg19':
        release_version = '19'
    elif genome_version == 'hg38':
        release_version = str(GENCODE_HG38_CURR)
    else:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': get_gtf was passed a non-compatible genome version.', logfile)
    gencode_file = ''.join([index_path, '/',
                            'gencode.v', release_version, '.annotation.gtf'])
    with open(''.join([gencode_file, '.gz']), 'w') as gff, \
          open(os.devnull, 'w') as dev_null:
        return_value = call(['curl', '-L',
                             ''.join(['ftp://ftp.sanger.ac.uk/pub/gencode/' +
                                      'Gencode_human/release_', release_version,
                                      '/gencode.v', release_version,
                                      '.annotation.gtf.gz'])], stdout=gff,
                            stderr=dev_null)
    if return_value != 0:
        raise pi_errors.MyRuntimeError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': Could not download the gtf file for the specified genome.',
            logfile)
    call(['gunzip', ''.join([gencode_file, '.gz'])])
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': Gencode ' +
          'annotation file downloaded.', file=logfile)
    return gencode_file


def py_mkdir(directory):
    '''
    This function accepts a directory path as a string and will attempt to
    create the directory, recursively if it's base directory does not exist.
    '''
    # Don't do anything if the directory exists
    if os.path.exists(directory):
        return None
    # Base Case
    if (os.path.exists(os.path.split(directory)[0]) or
            os.path.split(directory)[0] == ''):
        os.mkdir(directory)
    else:
        py_mkdir(os.path.split(directory)[0])
        os.mkdir(directory)
    return None


def py_where(prog):
    '''
    This function accepts a program as a string and returns the location of
    the file on the users PATH. It uses py_which  from the prepare.py file to
    get the path to the binary and then splits that to get the path to the
    program
    '''
    where = py_which(prog)
    if where == '404 - Not found':
        return '404 - Not found'
    else:
        return os.path.split(where)[0]


def move_output(params):
    '''
    This module is the last in the pipeline for each tool. It is used to move
    the output files from th working directory to the output directory.
    '''
    for out_file in os.listdir(params.working_dir):
        call(['mv', out_file, params.outdir])
    os.rmdir(params.working_dir)
    return None


def download_vcf(vcf_type, params):
    '''
    This module is used to download either a dbsnp or cosmic vcf for the given
    genome version.
    params contains the following variables:
        GENERAL
        logfile - the file to write logs to
        genome_version - the genome version to use [hg19, hg38]
        vcf_location - Where the vcf should be stored
        DBSNP SPECIFIC
        dbsnp_version - version of dbsnp to use
        dbsnp_genome_patch - patch of the genome to use
        COSMIC SPECIFIC
        cosmic_version - version of cosmic to use
        username - username to download from cosmic (sanger)
        password - password to download from cosmic (sanger)
    '''
    #  Only dbsnp and cosmic are allowed
    if vcf_type not in ['dbsnp', 'cosmic']:
        raise pi_errors.ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': download_vcf saw something weird and freaked out.',
            params.logfile)
    #  Process download dbsnp
    if vcf_type == 'dbsnp':
        #   Create the vcf location if it doesn't exist
        if not os.path.exists('/'.join([params.vcf_location,
                                        params.genome_version])):
            py_mkdir('/'.join([params.vcf_location, params.genome_version]))
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Downloading dbSNP vcf from NCBI.', file=params.logfile)
        dbsnp_file = '/'.join([params.vcf_location, params.genome_version,
                               '00-All.vcf.gz'])
        if params.genome_version == 'hg19':
            dbsnp_genome_version = 'GRCh37'
            if params.dbsnp_genome_patch == 'LATEST':
                params.dbsnp_genome_patch = 13
            else:
                params.dbsnp_genome_patch = int(params.dbsnp_genome_patch)
        else:
            dbsnp_genome_version = 'GRCh38'
            if params.dbsnp_genome_patch == 'LATEST':
                params.dbsnp_genome_patch = 2
            else:
                params.dbsnp_genome_patch = int(params.dbsnp_genome_patch)
        dbsnp_download_link = ''.join(['ftp://ftp.ncbi.nlm.nih.gov/snp/',
                                       'organisms/human_9606_b',
                                       str(params.dbsnp_version), '_',
                                       str(dbsnp_genome_version), 'p',
                                       str(params.dbsnp_genome_patch),
                                       '/VCF/00-All.vcf.gz'])
        with open(dbsnp_file, 'w') as dbsnp_gz:
            dbsnp_call = ['curl', '-L']
            dbsnp_call.append(dbsnp_download_link)
            return_value = call(dbsnp_call, stdout=dbsnp_gz)
            if return_value != 0:
                raise pi_errors.MyRuntimeError(
                    dt.now().strftime('%I:%M %p %b %d, %Y') + \
                    ': Could not download dbsnp vcf to ' + \
                    params.vcf_location + '.', params.logfile)
        try:
            dbsnp_file = correct_chromosomes_in_vcf(dbsnp_file, params.logfile)
        except:  # Catch any error
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Could not gunzip the downloaded dbsnp vcf.', params.logfile)
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Successfully downloaded dbsnp vcf.', file=params.logfile)
        return dbsnp_file.rstrip('.gz')
    elif vcf_type == 'cosmic':
        #   Create the vcf location if it doesn't exist
        if not os.path.exists('/'.join([params.vcf_location,
                                        params.genome_version])):
            py_mkdir('/'.join([params.vcf_location, params.genome_version]))
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Downloading cosmic vcf.', file=params.logfile)
        # Ensure user has provided credentials
        if params.username is None or params.password is None:
            raise pi_errors.ParameterError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Require ID and password to download from Sanger',
                params.logfile)
        # Ensure system has lftp installed
        lftp_executable = py_which('lftp')
        if lftp_executable == '404 - Not found':
            raise pi_errors.RequirementError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Require lftp to download from Sanger', params.logfile)
        cosmic_file = '/'.join([params.vcf_location, params.genome_version,
                                'CosmicCodingMuts.vcf.gz'])
        if params.genome_version == 'hg19':
            cosmic_genome_version = 'grch37'
        else:
            cosmic_genome_version = 'grch38'
        cosmic_download_link = ''.join(['get /files/',
                                        str(cosmic_genome_version), '/cosmic/v',
                                        str(params.cosmic_version),
                                        '/VCF/CosmicCodingMuts.vcf.gz'])
        with open('get_cosmic.sh', 'w') as get_cosmic_sh:
            print('lftp -u \"', params.username, '\",\"', params.password,
                  '\" sftp://sftp-cancer.sanger.ac.uk << EOF', sep='',
                  file=get_cosmic_sh)
            print(cosmic_download_link,
                  file=get_cosmic_sh)
            print('bye\nEOF', file=get_cosmic_sh)
        return_value = call(['sh', 'get_cosmic.sh'])
        if return_value != 0:
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Could not download cosmic vcf.', params.logfile)
        try:
            call(['mv', 'CosmicCodingMuts.vcf.gz', cosmic_file])
            call(['rm', 'get_cosmic.sh'])
            cosmic_file = correct_chromosomes_in_vcf(cosmic_file,
                                                     params.logfile)
        except:  # Catch any error
            raise pi_errors.MyRuntimeError(
                dt.now().strftime('%I:%M %p %b %d, %Y') + \
                ': Could not gunzip the downloaded cosmic vcf.', params.logfile)
        print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
              'Successfully downloaded cosmic vcf.', file=params.logfile)
        return cosmic_file.rstrip('.gz')
    else:
        if params.logfile.name != '<stderr>':
            params.logfile.close()
        raise SystemExit('Shouldn\'t see this ever')
    return None

def correct_chromosomes_in_vcf(vcfgz_file_in, logfile):
    '''
    This module takes in a gzipped vcf file (dbsnp/cosmic) and prefixes 'chr'
    to each snp record.  If the input is input.vcf.gz, the output file is
    input_chrfixed.vcf
    '''
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': ' +
          'Unzipping vcf file while fixing chromosome names.  This might take' +
          ' a couple of minutes', file=logfile)
    # 2 os.path.plisexts for .gz then .vcf
    vcf_file_out = '_'.join([os.path.splitext(os.path.splitext(vcfgz_file_in
                                                              )[0])[0],
                             'chrfixed.vcf'])
    with gzip.GzipFile(vcfgz_file_in, 'r') as infile, \
            open(vcf_file_out, 'w') as outfile:
        for line in infile:
            #  Leave header lines alone
            if line.startswith('#'):
                print(line, sep='', end='', file=outfile)
                continue
            #  Chromosome M should be chrM not chrMT
            elif line.startswith('MT'):
                print('chrM', line[2:], sep='', end='', file=outfile)
            #  Prefix chr to every other record
            else:
                print('chr', line, sep='', end='', file=outfile)
    #  Remove the input gzipped file
    call(['rm', vcfgz_file_in])
    print('PROGRESS ' + dt.now().strftime('%I:%M %p %b %d, %Y') + ': vcf ' +
          'unzipped.', file=logfile)
    return vcf_file_out

# Obtained and modified these functions from John's rna-seq pipeline
def docker_path(filepath):
    """
    Given a path, returns that files path inside the docker mount directory
    (/data).
    """
    return os.path.join('/data', os.path.basename(filepath))


def docker_call(tool, tool_parameters, work_dir):
    """
    Makes subprocess call of a command to a docker container.
    work_dir MUST BE AN ABSOLUTE PATH or the call will fail.
    """
    base_docker_call = 'sudo docker run -v {}:/data'.format(work_dir)
    call = base_docker_call.split() + [tool] + tool_parameters
    try:
        subprocess.check_call(call)
    except subprocess.CalledProcessError:
        raise RuntimeError('docker command returned a non-zero exit status for cmd {}'.format(call))
    except OSError:
        raise RuntimeError('docker not found on system. Install on all nodes.')