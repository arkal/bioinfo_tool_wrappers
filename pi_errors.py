#!/usr/bin/env python2.7
'''
Author : Arjun Arkal Rao
Affiliation : UCSC BME, UCSC Genomics Institute
File : docker_scripts/pi_errors.py

This file contains the assertions and generic error handlers used by all the
scripts in the precision immunology pipeline.
'''
from __future__ import division, print_function
from datetime import datetime as dt

import os

class MyCustomException(Exception):
    '''
    This is a custom version of the Exception Class and will be inherited by
    all custom Errors defined in this file. It accepts 2 arguments
    err_meg - The Error message passed to the Exception that will be printed on
              exit.
    logfile - The logfile where logs are being printed. This needs to be closed
              when an Error is raised. The default for logfile in the precision
              immuno pipeline is sys.stderr and that should NOT be closed as the
              Error needs to be printed there.
    '''
    def __init__(self, err_mesg, logfile):
        '''
        Redefining the __init__ to accept the error message AND the logfile such
        that the log file can be closed if it isn't STDERR.
        '''
        super(MyCustomException, self).__init__(err_mesg, logfile)
        self.err_mesg = err_mesg
        self.logfile = logfile
        if self.logfile.name != '<stderr>':
            self.logfile.close()

    def __str__(self):
        '''
        This module is used to print the error message passed to this error
        class
        '''
        return repr(self.err_mesg)

#  Exception for bad input files
class InputFileError(MyCustomException):
    '''
    This Error Class will be raised  in the case of a bad input provided.
    '''
    pass


#  Exception for bad parameters provided
class ParameterError(MyCustomException):
    '''
    This Error Class will be raised  in the case of a bad parameter provided.
    '''
    pass

#  Exception for failed requirements
class RequirementError(MyCustomException):
    '''
    This Error Class will be raised  in the case of a bad parameter provided.
    '''
    pass

#  Custom Exception for runtime error
class MyRuntimeError(MyCustomException):
    '''
    This Error Class will be raised  in the case of a bad parameter provided.
    '''
    pass


def test_param_value(params_param_val, param_name, param_flag, logfile):
    '''
    Does the provided value for PARAM_NAME point to a valid file?
        If not : raise erro and quit
        If so  : return the abspath for PARAMS_PARAM_VAL
    '''
    if not os.path.exists(params_param_val):
        raise ParameterError(
            dt.now().strftime('%I:%M %p %b %d, %Y') + \
            ': ' + param_name + ' was not found at specified or default ' + \
            'location. Please verify and retry with ' + param_flag + '.',
            logfile)
        
    else:
        return os.path.abspath(params_param_val)
