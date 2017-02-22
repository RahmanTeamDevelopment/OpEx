#!/usr/bin/env python

"""
Last updated: 15.06.2016
"""

from optparse import OptionParser
from subprocess import call
import os
import sys

##############################################################################################################

# Current directory
scriptdir = os.path.dirname(os.path.realpath(__file__))

# Version
ver = '1.0.0'

# Installation log file
logfile = open('install_log.txt', 'w')

# Command line argument parsing
descr = 'Installer script of OpEx ' + ver + '.'
parser = OptionParser(usage='python install.py <options>', version=ver, description=descr)
parser.add_option("-r", "--reference", default=None, dest='reference', action='store', help="Reference genome file")
(options, args) = parser.parse_args()

# Reference directory
if not options.reference is None:
    refdir = os.path.abspath(options.reference)
else:
    refdir = ''

# Print welcome message
print '\n'+'-'*50
print 'INSTALLING OPEX PIPELINE VERSION ' + ver
print '-'*50+'\n'
print '(More details in install_log.txt)\n'

# Call build_opex.sh script
sys.stdout.write('Downloading and building components (BWA, Stampy, CoverView, Platypus, CAVA) ...')
sys.stdout.flush()
call(['chmod', '+x', './build_opex.sh'])
call(['./build_opex.sh'], stdout = logfile, stderr = logfile)
print ' Done.'

# Create config file
config = open('config.txt', "wt")
config.write('ENSTDB = ' + scriptdir + '/exome_65_GRCh37.gz\n')
config.write('CAVA_CONFIG = ' + scriptdir + '/cava_config.txt\n')

# CoverView config files
call(['cp', 'templates/coverview_config_template', 'CoverView_default.json'])
call(['cp', 'templates/coverview_config_template', 'CoverView_full.json'])
defaultconfig = open('CoverView_default.json', "a")
defaultconfig.write('\t\"only_fail_profiles\": true,\n')
defaultconfig.write('\t\"transcript\":  {\"regions\": false, \"profiles\": false, \"poor\": true },\n')
defaultconfig.write('\t\"transcript_db\": \"' + scriptdir + '/exome_65_GRCh37.gz\"\n}\n')
fullconfig = open('CoverView_full.json', "a")
fullconfig.write('\t\"transcript\":  {\"regions\": true, \"profiles\": true, \"poor\": true },\n')
fullconfig.write('\t\"transcript_db\": \"' + scriptdir + '/exome_65_GRCh37.gz\"\n}\n')
defaultconfig.close()
fullconfig.close()

# CAVA config file
call(['cp', 'templates/cava_config_template', 'cava_config.txt'])
cavaconfig = open('cava_config.txt', "a")
cavaconfig.write('# Name of Ensembl transcript database file\n')
cavaconfig.write('# Possible values: string | Optional: yes (if not given, no transcript-based annotations are reported)\n')
cavaconfig.write('@ensembl = ' + scriptdir + '/exome_65_GRCh37.gz\n')
if not options.reference is None:
    cavaconfig.write('\n# Name of reference genome file\n')
    cavaconfig.write('# Possible values: string | Optional: no\n')
    cavaconfig.write('@reference = ' + refdir + '\n')
cavaconfig.close()

# Call index_genome.sh and add reference fields to config file
if not options.reference is None:
    sys.stdout.write('Adding default reference genome ...')
    sys.stdout.flush()
    call(['chmod', '+x', './index_genome.sh'])
    call(['./index_genome.sh', refdir], stdout=logfile, stderr=logfile)
    config.write('REFERENCE = ' + refdir + '\n')
    config.write('GENOME_INDEX = ' + scriptdir + '/index/ref\n')
    config.write('HASH = ' + scriptdir + '/index/ref\n')
    print ' Done.'
else:
    print '\n!!! Referemce genome must be added later.'

# Close config file and print goodbye message
config.close()
print '\n'+'-'*50
print 'OPEX v' + ver + ' INSTALLATION COMPLETED'
if not options.reference is None: print ' Test installation: python test_installation.py'
print '-'*50+'\n'

# Close log file
logfile.close()

