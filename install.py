#!/usr/bin/env python

from optparse import OptionParser
from subprocess import call
import os
import sys
import stat

# Make script executable
def makeExecutable(fn):
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IEXEC)

##############################################################################################################

# Current directory
scriptdir = os.path.dirname(os.path.realpath(__file__))

# Version
ver = '1.1.0'

# Installation log file
logfile = open('install_log.txt', 'w')

# Command line argument parsing
descr = 'Installer script of OpEx ' + ver + '.'
parser = OptionParser(usage='python install.py <options>', version=ver, description=descr)
parser.add_option("-r", "--reference", default=None, dest='reference', action='store', help="Reference genome file")
(options, args) = parser.parse_args()

if options.reference is None:
    print '\nReference genome must be specified\n'
    quit()

# Reference directory
refpath = os.path.abspath(options.reference)

# Print welcome message
print '\n'+'-'*50
print 'INSTALLING OPEX PIPELINE VERSION ' + ver
print '-'*50+'\n'
print '(More details in install_log.txt)\n'

# Call build_opex.sh script
sys.stdout.write('Building components of pipeline ... ')
sys.stdout.flush()
makeExecutable('build_opex.sh')
call(['./build_opex.sh'], stdout = logfile, stderr = logfile)
print 'Done.'

# Call index_genome.sh script
sys.stdout.write('Indexing reference genome ... ')
sys.stdout.flush()
makeExecutable('index_genome.sh')
call(['./index_genome.sh', refpath], stdout=logfile, stderr=logfile)
print 'Done.'

sys.stdout.write('Creating default configuration file ... ')
sys.stdout.flush()

# CAVA config file
cava_config = open('config/cava_config.txt', "w")
cava_config.write('@ensembl = ' + scriptdir + '/defaultdb/ensembl75s.gz\n')
cava_config.write('@reference = ' + refpath + '\n')
cava_config.close()

# CoverView config file
coverview_config = open('config/coverview_config.txt', "a")
coverview_config.write('\t\"transcript_db\": \"' + scriptdir + '/defaultdb/ensembl75s.gz\"\n}\n')
coverview_config.close()

# Default OpEx config file
default_config = open('config.txt', "w")
default_config.write('CAVA_CONFIG = ' + scriptdir + '/config/cava_config.txt\n')
default_config.write('COVERVIEW_CONFIG = ' + scriptdir + '/config/coverview_config.json\n')
default_config.write('REFERENCE = ' + refpath + '\n')
default_config.write('STAMPY_INDEX = ' + scriptdir + '/index/ref\n')
default_config.write('STAMPY_HASH = ' + scriptdir + '/index/ref\n')
default_config.close()

# Goodbye message
print 'Done.'
print '\n'+'-'*50
print 'OPEX v' + ver + ' INSTALLATION COMPLETED'
print 'Test installation: python test_installation.py'
print '-'*50+'\n'

# Close log file
logfile.close()

