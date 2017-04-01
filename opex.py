#!/usr/bin/env python

import os
import sys
import stat
from optparse import OptionParser
import subprocess
from collections import OrderedDict
import datetime

#########################################################################################################

# ...
def readConfigFile(scriptdir, fn):
    ret = OrderedDict()
    if fn is None: fn = scriptdir + '/config.txt'
    for line in open(fn):
        line = line.strip()
        if line.startswith('#') or line == '': continue
        if '=' not in line: continue
        [key, value] = line.split('=', 1)
        key, value = key.strip().upper(), value.strip()
        if value in ['', '.']: continue
        ret[key] = value

    # Check required options
    for k in ['REFERENCE', 'STAMPY_INDEX', 'STAMPY_HASH', 'CAVA_CONFIG', 'COVERVIEW_CONFIG']:
        if k not in ret: sys.exit('Required option %s not specified in configuration file' % k)

    # Check if file exists and convert to absolute path
    for k in ['REFERENCE', 'CAVA_CONFIG', 'COVERVIEW_CONFIG']:
        if os.path.isfile(ret[k]): ret[k] = os.path.abspath(ret[k])
        else: sys.exit('File specified by option %s (%s) does not exist.' % (k, ret[k]))

    return ret

# ...
def checkInputs(options):
    if options.fastq is None: sys.exit('\nInput files not specified.\n')
    x = options.fastq.split(',')
    if not len(x) == 2: sys.exit('\nIncorrect format for option --input.\n')
    if not x[0].endswith('.fastq.gz') or not x[1].endswith('.fastq.gz'): sys.exit('\nInput files must have .fastq.gz format.\n')
    if options.name is None: sys.exit('\nOutput file name not specified.\n')

# ...
def generateFile(params, fnin, fnout):
    with open(fnout, "wt") as fout:
        with open(fnin, "rt") as fin:
            for line in fin:
                for key, value in params.iteritems():
                    line = line.replace('@' + key, value)
                fout.write(line)

# Make script executable
def makeExecutable(fn):
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IEXEC)

# ...
def executeScript(cmd_str):
    cmd = cmd_str.split()
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""): yield stdout_line.strip()
    popen.stdout.close()
    return_code = popen.wait()
    if return_code: raise subprocess.CalledProcessError(return_code, cmd)

# ...
def report(info):
    msg = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': ' + info
    print msg
    sys.stdout.flush()

##############################################################################################################

scriptdir = os.path.dirname(os.path.realpath(__file__))
workingdir = os.getcwd()

# Version
ver = '1.1.0'

# Command line argument parsing
descr = 'OpEx (Optimised Exome) pipeline ' + ver + '.'
parser = OptionParser(usage='python opex.py <options>', version=ver, description=descr)
parser.add_option('-i', "--input", default=None, dest='fastq', action='store', help="fastq.gz files")
parser.add_option('-o', "--output", default=None, dest='name', action='store', help="Sample name (output prefix)")
parser.add_option('-b', "--bed", default=None, dest='bed', action='store', help="Bed file")
parser.add_option('-t', "--threads", default=1, dest='threads', action='store', help="Number of processes to use")
parser.add_option('-c', "--config", default=None, dest='config', action='store', help="Configuration file")
parser.add_option('-k', "--keep", default=False, dest='keep', action='store_true', help="Keep temporary files")
(options, args) = parser.parse_args()
checkInputs(options)

# Welcome message
print '\n' + '-' * 80
print 'OpEx pipeline version ' + ver
print '-' * 80 + '\n'

# Read configuration file
params = readConfigFile(scriptdir, options.config)

# Additional params
params['NAME'] = options.name
params['FASTQ1'], params['FASTQ2'] = options.fastq.split(',')
params['OPEXDIR'] = scriptdir

params['MORECV'] = '-c ' + params['COVERVIEW_CONFIG']
if options.bed is not None: params['MORECV'] = params['MORECV'] + ' -b ' + options.bed
if int(options.threads) > 1: params['MORECV'] = params['MORECV'] + ' -t ' + str(options.threads)

params['MORECAVA'] = ''
if int(options.threads) > 1: params['MORECAVA'] = params['MORECAVA'] + '-t ' + str(options.threads)

if options.keep: params['KEEPREMOVE'] = ''
else: params['KEEPREMOVE'] = 'rm -r ' + params['NAME'] + '_tmp'

# Genearate Bash script file
scriptfn = params['NAME'] + '_opex_pipeline.sh'
generateFile(params, scriptdir + '/templates/opex_pipeline_template', scriptfn)
makeExecutable(scriptfn)

# ...
logf = open(options.name + '_opex_log.txt', 'w')
for stdout_line in executeScript('./'+scriptfn):
    if stdout_line.startswith('OPEXMSG'):
        print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + ': ' + stdout_line[8:]; sys.stdout.flush()
        logf.write('\n'+'='*80+'\n'+stdout_line[8:]+'\n'+'='*80+'\n\n'); logf.flush()
    else:
        logf.write(stdout_line + '\n'); logf.flush()
logf.close()

# Goodbye message
print '\n' + '-' * 80
print 'OpEx pipeline finished.'
print '-' * 80 + '\n'
