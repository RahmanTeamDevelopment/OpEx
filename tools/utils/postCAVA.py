from __future__ import division
import sys
from subprocess import call

fn = sys.argv[1]
prefn = fn[:-4] + '_pre.vcf'
call(['mv', fn, prefn])
out = open(fn, 'w')
for line in open(prefn):
    if line.startswith('#CHROM'):
        cols = line.split('\t')
        x = cols[-1]
        cols[-1] = x[:x.find('_picard.bam')]
        out.write('\t'.join(cols) + '\n')
    else:
        out.write(line)
out.close()
call(['rm', prefn])

header = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'QUALFLAG', 'FILTER', 'TR', 'TC', 'SAMPLE', 'GT', 'TYPE', 'TRANSCRIPT',
          'GENE', 'GENEID', 'TRINFO', 'LOC', 'CSN', 'PROTPOS', 'PROTREF', 'PROTALT', 'CLASS', 'SO', 'IMPACT', 'ALTANN', 'ALTCLASS', 'ALTSO']
print '#' + '\t'.join(header)

for line in open(fn):
    line = line.strip()
    if line == '': continue
    if line.startswith('##'): continue

    cols = line.split('\t')

    if line.startswith('#'):
        sample = cols[9]
        continue

    chrom = cols[0]
    if chrom.startswith('chr'): chrom = chrom[3:]
    pos = cols[1]
    id = cols[2]
    ref = cols[3]
    alts = cols[4].split(",")
    qual = cols[5]
    filter = cols[6]
    info = cols[7]

    infobits = info.split(';')
    infodict = dict()
    for infobit in infobits:
        idx = infobit.find('=')
        if idx != -1:
            key = infobit[:idx].strip()
            value = infobit[idx + 1:].strip()
            infodict[key] = value

    TRANSCRIPT_byalt = infodict['TRANSCRIPT'].split(',')
    TYPE_byalt = infodict['TYPE'].split(',')
    GENE_byalt = infodict['GENE'].split(',')
    GENEID_byalt = infodict['GENEID'].split(',')
    TRINFO_byalt = infodict['TRINFO'].split(',')
    LOC_byalt = infodict['LOC'].split(',')
    CSN_byalt = infodict['CSN'].split(',')
    PROTPOS_byalt = infodict['PROTPOS'].split(',')
    PROTREF_byalt = infodict['PROTREF'].split(',')
    PROTALT_byalt = infodict['PROTALT'].split(',')
    CLASS_byalt = infodict['CLASS'].split(',')
    SO_byalt = infodict['SO'].split(',')
    IMPACT_byalt = infodict['IMPACT'].split(',')
    ALTANN_byalt = infodict['ALTANN'].split(',')
    ALTCLASS_byalt = infodict['ALTCLASS'].split(',')
    ALTSO_byalt = infodict['ALTSO'].split(',')

    TRs = infodict['TR'].split(',')
    TC = infodict['TC']

    if float(TC) == 0: continue

    GT = cols[9][:3]

    for i in range(len(alts)):
        alt = alts[i]
        transcripts = TRANSCRIPT_byalt[i].split(':')
        GENE_bytrans = GENE_byalt[i].split(':')
        GENEID_bytrans = GENEID_byalt[i].split(':')
        TRINFO_bytrans = TRINFO_byalt[i].split(':')
        LOC_bytrans = LOC_byalt[i].split(':')
        CSN_bytrans = CSN_byalt[i].split(':')
        PROTPOS_bytrans = PROTPOS_byalt[i].split(':')
        PROTREF_bytrans = PROTREF_byalt[i].split(':')
        PROTALT_bytrans = PROTALT_byalt[i].split(':')
        CLASS_bytrans = CLASS_byalt[i].split(':')
        SO_bytrans = SO_byalt[i].split(':')
        IMPACT_bytrans = IMPACT_byalt[i].split(':')
        ALTANN_bytrans = ALTANN_byalt[i].split(':')
        ALTCLASS_bytrans = ALTCLASS_byalt[i].split(':')
        ALTSO_bytrans = ALTSO_byalt[i].split(':')

        qualflag = ''
        if TYPE_byalt[i] == 'Substitution':
            if float(qual) >= 100: qualflag = 'high'
            else: qualflag = 'low'
        else:
            prop = float(TRs[i]) / float(TC)
            if prop > 0.2 and filter == 'PASS': qualflag = 'high'
            else: qualflag = 'low'

        for j in range(len(transcripts)):

            if CSN_bytrans[j] == '.': continue

            out = [chrom, pos, ref, alt, qual, qualflag, filter, TRs[i], TC, sample, GT, TYPE_byalt[i], transcripts[j],
                   GENE_bytrans[j], GENEID_bytrans[j], TRINFO_bytrans[j], LOC_bytrans[j], CSN_bytrans[j], PROTPOS_bytrans[j],
                   PROTREF_bytrans[j], PROTALT_bytrans[j], CLASS_bytrans[j], SO_bytrans[j], IMPACT_bytrans[j], ALTANN_bytrans[j],
                   ALTCLASS_bytrans[j], ALTSO_bytrans[j]]

            print '\t'.join(out)
