#!/usr/bin/env python
import sys
import os
import glob
import subprocess as sub
import time
import csv

"""
Re-implementation of:
https://github.com/rob123king/EMS_test_scripts/blob/master/Dragen_VCF_filtering2.pl

Helpful reading:
https://www.biostars.org/p/84951/
http://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
http://bedtools.readthedocs.io/en/latest/content/overview.html
https://stackoverflow.com/questions/5466451/how-can-i-print-literal-curly-brace-characters-in-python-string-and-also-use-fo
https://stackoverflow.com/questions/6997430/looping-over-input-fields-as-array

"""



SBATCH = 'sbatch -p {} --mem {} -c {} --wrap="{}"'
VCFMERGE = 'source vcftools-0.1.13; vcf-merge {} | grep -v \'#\' | cut -f 1,2 | awk -v OFS=\'\\t\' \'{{ print \$1,\$2-1,\$2; }}\' > {}; touch {};'
MULTICOV = 'source bedtools-2.26.0; bedtools multicov -q 1 -p -bams {} -bed {} | awk -v OFS=\'\\t\' -v mins={} \'{{ c=0; for(i = 4; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > {}; touch {};' 


BEDCOV = 'source bedtools-2.26.0; bedtools coverage -abam {} -b {} | cut -f 14 > {}; touch {};' 
COVMERGE = 'paste {} {}  | awk -v OFS=\'\\t\' -v mins={} \'{{ for(i = 3; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > {}; touch {};'
VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11 <(gzip -dc {0} | grep -v \'#\' | awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' | sort -k1,1) <(awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' | sort -k1,1) | tr " " "\\t" | cut -f 2- | gzip  >> {1}; touch {2};'
HHFILTER = 'gzip -dc {} | awk -v OFS=\'\\t\' \'{{ split($10, data, ":"); split(data[2], ad, ","); if (data[0] == "0/1" && (ad[0] > 0 || ad[1] > 0)) {{ if ( ad[1] / (ad[0] + ad[1]) > 0.1999 && ad[1] >= 5) print $0; }} else if (data[0] == "1/1" && ad[1] >= 3) print $0; }}\' | gzip > {}; touch {};'
"""

"""
def hh_filter(vcf_in, minhom=3, minhet=5, out=sys.stdout):
    """
    chr1A_part1	800103	.	G	T	41.74	.	AC=2;AF=1.000;AN=2;DP=2;FS=0.000;MQ=60.00;QD=20.87;SOR=0.693;FractionInformativeReads=1.000	GT:AD:DP:GQ:PL:SB	1/1:0,2:2:6:69,6,0:0,0,1,1
    """
    for line in vcf_in:
        data = line.strip().split('\t')[9].split(':')
        AD = list(map(int, data[1].split(',')))
        if data[0] == '0/1':
            if AD[0] == AD[1] == 0:
                pass
            else:
                if AD[1] / sum(AD) > 0.1999 and AD[1] >= minhet:
                    print(line.strip(), file=out)
        elif data[0] == '1/1' and AD[1] > minhom:
            print(line.strip(), file=out)
        



if __name__ == '__main__':

    with open(sys.argv[1]) as vcf_in:
        directories = list(line.strip() for line in vcf_in)

    vcffiles = list(os.path.join('/'.join(path.strip().split('/')), path.strip().split('/')[-1] + '.iwgsc10.vcf.nod.gz') for path in directories)
    bamfiles = list(os.path.join('/'.join(path.strip().split('/')), path.strip().split('/')[-1] + '.iwgsc10.bam') for path in directories)
    print(bamfiles, vcffiles) 
    if not os.path.exists('vcfmerge.done'): 
        cmd = VCFMERGE.format(' '.join(vcffiles), 'merged.bed', 'vcfmerge.done')
        sbatch = SBATCH.format('ei-medium', '64GB', 1, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        err, out = pr.communicate()
        tstart = time.time()
        while not os.path.exists('vcfmerge.done') and time.time() - tstart < 24 * 3600:
            time.sleep(300)

    mbedfiles = list()
    with open('merged.bed') as bed_in:
        currentChr = ''
        outf = None
        for row in csv.reader(bed_in, delimiter='\t'):
            if row[0] != currentChr:
                if outf is not None:
                   outf.close()
                mbedfiles.append('merged.bed.' + row[0])
                outf = open(mbedfiles[-1], 'w') 
                currentChr = row[0]
            print(*row, sep='\t', file=outf)
        if outf is not None:
            outf.close()

    donefiles = list()
    mincovfiles = list()
    for mbed in mbedfiles:
        donefiles.append(mbed + '.mincov.done')
        mincovfiles.append(mbed + '.mincov')
        cmd = MULTICOV.format(' '.join(bamfiles), mbed, 24, mincovfiles[-1], donefiles[-1])
        sbatch = SBATCH.format('ei-medium', '20GB', 1, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        err, out = pr.communicate()
    tstart = time.time()   

    while True:
        donefiles = list(f for f in donefiles if not os.path.exists(f))
        if not donefiles or time.time() - tstart > 24 * 3600:
            break
        time.sleep(300)

    cmd = 'cat {} > {}; touch {};'.format(' '.join(mincovfiles), 'mincov.bed', 'multicov.done')
    
    pr = sub.Popen(cmd, shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = pr.communicate()

    """ this is still needed for scaffold-based calling """
    """
    if not os.path.exists('multicov.done'):
        cmd = MULTICOV.format(' '.join(bamfiles), 'merged.bed', 24, 'mincov.bed', 'multicov.done')
        sbatch = SBATCH.format('ei-medium', '20GB', 1, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        out, err = pr.communicate()
        tstart = time.time()

        while not os.path.exists('multicov.done') and time.time() - tstart < 24 * 3600:
            time.sleep(300) 
    """

    """
    donefiles = list()
    covfiles = list()
    for bf in bamfiles:
        prefix = os.path.basename(bf)
        donefiles.append(prefix + '.done')
        covfiles.append(prefix + '.cov')
        if os.path.exists('bedcov.done'):
            continue
        cmd = BEDCOV.format(bf, 'merged.bed', covfiles[-1], donefiles[-1])
        sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        err, out = pr.communicate()
    tstart = time.time()
    
    while True:
        donefiles = list(f for f in donefiles if not os.path.exists(f))
        if not donefiles or time.time() - tstart > 24 * 3600:
            break
        time.sleep(300)
	
    pr = sub.Popen('touch bedcov.done;', shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = pr.communicate()

    if not os.path.exists('covmerge.done'):
        cmd = COVMERGE.format('merged.bed', ' '.join(covfiles), 24, 'mincov.bed', 'covmerge.done')
        sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        out, err = pr.communicate()
        tstart = time.time()

        while not os.path.exists('covmerge.done') and time.time() - tstart < 24 * 3600:
            time.sleep(300)
    """

    if not os.path.exists('vcffilter.done'):
        donefiles = list()
        fvcffiles = list()
        for vf in vcffiles:
            prefix = os.path.basename(vf)
            donefiles.append(prefix + '.done')
            fvcffiles.append(prefix.replace('.vcf.gz', '.filtered.vcf.gz'))
            cmd = VCFFILTER.format(vf, fvcffiles[-1], donefiles[-1])
            sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
            pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
            out, err = pr.communicate()
        tstart = time.time()

        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)
	
        pr = sub.Popen('touch vcffilter.done;', shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()

    if not os.path.exists('hhfilter.done'):
        donefiles = list()
        hhfiles = list()
        for fvf in fvcffiles:
            prefix = os.path.basename(fvf)
            donefiles.append(prefix + '.done')
            hhfiles.append(prefix.replace('.vcf.gz', '.hh.vcf.gz'))
            cmd = HHFILTER.format(fvf, hhfiles[-1], donefiles[-1])
            sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
            pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
            out, err = pr.communicate()
        tstart = time.time()
        
        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)

        pr = sub.Popen('touch hhfilter.done;', shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()
"""

        
"""
        
    


