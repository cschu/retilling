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
https://stackoverflow.com/questions/32481877/what-is-nr-fnr-in-awk
"""

SBATCH = 'sbatch -p {} --mem {} -c {} --wrap="{}"'
SBATCH_ARRAY = 'sbatch -p {} --mem {} -c {} --wrap="{}" --array=1-{}'

VCFMERGE = 'source vcftools-0.1.13; vcf-merge {} | grep -v \'#\' | cut -f 1,2 | awk -v OFS=\'\\t\' \'{{ print \$1,\$2-1,\$2; }}\' > {}; touch {};'
MULTICOV = 'source bedtools-2.26.0; bedtools multicov -q 1 -p -bams {} -bed {} | awk -v OFS=\'\\t\' -v mins={} \'{{ c=0; for(i = 4; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > {}; touch {};' 
MULTICOV_ARRAY = 'source bedtools-2.26.0; bed=\$(head -n \$SLURM_ARRAY_TASK_ID {} | tail -n 1); bedtools multicov -q 1 -p -bams {} -bed \$bed | awk -v OFS=\'\\t\' -v mins={} \'{{ c=0; for(i = 4; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > \$bed.mincov; touch \$bed.done;'


#BEDCOV = 'source bedtools-2.26.0; bedtools coverage -abam {} -b {} | cut -f 14 > {}; touch {};' 
#COVMERGE = 'paste {} {}  | awk -v OFS=\'\\t\' -v mins={} \'{{ for(i = 3; i <= NF; i++) {{ if (\$i > 0) c++; }}; if (c >= mins) print \$0; }}\' > {}; touch {};'
#VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11 <(gzip -dc {0} | grep -v \'#\' | awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' | sort -k1,1) <(awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' | sort -k1,1) | tr " " "\\t" | cut -f 2- | gzip  >> {1}; touch {2};'
#VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11 <(gzip -dc {0} | grep -v \'#\' | awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' | sort -k1,1) <(awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' {3} | sort -k1,1) | tr " " "\\t" | gzip  >> {1}; touch {2};'
VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; gzip -dc {0} | grep -v \'#\' | awk -F\'\\t\' \'NR==FNR{{c[\$1\$3]++;next}};c[\$1\$2] > 0\' {3} - | gzip >> {1}; touch {2};'

HHFILTER = 'gzip -dc {} | awk -v OFS=\'\\t\' \'{{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= 5) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= 3) print \$0; }}\' | gzip > {}; touch {};'

DIFILTER = 'source vcftools-0.1.13; gzip -dc {0} | awk -F\'\\t\' \'/^#/{{ print \$0; next; }}; NR==FNR{{c[\$1\$3]++; next; }}; c[\$1\$2] > 0\' {1} - | awk -v OFS=\'\\t\' \'/^#/{{print \$0; next;}} {{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= 5) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= 3) print \$0; }}\' | bgzip -c > {2} && tabix -p vcf {2}; touch {3};'

INTERSECT = 'source vcftools-0.1.13; vcf-isec -f -c {} {} > {}; touch {};'

#vcfout=$(dirname $vcfin)/$(basename $vcfin .gz).nod.gz
#gzip -dc $vcfin | grep -v '##DRAGEN' | bgzip -c > $vcfout && tabix -p vcf $vcfout;

#system("vcf-isec -f -c $data2[$string0] $data[$string1] $data[$string2] $data[$string3] $data[$string4] $data[$string5] $data[$string6] $data[$string7] $data[$string8] $data[$string9] $data[$string10] $data[$string11] $data[$string12] $data[$string13] $data[$string14] $data[$string14] $data[$string15] $data[$string16] $data[$string17] $data[$string18] $data[$string19] $data[$string20] $data[$string21] $data[$string22] $data[$string23] $data[$string24] $data[$string25] $data[$string26] $data[$string27] /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1116_LIB10020_LDI8275.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1128_LIB10021_LDI8276.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1128_LIB10022_LDI8277.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1132_LIB10313_LDI8531.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1132_LIB10314_LDI8532.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1132_LIB10315_LDI8533.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1132_LIB10316_LDI8534.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1133_LIB10317_LDI8535.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1133_LIB10318_LDI8536.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1134_LIB10495_LDI8690.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1147_LIB10700_LDI8902.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1147_LIB10701_LDI8903.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1149_LIB10703_LDI8905.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1158_LIB11007_LDI9084.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1162_LIB11008_LDI9085.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1162_LIB11009_LDI9086.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1170_LIB11492_LDI9431.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1170_LIB11493_LDI9432.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1171_LIB11494_LDI9433.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1171_LIB11495_LDI9434.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1176_LIB11490_LDI9429.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1176_LIB11491_LDI9430.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1238_LIB10702_LDI8904.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1238_LIB10704_LDI8906.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1238_LIB10705_LDI8907.vcf.2.gz /home/data/wheat_ems_mutation/Dragen_cadenza_controls/Sample_1238_LIB10706_LDI8908.vcf.2.gz > $data[$string0].filtered.unique.vcf\n") == 0  or die "system vcf isec filter job failed: $?"; 


"""

"""
def splitInputFile(inputfile):
    nlines = sum(1 for line in open(inputfile))        
    with open(inputfile) as _in:
        part = -1
        outf = None
        for i, line in enumerate(_in):
            current = i // 20000
            if current != part:
                part = current
                if outf is not None:
                    outf.close()
                f = inputfile + '.part{}'.format(current)
                yield f
                outf = open(f, 'w')
            print(line.strip(), file=outf)
        if outf is not None:
            outf.close()
    """
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
    """
    pass

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

    workdir = 'test'
    VCFMERGE_DONE = os.path.join(workdir, 'vcfmerge.done')
    VCFMERGE_OUTPUT = os.path.join(workdir, 'merged.bed')
    MULTICOV_DONE = os.path.join(workdir, 'multicov.done')
    MULTICOV_OUTPUT = os.path.join(workdir, 'mincov.bed')
    MULTICOV_INPUT = os.path.join(workdir, 'split_mbedfiles_in')
    VCFFILTER_DONE = os.path.join(workdir, 'vcffilter.done')
    HHFILTER_DONE = os.path.join(workdir, 'hhfilter.done')
    DIFILTER_DONE = os.path.join(workdir, 'difilter.done')
    INTERSECT_DONE = os.path.join(workdir, 'intersect.done')
    ALL_DONE = os.path.join(workdir, 'all.done')
    FINAL_OUTPUT = 'unique_snps.vcf'

    with open(sys.argv[1]) as vcf_in:
        directories = list(line.strip() for line in vcf_in)

    vcffiles = list(os.path.join('/'.join(path.strip().split('/')), path.strip().split('/')[-1] + '.iwgsc10.vcf.nod.gz') for path in directories)
    bamfiles = list(os.path.join('/'.join(path.strip().split('/')), path.strip().split('/')[-1] + '.iwgsc10.bam') for path in directories)
    print(bamfiles, vcffiles) 
    if not os.path.exists(VCFMERGE_DONE): 
        cmd = VCFMERGE.format(' '.join(vcffiles), VCFMERGE_OUTPUT, VCFMERGE_DONE)
        sbatch = SBATCH.format('ei-medium', '64GB', 1, cmd)
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        err, out = pr.communicate()
        tstart = time.time()
        while not os.path.exists(VCFMERGE_DONE) and time.time() - tstart < 24 * 3600:
            time.sleep(300)


    if not os.path.exists(MULTICOV_DONE):
        mbedfiles = list(splitInputFile(VCFMERGE_OUTPUT))
        with open(MULTICOV_INPUT, 'w') as mbedfiles_f:
            print(*mbedfiles, sep='\n', file=mbedfiles_f)
        donefiles = list(f + '.done' for f in mbedfiles)
        mincovfiles = list(f + '.mincov' for f in mbedfiles)
        cmd = MULTICOV_ARRAY.format(MULTICOV_INPUT, ' '.join(bamfiles), 24)
        sbatch = SBATCH_ARRAY.format('ei-medium', '8GB', 1, cmd, len(mbedfiles))
        pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
        err, out = pr.communicate()
        tstart = time.time()
    
        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)
    
        cmd = 'cat {} > {}; touch {};'.format(' '.join(mincovfiles), MULTICOV_OUTPUT, MULTICOV_DONE)
        
        pr = sub.Popen(cmd, shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()

    """
VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11 <(gzip -dc {0} | grep -v \'#\' | awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' | sort -k1,1) <(awk -v OFS=\'\\t\' \'{{ print \$1":"\$2,\$0; }}\' {3} | sort -k1,1) | tr " " "\\t" | gzip  >> {1}; touch {2};'

VCFFILTER = 'gzip -dc {0} | grep \'^#\' | gzip > {1}; gzip -dc {0} | grep -v \'#\' | awk -F'\t' 'NR==FNR{c[\$1\$3]++;next};c[\$1\$2] > 0' {3} - | gzip >> {1}; touch {2};'
< tmp2.tsv awk -F '\t' 'FNR==NR{ a[$1] = $2; next }{ print $1 FS a[$1] }' tmp1.tsv -

DIFILTER = 'gzip -dc {} | awk -F\'\\t\' \'/^#/{{ print \$0; next; }}; NR==FNR{{c[\$1\$3]++; next; }}; c[\$1\$2] > 0\' {} - | awk -v OFS=\'\\t\' \'/^#/{{print \$0; next;}} {{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= 5) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= 3) print \$0; }}\' > {}; touch {};'

DIFILTER = 'source vcftools-0.1.13; gzip -dc {0} | awk -F\'\\t\' \'/^#/{{ print \$0; next; }}; NR==FNR{{c[\$1\$3]++; next; }}; c[\$1\$2] > 0\' {1} - | awk -v OFS=\'\\t\' \'/^#/{{print \$0; next;}} {{ split(\$10, data, \\":\\"); split(data[2], ad, \\",\\"); if (data[1] == \\"0/1\\" && (ad[1] > 0 || ad[2] > 0)) {{ if ( ad[2] / (ad[1] + ad[2]) > 0.1999 && ad[2] >= 5) print \$0; }} else if (data[1] == \\"1/1\\" && ad[2] >= 3) print \$0; }}\' | bgzip -c > {2} && tabix -p vcf {2}; touch {3};'

    """
    if not os.path.exists(DIFILTER_DONE):
        donefiles, fvcffiles = list(), list()
        for vf in vcffiles:
            prefix = os.path.join(workdir, os.path.basename(vf))
            donefiles.append(prefix + '.done')
            fvcffiles.append(prefix.replace('.vcf.nod.gz', '.filtered.vcf.gz'))
            cmd = DIFILTER.format(vf, MULTICOV_OUTPUT, fvcffiles[-1], donefiles[-1])
            sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
            pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
            out, err = pr.communicate()
        tstart = time.time()
        print('FVCFFILES', fvcffiles)
        
        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)
	
        pr = sub.Popen('touch {};'.format(DIFILTER_DONE), shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()

    else:
        fvcffiles = list(os.path.join(workdir, os.path.basename(vf)).replace('.vcf.nod.gz', '.filtered.vcf.gz') for vf in vcffiles) 
        print('FVCFFILES', fvcffiles)

    """
    INTERSECT = 'source vcftools-0.1.13; vcf-isec -f -c {} {} > {}; touch {};'
    """
    if not os.path.exists(INTERSECT_DONE):
        donefiles, uniqfiles = list(), list()
        isect_sets = set(fvcffiles)
        for vf in fvcffiles:
            donefiles.append(vf + '.done')
            uniqfiles.append(vf.replace('.filtered.vcf.gz', '.uniq.vcf'))
            cmd = INTERSECT.format(vf, ' '.join(isect_sets.difference([vf])), uniqfiles[-1], donefiles[-1])
            sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
            pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
            out, err = pr.communicate()
        tstart = time.time()

        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)

        pr = sub.Popen('touch {};'.format(INTERSECT_DONE), shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()
    else:
        uniqfiles = list(vf.replace('.filtered.vcf.gz', '.uniq.vcf') for vf in fvcffiles)

    vcfheader = []
    snps = list()
    for uf in uniqfiles: 
        sampleID = os.path.basename(uf).split('.')[0]
        with open(uf) as _in:
            for row in csv.reader(_in, delimiter='\t'):
                if not row[0].startswith('#'):
                    row[1], row[2] = int(row[1]), sampleID
                    snps.append(row)
    with open(FINAL_OUTPUT, 'w') as final_out:
        for line in open(uniqfiles[0]):
            if line.startswith('#') and not line.startswith('##reference') and not line.startswith('##source'):
                print(line.strip(), file=final_out)
            else:
                break
        for snp in sorted(snps, key=lambda x:(x[0], x[1])):
            print(*snp, sep='\t', file=final_out)

    pr = sub.Popen('touch {};'.format(ALL_DONE), shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
    out, err = pr.communicate() 

    sys.exit(1)

    """
    if not os.path.exists(VCFFILTER_DONE):
        donefiles = list()
        fvcffiles = list()
        for vf in vcffiles:
            prefix = os.path.join(workdir, os.path.basename(vf))
            donefiles.append(prefix + '.done')
            fvcffiles.append(prefix.replace('.vcf.nod.gz', '.filtered.vcf.nod.gz'))
            cmd = VCFFILTER.format(vf, fvcffiles[-1], donefiles[-1], MULTICOV_OUTPUT)
            sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
            pr = sub.Popen(sbatch, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE, shell=True)
            out, err = pr.communicate()
        tstart = time.time()
        print('FVCFFILES', fvcffiles)
        
        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)
	
        pr = sub.Popen('touch {};'.format(VCFFILTER_DONE), shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()

    else:
        fvcffiles = list(os.path.join(workdir, os.path.basename(vf)).replace('.vcf.nod.gz', '.filtered.vcf.nod.gz') for vf in vcffiles) 
        print('FVCFFILES', fvcffiles)

    print(fvcffiles)
    if not os.path.exists(HHFILTER_DONE):
        donefiles = list()
        hhfiles = list()
        for fvf in fvcffiles:
            prefix = os.path.join(workdir, os.path.basename(fvf))
            donefiles.append(prefix + '.done')
            hhfiles.append(prefix.replace('.vcf.nod.gz', '.hh.vcf.nod.gz'))
            cmd = HHFILTER.format(fvf, hhfiles[-1], donefiles[-1])
            sbatch = SBATCH.format('ei-medium', '8GB', 1, cmd)
            print('HH-SBATCH:', sbatch)
            pr = sub.Popen(sbatch, shell=True, stdout=sub.PIPE, stderr=sub.PIPE, stdin=sub.PIPE)
            out, err = pr.communicate()
        tstart = time.time()
        
        while True:
            donefiles = list(f for f in donefiles if not os.path.exists(f))
            if not donefiles or time.time() - tstart > 24 * 3600:
                break
            time.sleep(300)

        pr = sub.Popen('touch {};'.format(HHFILTER_DONE), shell=True, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.PIPE)
        out, err = pr.communicate()
    """
"""

        
"""
        
    



