#/usr/bin/env python3

import pandas as pd
import subprocess
import tempfile
import requests
import gzip
import os
import sys
import argparse

def directory_maker(directory):
    if not os.path.exists(directory):
        print('Creating output directory: {}'.format(directory))
        os.makedirs(directory, exist_ok=True)
    else:
        print('Output directory already exists: {}'.format(directory))

class genomeDownload:
    def __init__(self, args):
        if args.genome == 'hg38':
            self.refGenome = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
            self.genomeIndex = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
            self.gtRNAdb = 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz'
            self.genesGtf = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz'
            self.trnasGtf = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.tRNAs.gtf.gz'
        if args.genome == 'mm10':
            self.refGenome = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz'
            self.genomeIndex = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
            self.gtRNAdb = 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz'
            self.genesGtf = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.basic.annotation.gtf.gz'
            self.trnasGtf = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.tRNAs.gtf.gz'
    
    def run(self):
        print('Downloading genome files...')
        self.download(self.refGenome, 'genomes')
        self.download(self.genomeIndex, 'genomes', gz=False)
        print('Building Blast database...')
        subprocess.run(f'makeblastdb -in genomes/{self.refGenome.split("/")[-1].replace(".gz","")} -parse_seqids -dbtype nucl', shell=True)
        print('Downloading gtRNAdb files...')
        self.download(self.gtRNAdb, 'gtRNAdb', tar=True)
        print('Downloading GTF files...')
        self.download(self.genesGtf, 'genes')
        self.download(self.trnasGtf, 'genes')
        print('Creating Bed from GTF files...')
        subprocess.run(f'awk \'OFS="\t" {{print $1,$4-1,$5,$10,$14,$7}}\' genes/{self.genesGtf.split("/")[-1].replace(".gz","").replace(".gtf","")}.gtf | tr -d \'\\"\;\' > genes/{self.genesGtf.split("/")[-1].replace(".gz","").replace(".gtf","")}.bed', shell=True)
        subprocess.run(f'awk \'OFS="\t" {{print $1,$4-1,$5,$10,$14,$7}}\' genes/{self.trnasGtf.split("/")[-1].replace(".gz","").replace(".gtf","")}.gtf | tr -d \'\\"\;\' > genes/{self.trnasGtf.split("/")[-1].replace(".gz","").replace(".gtf","")}.bed', shell=True)

    def download(self, url, directory, gz=True, tar=False):
        r = requests.get(url)
        name = url.split('/')[-1]
        if gz == True:
            name = name.replace('.gz','')
        with open(f'{directory}/{name}','wb') as file:
            if gz == True:
                file.write(gzip.decompress(r.content))
                if tar == True:
                    tarName = name.replace('.tar','')
                    directory_maker(f'{directory}/{tarName}')
                    subprocess.run(f'tar -xf {directory}/{name} -C {directory}/{tarName}', shell=True)
                    os.remove(f'{directory}/{name}')
            else:
                file.write(r.content)

class nucleicAcidTools:
    def compliment(seq):
        return ''.join([{'A':'T','T':'A','G':'C','C':'G','U':'A','-':'-'}[i] for i in seq])

    def reverse_compliment(seq):
        return nucleicAcidTools.compliment(seq)[::-1]

    def to_rna(seq):
        return seq.replace('T','U')

    def to_dna(seq):
        return seq.replace('U','T')
    
    def is_dna(seq):
        return all(i in 'AGTC' for i in seq)
    
    def is_rna(seq):
        return all(i in 'AGUC' for i in seq)

    def fasta_to_single_line(fasta):
        with open(fasta) as file:
            for line in file:
                l = line.strip()
                if not l:
                    continue
                if l.startswith('>'):
                    if not l.endswith('_'):
                        print('Fasta file must be single line per sequence')
                        sys.exit()
                    continue
                yield l
    
class oligoGen:
    def __init__(self, args):
        self.fasta = args.fasta
        self.bed = args.bed
        self.genomeIndex = args.genomeindex
        self.output = args.output
        self.kmerRange = args.kmerrange
        self.fastaDict = {}
        self.bedDict = {}
        self.kmerDict = {}
        self.fivemerDict = {}

    def run(self):
        # Generate a fasta dictionary from the fasta file found on gtRNAdb
        print('Reading input files...')
        with open(self.fasta) as file:
            for line in file:
                l = line.strip()
                if not l:
                    continue
                if l.startswith('>'):
                    asn = l[1:].split(' ')[0].split('_')[-1]
                    if asn not in self.fastaDict:
                        self.fastaDict[asn] = ''
                    continue
                self.fastaDict[asn] += l
        # Create a list of unique tRNA sequences
        unique_trnas = pd.Series(list(self.fastaDict.values())).unique()
        # Create a bed dictionary from the bed file found on gtRNAdb
        bedDf = pd.read_csv(self.bed, delimiter='\t', header=None)
        self.bedDict = {i[1][3]:(i[1][0],i[1][1],i[1][2],i[1][5]) for i in bedDf.iterrows()}

        # Create kmer dict
        print('Creating oligo sequences...')
        for trna in unique_trnas:
            self.kmerDict = self.kmer_counter(trna,self.kmerDict,self.kmerRange)

        # Create fivemer dict for scoring
        for trna in unique_trnas:
            self.fivemerDict = self.kmer_counter(trna,self.fivemerDict,5)

        scoreDf = None
        for i in self.fastaDict.keys():
            if scoreDf is None:
                scoreDf = self.bed_score(i,bedDf)
            else:
                scoreDf = pd.concat([self.bed_score(i,bedDf),scoreDf])
        scoreDf = scoreDf.reset_index(drop=True)

        # Fix bed coordinates by moveing start and end based on direction
        scoreDf['kmer'] = [scoreDf['kmer'][i][1:] if scoreDf['direction'][i] == '+' else scoreDf['kmer'][i][:-1] for i in range(len(scoreDf))]
        scoreDf['start'] = [scoreDf['start'][i] + 1 if scoreDf['direction'][i] == '-' else scoreDf['start'][i] for i in range(len(scoreDf))]
        scoreDf['end'] = [scoreDf['end'][i] - 1 if scoreDf['direction'][i] == '+' else scoreDf['end'][i] for i in range(len(scoreDf))]

        # Compliment the kmer sequence so that the oligo will bind to the target
        scoreDf['kmer'] = [nucleicAcidTools.reverse_compliment(i) for i in scoreDf['kmer']]

        # Save the oligo sequences to a csv file
        print(f'Saving oligo sequences to: {self.output}/all_oligos.csv')
        scoreDf.to_csv(f'{self.output}/all_oligos.csv',index=None)

        # Create a bed file of the oligo sequences for genome browser
        bscoreDf = scoreDf[['chr','start','end','name']]
        print(f'Saving oligo bed file to: {self.output}/all_oligos.bed')

        with open(f'{self.output}/all_oligos.bed', 'w') as file:
            bscoreDf.to_csv(file,sep='\t',header=None,index=None)
        
        subprocess.run(['sort','-k','1,1','-k','2,2n',f'{self.output}/all_oligos.bed'],stdout=open(f'{self.output}/all_oligos.sorted.bed','w'))
        bscoreDf = pd.read_csv(f'{self.output}/all_oligos.sorted.bed', delimiter='\t', header=None)
        print(f'Saving sorted oligo bed file to: {self.output}/all_oligos.sorted.bed')

        with open(f'{self.output}/all_oligos.sorted.bed', 'w') as file:
            #file.write('track name="oligo_targets"\n')
            bscoreDf.to_csv(file,sep='\t',header=None,index=None)

        # bedToBigBed all_top25.sorted.bed ../..g/data_large/rnadb/rnadb-hg38/genome.fa.fai all_top25.sorted.bb
        if self.genomeIndex is not None:
            print(f'Creating bigBed file: {self.output}/all_oligos.sorted.bb')
            subprocess.run(['bedToBigBed',f'{self.output}/all_oligos.sorted.bed',self.genomeIndex,f'{self.output}/all_oligos.sorted.bb'])
        else:
            print('No genome index provided, skipping bigBed creation.')

    def kmer_counter(self,seq,tdict,trange):
        if type(trange) is not int:
            for klen in range(trange[0],trange[1]):
                for i in range(len(seq)+1-klen):
                    kmer = seq[i:klen+i]
                    rkmer = kmer[::-1]
                    if kmer not in tdict:
                        tdict[kmer] = 0
                    if rkmer not in tdict:
                        tdict[rkmer] = 0
                    tdict[kmer] += 1
                    if rkmer != kmer:
                        tdict[rkmer] += 1
            return tdict
        else:
            for i in range(len(seq)+1-trange):
                kmer = seq[i:trange+i]
                rkmer = kmer[::-1]
                if kmer not in tdict:
                    tdict[kmer] = 0
                if rkmer not in tdict:
                    tdict[rkmer] = 0
                tdict[kmer] += 1
                if rkmer != kmer:
                    tdict[rkmer] += 1
            return tdict
        
    def kmer_score(self,seq,score_dict):
        tscore = 0
        klen = len(list(score_dict.keys())[0])
        for i in range(len(seq)+1-klen):
            kmer = seq[i:klen+i]
            tscore += score_dict[kmer]
        return tscore/len(seq)

    def bed_score(self,target,df):
        tdict = {}
        seq = self.fastaDict[target]
        bedrange = self.bedDict[target]

        for klen in range(self.kmerRange[0],self.kmerRange[1]):
            for i in range(len(seq)+1-klen):
                kmer = seq[i:klen+i]
                if kmer not in tdict:
                    tdict[kmer] = []
                tdict[kmer].append(self.kmerDict[kmer])
                tdict[kmer].append(bedrange[0])
                tdict[kmer].append(bedrange[3])
                # Add bed range based on direction
                if bedrange[3] == '+':
                    tdict[kmer].append(bedrange[1]+i)
                    tdict[kmer].append(bedrange[1]+i+klen)
                else:
                    tdict[kmer].append(bedrange[2]-i-klen)
                    tdict[kmer].append(bedrange[2]-i)
                tdict[kmer].append(target)

        df = pd.DataFrame(tdict,index=['unique_count','chr','direction','start','end','trna']).T
        df['kmer'] = df.index
        df['kmer_length'] = [len(i)-1 for i in df['kmer']] # Subtract 1 to fix a misscount earlier and I don't want to regenerate the wiggle track cause that will alter all subsequent results
        df = df[df['unique_count']==1]
        df = df.reset_index(drop=True)

        kscore_list = []
        for i in df['kmer']:
            kscore_list.append(self.kmer_score(i,self.fivemerDict))

        df['kscore'] = kscore_list
        df = df.sort_values('kscore').reset_index(drop=True)
        df['target'] = df.index.astype('str').values
        df['name'] = df['trna'] + '_target_' + df.index.astype('str').values

        df['start'] = [df['start'][i] + 1 for i in range(len(df))]

        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)

        df = df.drop(['unique_count'],axis=1)

        return df
    
class oligoValidate:
    def __init__(self, args):
        self.fasta = args.fasta
        self.oligoDf = pd.read_csv(args.oligos,header=0)
        self.targets = args.targets
        self.output = args.output
        self.idtConfig = args.idtconfig
        self.blastdb = args.blastdb
        self.bedgtf = args.bedgtf

    def run(self):
        if os.path.exists('.'.join(self.fasta.split('.')[:-1]) + '.singleline.fa'):
            self.fasta = '.'.join(self.fasta.split('.')[:-1]) + '.singleline.fa'
        else:
            print('.'.join(self.fasta.split('.')[:-1]) + '.singleline.fa does not exist, creating now.')
            with open(self.fasta) as file:
                with open('.'.join(self.fasta.split('.')[:-1]) + '.singleline.fa','w') as outfile:
                    # if first line dont add newline
                    outfile.write(file.readline())
                    for l in file.readlines():
                        if l[0] == '>':
                            outfile.write('\n')
                            outfile.write(l)
                        else:
                            outfile.write(l.strip())
            self.fasta = '.'.join(self.fasta.split('.')[:-1]) + '.singleline.fa'

        with open(self.targets) as file:
            with open(f'{self.output}/combine_IDT.csv','w') as combineFile:
                for l in file.readlines():
                    oligo = l.strip().split()
                    self.validate(oligo,combineFile)

    def validate(self,oligo,combineFile):
        oligo_name = oligo[0]
        if len(self.oligoDf[self.oligoDf['name']==oligo[1]]) > 0:
            print(f'{oligo[0]} found by oligo target name: {oligo[1]}')
            oligo_seq = self.oligoDf[self.oligoDf['name']==oligo[1]]['kmer'].values[0]
        elif len(self.oligoDf[self.oligoDf['kmer']==oligo[1]]) > 0:
            print(f'{oligo[0]} found by oligo target sequence: {oligo[1]}')
            oligo_seq = oligo[1]
        else:
            oligo[1] = oligo[1].upper()
            if nucleicAcidTools.is_dna(oligo[1]):
                print(f'{oligo[0]} is not found in oligo target list by name or seuqence but is a valid DNA sequence that will be used.')
                oligo_seq = oligo[1]
            elif nucleicAcidTools.is_rna(oligo[1]):
                print(f'{oligo[0]} is not found in oligo target list by name or seuqence but is a valid RNA sequence that will be converted to DNA and used.')
                oligo_seq = nucleicAcidTools.to_dna(oligo[1])
            else:
                print(f'{oligo[0]} is not found in oligo target list by name or seuqence and is not a valid DNA or RNA sequence and will be skipped.')
                return

        if self.alignment_check(oligo_name,oligo_seq):
            print(f'{oligo_name} is a valid oligo and will be added to the combine_IDT.csv file.')
            combineFile.write(f'{oligo_name},{self.lna(oligo_seq)}{self.idtConfig}\n')

    def lna(self,seq):
        return ''.join(['+' + seq[i] if i%2 == 0 else seq[i] for i in range(len(seq))])

    def smith_waterman(self, seq1, seq2, match_score=1, mismatch_score=0, gap_penalty=-1):
        if len(seq1) > len(seq2):
            seq1, seq2 = seq2, seq1
        # Initialize the scoring matrix
        n = len(seq1)
        m = len(seq2)
        score_matrix = [[0 for j in range(m+1)] for i in range(n+1)]
        max_score = 0
        max_pos = (0,0)
        for i in range(1, n+1):
            for j in range(1, m+1):
                match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
                delete = score_matrix[i-1][j] + gap_penalty
                insert = score_matrix[i][j-1] + gap_penalty
                score_matrix[i][j] = max(0, match, delete, insert)
                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)
        
        # Traceback to find the optimal alignment
        aligned_seq1 = ""
        aligned_seq2 = ""
        alignment = ""
        i = max_pos[0]
        j = max_pos[1]
        while i > 0 or j > 0:
            if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                if seq1[i-1] == seq2[j-1]:
                    alignment = "|" + alignment
                else:
                    alignment = " " + alignment
                i -= 1
                j -= 1
            elif i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap_penalty:
                aligned_seq1 = seq1[i-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                alignment = " " + alignment
                i -= 1
            else:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[j-1] + aligned_seq2
                alignment = " " + alignment
                j -= 1
        
        aligned_seq2_end = aligned_seq2 + seq2[max_pos[1]:]
        aligned_seq1_end = aligned_seq1 + "-" * len(seq2[max_pos[1]:])
        aligned_row = ""
        for i in range(len(aligned_seq1_end)):
            if aligned_seq1_end[i] == aligned_seq2_end[i]:
                aligned_row += "|"
            else:
                aligned_row += " "

        return [aligned_seq2_end,aligned_row,nucleicAcidTools.compliment(aligned_seq1_end)]

    def alignment_check(self, name, seq):
        dashed_l = '\n' + '-'*120 + '\n'
        try:
            out = subprocess.check_output(['grep', nucleicAcidTools.reverse_compliment(seq), self.fasta, '-B', '1'], text=True).strip()
        except:
            print(f'{name} not found in {self.fasta}')
            return None
        
        name_out = f'{self.output}/{name}.txt'
        with open(name_out,'w') as f:
            f.write(name+'\noligo: '+seq+'\nlength: '+str(len(seq))+'\n')
            for i in out.split('\n'):
                if i[0] == '>': 
                    f.write(i+'\n')
            for i in self.smith_waterman(out.split('\n')[1],nucleicAcidTools.reverse_compliment(seq)):
                f.write(i+'\n')

            if self.blastdb:
                temp = tempfile.NamedTemporaryFile(mode='w', delete=False)
                temp.write(f'>{name}\n'+nucleicAcidTools.reverse_compliment(seq))
                temp.close()
                f.write(dashed_l)
                blastout = subprocess.check_output(['blastn', '-db', self.blastdb, '-query', temp.name, '-task', 'blastn-short', \
                                                    '-outfmt', '6 qseqid sseqid pident length mismatches gaps qstart qend sstart send evalue bitscore sstrand', \
                                                    '-perc_identity', '95', '-qcov_hsp_perc', '95'], text=True).strip()
                if blastout:
                    f.write('\nBLAST Results:\n')
                    f.write('qseqid\tsseqid\tpident\tlength\tmismatches\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsstrand\n')
                    for i in blastout.split('\n'):
                        f.write(i+'\n')
                    if self.bedgtf:
                        f.write(dashed_l)
                        # Convert blastout to bed format
                        blastoutbed = ''
                        for i in blastout.split('\n'):
                            if i:
                                i = i.split('\t')
                                if i[11] == 'plus':
                                    blastoutbed += f'{i[1]}\t{i[7]}\t{i[8]}\t{i[0]}\t+\n'
                                else:
                                    blastoutbed += f'{i[1]}\t{i[8]}\t{i[7]}\t{i[0]}\t-\n'
                                
                        temp = tempfile.NamedTemporaryFile(mode='w', delete=False)      
                        temp.write(blastoutbed)
                        temp.close()

                        for i in self.bedgtf:
                            f.write('\n'+i+' Results:\n')
                            bedout = subprocess.check_output(['bedtools', 'intersect', '-a', i, '-b', temp.name, '-wa'], text=True).strip()
                            if bedout:
                                f.write('chr\tstart\tend\tname\tstrand\n')
                                j = ''
                                sublist = []
                                for i in bedout.split('\n'):
                                    if i != j:
                                        bsize = int(i.split('\t')[2]) - int(i.split('\t')[1])
                                        if bsize < 500:
                                            sublist.append((i,bsize))
                                        else:
                                            f.write(i+'\n')
                                    j = i
                                if len(sublist) > 0:
                                    f.write('\n!WARNING! Sub 500bp:\n')
                                    f.write('chr\tstart\tend\tname\tstrand\n')
                                    for i,bsize in sublist:
                                        f.write(i+'\t'+str(bsize)+'\n')
                            else:
                                f.write('No BED hits found\n')
                else:
                    f.write('\nBLAST Results:\nNo BLAST hits found\n')
        return out

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='oligoTools.py',
        description='Generate oligos for tRNA targets')
    
    subparsers = parser.add_subparsers(
        title='Operating modes',
        description='Choose between generating oligos or validating oligos',
        dest='mode',
        required=True,
    )

    parser_download = subparsers.add_parser("download", help="Download required files for oligoTools for hg38 or mm10")
    parser_download.add_argument('-g', '--genome', help='Genome to download (hg38 or mm10)', required=True, choices=['hg38', 'mm10'])
    parser_download.add_argument('--log', help='Log output to file (optional)', default=None)

    parser_gen = subparsers.add_parser("generate", help="Generate oligos for tRNA targets")
    parser_gen.add_argument('-f', '--fasta', help='Fasta file of tRNA sequences', required=True)
    parser_gen.add_argument('-b', '--bed', help='Bed file of tRNA coordinates', required=True)
    parser_gen.add_argument('-g', '--genomeindex', help='Genome fasta index genome.fa.fai (optional)', required=False)
    parser_gen.add_argument('-k', '--kmerrange', help='Range of kmer lengths to use (Default: [20, 30])', nargs=2, default=[20, 30])
    parser_gen.add_argument('-o', '--output', help='Output directory (Default: oligos)', default='oligos')
    parser_gen.add_argument('--log', help='Log output to file (optional)', default=None)

    parser_val = subparsers.add_parser("analyze", help="Validate oligos for tRNA targets")
    parser_val.add_argument('-f', '--fasta', help='Fasta file of tRNA sequences', required=True)
    parser_val.add_argument('-l', '--oligos', help='Oligo targets file that matches generated bw and bb files (all_oligos.csv)', required=True)
    parser_val.add_argument('-t', '--targets', help='Oligo targets or oligo sequences to validate in .tsv format', required=True)
    parser_val.add_argument('-o', '--output', help='Output directory (Default: targets)', default='targets')
    parser_val.add_argument('-d', '--idtconfig', help='IDT config file default: "/3Bio/,100nm"', default="/3Bio/,100nm")
    parser_val.add_argument('-b', '--blastdb', help='Blast database to use (optional)', default=None)
    parser_val.add_argument('-g', '--bedgtf', help='Bed file/s that are converted from GTF format (optional)', nargs='+', default=None)
    parser_val.add_argument('--log', help='Log output to file (optional)', default=None)

    args = parser.parse_args()

    sys.stdout = open(args.log, 'w') if args.log else sys.stdout

    if args.mode == 'download':
        directory_maker('genes')
        directory_maker('genomes')
        directory_maker('gtRNAdb')
        genomeDownload(args).run()

    if args.mode == 'generate':
        directory_maker(args.output) 
        oligoGen(args).run()

    if args.mode == 'analyze':
        directory_maker(args.output)
        oligoValidate(args).run()