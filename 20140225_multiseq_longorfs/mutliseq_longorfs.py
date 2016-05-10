""" small pipeline generating an fasta file with long orfs using glimmer from a fasta files containing multiple contigs

Usage:
  multiseq_longorfs.py -o <outfile> -i <infile>

Options:
   -h --help       Show this screen.
   -i <infile>     input fasta file
   -o <outfile>    output fasta file
"""
from docopt import docopt
import os
import sh
from gefes.fasta.single import FASTA
from sh import ErrorReturnCode

updist=25-1

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)

seqs=FASTA(arguments['-i'])

motif_count = sh.Command('get-motif-counts.awk')
sh.touch('temp.genes.fasta')

for s in seqs:
    temp=FASTA('./temp.fasta')
    s.description=s.name
    temp.write(s)
    try:
        sh.long_orfs('-l','-g', 30 ,'-n', '-t', 1.15, './temp.fasta', './temp.longorfs')
        sh.extract('-t', './temp.fasta', './temp.longorfs', _out="temp2.fasta")
        sh.sed( '-n', 'wtemp.fasta', 'temp.genes.fasta', 'temp2.fasta')
        sh.mv('temp.fasta', 'temp.genes.fasta')
    except ErrorReturnCode:
        sh.rm('temp.fasta')
        print "kicked one out"
    


seqs.close()
sh.rm('temp2.fasta')
sh.rm('temp.longorfs')

with open('temp.genes.fasta') as inpo:
    sh.build_icm('-r','temp.icm',_in=inpo)
    
sh.glimmer3('-l','-o50', '-g110', '-t30', arguments['-i'], 'temp.icm', 'temp')

with open('temp.predict') as f_in:
    content = f_in.readlines()
    with open('temp.gene.coord','w') as g_out:
        with open('temp.upstream.coord','w') as up_out:
            for line in content:
                if '>' in line:
                    iD=line[1:-1].split(" ")[0]
                else:
                    line=[w for w in line.split(" ") if w!='']
                    line_gene= "\t".join([line[0],iD,"\t".join(line[1:])])
                    g_out.write(line_gene)
                    pos = int(line[3]) > 0
                    start=int(line[1])
                    end=int(line[2])
                    if start > end:
                        start = str(int(start)+1)
                        end = str(int(start)+updist)
                    else:
                        start = str(int(start)-1)
                        end = str(int(start)-updist)
                    line[1]=start
                    line[2]=end
                    line= "\t".join([line[0],iD,"\t".join(line[1:3])])+"\n"
                    up_out.write(line)
                
sh.multi_extract(arguments['-i'], 'temp.upstream.coord', _out='temp.upstream')
motif_count(sh.elph('temp.upstream', 'LEN=6'), _out='temp.motif')
#add start codon use
sh.glimmer3('-l', '-o50', '-g110', '-t30','-b','temp.motif', arguments['-i'], 'temp.icm', 'final')

with open('final.predict') as f_in:
    content = f_in.readlines()
    with open('gene.coord','w') as g_out:
        for line in content:
            if '>' in line:
                iD=line[1:-1].split(" ")[0]
            else:
                line=[w for w in line.split(" ") if w!='']
                line_gene= "\t".join([line[0],iD,"\t".join(line[1:])])
                g_out.write(line_gene)
sh.multi_extract('-d', '-t',arguments['-i'], 'gene.coord', _out= 'temp.genes')

in_file=FASTA('./temp.genes')
sh.touch('./'+arguments['-o'])
with FASTA('./'+arguments['-o']) as out_file:
    for seq in in_file:
        seq.description=seq.description.replace("  ",";",1)
        seq.name=seq.description.split("  ")[0]
        seq.id=seq.name
        out_file.add_seq(seq)

        
sh.rm('gene.coord')
sh.rm('temp.motif')
sh.rm('temp.genes')
sh.rm('temp.genes.fasta')
sh.rm('temp.icm')
sh.rm('temp.predict')
sh.rm('temp.gene.coord')
sh.rm('temp.upstream.coord')
sh.rm('temp.upstream')
sh.rm('final.predict')
sh.rm("final.detail")
sh.rm("temp.detail")
#c(228,300,304,268,360,300,264,140)
# blastall -p tblastx -d /home/moritz/glob/data/blast.dbs/OD1_db.fasta -i bin_0.fasta.genes -o test.blast -a 16 -K1 -m8 -P1


