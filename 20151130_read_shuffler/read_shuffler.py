#!/home/moritz/.pyenv/shims/python -u
#SBATCH -D /home/moritz/repos/Pyscratches/20151130_read_shuffler/
#SBATCH -J gambit
#SBATCH -o /home/moritz/repos/Pyscratches/20151130_read_shuffler/gambit_2.out
#SBATCH -e /home/moritz/repos/Pyscratches/20151130_read_shuffler/gambit_2.err
#SBATCH -A b2013127
#SBATCH -t 7-00:00:00
#SBATCH -n 1
#SBATCH -p core
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

""" shuffles a bunch of libraries

Usage:
  read_shuffler.py [-g -n <bps>] -i <infastq> -o <outhead>
  read_shuffler.py -h

Options:
    -o <outfastq>   output fastq
    -i <infastq>    input folder in which to find fastq
    -g              input is gzipped
    -n <bps>        number of basepairs to trim if you want trimming
    -h              this help
"""

from docopt import docopt
import sh
import os
from random import randint
from tqdm import tqdm
import sys
import gzip
import os
from os.path import join


class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)




if __name__ == '__main__':
    sys.stdout = Unbuffered(sys.stdout)
    arguments = docopt(__doc__)
    input_folder = arguments['-i']
    output_head = arguments['-o']
    bps_fwd = int(arguments['-n']) if arguments['-n'] else 0
    bps_rev = bps_fwd + 1
    file_opener = gzip.open if arguments['-g'] else open
    lib_list_fwd = [f[:-1] for f in sh.find(input_folder) if "1P_fastq.gz" in f ]
    lib_list_rev = [f.replace("1P","2P") for f in lib_list_rev]

#print lib_list
handles = [open(l) for l in lib_list]
out_file = output_head + "_%s.fasta"
entry = ""
ll = len(handles)-1
i=0
handle = None
fileID = 0
min_seq_len = 40
g_counter = 0
trash_count = 0

while ll > -1:
    if not handle:
        fileID += 1
        handle = open(out_file % str(fileID).zfill(5), "w")

    which = randint(0,ll)
    tt=handles[which]
    temp =  next(tt, None)
    if temp:
        temp2 = next(tt, None)
        if len(temp2) > min_seq_len:
            entry += temp
            entry += temp2[bps_fwd:-bps_rev] + "\n"
            entry += next(tt, None)
            entry += next(tt, None)[bps_fwd:-bps_rev] + "\n"
            i += 1
        else :
            trash = next(tt, None)
            trash = next(tt, None)
            trash_count += 1

        if len(entry) > 1000000000:
            handle.write(entry)
            entry = ""
            print i, "reads", "of which ", str(float(trash_count) / i) + "% dirties"
            g_counter += 1
            if g_counter > 100:
                print "Closing ", handle.name
                handle.close()
                handle = None
                g_counter = 0
    else:
        print "Closing", tt.name
        handles.remove(tt)
        ll -= 1

if handle:
    handle.write(entry)
    handle.close()
