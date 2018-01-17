from ftplib import FTP
from os.path import join as pjoin
import os
from tqdm import tqdm
from pandas import DataFrame
from joblib import Parallel, delayed
import multiprocessing
import shutil
from subprocess import call, check_output
import signal
from io import StringIO,BytesIO


ncbi = "ftp.ncbi.nlm.nih.gov"
num_cores = 4

def download(info):

    def handler(signum, frame):
        print("FTP download timeout")
        raise(Exception("FTP Downlaod timeout"))

    ftp_head, genome_path, genome_file = info

    dpat = ftp_head.split("/")
    dir = "/".join(dpat [3:])

    fhead = dpat[-1]
    signal.signal(signal.SIGALRM, handler)
    signal.alarm(60)

    try :
        ftp = FTP(ncbi)
        ftp.login()
        ftp.cwd(dir)
        ftp.retrbinary("RETR " + fhead + "_genomic.gff.gz", open(pjoin(genome_path, fhead +".gff.gz"),"wb").write)
        ftp.close()

#        gffs = [f for f in os.listdir(genome_path) if f.endswith(".gff.gz")]
#        with open(genome_file, "w") as handle:
#            call(["zcat", pjoin(genome_path, *gffs)], stdout = handle)
    except Exception as exc:
        print(exc)

    signal.alarm(0)

if __name__ == '__main__':
    ftp =  FTP(ncbi)
    data_path = "genomes"
    if not os.path.exists(data_path) :
        os.makedirs(data_path)
    print("Getting metadata from ncbi")
    FNULL = open(os.devnull, 'w')
    ftp.login()
    ftp.cwd('genomes/refseq/bacteria/')
    info = BytesIO()
    ftp.retrbinary("RETR " + "assembly_summary.txt", info.write)
    info.seek(0)
    metadata = DataFrame.from_csv(info, sep="\t", header=1)
    ftp.close()
    metadata['assembly_level'] = metadata['assembly_level'].apply(lambda x: x.replace(" ","_"))
    metadata = metadata.transpose().to_dict()

    to_dl = []
    for k,v in tqdm(metadata.items()):
        genome_path = pjoin(data_path, v['assembly_level'].replace(" ","_"), k)
        genome_file = pjoin(genome_path, k + ".gff.gz")
        if not os.path.exists(genome_path) :
            os.makedirs(genome_path)
        if not os.path.exists(genome_file):
            to_dl += [ (v['ftp_path'],genome_path, genome_file)]

    dlstuff= Parallel(n_jobs=num_cores)(delayed(download)(i) for i in tqdm(to_dl))
