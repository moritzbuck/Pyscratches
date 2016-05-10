""" merges and nests gff files 


Usage:
  gffmerger.py <file1.gff> <file2.gff>
  fastasubstracter.py -h
  
Options:
   -h --help       Show this screen.
"""
from docopt import docopt
import gffutils
import os

class GFFmerger(object):

    def __init__(self, prokka_gff=None, interpro_gff=None, db_file = None):
        self.prokka_gff = prokka_gff
        self.interpro_gff = interpro_gff
        self.db_file = db_file
        self.prokka_db = self.correct_prokka_gff(self.prokka_gff)
        self.interpro_db = self.attach_interpro_gff(self.prokka_db,self.interpro_gff)
        
    def id_transform(self,dodo):
        if dodo.featuretype == 'repeat_region' : return  None
        if dodo.featuretype == 'gene':
            return dodo['locus_tag'][0]
        else:
            return dodo.id
    
    def entry_transform(self,dodo):
        if dodo.featuretype == 'repeat_region' :
            return  dodo
        if dodo.featuretype  == 'gene':
            dodo['locus_tag'] = dodo['locus_tag ']
            del dodo.attributes['locus_tag ']
            return dodo
        else:
            dodo.id = dodo['locus_tag'][0] + "_" + dodo.featuretype
            dodo['ID'] = dodo['locus_tag'][0] + "_" + dodo.featuretype
            dodo['Parent'] = dodo['locus_tag']
            dodo['parent_type'] = 'gene'
            return dodo       

    def attach_transform(self,dodo):
        dodo.seqid = dodo.seqid=self.prokka_db[dodo.seqid + "_CDS"].seqid 
        if dodo.featuretype  == 'polypeptide':
            cds =  self.prokka_db[dodo['ID'][0]]
            dodo['Parent'] = dodo['ID'][0] + "_CDS"
            dodo['parent_type'] = 'CDS'
            dodo['ID'] = dodo['ID'][0] + "_PP"
            dodo.id = dodo['ID'][0] + "_PP"
        else:
            dodo['Parent'] = dodo['Target'][0].split(" ")[0] +"_PP"
            dodo['parent_type'] = 'polypeptide'
            dodo['ID'] = dodo['Target'][0].split(" ")[0] + "_" + dodo.source
            cds = self.prokka_db[dodo['Target'][0].split(" ")[0] + "_CDS"]
        if cds.strand == '+' :
            dodo.start = cds.start + dodo.start -1
            dodo.end = cds.start + dodo.end*3 +2
        else :
            dodo.strand = "-"
            end = cds.end - dodo.end*3 -2
            start = cds.end - dodo.start +1
            dodo.start = end
            dodo.end = start
        return dodo       

    
    def correct_prokka_gff(self,file_name):
        db = gffutils.create_db(file_name, self.db_file, force = True, transform=self.entry_transform, id_spec=self.id_transform)
        return db

    def attach_interpro_gff(self,db,file_name):
        ip = gffutils.create_db(file_name, "tmp.db", force = True, transform = self.attach_transform,merge_strategy="create_unique")
        db.update(ip)
        os.remove("tmp.db")

    def write_db(self, file_name):
        with open(file_name,"w") as stream:
            for f in self.prokka_db.all_features():
                stream.write(str(f) + "\n")


    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
