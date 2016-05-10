import sh
import os
from tqdm import tqdm
import os.path

class orthoMCL(object):
    orthoMCL_path = "/home/moritz/share/orthomclSoftware-v2.0.9/bin/"
    db_cfg = "/home/moritz/repos/Pyscratches/20150610_orthomcl_tools/orthmcl_tools/orthomcl.config.template"
    mysql_line = "mysql --defaults-file=~/sandboxes/msb_5_1_73/mysql.cnf -u orthomcl --protocol=TCP -P 5173"

    def __init__(self, out_dir, proteoms, name = "MCL"):
        self.out_dir = out_dir
        # list of prokka faa file
        self.proteoms = proteoms
        self.name = name
        self.out_mcl = out_dir+ "final_clusters.tsv"
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
   
    def full_pipe(self):
        self.empty_all_tables()
        print "Installing:"
  #      self.install_schemes()
        self.make_compliant_fastaa()
        self.make_good_fasta()
        self.do_the_blast()
        self.parse_blastout()
        self.load_into_db()
        self.make_pairs()
        self.dump_pairs()
        self.run_mcl()
        self.parse_mcl()

    def pre_blast(self):
        self.make_compliant_fastaa()
        self.make_good_fasta()


    def post_blast(self):
#        self.start_server()
#        self.empty_all_tables()
#        self.install_schemes()
        self.load_into_db()
        self.make_pairs()
        self.dump_pairs()
        self.stop_server()
        self.run_mcl()
        self.parse_mcl()

    def start_server(self):
        print "Start MySQL server"
        start_cmd = sh.Command('/home/moritz/sandboxes/msb_5_1_73/start')
        start_cmd("--myisam_max_sort_file_size=6400000000", "--myisam_sort_buffer_size=7000000000","--read_buffer_size=250000000")

    def stop_server(self):
        print "Stopping MySQL server"
        stop_cmd = sh.Command('/home/moritz/sandboxes/msb_5_1_73/stop')
        stop_cmd()
        
        
    def install_schemes(self):
        print "installing orthoMCL DB-schemes"
        script = sh.Command(self.orthoMCL_path + "orthomclInstallSchema")
        script(self.db_cfg)

    def make_compliant_fastaa(self):
        print "making compliant fastaas"
        script = sh.Command(self.orthoMCL_path + "orthomclAdjustFasta")
        if not os.path.exists(self.out_dir + "temp"):
            os.makedirs(self.out_dir + "temp")
        for f in tqdm(self.proteoms):
            script(f.split("/")[-1].split(".")[0], f, 1)
            sh.mv(f.split("/")[-1].split(".")[0] + ".fasta",self.out_dir + "temp/" )

    def make_good_fasta(self, min_len = 10, max_stop_freq=20):
        print "making filtered fastaa"
        script = sh.Command(self.orthoMCL_path + "orthomclFilterFasta")
        script(self.out_dir + "temp/", min_len, max_stop_freq)
        sh.mv("goodProteins.fasta", self.out_dir)

    def do_the_blast(self, processors = 16):
        print "running blast"
        sh.formatdb( "-i", self.out_dir + "goodProteins.fasta",  "-p", "T")
        blastall_cmd = ["blastall", "-p", "blastp", "-F", "'m S'", "-v", 100000, "-b", 100000, "-z", 0,  "-e", 1e-5, "-m", 8, "-a", processors, "-d", self.out_dir + "goodProteins.fasta", "-i", self.out_dir + "goodProteins.fasta", ">", self.out_dir + "blast_all_v_all.tsv"]
        blastall_cmd = [str(c) for c in blastall_cmd]
        os.system(" ".join(blastall_cmd))

    def parse_blastout(self):
        print "parsing the blastout"
        script = sh.Command(self.orthoMCL_path + "orthomclBlastParser")
        script(self.out_dir + "blast_all_v_all.tsv", self.out_dir + "temp/",  _out = self.out_dir + "blast_all_v_all.parsed.tsv")

        
    def load_into_db(self):
        print "loading into the DB"
        script = sh.Command(self.orthoMCL_path + "orthomclLoadBlast")
        script(self.db_cfg, self.out_dir + "blast_all_v_all.parsed.tsv")


    def empty_all_tables(self):
        print "clean db"
        os.system('echo "TRUNCATE orthomcl.SimilarSequences;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.CoOrtholog;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.InParalog;" | ' + self.mysql_line)
#        os.system('echo "TRUNCATE orthomcl.InterTaxonMatch;" | ' + self.mysql_line)
        os.system('echo "TRUNCATE orthomcl.Ortholog;" | ' + self.mysql_line)
        

    def make_pairs(self):
        print "making pairs with db-magik"
        script = sh.Command(self.orthoMCL_path + "orthomclPairs")
        script(self.db_cfg, self.out_dir + "pairs.log", "cleanup=yes")

    def dump_pairs(self):
        if os.path.exists( "pairs/"):
            sh.rm("-r", "pairs/")
        if os.path.exists(self.out_dir + "pairs/"):
            sh.rm("-r", self.out_dir + "pairs/")
        print "dumping pairs"
        script = sh.Command(self.orthoMCL_path + "orthomclDumpPairsFiles")
        script(self.db_cfg)
        sh.mv("mclInput", self.out_dir)
        sh.mv( "pairs", self.out_dir)

    def run_mcl(self):
        print "Running MCL"
        sh.mcl(self.out_dir + "mclInput",  "--abc",  "-I", 1.5, "-o", self.out_dir + "mclOutput")

    def parse_mcl(self):
        print "parsing MCL"
        os.system(" ".join([self.orthoMCL_path + "orthomclMclToGroups", self.name + "_", "00000000", "<", self.out_dir + "mclOutput", ">", self.out_mcl]) )
