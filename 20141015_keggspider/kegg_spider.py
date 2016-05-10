import collections
from restful_lib import Connection
from tqdm import tqdm

kegg_url = "http://rest.kegg.jp"
conn = Connection(kegg_url)

all_paths = conn.request_get('list/path', headers={'Accept':'text/json'})

path_list = [str(p.split("\t")[0]).split(":")[1] for p in  all_paths['body'].split('\n') if p!=""]

all_mds = conn.request_get('list/md', headers={'Accept':'text/json'})
mds_list = [str(p.split("\t")[0]).split(":")[1] for p in  all_mds['body'].split('\n') if p!=""]

pathways2ko = { p : [str(l.split(" ")[-1]) for l in conn.request_get('get/' + p, headers={'Accept':'text/json'})['body'].split("\n") if "KO_PATHWAY" in l] for p in tqdm(path_list)}
print len([p[0] for p in pathways2ko.values()  if len(p) > 1 ]), "more than one ko"
all_ko_ids = [p[0] for p in pathways2ko.values()  if len(p) > 0 ]

all_kos = { p : .split(" ")[-1]) for l in conn.request_get('get/' + p, headers={'Accept':'text/json'})['body'].split("\n") if "KO_PATHWAY" in l] for p in tqdm(path_list)}

