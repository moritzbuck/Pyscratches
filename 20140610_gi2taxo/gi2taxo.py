gi_taxo = "/home/moritz/glob/data/kraken_large/taxonomy/gi_taxid_nucl.dmp"

gi = "379009891"

with open(gi_taxo) as f:
    for line in f:
        if gi in line: break
