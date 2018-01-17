metacyc_path = "/home/moritzbuck/data/metacyc/20.0/data/"

        
class Reaction(object):
    def __repr__(self):
        return "The " + self.name + " reaction from metacyc with " + (("EC number " + self.ec ) if self.ec else "no EC number" )

    def __init__(self, block) :
        self.__attribute__ = {k : [] for k in set([l[0] for l in block])}
        for l in block:
            self.__attribute__[l[0]] += [l[1]]
        self.id = self.__attribute__['UNIQUE-ID'][0]
        self.name = self.__attribute__['COMMON-NAME'][0] if  self.__attribute__.has_key('COMMON-NAME') else "Unnamed"
        self.ec = self.__attribute__['EC-NUMBER'][0].replace("|","") if  self.__attribute__.has_key('EC-NUMBER') else None

    def get_reactions(self):
        return [self]
        
class Pathway(object):
    def __repr__(self):
        return "The " + self.name + " pathway from metacyc with " + str(len(self.reactions)) + " reactions"

    def __init__(self, block, reaction_dict) :
        self.__attribute__ = {k : [] for k in set([l[0] for l in block])}
        for l in block:
            self.__attribute__[l[0]] += [l[1]]
        self.id = self.__attribute__['UNIQUE-ID'][0]
        self.name = self.__attribute__['COMMON-NAME'][0] if  self.__attribute__.has_key('COMMON-NAME') else "Unnamed"
        self.reactions = None
        
    def update(self, reaction_dict, pathway_dict):
        self.reactions = [ reaction_dict[r] if reaction_dict.has_key(r) else pathway_dict[r]  for r in self.__attribute__['REACTION-LIST']]

    def get_reactions(self):
        return sum([r.get_reactions() for r in self.reactions],[])

    def get_ecs(self):
        return set([r.ec.split("EC-")[1] for r in self.get_reactions() if r.ec])

    def eced_ratio(self):
        return float(len(self.get_ecs()))/len(self.get_reactions())

    def completness(self, ec_set):
        pat_ecs = self.get_ecs()
        overlap_len = len(ec_set.intersection(pat_ecs))
        return float(overlap_len)/len(pat_ecs)
    
def parse_metacyc_file(file):
    with open(file) as handle:
        lines = [l[:-1] for l in  handle.readlines() if l[0] != "#"]
        
    blocks = []
    block = []
    for l in lines:
        if  l=="//":
            blocks += [ block ]
            block = []
        else :
            if l[0] == "/":
                block[-1][1] +=  "\n" + l[1:]

            else :
                block += [ l.split(" - ") ]
    return blocks
                
    

blocks = parse_metacyc_file(metacyc_path + "/reactions.dat")
reactions = [Reaction(b) for b in blocks]
reactions = {b.id : b for b in reactions}

blocks = parse_metacyc_file(metacyc_path + "/pathways.dat")
pathways = [Pathway(b, reactions) for b in blocks]
pathways = {b.id : b for b in pathways}
for p in pathways.values():
    p.update(reactions, pathways)

pathways = {k : v for k,v in pathways.iteritems() if len(v.get_ecs()) >0 }
