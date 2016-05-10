""" small script doing some stats and plots from a spades generated fasta


Usage:
  spades_fasta_stats.py [-a -l <l_min> -L <L_max> -c <c_min> -C <C_max> ] <fasta> <plot_header>
  spades_fasta_stats.py -h
  
Options:
   -h --help       Show this screen.
   -a              Auto-filtering (2 independent bimodal guassians)
   -l <l_min>      Length min
   -c <c_min>      Coverage min
   -L <L_max>      Length max
   -C <C_max>      Coverage max
   
"""
from docopt import docopt
from Bio import SeqIO
from tqdm import tqdm
from pandas import DataFrame
from numpy import array
import numpy
from matplotlib import pyplot
from sklearn import mixture
from matplotlib import cm
import numpy as np
import pylab as pl
from sklearn.base import BaseEstimator
from sklearn.utils import check_random_state
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances, manhattan_distances
from sklearn.datasets.samples_generator import make_blobs
from weighted import median as wmedian
from numpy import array
import itertools


np.random.seed(0)

class WKMedians(BaseEstimator):

    def __init__(self, k, max_iter=100, random_state=0, tol=1e-4):
        self.k = k
        self.max_iter = max_iter
        self.random_state = random_state
        self.tol = tol

    def _e_step(self, X,w):
        self.labels_ = array([v[1]*v[0] for v in zip(w,manhattan_distances(X,self.cluster_centers_))]).argmin(axis=1)
        
    def _average(self, X,w):
        return wmedian(array(list(itertools.chain(*X))), w)

    def _m_step(self, X,w):
        X_center = None
        for center_id in range(self.k):
            center_mask = self.labels_ == center_id
            if not np.any(center_mask):
                # The centroid of empty clusters is set to the center of
                # everything
                if X_center is None:
                    X_center = self._average(X,w)
                self.cluster_centers_[center_id] = X_center
            else:
                self.cluster_centers_[center_id] = \
                    self._average(X[center_mask],w[center_mask])

    def fit(self, X, w):
        if len(X.shape) == 1: X = array([[v] for v in X])
            
        n_samples = X.shape[0]
        vdata = np.mean(np.var(X, 0))

        random_state = check_random_state(self.random_state)
        self.labels_ = random_state.permutation(n_samples)[:self.k]
        self.cluster_centers_ = X[self.labels_]


        for i in xrange(self.max_iter):
            centers_old = self.cluster_centers_.copy()

            self._e_step(X,w)
            self._m_step(X,w)

            if np.sum((centers_old - self.cluster_centers_) ** 2) < self.tol * vdata:
                break

        return self

    def predict(self, X):
        return manhattan_distances(X, self.cluster_centers_).argmin(axis=1)





def filter_assembly(fasta, out_head, auto=True,  l_min = 0 ,c_min = 0 , L_max = None , C_max = None) : 

    with open(fasta,"r") as seqs:
        data = DataFrame.from_dict({s.id : array(s.id.split("_"))[[3,5]] for s in SeqIO.parse(seqs,"fasta")}, orient='index',dtype = numpy.float)

    L_max = L_max if L_max is not None else data[0].max()+1
    C_max = C_max if C_max is not None else data[1].max()+1

    good_contig = lambda x :  x[0] > l_min and x[0] < L_max and x[1] > c_min and x[1] < C_max 
    
    print "median length:\t"+str(numpy.median(data[0]))
    print "median coverage:\t"+str(numpy.median(data[1]))


    if(auto):
        classes = array([0]*data.shape[0])
        sub_data = data.loc[data.apply(good_contig,axis=1)]
        kmeds = WKMedians(2)
        kmeds.fit(array(numpy.log10(sub_data[1])),array(sub_data[0]))
        clust_y = kmeds.labels_
        kmeds = WKMedians(2)
#        gmm=mixture.GMM(n_components=2)
#        gmm.fit(numpy.log10(data[0]))
#        clust_x = gmm.predict(numpy.log10(data[0]))
        kmeds.fit(array(numpy.log10(sub_data[0])),array([1]* sub_data.shape[0]))
        clust_x = kmeds.labels_
        classes[array(data.apply(good_contig,axis=1))] = 1 + clust_x + 2*clust_y
    else:
        classes = array(data.apply(good_contig,axis=1).astype(numpy.int))
        print "min length is " + str(l_min) + ", max length is " + str(L_max)
        print "min coverage is " + str(c_min) + ", max coverage is " + str(C_max)
    
 


    lens = data[0].groupby(classes).apply(sum)
    plens = data[0].groupby(classes).apply(lambda x: str(numpy.round(sum(x)/1000))+"kb")
    medians_lens = data[0].groupby(classes).apply(lambda x: str(numpy.round(numpy.median(x)/1000))+"kb")
    medians_coverage = data[1].groupby(classes).apply(lambda x: str(numpy.round(numpy.median(x)))+"X")
    
    print "total lengths of bins: " + "/".join(plens)
    print "median length of contigs in bins: " + "/".join(medians_lens)
    print "median coverage of contigs in bins: " + "/".join(medians_coverage)

    out = {'lengths' : lens, 'median_lengths' : medians_lens, 'median_coverages' : medians_lens}
    data[classes == lens.idxmax()]
    keepers = data[classes == lens.idxmax()].index

    fig = pyplot.figure()
    axes = fig.add_subplot(111)
    axes.set_xscale('log')
    axes.set_yscale('log')
    pyplot.xlabel("length")
    pyplot.ylabel("coverage")
    pyplot.scatter(data[0],data[1],c=classes,cmap=cm.rainbow)
    fig.savefig(out_head+"_scatter.pdf")        
    pyplot.close()
    
    with open(fasta,"r") as seqs:
        kept_sequences = [s for s in SeqIO.parse(seqs,"fasta") if s.id in keepers]
    with open(out_head+".filtered.fasta","w") as filtered:
        SeqIO.write(kept_sequences,filtered,"fasta")

    return out


if __name__ == '__main__':
    arguments = docopt(__doc__)
    print arguments
    
    fasta = arguments['<fasta>']
    out_head = arguments['<plot_header>']
    l_min = int(arguments['-l']) if arguments['-l'] is not None else 0
    c_min = int(arguments['-c']) if arguments['-c'] is not None else 0
    L_max = int(arguments['-L']) if arguments['-L'] is not None else None
    C_max = int(arguments['-C']) if arguments['-C'] is not None else None
    auto = arguments['-a']
    filter_assembly(fasta,out_head,auto,l_min,c_min,L_max,C_max)
