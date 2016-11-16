#!/bin/python

class SomAnalysis:
    '''Analysis framework for R's Self-Organized-Maps
    '''

    def __init__(self, SOM_clusterfile, SOM_Mapfile):
        self.bins     = []
        self.clusters = []







    def __iter__(self):
        return self.bins

#    def __next__(self): #Attention: __next__() in python 3, just set __next__ = next for now
#        pass
#        if True:
#            raise StopIteration
#
#    next = __next__ #for python 3 compatibility

class SOM:
    '''Data container'''

    #Attributes:
    data
    clusters
    bins = []

    def  __init__(self, filename=None, mapfile=None, codes=None):
        if filename is None or mapfile is None:
            raise ValueError('Not all files provided')
        self.hits, self.cms = _ParseSom(filename)
        self.uid_map = _ParseMapfile(mapfile)

        if codes not is None:
            self.code_vectors = _ParseCodes(self, codes)

    def _ParseSOM(self, filename):
        hits = []
        with open(filename, 'r') as SOM_OUT:
            header = SOM_OUT.readline().split()
            cms = header[3:]
            for line in SOM_OUT:
                line_tmp = line.split('\t')
                hits.append((line_tmp[0], line_tmp[1], line_tmp[3]))

        return hits, cms

    def _ParseMapfile(self, mapfile):
        with open(mapfile, 'r') as MAP:
            pass

    def _ParseCodes(self, codefile):
        pass
