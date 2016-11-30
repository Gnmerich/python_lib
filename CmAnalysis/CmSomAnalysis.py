#!/bin/python

def main():
    pass


class SOM:
    '''Data container'''

    def  __init__(self, filename=None, mapfile=None, codes=None):
        if filename is None or mapfile is None:
            raise ValueError('Not all files provided')
        if codes not is None: #Optional
            self.code_vectors = _ParseCodes(self, codes)

        #Parse Data and auxiliary Mapping
        self.hits, self.cms = _ParseSom(filename)
        self.uid_map = _ParseMapfile(mapfile)
        self.bins = []
        self._binmap = {}
        #Build internal Data Structure
        for hit in self.hits:
            hit.seqname = self.uid_map[hit.uid]['seqname']
            hit.range = (self.uid_map[hit.uid]['from'], self.uid_map[hit.uid]['to'])
            hit.strand = self.uid_map[hit.id]['strand']
            #Build bins
            try:
                self._binmap[hit.bin].append(hit)
            except KeyError:
                self._binmap[hit.bin] = []
                self._binmap[hit.bin].append(hit)
        #Build SomBin objects
        for bin_id in self._binmap
            self.bins.append(SomBin(bin_id, self._binmap[bin_id]))

    def _ParseSOM(self, filename):
        hits = []
        with open(filename, 'r') as SOM_OUT:
            header = SOM_OUT.readline().split()
            cms = header[3:]
            for l in SOM_OUT:
                l_split = l.split('\t')
                hits.append(SomHit(l_split[:3])) #CLUSTER BIN UID
        return hits, cms

    def _ParseMapfile(self, mapfile):
        uid_map = {}
        with open(mapfile, 'r') as MAP:
            for line in MAP:
                tmp = line.split('\t')
                uid_map[tmp[0]] = {'from' : tmp[1], 'to': tmp[2], 'strand': tmp[3], 'seqname': tmp[4]}
        return uid_map

    def _ParseCodes(self, codefile):
        pass

class SomBin:
    '''Data Container for Bins, since i hate dictionaries'''
    def __init(self, bin_id, hit_list):
        self.bin_id = bin_id
        self.hits = hit_list


class SomHit:
    '''A unique locus with information regarding '''
    def __init__(self, content=None):
        self.uid = content[2]
        self.seqname = None
        self.range = None
        self.strand = None
        self.cluster = content[0]
        self.bin = content[1]

if __name__ == '__main__':
    main()
