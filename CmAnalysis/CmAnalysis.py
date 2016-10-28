__author__ = "romanoch"
__date__ = "$Oct 24, 2016 2:04:42 PM$"

if __name__ == "__main__":
    print "Not executable"
    sys.exit()


class CmAnalysis:
    '''CmAnalysis - Store results of cmsearch and analyze several criteria
    '''
    import Parsing_Utils as PU
    from Bio import SeqIO
    from CmHit import CmHit, CompoundCmHit

    def __init__(self, tblout_fi=None, seq_fi=None):
        self.cm_hits = None
        self.seq_records = None
        self.seqids_general = None
        self.seqids_withhits = None
        self.seen_cms = None
        self.unique_loci = None
        self.seqname_to_hits = {} #TODO rewrite __init__, is too bloated

        #Parse CMsearch/Sequence Data if present
        if tblout_fi:
            self.cm_hits = _GetCmHits(tblout_fi)
            self.seen_cms = set()
            self.seqids_withhits = set()
        if seq_fi:
            self.seq_records = _GetSeqRecords(seq_fi)
            self.seqids_general = set()

        #Gather information about how many hits we have and how they are mapped
        for hit in self.cm_hits:
            self.seen_cms.add(hit.cm)
            self.seqids_withhits.add(hit.seqname)
            if hit.seqname not in self.seqname_to_hits:
                self.seqname_to_hits[hit.seqname] = []
                self.seqname_to_hits[hit.seqname].append(hit)
            else:
                self.seqname_to_hits[hit.seqname].append(hit)

        for rec in self.seq_records:
            self.seqids_general.add(rec.id)

        _GetUniqueLoci(self)

    def _GetCmHits(self, tblout):
        hit_objects = []
        for dic in PU.ParseCmscanTblout:
            hit_objects.append(CmHit(dic))
        return hit_objects

    def _GetSeqRecords(self, seq_fi):
        seq_records = []
        for rec in PU.ParseSequenceFile(seq_fi):
            seq_records.append(rec)
        return seq_records

    def _GetUniqueLoci(self):
        #Not finished TODO
        self.unique_loci = []

        for seqname in self.seqname_to_hits:
            for hit in seqname_to_hits[seqname]:
                loci = []
                for l in loci:
                    if l.add_hit(hit): #if adding successful then break
                        break
                    else: #if not - continue checking all other loci
                        continue
                loci.append(CompoundCmHit(hit)) #if the hit couldn't be added to any existing hit -> create new locus
            self.unique_loci += loci


    def _MapHits(self):
        #Move sequence name -> hits mapping from __init__ to here
        pass

    #Methods for post-intit data structure manipulation
    def add_seq(self, seq_records):
        '''Add a new SeqRecord to the pool and update dependent data strucs
        Args:
            seq_records - List of Biopython Sequence Record objects
        '''
        if self.seqids_general is None: self.seqids_general = set()
        if self.seq_records is None: self.seq_records = []

        for rec in seq_records:
            self.seqids_general.add(rec.id)
            self.seq_records.append(rec)

    def add_cmhit(self, cm_hits):
        '''Add a new CMsearch hit to the pool and update dependent data strucs
        Args:
            cm_hits - List of Dictionaries (see ParseCmScanTblout dict)
        '''
        if self.cm_hits is None: self.cm_hits = []
        if self.seen_cms is None: self.seen_cms = set()
        if self.seqids_withhits is None: self.seqids_withhits = set()

        for hit in self.cm_hits:
            new_hit = CmHit(hit)
            self.cm_hits.add(new_hit)
            self.seen_cms.add(new_hit.cm)
            self.seqids_withhits.add(new_hit.seqname)
