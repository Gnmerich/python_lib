__author__ = "romanoch"
__date__ = "$Oct 24, 2016 2:04:42 PM$"

if __name__ == "__main__":
    print "Hello World"


class CmAnalysis:
    '''CmAnalysis - Store results of cmsearch and analyze several criteria
    '''
    import Parsing_Utils as PU
    from Bio import SeqIO

    def __init__(self, tblout_fi=None, seq_fi=None):
        self.cm_hits = None
        self.seq_records = None
        self.seqids_general = None
        self.seqids_withhits = None
        self.seen_cms = None

        #Parse CMsearch/Sequence Data if present
        if tblout:
            self.cm_hits = PU.ParseCmScanTblout(tblout)
            self.seen_cms = set()
            self.seqids_withhits = set()
        if fasta:
            self.seq_records = PU.ParseSequenceFile(seq_fi)
            self.seqids_general = set()

        #Gather information about how many hits we have
        for hit in self.cm_hits:
            self.seen_cms.add(hit['query_name'])
            self.seqids_withhits.add(hit['name'])
        for rec in self.seq_records:
            self.seqids_general.add(rec.id)


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
            self.cm_hits.add(hit)
            self.seen_cms.add(hit['query_name'])
            self.seqids_withhits.add(hit['name'])

    def UniqueLoci(self):
        from util import Overlap
        pass



    def _members(self, foo):

        pass
