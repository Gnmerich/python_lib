__author__ = "romanoch"
__date__ = "$Oct 24, 2016 2:04:42 PM$"

if __name__ == "__main__":
    print "Hello World"


class CmAnalysis:
    '''CmAnalysis - Store results of cmsearch and analyze several criteria
    '''
    import Parsing_Utils as PU
    from Bio import SeqIO
    from CmHit import CmHit

    def __init__(self, tblout_fi=None, seq_fi=None):
        self.cm_hits = None
        self.seq_records = None
        self.seqids_general = None
        self.seqids_withhits = None
        self.seen_cms = None
        self.unique_loci = None

        #Parse CMsearch/Sequence Data if present
        if tblout:
            self.cm_hits = _GetCmHits(tblout)
            self.seen_cms = set()
            self.seqids_withhits = set()
        if fasta:
            self.seq_records = _GetSeqRecords(fasta)
            self.seqids_general = set()

        #Gather information about how many hits we have
        for hit in self.cm_hits:
            self.seen_cms.add(hit.cm)
            self.seqids_withhits.add(hit.seqname)
        for rec in self.seq_records:
            self.seqids_general.add(rec.id)


    def _GetCmHits(tblout):
        hit_objects = []
        for dic in PU.ParseCmscanTblout:
            hit_objects.append(CmHit(dic))
        return hit_objects

    def _GetSeqRecords(self):
        seq_records = []
        for rec in PU.ParseSequenceFile(seq_fi):
            seq_records.append(rec)
        return seq_records

    def _GetUniqueLoci(self):
        #Not finished TODO
        # from util import Overlap
        uniq_loci = []

        for hit in self.cm_hits:
            name = hit.seqname
            start = hit.seq_from
            end = hit.seq_to
            name_mdl = hit.cm
            e_val = abs(log10(hit.evalue))

            if len(uniq_loci) == 0:  #Start, take first hit as locus
                locus = {'name' : name, 'start' : start, 'end' : end, name_mdl : e_val}
                uniq_loci.append(locus)
                continue
            else:
                for locus in uniq_loci:    #Check if hit overlaps with existing hits
                    if Overlap(start, end, int(locus['start']), int(locus['end'])) > 0.90:  #Overlapping hits found
                        locus[name_mdl] = e_val
                        break
                else: #New Locus found
                    #print 'no overlap', overlap(start, end, locus['start'], locus['end'])
                    locus = {'name' : name, 'start' : start, 'end' : end, name_mdl : e_val}
                    uniq_loci.append(locus)

        return uniq_loci

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
