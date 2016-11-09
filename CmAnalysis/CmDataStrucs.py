#Data Structures for Storage of CmScan/CmSearch Output

class CmHit:
    '''CmHit - Data Container for CMsearch output
    '''

    def __init__(self, dic=None):
        if dic is None:
            raise ValueError('No dictionary given')

        self.seqname = dic['name']
        self.seq_accession = dic['accession']
        self.cm = dic['query_name']
        self.cm_accession = dic['accession_mdl']

        self.cm_from  = int(dic['mdl_from'])
        self.cm_to    = int(dic['mdl_to'])
        self.cm_range = (self.cm_from, self.cm_to)

        self.seq_from  = int(dic['from'])
        self.seq_to    = int(dic['to'])
        self.seq_range = (self.seq_from, self.seq_to)
        #TODO Hits on - strand have seq_to > seq_from!!
        self.strand = dic['strand']
        self.trunc = dic['trunc']

        self.gc = dic['gc']
        self.bias = dic['bias']
        self.bitscore = dic['score']
        self.evalue = float(dic['e-value'])
        self.inc = dic['inc']

        self.struc = None
        self.seq = None #Just gets filled when parsing from stdout

        # Omitted, serve no practical purpose
        # self.pass = dic['pass']
        # self.mdl = dic['mdl']
        # self.description = dic['description']

    def _validate(self):
        '''Verify data structure integrity'''
        pass

    def seq_record(self):
        '''Return SeqRecord Object when self.seq is set
        '''
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC

        if self.seq:
            rec = SeqRecord(Seq(self.seq, IUPAC.unambiguous_rna),
                id=self.seqname,
                name=self.seqname,
                description='')
            #TODO Add annotation to SeqRecord
            return rec
        else:
            raise ValueError('No Sequence found')


class CompoundCmHit:
    '''CompoundCmHit - Store overlapping hits of cmsearch
    Args:
        seed_hit(CmHit object)
    Raises:
        ValueError - No Seed was given for initialization
        TypeError - Seed was not CmHit object
    '''
    from utils import Overlap

    def __init__(self, seed_hit=None):
        #Check input
        if seed_hit is None:
            raise ValueError('Needs Seed for initialization')
        elif not isinstance(seed_hit, CmHit):
            raise TypeError('Wrong Seed Hit format (must be CmHit)')
        else:
            self.overlap_threshold = 0.90
            self.seq      = seed_hit.seq
            self.seqname  = seed_hit.seqname
            self.seqID    = seed_hit.accession
            self.seq_from = seed_hit.seq_from
            self.seq_to   = seed_hit.seq_to
            self.range    = seed_hit.range
            self.cms      = [seed_hit.cm]
            self.evalues  = {seed_hit.cm : seed_hit.evalue}
            self.ranges   = {seed_hit.cm : seed_hit.seq_range}

    def add_hit(self, hit):
        '''Add a hit to an existing CompountCmHit object
        Args:
            hit(CmHit)
        Returns:
            True - Hit added successfully
            False - Hit is not overlapping with seed hit and was not added
        Raises:
            TypeError - Hit is not CmHit object
            ValueError - Hit is on different sequence than seed hit
        '''
        if not isinstance(hit, CmHit):
            raise TypeError('Wrong Seed Hit format (must be CmHit)')
        elif not hit.seqname == self.seqname:
            raise ValueError('Hit belongs to different seq, seed: ' + self.seqname)
        elif Overlap(self.seq_from, self.seq_to, hit.seq_from, hit.seq_to) > self.overlap_threshold:
            self.cms.append(hit.cm)
            self.evalues[hit.cm] = hit.evalue
            self.ranges[hit.cm] = hit.seq_range
            return True
        else:
            return False

    def seq_record(self):
        '''Return SeqRecord Object when self.seq is set
        '''
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC

        if self.seq:
            rec = SeqRecord(Seq(self.seq, IUPAC.unambiguous_rna),
                id=self.seqname,
                name=self.seqname,
                description='')
            #TODO Add cm annotation to SeqRecord
            return rec
        else:
            raise ValueError('No Sequence found')
