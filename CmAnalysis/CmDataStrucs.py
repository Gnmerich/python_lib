#Data Structures for Storage of CmScan/CmSearch Output

class CmHit:
    '''CmHit - Store one result of cmsearch in easy and accessible manner
    '''

    def __init__(self, dic):
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

        self.strand = dic['strand']
        self.trunc = dic['trunc']

        self.gc = dic['gc']
        self.bias = dic['bias']
        self.bitscore = dic['score']
        self.evalue = float(dic['e-value'])
        self.inc = dic['inc']

        self.struc = None
        self.seq = None

        # Omitted, serve no practical purpose
        # self.pass = dic['pass']
        # self.mdl = dic['mdl']
        # self.description = dic['description']

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
    '''CompoundCmHit - Store one result of cmsearch in easy and accessible manner
    '''
    pass
