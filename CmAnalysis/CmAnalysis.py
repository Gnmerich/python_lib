__author__ = "romanoch"
__date__ = "$Oct 24, 2016 2:04:42 PM$"

if __name__ == "__main__":
    print "Not executable"
    sys.exit()

import Parsing_Utils as PU
from Bio import SeqIO
from CmDataStrucs import CmHit, CompoundCmHit
from math import log10

class CmAnalysis:
    '''CmAnalysis - Store results of cmsearch and analyze several criteria
    '''

    def __init__(self, tblout_fi=None, seq_fi=None):
        if tblout_fi is None or seq_fi is None:
            raise ValueError('No arguments provided!')
        #General data
        self.cm_hits, self.cm_hits_UIDs  = self._GetCmHits(tblout_fi)
        self.seq_records    = self._GetSeqRecords(seq_fi)
        self.hit_map        = {}
        self.loci_map       = {}
        #Sets
        self.seqids_general  = set()
        self.seqids_withhits = set()
        self.seen_cms        = set()
        self._FillSets()
        #unique loci
        self.unique_loci, self.unique_loci_uids = self._GetUniqueLoci()
        #extract rna sequence for each hit from fasta file
        self._GetHitSequences()


    def MapHits(self, attribute):
        '''Generate a Mapping for your desired attribute
            eg. mapping['attribute'] -> [matching hits]
            Seqname: mapping['Dengue Virus 1'] -> [CmHits in Dengue Virus 1]
            CM: mapping['SL.KOKV.1'] -> [CmHits with SL.KOKV.1]
            Evalue (stupid but possible): mapping[evalue=0.00045] -> [CmHits with evalue=0.00045]
        '''
        if attribute in self.hit_map:                  #Dictionary already exists
            return self.hit_map[attribute]
        elif attribute not in self.cm_hits[0].__dict__: #Stop if attribute does not exist
            return False
        #TODO
        #Do this via setattr(), getattr() and hasattr() !
        else:
            tmp_dic = {}
            for hit in self.cm_hits:
                key = hit.__dict__[attribute]
                try:
                    tmp_dic[key].append(hit)
                except KeyError:
                    tmp_dic[key] = []
                    tmp_dic[key].append(hit)

            self.hit_map[attribute] = tmp_dic
            return self.hit_map[attribute]

    def MapLoci(self, attribute):
        '''Same as above just with unique loci and not hits (caution: just possible for seqname and accession numbers)
        '''
        if attribute in self.loci_map:
            return self.loci_map[attribute]
        elif attribute not in ['seqname', 'accession']:
            return False
        else:
            tmp_dic = {}
            for locus in self.unique_loci:
                key = locus.__dict__[attribute]
                try:
                    tmp_dic[key].append(locus)
                except KeyError:
                    tmp_dic[key] = []
                    tmp_dic[key].append(locus)

            self.loci_map[attribute] = tmp_dic
            return self.loci_map[attribute]


    def _UpdateMap(self, CmHit):
        '''Update the mappings after hits are added manually
        '''
        pass

    def _GetCmHits(self, tblout):
        hit_objects = []
        cm_hits_UIDs = {}
        uid = 0
        for dic in PU.ParseCmscanTblout(tblout):
            tmp = CmHit(dic)
            tmp.uid = uid
            hit_objects.append(tmp)
            cm_hits_UIDs[uid] = tmp
            uid += 1
        return hit_objects, cm_hits_UIDs

    def _GetSeqRecords(self, seq_fi):
        seq_records = {}
        seq_records = PU.ParseSequenceFile(seq_fi)
        return seq_records

    def _GetHitSequences(self):
        for hit in self.cm_hits:
            hit.seq = self.seq_records[hit.seqname].seq[hit.seq_from:hit.seq_to]

    def _FillSets(self):
        '''Iterate over hits and fill sets'''
        for hit in self.cm_hits:
            self.seen_cms.add(hit.cm)
            self.seqids_withhits.add(hit.seqname)
        #From sequence (FASTA) file
        for rec_name in self.seq_records:
            self.seqids_general.add(self.seq_records[rec_name].id)

    def _GetUniqueLoci(self):
        self.MapHits('seqname') #Generate dictionary where seqname->[list of cm hits]
        unique_loci = []
        unique_loci_uids = {}

        for seqname in self.hit_map['seqname']:
            loci = []
            for hit in self.hit_map['seqname'][seqname]:
                for l in loci:
                    if l.add_hit(hit): #if adding successful then break
                        break
                else: #if not added to any existing locus - creat new one
                    loci.append(CompoundCmHit(hit)) #if the hit couldn't be added to any existing hit -> create new locus
            unique_loci += loci
        #Assign a Unique ID to each locus
        uid = 0
        for locus in unique_loci:
            locus.uid = uid
            unique_loci_uids[uid] = locus
            uid += 1
        return unique_loci, unique_loci_uids

    # def HitsInRange(self, seqname, start, end):
    #     '''Return all unique loci which overlap to some degree
    #     TODO - handle minus strand hits
    #     '''
    #     self.MapLoci('seqname') #Generate dictionary where seqname->[list of unique loci]
    #     overlaps = []
    #     for locus in self.loci_map[seqname]:
    #         if locus.seq_from > start and locus.seq_from < end:
    #             overlaps.append(locus)
    #         elif locus.seq_to > start and locus.to < end:
    #             overlaps.append(locus)
    #         else:
    #             continue
    #     return overlaps

    #Methods for post-intit data structure manipulation (TODO: get rid of those methods)
    def add_seq(self, seq_records):
        '''Add a new SeqRecord to the pool and update dependent data strucs
        Args:
            seq_records - List of Biopython Sequence Record objects
        '''
        for rec in seq_records:
            self.seqids_general.add(rec.id)
            self.seq_records[rec.id] = rec

    def add_cmhit(self, cm_hits):
        '''Add a new CMsearch hit to the pool and update dependent data strucs
        Args:
            cm_hits - List of Dictionaries (see ParseCmScanTblout dict)
        '''
        for hit in self.cm_hits:
            new_hit = CmHit(hit)
            self.cm_hits.add(new_hit)
            self.seen_cms.add(new_hit.cm)
            self.seqids_withhits.add(new_hit.seqname)

    def cmVectors(self, filename=None):
        '''Compose E-value vectors of all unique loci, writes Vectorfile and Mapfile
        Args:
            filename (string): Name of Output File
        '''
        cm_sorted = sorted(self.seen_cms)
        if filename is None:
            filename = 'vectors'

        with open(filename+'.vec', 'w') as VEC_OUT, open(filename+'.map', 'w') as MAPFILE:
            #Headers
            VEC_OUT.write('\t'.join(cm_sorted) + '\n')
            MAPFILE.write('#UID\tFROM\tTO\tSTRAND\tSEQNAME\n')
            #Data
            for locus in self.unique_loci:
                evals = []
                for cm in cm_sorted:
                    try:
                        evals.append(str(abs(log10(locus.evalues[cm]))))
                    except KeyError:
                        evals.append(str(0))
                VEC_OUT.write(str(locus.uid) + '\t')
                VEC_OUT.write('\t'.join(evals) + '\n')
                mapstr = [str(locus.uid), str(locus.seq_from), str(locus.seq_to), locus.strand, locus.seqname]
                MAPFILE.write('\t'.join(mapstr) + '\n')
