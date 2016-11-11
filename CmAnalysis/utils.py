from math import log10

def _UniqueLoci(hit_list): #DEPRECATED
    unique_positions = []
    for hit in hit_list:
        name = hit['name']
        start = int(hit['from'])
        end = int(hit['to'])
        name_mdl = hit['query_name']
        e_val = abs(log10(float(hit['e-value'])))  #try log10

        if len(unique_positions) == 0:  #Start, take first hit as locus
            locus = {'name' : name, 'start' : start, 'end' : end, name_mdl : e_val}
            unique_positions.append(locus)
            continue
        else:
            for pos in unique_positions:    #Check if hit overlaps with existing hits
                if Overlap(start, end, int(pos['start']), int(pos['end'])) > 0.90:  #Overlapping hits found
                    pos[name_mdl] = e_val
                    break
            else: #New Locus found
                #print 'no overlap', overlap(start, end, pos['start'], pos['end'])
                locus = {'name' : name, 'start' : start, 'end' : end, name_mdl : e_val}
                unique_positions.append(locus)

    return unique_positions


def Overlap(start1, end1, start2, end2):
    '''Ad-hoc measure for quantifying the amount of overlap between 2 regions on 1 genome'''
    len1 = end1-start1
    len2 = end2-start2

    startdiff = start2-start1
    enddiff = end2-end1

    diff_sum = abs(startdiff)+abs(enddiff)
    shared1 = 1-(diff_sum/len1)
    shared2 = 1-(diff_sum/len2)

    frac_overlap = (shared1 + shared2) / 2

    return frac_overlap


def Evalue4Bed(evalue, scale_factor=50, max_score=1000):
    '''Convert cmsearch evalue (float) to an integer for the BED-file score field
    Args:
        evalue (float)
        scale_factor (int): scale log10(evalue) by this factor, Default = 50
        max_score (int): cutoff, values > max_score are set to max_score, Default = 1000
    Returns:
        scaled evalue (int)
    '''
    ev_log = log10(evalue)
    ev_int = int(ev_log)

    if ev_int > 0:
        return 0
    elif abs(ev_int) * scale_factor > max_score:
        return max_score
    else:
        return abs(ev_int) * scale_factor


def AlignLoci(outname, fasta, loci):
    '''AlignLoci - build a sequence alignment from several loci
    Args:
        outname - name of the outfile
        fasta(string) - File containing Sequences in FASTA-format
        loci - list of tuples (FASTA-header, from, to, inc)
    Returns:

    '''
    from Bio import SeqIO

    records = {}
    FAin = open(fasta, "rU")
    for rec in SeqIO.parse(FAin, 'fasta'):
        records[rec.id] = rec

    FAout = open(outname, 'w')
    for i in loci:
        id, fr, to, inc = i
        seq = records[id]
        SeqIO.write(seq[int(fr):int(to)], FAout, 'fasta')

    #Then perform alifold, mlocarna, (or AliDot) on it
    return 0


__author__ = "romanoch"
__date__ = "$Oct 25, 2016 5:13:02 PM$"
