def ParseCmscanTblout(filename):
    '''Parse the Output of cmscan --tblout
    Args:
        filename(string)
    Returns:
        query_summmary(dic) - Dictionary{QueryName : [List of CMscan Hits as dics]}
    '''
    hits = []

    #Handle non-existing file
    try:
        TBLOUT = open(filename)
    except IOError:
        print 'No such file!'
        sys.exit() #remains here until i figure out how to solve properly TODO

    for line in TBLOUT:
        line = line.rstrip('\n')

        if line[0] == '#':
            continue
        else:
            data = line.split()
            #Column 17 can contain whitespaces - fix before error-checking
            #if len(data) != 17:
            #    raise ValueError('Could not parse ' + filename)
            dic = {'name': data[0], 'accession': data[1],
            'query_name': data[2], 'accesion_mdl': data[3],
            'mdl': data[4],
            'mdl_from': data[5], 'mdl_to': data[6],
            'from': data[7], 'to': data[8],
            'strand': data[9], 'trunc': data[10],
            'pass': data[11], 'gc': data[12],
            'bias': data[13],'score': data[14],
            'e-value': data[15], 'inc': data[16]}

            hits.append(dic)

    TBLOUT.close()
    return hits


def ParseCmscanStdOut(filename):
    print 'Not implemented yet'
    pass


def ParseSequenceFile(filename):
    formats = ['fasta', 'genbank', 'gb', 'embl', 'seqxml', 'tab']
    recs = []

    for fmt in formats:
        try:
            for r in SeqIO.parse(filename, fmt):
                recs.append(r)
        except IOError:
            print('Could not open ' + filename)
            #sys.exit() #TODO
        except ValueError:
            print('Could not use format: ' + fmt)
            continue
    return recs


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

def WussToVienna(wuss):
    '''wuss2vienna - convert a WUSS Secondary Structure to Dot-Bracket Format
    Args:
        wuss(str): Secondary structure in WUSS format
    Returns:
        vienna(str): secondary structure in dot-bracket format
    Raises:
        atm - Nothing (TypeError: if wuss is not a string)
    '''
    import re

    vienna = wuss
    #Base pairs
    vienna = re.sub('[\[{<]', '(', vienna)
    vienna = re.sub('[\]}>]', ')', vienna)
    #Unpaired Bases
    vienna = re.sub('[:_\-,~]', '.', vienna)
    #Pseudoknots
    vienna = re.sub('[aA-zZ]', '.', vienna)

    return vienna


if __name__ == "__main__":
    print "Just a module.."

__author__ = "romanoch"
__date__ = "$Oct 24, 2016 3:49:21 PM$"
