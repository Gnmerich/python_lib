def ParseCmscanTblout(filename):
    '''Parse the Output of cmscan --tblout
    For detailed description of the .tblout format, please refer to
    http://eddylab.org/infernal/Userguide.pdf page 59.
    Args:
        filename(string)
    Returns:
        hits(List) - List[{Column_name : Column_content}, ..]
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
            'query_name': data[2], 'accession_mdl': data[3],
            'mdl': data[4],
            'mdl_from': data[5], 'mdl_to': data[6],
            'from': data[7], 'to': data[8],
            'strand': data[9], 'trunc': data[10],
            'pass': data[11], 'gc': data[12],
            'bias': data[13],'score': data[14],
            'e-value': data[15], 'inc': data[16],
            'description' : ' '.join(data[17:])}

            hits.append(dic)

    TBLOUT.close()
    return hits


def ParseCmscanStdOut(filename):
    # TODO - Compared to the tblout-version, we also parse the seed/query alignment along with the appropriate secondary structure
    '''Parse the STDOUT of cmscan
    For detailed description of the format, please refer to
    http://eddylab.org/infernal/Userguide.pdf page 97.
    Args:
        filename(string)
    Returns:
        hits(List) - List[{Column_name : Column_content}, ..]
    '''
    print 'Not implemented yet'
    sys.exit()
    hits = []

    #Handle non-existing file TODO - use with keyword (much more elegant)
    try:
        STDOUTFILE = open(filename)
    except IOError:
        print 'No such file!'
        sys.exit() #remains here until i figure out how to solve properly TODO

    for l in STDOUTFILE:
        break

    STDOUTFILE.close()
    return


def ParseSequenceFile(filename):
    formats = ['fasta', 'genbank', 'gb', 'embl', 'seqxml', 'tab']
    recs = []

    for fmt in formats:
        try:
            for r in SeqIO.parse(filename, fmt):
                recs.append(r)
        except IOError:
            print('Could not open ' + filename)
            #sys.exit() #TODO Don't know whether i need those excepts
        except ValueError:
            print('Could not use format: ' + fmt)
            continue
    return recs


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


__author__ = "romanoch"
__date__ = "$Oct 24, 2016 3:49:21 PM$"
