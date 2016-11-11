#!/bin/python

import sys, os
import argparse
import RNA
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC

def main():
    ###Argument Parsing###
    p=argparse.ArgumentParser(description='Do various things with alignments')
    #What to do
        #default: look at whole alignment    #just look a slice    #do a sliding window     #build a CM?
    p.add_argument('--mode', default='whole', choices=['whole', 'slice', 'window'], help='Select mode of action')

    p.add_argument('--start', type=int, default=1, help='Slice Start,(1-indexed), Default: 1')
    p.add_argument('--end',   type=int, default=None, help='Slice End, Default: len(aln[0])')
    p.add_argument('--win',   type=int, default=100,  help='Set window size')
    #Metrics:
        #Calculate an SCI    #Calculate the MPI    #Provide Consensus structures/MFEs
    p.add_argument('--sci',  default=None, action='store_true', help='Compute SCI for Alignment')                 #TODO
    p.add_argument('--mpi',  default=None, action='store_true', help='Compute MPI for Alignment')                 #TODO
    p.add_argument('--cons', default=None, action='store_true', help='Print consensus Structure to Alignment')    #TODO
    p.add_argument('--ribo', default=None, action='store_true', help='Use RIBOSUM-Scoring')                       #TODO

    #Output control
        #dealigned    #split into files    #extracted sequences    #desired format
    p.add_argument('--dealign', default=None, action='store_true', help='Remove gaps from sequences')
    p.add_argument('--split',   default=None, action='store_true', help='Save sequences in separate Files')
    p.add_argument('--extract', type=str, default=None, help='Extract Sequences with corresponding IDs (provide list whitespace-separated)')
    p.add_argument('--format',  type=str, default='clustal', help='Input alignment format (Default: clustal)')
    p.add_argument('--v', type=int, default=0, help='verbosity flag')

    #Alignment file
    p.add_argument('FILE', help='Alignment file')
    ARGS=p.parse_args()

    #Preprocess Arguments and perform idiot check
    aln = ReadAln(ARGS.FILE, ARGS.format)
    index_last = len(aln[0]) - 1
    ARGS.start -= 1
    if ARGS.end != None: ARGS.end -= 1

    if ARGS.start < 0:
        ARGS.start = 0
    elif ARGS.start > index_last:
        print('--start is greater than alignment length!')
        sys.exit()
    else:
        pass

    if ARGS.end == None:
        ARGS.end = index_last
    elif ARGS.end > index_last:
        ARGS.end = index_last
    elif ARGS.end < 0:
        ARGS.end = index_last
    elif ARGS.end < ARGS.start:
        print('--end must be greater than --start')
        sys.exit()

    output_options= {
        'dealign' : ARGS.dealign,   'split' : ARGS.split,
        'extract' : ARGS.extract,   'format' : ARGS.format,
        'sci' : ARGS.sci,           'mpi' : ARGS.mpi,
        'ribo': ARGS.ribo,          'cons': ARGS.cons,
        'name': ARGS.FILE
    }

    if ARGS.mode == 'whole':
        WriteOutput(aln, output_options)
    elif ARGS.mode == 'slice':
        aln_slice = Slice(aln, ARGS.start, ARGS.end)
        WriteOutput(aln_slice, output_options)
    elif ARGS.mode == 'window':
        WindowMode(aln, ARGS.start, ARGS.end, ARGS.win)
    else:
        print('No valid option selected!')

##### END MAIN #####

def WriteOutput(aln, output_options):
    #Alignment Postprocessing
    if output_options['extract']:
        #extract_ids = map(int, output_options['extract'].split())
        extract_ids = output_options['extract'].split()
        aln = ExtractSeqs(aln, extract_ids)
    if output_options['dealign']:
        aln = Dealign(aln)

    #Calculations
    if output_options['mpi']:
        mpi = MPI(aln)
    if output_options['sci']:
        sci = SCI(aln, output_options['ribo'])
    if output_options['cons']:
        consensus = Consensus(aln, output_options['ribo'])

    #Summary and Print
    if output_options['mpi']: print('MPI: '+str(mpi))
    if output_options['sci']: print('SCI: '+str(sci))
    if output_options['cons']: print('Consensus' + consensus)

    if output_options['dealign']:   #For dealigned seqs
        if output_options['split']: #In case you want it split up
            for rec in aln:
                SeqIO.write(rec, output_options['name']+'_'+rec.id, 'fasta')
        else:
            SeqIO.write(aln, output_options['name']+'_slice', 'fasta')

    else:
        #AlignIO.write(aln, 'NAME GOES HERE', 'stockholm')
        if output_options['split']:
            for rec in aln:
                SeqIO.write(rec, output_options['name']+'_'+rec.id, 'fasta')
        else:
            AlignIO.write(aln, output_options['name']+'_slice', output_options['format'])


def ReadAln(name, aln_format):
    '''Reads alignment file'''
    try:
        aln = AlignIO.read(name, aln_format)
    except ValueError as val_err:
        print 'Couldn\'t read alignment - wrong format (check spelling?)', val_err
        sys.exit()
    except IOError as ioerr:
        print ioerr
        sys.exit()
    return aln


def Slice(aln, start, end):
    '''Return Slice of an alignment''' #Exception handling?
    aln_slice = aln[:, start:end]
    return aln_slice


def Dealign(aln):
    '''Dealign Sequences (coded ugly as fuck)'''
    seqs = []
    for rec in aln:
        new_rec = SeqRecord(Seq(str(rec.seq).replace('-',''), IUPAC.ambiguous_rna), id=rec.id, name=rec.name, description='')
        seqs.append(new_rec)
    return seqs


def ExtractSeqs(aln, extract_ids):
    '''Extract Sequences from alignment and build new alignment'''
    aln_ids = {}
    extracted_recs = []

    for rec in aln:     #Sort into Dictionary
        aln_ids[rec.id] = rec
    for ID in extract_ids:
        if ID in aln_ids:
            extracted_recs.append(aln_ids[ID])
        else:
            print('ID '+str(ID)+' not found in alignment!')

    extracted_aln = MultipleSeqAlignment(extracted_recs)
    return extracted_aln


def WindowMode(aln, start, end, win):
    #Create a Dir and save all slices to this dir
    try:
        os.stat('SLICE')
    except:
        os.mkdir('SLICE')

    for i in range(start, end, win):
        aln_slice = Slice(aln, i, i+win)
        WriteOutput(aln_slice, output_options)
        AlignIO.write(aln_slice, 'SLICE/'+str(i)+'_'+str(i+win), 'clustal')


def SCI(aln, ribosum):
    '''Calculate the Structure Conservation Index as done in RNAz
    Args:
        aln - AlignIO alignment object
    Returns:
        sci - structure conservation index
    '''
    aln_seqs = []
    for rec in aln:
        aln_seqs.append(str(rec.seq))
        #if ARGS.v >= 1: print rec.seq

    if ribosum == True:
        RNA.cvar.ribo = 1
    else:
        RNA.cvar.ribo = 0

    #Call AliFold
    cons_str, cons_en = RNA.alifold(aln_seqs)
    #if ARGS.v >= 1: print(cons_str, cons_en)

    aln_e = 0
    #Call RNAfold for single Sequences
    for s in aln_seqs:
        s = s.replace('-','') #might be a problem for other gap chars, use regex :(
        if s == '': continue #Exclude pure gaps
        stru, en = RNA.fold(s)
        aln_e += en

    #sci = E_consensus / (1/N) * Sum(E_single seq)
    try:
        sci = float(cons_en) / ( float(aln_e) / len(aln_seqs) )
        #if ARGS.v >= 1:
        #    print 'Mean E: ',( float(aln_e) / len(aln_seqs) ), 'Consensus E: ', cons_en #For Debug
    except ZeroDivisionError:
        sci = 0
        #if ARGS.v >= 2: print 'Sum of all free energies = 0 (Caught ZeroDivisionError)'
    return sci


def MPI(aln, maxgap=0.25):
    '''Calculate the Mean Pairwise Identity of an Alignment
    Args:
        aln -  Alignment as Biopython AlignIO object
        maxgap - maximum fraction of gaps in a sequence
    Returns:
        mpi - mean Pairwise Identity
    '''
    #TODO - Exclude heavily gapped sequences, they fuck up everything
    mpi = 0
    mpi_list = []
    rows=len(aln)
    cols=len(aln[0])
    bad_columns = 0

    for col in range(cols):
        here = {}
        no_res = rows
        no_gaps = 0
        mpi_here = 0

        #Count Occurrence of each residue
        for row in range(rows):
            res = aln[row][col]
            if res == '-':
                no_gaps += 1
                no_res -= 1
            elif res not in here:
                here[res] = 1
            else:
                here[res] += 1

        if no_gaps >= no_res:
            bad_columns += 1

        #Sum MPI for each residue
        for i in here:
            mpi_here += ((float(here[i]) ** 2) / float(no_res))

        mpi_here /= float(no_res)
        mpi += mpi_here
        mpi_list.append(mpi_here)

    #print mpi_list
    mpi /= cols
    return mpi


def Consensus(aln, ribosum):
    '''Calculate the consensus MFE from an Biopython Alignment object'''
    if ribosum == True:
        RNA.cvar.ribo = 1
    else:
        RNA.cvar.ribo = 0

    aln_seqs = []
    for rec in aln:
        aln_seqs.append(str(rec.seq))

    cons_str, cons_en = RNA.alifold(aln_seqs)
    return (cons_str, cons_en)


def Realign(aln):
    import subprocess, tempfile
    #Remove aln characters
    recs=Dealign(aln)
    #Convert aln into string
    aln_str = ''
    for r in recs:
        aln_str += r.format('fasta')
    #Create Temporary File, create handle to virtual file
    with tempfile.TemporaryFile(mode='r+b') as TMP_handle:
        TMP_handle.write(aln_str)
        #Reset file index
        TMP_handle.seek(0)
        #Execute Process, Read stdout and stderr
        cmd = ['clustalo', '--in', '-', '--outfmt', 'clu', '--threads', '6']
        p = subprocess.Popen(cmd, stdin=TMP_handle, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        aln = AlignIO.read(p.stdout, 'clustal')
    return aln

if __name__ == '__main__':
    main()
