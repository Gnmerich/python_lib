import re

class DotPlot:
    """DotPlot - Python Class for Vienna RNApackage Dotplot-related stuff
    Attributes:
        name(str):  Name, typically the FASTA-header
        seq(str):   RNA Sequence (A,C,G,U)
        type(str):  Program used for Calculation, i.e. [RNAfold/RNAplfold]
        lbox(list): lbox-part of the Dotplot (...)
        ubox(list): ubox-part of the Dotplot (...)
    Methods:
        ParseDP: Parse the Dotplot-file given the provided file name
        Reparse: Parse the file again, using the name of the initial parsing step
        MfePairs: returns a list of all pairs in the MFE structure
        PfPairs: returns a list of all possible pairs in the ensemble given a minimum probability cutoff

        compare
        parse
    """

    def __init__(self, filename, type):
        self.__file = filename
        self.seq = ''
        self.type = type
        self.lbox = []
        self.ubox = []
        self.cutoff = 0.5    #dirty hack
        self.ParseDP(filename)

    def ParseDP(self, filename):
	'''ParseDP - parses the output of RNAplfold/RNAfold -p
        Args:
            filename
	Returns:
            ubox, lbox; Both list of tuples (i,j,prob)
        '''
        ubox = []
        lbox = []
        sequence = ''

  	with open(filename) as DP:
		switch = False

		for line in DP:
			line = line.rstrip('\n')

			if line == '%start of base pair probability data':
                            switch = 'DP'
                        elif line == '/sequence { (\\':
                            switch = 'SEQ'
                        elif line == ') } def':
                            switch = False
                        elif line == 'showpage':
				break

                        elif switch == 'SEQ': #Parse Sequence
                            line = line.rstrip('\\')
                            sequence = sequence + line

                        elif switch == 'DP':
                            i,j,p,box = line.split() #TODO catch ValueError

                            if box == 'ubox':
                                ubox.append((i,j,p))
                            elif box == 'lbox':
                                lbox.append((i,j,p))
                            else:
                                raise ValueError('Entry not lbox or ubox!')

        self.seq=sequence
        self.ubox=ubox
        self.lbox=lbox
        #TODO parse Header, Dotplot-type (local, global..)

    def Reparse(self):
        self.ParseDP(self.__file)

    def MfePairs(self):
        '''Return basepairs in mfe, error if Lfold'''
        if self.type == 'RNAplfold':
            raise ValueError
        return self.lbox

    def DotBracket(self):
        '''Return DotBracket notation of MFE structure'''
        db = ['.'] * len(self.seq)
        for bp in self.lbox:
            db[bp[0]-1] = '('
            db[bp[1]-1] = ')'
        return ''.join(db)

    def PfPairs(self):
        '''PfPairs - return most probable base pairs
        Args:
            cutoff(float) - minimum basepair probability
        Returns:
            pairs_filtered
        '''
        pairs ={}

        for bp in self.ubox:
            i,j,p=bp
            if i not in pairs:
                pairs[i]=[]
                pairs[i].append((j,p))
            else:
                pairs[i].append((j,p))

        #Filter for most probable bp > cutoff
        pairs_filtered = []
        for i in pairs:
            n,p = 0,0
            for bp in pairs[i]:
                if bp[1] > p:
                    n,p = bp
            #Remember this pair if p < cutoff
            if float(p) >= self.cutoff:
                pairs_filtered.append((i,n,p))

        return pairs_filtered


    def PrintSeq(self):
        print self.seq

    def PairingDistances(self, type):
        '''PairingDistances - Returns the span of all base pairs
        Returns:
            dist(list): list with N entries, for each bp(i,j) list[i] = abs(j-i)
        '''
        dist = [0]*len(self.seq)

        if type == 'rnafold':
            for bp in self.lbox:
                i, j, prob = bp
                i = int(i) - 1
                j = int(j) - 1

                dist[i] = abs(j-i)
                dist[j] = abs(j-i)
            return dist


        elif type == 'plfold':
            #Get the most probable Pair for each position
            lbox_mostprob = self.PfPairs()

            for bp in lbox_mostprob:
                i, j, prob = bp
                i = int(i) - 1
                j = int(j) - 1

                dist[i] = abs(j-i)
                dist[j] = abs(j-i)

            return dist

        else:
            pass

    def __sub__(self, rightDP):
        '''Compare the MFE structures of 2 Dotplot Objects
        Args:
            self - left DotPlot object
            rightDP - right DotPlot object
        Returns:
            shared - Base Pairs in both structures
            leftonly - Base Pairs just in the 'left' structure
            rightonly - Base Pairs just in the right structure
        '''
        bps_left = []
        bps_right = []

        for bp in self.lbox:
            bps_left.append((bp[0],bp[1]))

        for bp in rightDP.lbox:
            bps_right.append((bp[0],bp[1]))

        #Compare both sets
        leftonly = []
        shared = []
        rightonly = []

        for bp in bps_left:
            if bp in bps_right:
                shared.append(bp)
                bps_right.remove(bp) #remove so in the end just non-shared bps are left
            else:
                leftonly.append(bp)

        rightonly = bps_right
        return leftonly, shared, rightonly


def ComparePF(dp1, dp2):
    '''ComparePF - check the number of basepairs shared by 2 PF
    Args:
        dp1, dp2 - DotPlot Objects
        cutoff(float) - matching basepair has to be the most probable
    Returns:
        shared_bps
        missed_bps
    '''
    bps_left = []
    bps_right = []

    for bp in dp1.PfPairs():
        bps_left.append((bp[0],bp[1]))

    for bp in dp2.PfPairs():
        bps_right.append((bp[0],bp[1]))

    #Compare both sets
    leftonly = []
    shared = []
    rightonly = []

    for bp in bps_left:
        if bp in bps_right:
            shared.append(bp)
            bps_right.remove(bp) #remove so in the end just non-shared bps are left
        else:
            leftonly.append(bp)

    rightonly = bps_right
    return leftonly, shared, rightonly
