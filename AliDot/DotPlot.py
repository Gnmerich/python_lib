import re

class DotPlot:

    def __init__(self, filename=None):
        self.name   = None
        self.file   = None
        self.seq    = None
        self.ubox   = None
        self.lbox   = None


    def _parse(self):
        with open(self.file, 'r') as FI:
            for line in FI:
                line = line.rstrip('\n')

        pass

    def _parseFromMatrix(self):
        pass
