# make top from tpr
import sys

try:
    from TprReader import TprReader
except ImportError:
    sys.exit('Can not import TprReader')

def make_top_from_tpr(fname:str = 'md.tpr', topfile:str='md.top'):
    # get tpr handle
    rd = TprReader(fname)

    # get atomtype


    del rd

if __name__ == '__main__':
    make_top_from_tpr('../test/md.tpr', 'md.top')
