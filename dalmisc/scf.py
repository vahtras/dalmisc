import itertools
import sys

from .scf_iter import DaltonOutputIterator

files = sys.argv[1:]

iterators = [DaltonOutputIterator(dalton_output=open(f)) for f in files]

for i, egs in enumerate(itertools.zip_longest(*iterators, fillvalue=(0.0, 0.0))):
    print(f'{i:2d}:' + "  ".join(f"{e:14.9f} {g:.5e}" for e, g in egs))
    # print(i, egs)
