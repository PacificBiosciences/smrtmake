#!/usr/bin/env python

import sys
from pbcore.io import BaxH5Reader

def main():
    """Sorts and prints a list of subread lengths from a base file."""
    basfofn = sys.argv[1]
    rgnfofn = sys.argv[2]

    basfiles = open(basfofn).read().splitlines()
    rgnfiles = open(rgnfofn).read().splitlines()

    for basfile, rgnfile in zip(basfiles, rgnfiles):
        reader = BaxH5Reader(basfile, rgnfile)
        for l in sorted([len(x) for x in reader.subreads()], reverse=True):
            print l

if __name__ == '__main__':
    sys.exit(main())
