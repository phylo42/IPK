#!/usr/bin/env python

# This is an example of usage pyxpas.
# An input database is supposed to be built by rappas2.


__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import sys
# Make sure that xpas.so is visible by this script
import xpas


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python xpas.py DATABASE_FILE")
    else:
        db_filename = sys.argv[1]
        db = xpas.load(db_filename)
        print(db.size())

        # See core::encode_kmer
        kmer_value = 45

        found = db.search(kmer_value)
        for entry in found.values:
            print(entry.branch, entry.score)
