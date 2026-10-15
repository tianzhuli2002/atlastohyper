#!/usr/bin/env python3
"""Established Run-2 dilepton VyPER H5 converter entry point.

The implementation lives in ``TtbarDileptonvyper.py`` so its older lowercase
entry point remains available to existing scripts.
"""

from TtbarDileptonvyper import MakeDatasetTtbarDileptonvyper
from Common import argparser, main


if __name__ == "__main__":
    parser = argparser()
    parser.add_argument("-ds", "--dataset-size", type=int, help="Number of events to save")
    parser.add_argument("-fm", "--fully-matched", action="store_true", help="Keep fully matched events")
    parser.add_argument("-b", "--btag", action="store_true", help="Include btag features")
    args = parser.parse_args()
    main(MakeDatasetTtbarDileptonvyper, args)
