#!/usr/bin/env python
# -*- coding: utf-8 -*-

# this small helper script generates test matrices of the given dimension and
# outputs them to stdout so that you can pipe it into files. The elements are
# separated by '\t'.
#
# Synopsis:
#   ./gen_matrices.py 4 4
#


import sys
import random

try:
    dim_x, dim_y = int(sys.argv[1]), int(sys.argv[2])
except Exception as e:
    sys.stderr.write("error parsing matrix dimensions ...\n")
    raise


for row in range(0, dim_x):
    print("\t".join([ str(random.uniform(0, 9999)) for x in range(0, dim_y) ]))
