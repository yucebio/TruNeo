#!/usr/bin/env python3

import sys
import os

file_1=sys.argv[1]
file_2=sys.argv[2]

count1 = len(open(file_1,'rU').readlines())

count2 = len(open(file_2,'rU').readlines())

if count1 != count2:
    print("not the same line number")
    sys.exit(1)

