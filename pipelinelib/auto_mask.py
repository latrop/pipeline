#! /usr/bin/env python

from math import pi
import subprocess
import sys

galFileName = sys.argv[1]
galName = galFileName.split(".")[0]
callString = "sex %s " % galFileName
callString += "-c ./lib/default.sex"

subprocess.call(callString, shell=True)

objCat = []
for line in open("field.cat"):
    if line.startswith("#"):
        continue
    params = line.split()
    xCen = float(params[1])
    yCen = float(params[2])
    kron = float(params[8])
    ellA = kron * float(params[4])
    ellB = kron * float(params[5])
    PA = float(params[6])
    ellArea = pi * ellA * ellB
    objCat.append((xCen, yCen, ellA, ellB, PA, ellArea))

# we interested in all objects except the main galaxy.
# let us find it by its area (it supposed to be biggest
# object on image)
objCat.sort(key=lambda x: x[-1])
objCat.pop()  # pop entry with galaxy


fout = open("%s.reg" % galName, "w")
fout.truncate(0)
fout.write("# Region file format: DS9 version 4.1\n")
fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
fout.write('image\n')
for params in objCat:
    xCen = params[0]
    yCen = params[1]
    ellA = params[2]
    ellB = params[3]
    PA = params[4]
    fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (xCen, yCen,
                                                             ellA, ellB,
                                                             PA))
fout.close()
