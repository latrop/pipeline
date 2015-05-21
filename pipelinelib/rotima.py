#! /usr/bin/env python

from math import *
import sys
import os

from ctypes import *

def rotPoint(x, y, x0, y0, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(x-x0) - sin(angle)*(y-y0) + x0
    y1 = sin(angle)*(x-x0) + cos(angle)*(y-y0) + y0
    return x1, y1


def crotima(name, outName, xOrig, yOrig, angle, layer=1):
    dllpath = "%s/crotima.so" % (os.path.dirname(__file__))
    crotlib = CDLL(dllpath)
    crotlib.crotima_crop.restype = c_int
    newCoordsOfOrigin = (c_double * 2)(0.0, 0.0)
    r = crotlib.crotima(name, outName, c_double(xOrig), c_double(yOrig),
                        c_double(-angle), c_int(layer), newCoordsOfOrigin)
    xNew = newCoordsOfOrigin[0]
    yNew = newCoordsOfOrigin[1]
    return xNew, yNew

def cstretchima(fitsName, outName, scale, workHDU=1):
    dllpath = "%s/crotima.so" % (os.path.dirname(__file__))
    crotlib = CDLL(dllpath)
    crotlib.stretchima.restype = c_int
    r = crotlib.stretchima(fitsName, outName, c_double(scale), c_int(workHDU))
    return r

def cdeproject(infile, outfile, xCen, yCen, incl, posang, workHDU):
    dllpath = "%s/crotima.so" % (os.path.dirname(__file__))
    crotlib = CDLL(dllpath)
    crotlib.stretchima.deproject = c_int
    r = crotlib.deproject(infile, outfile, c_double(xCen), c_double(yCen), c_double(incl), c_double(posang), c_int(workHDU))
    return r



if __name__ == "__main__":
    if len(sys.argv) < 6:
        print "Usage:"
        print "./rotima.py in out xCen yCen angle"
        exit()

    fileToRotate = sys.argv[1]
    outName = sys.argv[2]
    xOrig = float(sys.argv[3])
    yOrig = float(sys.argv[4])
    angle = float(sys.argv[5])
    coords = crotima(name=fileToRotate,
                     outName=outName,
                     xOrig=33.6,
                     yOrig=32.9,
                     angle=-62.0)

    print "New coordinates:", coords

