#! /usr/bin/env python

import os
import sys
import pyfits


def crop2box(regFile, fitsName, outFitsName=None):
    for line in open(regFile):
        if "box" in line:
            params = line[4:-2].split(",")
            xCen = float(params[0])
            yCen = float(params[1])
            width = float(params[2])
            height = float(params[3])
            PA = float(params[4])
            if PA != 0.0:
                print "Warning (crop2box) ",
                print "position angle of box is not equal to zero"
            break

    xMin = int(round(xCen - width/2.0))
    xMax = int(round(xCen + width/2.0))
    yMin = int(round(yCen - height/2.0))
    yMax = int(round(yCen + height/2.0))

    if outFitsName is None:
        outFitsName = fitsName.split(".")[0] + "_crop.fits"
    hdu = pyfits.open(fitsName)
    data = hdu[0].data
    header = hdu[0].header

    outData = data[yMin:yMax, xMin:xMax]
    outHDU = pyfits.PrimaryHDU(data=outData, header=header)
    if os.path.exists(outFitsName):
        os.remove(outFitsName)
    outHDU.writeto(outFitsName)


if __name__ == "__main__":
    regFile = sys.argv[1]
    fitsName = sys.argv[2]
    if len(sys.argv) > 3:
        outFitsName = sys.argv[3]
    else:
        outFitsName = None
    crop2box(regFile, fitsName, outFitsName)
