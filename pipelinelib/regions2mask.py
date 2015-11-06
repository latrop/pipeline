#! /usr/bin/env python

import sys
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import pyfits
from os.path import exists
from os import remove



class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)


def triangle_area(a, b, c):
    return ((c.x*b.y - b.x*c.y) - (c.x*a.y - a.x*c.y) + (b.x*a.y - a.x*b.y))/2


def circle_mask(cen, r, maskFile, xSize, ySize):
    xMin = max(0, cen.x-r)
    xMax = min(xSize, cen.x+r)
    yMin = max(0, cen.y-r)
    yMax = min(ySize, cen.y+r)
    for x in xrange(xMin, xMax):
        for y in xrange(yMin, yMax):
            d = hypot(x-cen.x, y-cen.y)
            if (d < r):
                maskFile.write("%i %i\n" % (x, y))


def ellipse_mask(cen, ellA, ellB, ellPA, maskFile, xSize, ySize):
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
        focus20 = Point(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        for x in xrange(xMin, xMax+1):
            for y in xrange(yMin, yMax+1):
                dFocus1 = hypot(x-focus1.x, y-focus1.y)
                dFocus2 = hypot(x-focus2.x, y-focus2.y)
                dPoint = dFocus1 + dFocus2
                if dPoint < dEll:
                    maskFile.write("%i %i\n" % (x, y))


def box_mask(cen, width, height, PA, maskFile, xSize, ySize):
    a = Point(cen.x - width/2.0, cen.y - height/2.0)
    b = Point(cen.x - width/2.0, cen.y + height/2.0)
    c = Point(cen.x + width/2.0, cen.y + height/2.0)
    d = Point(cen.x + width/2.0, cen.y - height/2.0)
    # Rotate points
    a = rot_point(a, cen, radians(PA))
    b = rot_point(b, cen, radians(PA))
    c = rot_point(c, cen, radians(PA))
    d = rot_point(d, cen, radians(PA))
    xmin = max(0, int(min(a.x, b.x, c.x, d.x)))
    xmax = min(xSize, int(max(a.x, b.x, c.x, d.x)))
    ymin = max(0, int(min(a.y, b.y, c.y, d.y)))
    ymax = min(ySize, int(max(a.y, b.y, c.y, d.y)))
    for x in xrange(xmin, xmax+1):
        for y in xrange(ymin, ymax+1):
            p = Point(x, y)
            t1Area = triangle_area(p, a, b)
            t2Area = triangle_area(p, b, c)
            t3Area = triangle_area(p, c, d)
            t4Area = triangle_area(p, d, a)
            sumOfSigns = (sign(t1Area) + sign(t2Area) +
                          sign(t3Area) + sign(t4Area))
            if abs(sumOfSigns) == 4:
                maskFile.write("%i %i\n" % (x, y))


def regions2mask(regFileName, maskFileName, fitsFileName):
    regFile = open(regFileName)
    maskFile = open(maskFileName, "w")
    maskFile.truncate(0)
    HDU = pyfits.open(fitsFileName)
    ySize, xSize = HDU[0].data.shape
    for line in regFile:
        if "circle" in line:
            params = line[7:-2].split(",")
            cen = Point(int(round(float(params[0]))),
                        int(round(float(params[1]))))
            r = int(round(float(params[2])))
            circle_mask(cen, r, maskFile, xSize, ySize)

        if "ellipse" in line:
            params = line[8:-2].split(",")
            cen = Point(float(params[0]),
                        float(params[1]))
            ellA = float(params[2])
            ellB = float(params[3])
            ellPA = float(params[4])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90
            ellipse_mask(cen, ellA, ellB, ellPA, maskFile, xSize, ySize)

        if "box" in line:
            params = line[4:-2].split(",")
            cen = Point(float(params[0]),
                        float(params[1]))
            width = float(params[2])
            height = float(params[3])
            PA = float(params[4])
            box_mask(cen, width, height, PA, maskFile, xSize, ySize)

    maskFile.close()


def fits_mask(maskASCIIFile, fitsDataFile, fitsMaskFile):
    HDU = pyfits.open(fitsDataFile)
    data = HDU[0].data
    mask = zeros_like(data)
    for line in open(maskASCIIFile):
        x = int(line.split()[0])
        y = int(line.split()[1])
        mask[y-1, x-1] = 1
    maskHDU = pyfits.PrimaryHDU(data=mask)
    if exists(fitsMaskFile):
        remove(fitsMaskFile)
    maskHDU.writeto(fitsMaskFile)
    

if __name__ == "__main__":
    regFileName = "dust.reg"   # "./results/mask.reg"  
    maskFileName = "dust.dat"  # "./results/mask.dat" 
    fitsFileName = sys.argv[1]  #"./results/cropped.fits" 
    regions2mask(regFileName, maskFileName, fitsFileName)
    fits_mask(maskFileName, fitsFileName, "dust_mask.fits")
