#! /usr/bin/env python

import sys
import subprocess
import pyfits
from pylab import *
import itertools
import os
from os.path import exists
from os import remove
from scipy.spatial import cKDTree
from scipy.optimize import fmin_tnc, fmin

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class ObjParams:
    def __init__(self, xCen, yCen, ellA, ellB, PA):
        self.cen = Point(xCen, yCen)
        self.ellA = ellA
        self.ellB = ellB
        self.PA = PA


def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)


def call_SE(fitsFile):
    # detect the name of SExtractor
    if subprocess.call("which sex >/dev/null", shell=True) == 0:
        callString = "sex %s " % fitsFile
    elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:
        callString = "sextractor %s " % fitsFile
    else:
        print "SExtractor was not found. Exiting..."
        exit(1)
    callString += "-c %s/default.sex" % (os.path.dirname(__file__))
    callString += " -VERBOSE_TYPE=QUIET"
    subprocess.call(callString, shell=True)


def get_SE_results():
    objects = []
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        n = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = 2*kron * float(params[4])
        ellB = 2*kron * float(params[5])
        PA = float(params[6])
        objects.append(ObjParams(xCen, yCen, ellA, ellB, PA))
    return objects


def create_mask(objects, xSize, ySize):
    mask = zeros((ySize, xSize))
    for obj in objects:
        cospa = cos(radians(obj.PA))
        sinpa = sin(radians(obj.PA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = obj.cen.x + obj.ellA * cose * cospa - obj.ellB * sine * sinpa
            y = obj.cen.y + obj.ellB * sine * cospa + obj.ellA * cose * sinpa
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
        focusR = (obj.ellA ** 2.0 - obj.ellB ** 2.0) ** 0.5
        focus10 = Point(obj.cen.x + focusR, obj.cen.y)  # Unrotated
        focus20 = Point(obj.cen.x - focusR, obj.cen.y)  #
        focus1 = rot_point(focus10, obj.cen, radians(obj.PA))
        focus2 = rot_point(focus20, obj.cen, radians(obj.PA))
        # Find pixels inside of the ellipse
        dEll = 2 * obj.ellA
        for x in xrange(xMin, xMax):
            for y in xrange(yMin, yMax):
                dFocus1 = hypot(x-focus1.x, y-focus1.y)
                dFocus2 = hypot(x-focus2.x, y-focus2.y)
                dPoint = dFocus1 + dFocus2
                if dPoint < dEll:
                    mask[y, x] = 1
    return mask


def buildDiagram(mask):
    maskedPixels = where(mask > 0)
    maskKDTree = cKDTree(data=zip(maskedPixels[0], maskedPixels[1]),
                         leafsize=100)
    diagram = zeros_like(mask)
    listOfAllCoords = list(itertools.product(range(mask.shape[0]), range(mask.shape[1])))
    distances, nearestPoints = maskKDTree.query(listOfAllCoords)
    diagram = reshape(distances, (mask.shape[0], mask.shape[1]))
    diagram = diagram / sum(diagram)
    diagramHDU = pyfits.PrimaryHDU(diagram)
    if exists("diagram.fits"):
        remove("diagram.fits")
    diagramHDU.writeto("diagram.fits")
    return diagram


def create_2d_poly_of_n_degree(n):
    """Function returns function that represents 2d polynomial
    of n degree"""
    def pk(coeffs, xy):
        r = coeffs[0]
        for k in xrange(1, n+1):
            r += coeffs[k] * xy[0]**k + coeffs[k+n] * xy[1]**k
        return r
    return pk


def fitness_function(data, poly, coeffs, diagram):
    xSize, ySize = data.shape
    yCen = int(ySize/2.0)
    xCen = int(xSize/2.0)
    coords = meshgrid(xrange(ySize), xrange(xSize))
    coords[0] = coords[0] - yCen
    coords[1] = coords[1] - xCen
    difference = sum(diagram * (data - poly(coeffs, coords))**2) / (xSize*ySize)
    return difference


def fit_sky_by_nd_poly(data, diagram, degree):
    meanSky = average(data, weights=diagram)
    print "Free coeff:", meanSky
    # if the degree of the polynomial is equal to zero
    # than return just weighted average of the sky level
    if degree == 0:
        return meanSky
    coeffs = [meanSky]
    # if the degree is not equal to zero.
    # this cycle fitts the data with polynomials of
    # increasing degree (from 1 to degree=degree). At each
    # inerration the results of previous iteration (k-1 degree poly)
    # are being used as initial guess. Coefficients for higher degrees
    # of x and y are set to be zero as initial guesses
    for k in xrange(1, degree+1):
        print "Fitting by %i degree polynomial" % k
        kdPoly = create_2d_poly_of_n_degree(k)
        fitFunc = lambda c: fitness_function(data, kdPoly, c, diagram)
        # unpack results of previous iteration
        freeCoeff = coeffs[0]
        xCoeffs = list(coeffs[1:k])
        yCoeffs = list(coeffs[k:])
        coeffs = None
        # new highest degree coefficients are zeros
        xCoeffs.append(0.0)
        yCoeffs.append(0.0)
        # assemble all coefficients are new initial guesses
        coeffs = [freeCoeff]
        coeffs.extend(xCoeffs)
        coeffs.extend(yCoeffs)
        # run gradient descent
        print "coeffs before fmin", coeffs, xCoeffs, yCoeffs
        res = fmin_tnc(fitFunc, x0=coeffs, approx_grad=True, maxfun=1000)
        coeffs = res[0]
        print "res:", coeffs
    return coeffs


def get_mean_sky(data, diagram):
    meanSky = average(data, weights=diagram)
    print meanSky
    return meanSky


def subtract_back(inFitsName, outFitsName, degree, coeffs):
    inHDU = pyfits.open(inFitsName)
    data = inHDU[0].data
    header = inHDU[0].header
    xSize, ySize = data.shape
    yCen = int(ySize/2.0)
    xCen = int(xSize/2.0)
    if degree == 0:
        outData = data - coeffs
    else:
        poly = create_2d_poly_of_n_degree(degree)
        coords = meshgrid(xrange(ySize), xrange(xSize))
        coords[0] = coords[0] - yCen
        coords[1] = coords[1] - xCen
        outData = data - poly(coeffs, coords)
    outHDU = pyfits.PrimaryHDU(outData, header=header)
    if exists(outFitsName):
        remove(outFitsName)
    outHDU.writeto(outFitsName)


def estimate_sky(fitsFileName, cleanName, degree):
    HDU = pyfits.open(fitsFileName)
    data = HDU[0].data
    ySize, xSize = data.shape
    call_SE(fitsFileName)
    objects = get_SE_results()
    mask = create_mask(objects, xSize, ySize)
    diagram = buildDiagram(mask)
    if degree == 0:
        sky = get_mean_sky(data, diagram)
    else:
        sky = fit_sky_by_nd_poly(data, diagram, degree)

    subtract_back(fitsFileName, cleanName, degree, sky)
    return sky


if __name__ == "__main__":
    sky = estimate_sky(sys.argv[1], "back_clean.fits", 1)
    print sky
