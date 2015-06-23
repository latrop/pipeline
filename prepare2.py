#! /usr/bin/env python

import sys
import subprocess
import os
import shutil
import argparse
from os.path import exists
from math import hypot, pi, log10
from numpy import ma, mean, std, zeros_like, where, arange, exp, zeros
from numpy import copy as npcopy
from numpy import sum as npsum
from pylab import plot, show, savefig, xlabel, ylabel, vlines, clf
import pyfits
from scipy.odr.odrpack import *

from pipelinelib.rotima import crotima
from pipelinelib.regions2mask import regions2mask, fits_mask
from pipelinelib.backEst import estimate_sky


def move(src, dst):
    if not os.path.isfile(src):
        print "File %s not found and cannot be moved" % (src)
        return
    shutil.copy(src, dst)
    os.remove(src)


def gauss(B, x):
    a = B[0]
    x0 = B[1]
    sigma = B[2]
    return a * exp(-(x-x0)**2 / (2*sigma**2))


def fit_by_gauss(xList, yList):
    aGuess = max(yList)
    x0Guess = (xList[0] + xList[-1]) / 2
    sigmaGuess = x0Guess - xList[0]  # FIXME
    fitting = ODR(RealData(xList, yList), 
                  Model(gauss),
                  [aGuess, x0Guess, sigmaGuess])
    fitting.set_job()
    result = fitting.run()
    B1 = (result.beta[0], result.beta[1], result.beta[2])
    # second interation
    dyLust = abs(yList-gauss(B1, xList)) ** 0.5
    fitting = ODR(RealData(xList, yList, sy=dyLust), 
                  Model(gauss),
                  B1)
    fitting.set_job()
    result = fitting.run()
    B2 = (result.beta[0], result.beta[1], result.beta[2])
    return B2


def call_SE(fitsFile):
    # detect the name of SExtractor
    if subprocess.call("which sex >/dev/null", shell=True) == 0:
        callString = "sex %s " % fitsFile
    elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:
        callString = "sextractor %s " % fitsFile
    else:
        print "SExtractor was not found. Exiting..."
        exit(1)
    callString += "-c ./pipelinelib/default.sex"
    callString += " -VERBOSE_TYPE=QUIET"
    subprocess.call(callString, shell=True)


def get_galaxy_params(fitsFile):
    """ Function finds object nearest to the center
    of the field."""
    hdu = pyfits.open(fitsFile)
    ySize, xSize = hdu[0].data.shape
    xCenField = xSize / 2.0
    yCenField = ySize / 2.0
    hdu.close()
    minArea = 100  # Minimal area of interested objects [pix^2]
    minCenterDist = 1e10
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        N = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = kron * float(params[4])
        ellB = kron * float(params[5])
        PA = float(params[6])
        ellArea = pi * ellA * ellB
        if ellArea > minArea:
            centerDist = hypot(xCen-xCenField, yCen-yCenField)
            if centerDist < minCenterDist:
                minCenterDist = centerDist
                galN = N
                galXCen = xCen
                galYCen = yCen
                galEllA = ellA
                galEllB = ellB
                galPA = PA
    return galN, galXCen, galYCen, galEllA, galEllB, galPA


def get_backgr_params(fitsFile):
    """Function finds the mean value and the sigma 
    of the backgound. Only pixels outside of all objects
    are being used (based on SE segmentation map)"""
    imageHDU = pyfits.open(fitsFile)
    imageData = imageHDU[0].data
    segmHDU = pyfits.open("segm.fits")
    segmData = segmHDU[0].data
    imageDataMasked = ma.masked_array(imageData, segmData)
    backMean = mean(imageDataMasked)
    backSTD = std(imageDataMasked)
    imageHDU.close()
    segmHDU.close()
    return backMean, backSTD


def mask_all_except_galaxy(fitsFile, nObj, backMean, backSTD):
    """ Function sets values of all background pixels
    and pixels of all objects except of galaxy equal to zero"""
    imageHDU = pyfits.open(fitsFile)
    imageData = imageHDU[0].data
    segmHDU = pyfits.open("segm.fits")
    segmData = segmHDU[0].data
    maskedData = npcopy(imageData)
    maskedData[where(maskedData < backMean + 3*backSTD)] = 0.0
    maskedData[where((segmData !=0) & (segmData != nObj))] = 0.0
    maskedHDU = pyfits.PrimaryHDU(data=maskedData)
    if exists("galaxy_only.fits"):
        os.remove("galaxy_only.fits")
    maskedHDU.writeto("galaxy_only.fits")


def get_pa(xCen, yCen, sePA, ellA, ellB):
    # Finds PA as an angle where biggest number of galaxy pixels are located at
    # a narrow horisontal line
    sliceWidth = min(ellA, ellB) / 2.0
    paList = []
    nPixList = []
    for curPA in arange(sePA-10.0, sePA+10.0, 0.25):
        xCenRot, yCenRot = crotima("galaxy_only.fits", "rot_tmp.fits", xCen, yCen, curPA)
        rotHDU = pyfits.open("rot_tmp.fits")
        rotData = rotHDU[0].data
        xSizeRot, ySizeRot = rotData.shape
        # Find number of galaxy pixels located inside of horisontal slice of given width
        galPixelsInSlice = 0
        upperEdge = int(xCenRot + sliceWidth)
        lowerEdge = int(xCenRot - sliceWidth)
        for sliceLoc in xrange(lowerEdge, upperEdge+1):
            row = rotData[sliceLoc]
            galPixelsInSlice += len(row[where(row != 0)])
        paList.append(curPA)
        nPixList.append(galPixelsInSlice)
        os.remove("rot_tmp.fits")
    params = fit_by_gauss(paList, nPixList)
    optPosAng = params[1]
    plot(paList, nPixList)
    plot(paList, gauss(params, paList))
    vlines(optPosAng, min(nPixList), max(nPixList)*1.1, linestyles=":")
    xlabel("posang")
    ylabel("pixels in slice")
    savefig("pa_gauss.png")
    clf()
    return optPosAng


def crop_fits(inFitsName, outFitsName, xCen, yCen, ellA, ellB):
    """ Functions crops image to size twice bigger than SE ellipse"""
    inHDU = pyfits.open(inFitsName)
    inData = inHDU[0].data
    ySize, xSize = inData.shape
    inHeader = inHDU[0].header
    if ellA < ellB:
        ellA, ellB = ellB, ellA
    xMin = max(1, xCen - ellA)
    xMax = min(xCen + ellA, xSize)
    yMin = max(1, yCen - ellB)
    yMax = min(yCen + ellB, ySize)
    outData = inData[yMin:yMax, xMin:xMax]
    ySize, xSize = outData.shape
    print "Center:", xSize/2.0, ySize/2.0
    outHDU = pyfits.PrimaryHDU(data=outData, header=inHeader)
    if exists(outFitsName):
        os.remove(outFitsName)
    outHDU.writeto(outFitsName)


def norm_fits(inFitsName, outFitsName):
    """ Normalise flux of a given fits file"""
    inHDU = pyfits.open(inFitsName)
    inData = inHDU[0].data
    inHeader = inHDU[0].header
    totalFlux = npsum(inData)
    outData = inData / totalFlux
    m0 = float(inHeader["m0"])
    m0Norm = m0 - 2.5 * log10(totalFlux)
    inHeader["m0"] = m0Norm
    outHDU = pyfits.PrimaryHDU(data=outData, header=inHeader)
    if exists(outFitsName):
        os.remove(outFitsName)
    outHDU.writeto(outFitsName)


def make_combined_fits(inFitsName, maskName, outFitsName):
    """ Creates fits file, where all masked pixels are equal to zero"""
    inHDU = pyfits.open(inFitsName)
    inData = inHDU[0].data
    inHeader = inHDU[0].header
    maskHDU = pyfits.open(maskName)
    maskData = maskHDU[0].data
    outData = npcopy(inData)
    outData[maskData>0] = 0.0
    outHDU = pyfits.PrimaryHDU(data=outData, header=inHeader)
    if exists(outFitsName):
        os.remove(outFitsName)
    outHDU.writeto(outFitsName)


def auto_mask(band):
    """Function creates ds9 *.reg file with regions that cover
    all objects except the galaxy"""
    call_SE("cropped_%s.fits" % (band))
    # find SE number of our galaxy (it is in the centre of the image now
    # so we can just get the intensity of the central pixel)
    segmHDU = pyfits.open("segm.fits")
    segmData = segmHDU[0].data
    ySize, xSize = segmData.shape
    xCen = int(xSize / 2.0)
    yCen = int(ySize / 2.0)
    galN = int(segmData[yCen, xCen])
    fout = open("mask_%s.reg" % (band), "w")
    fout.truncate(0)
    fout.write("# Region file format: DS9 version 4.1\n")
    fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
    fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
    fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
    fout.write('image\n')
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        n = int(params[0])
        if n == galN:
            continue
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = kron * float(params[4])
        ellB = kron * float(params[5])
        PA = float(params[6])
        fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (xCen, yCen,
                                                                 ellA, ellB,
                                                                 PA))
    fout.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare images for decomposition")
    parser.add_argument("objName", help="Object name")
    parser.add_argument("filters", help="List of filters (for example gri or uiz)")
    parser.add_argument("--combined", action="store_true", default=False,
                        help="Create combined images with mask and data")
    parser.add_argument("--norm", action="store_true", default=False,
                        help="Create normalised images")
    args = parser.parse_args()

    objName = args.objName
    bandList = args.filters
    # First step: background estimation for every passband
    for band in bandList:
        inFile = "%s_%s_trim.fits" % (objName, band)
        outFile = "back_clean_%s.fits" % (band)
        estimate_sky(inFile, outFile, 1)
        os.remove("segm.fits")
        os.remove("field.cat")
        os.remove("diagram.fits")
    # Second step: launching of the SE for background subtracted images
    # to obtain some geometric parameters of the galaxy
    galNAll = {}
    xCenAll = {}
    yCenAll = {}
    ellAAll = {}
    ellBAll = {}
    ellPAAll = {}
    for band in "irgzu": # i-band is more sutable to find PA, r is the next and so on
        if band in bandList:
            print "Position angle will be determined in %s band" % (band)
            fitsFile = "back_clean_%s.fits" % (band)
            call_SE(fitsFile)
            galN, xCen, yCen, ellA, ellB, ellPA = get_galaxy_params(fitsFile)
            galNAll[band] = galN
            xCenAll[band] = xCen
            yCenAll[band] = yCen
            ellAAll[band] = ellA
            ellBAll[band] = ellB
            ellPAAll[band] = ellPA
            # Now we want to estimate position abgle in chosen passband
            backMean, backSTD = get_backgr_params(fitsFile)
            mask_all_except_galaxy(fitsFile, galN, backMean, backSTD)
            optPosAng = get_pa(xCen, yCen, ellPA, ellA, ellB)
            print "Position angle = %1.1f" % (optPosAng)
            break
    # Rotate every image using the mean position angle
    xCenRot = {}
    yCenRot = {}
    meanXCen = mean(xCenAll.values())
    meanYCen = mean(yCenAll.values())
    for band in bandList:
        fitsName = "back_clean_%s.fits" % (band)
        outFitsName = "rotated_%s.fits" % (band)
        xCenRot[band], yCenRot[band] = crotima(fitsName, outFitsName,
                                               meanXCen, meanYCen,
                                               optPosAng)
    # Now lets show one of rotated images, so user can see if rotation was ok
    # and make some corrections if it is nessesary
    if "r" in bandList:
        # r-band file is usually more deep, so it is better to show it
        fileToShow = "rotated_r.fits"
    elif "g" in bandList:
        # g-band image is the second in the row
        fileToShow = "rotated_g.fits"
    elif "i" in bandList:
        fileToShow = "rotated_i.fits"
    else:
        # if there are not g, i, or r bands in processing list,
        # pick up just first band in the list
        fileToShow = "rotated_%s.fits" % (bandList[0])
    ds9Proc = subprocess.Popen(["ds9", str(fileToShow),
                                "-scale", "log"])
    ds9Proc.wait()
    inputString = raw_input("Correction to position angle: ").strip()
    if inputString:
        deltaPosAng = float(inputString)
    else:
        deltaPosAng = 0.0
    # if the correction to position angle is not equal to zero
    # perform rotation once again with new position angle value
    for band in bandList:
        fitsName = "back_clean_%s.fits" % (band)
        outFitsName = "rotated_%s.fits" % (band)
        os.remove(outFitsName)
        xCenRot[band], yCenRot[band] = crotima(fitsName, outFitsName,
                                               meanXCen, meanYCen,
                                               optPosAng-deltaPosAng)

    # Croping fields goes as follows: using ellipse parameters of galaxy
    # script creates ds9 reg-file with box to show region to crop.
    # This box is being shown to user, so he can change it if nessesary
    # Parameters can be somewhat different in different filters
    # so we will use maximal ones as box values.
    # Box will be shown overlapped on sum of normalised images in all filters
    # So let us begin:
    # 1) Creation of the reg-file
    maxEllA = max(ellAAll.values())
    maxEllB = max(ellBAll.values())
    meanXCenRot = mean(xCenRot.values())
    meanYCenRot = mean(yCenRot.values())
    fout = open("crop_box.reg", "w")
    fout.write("box(%1.1f,%1.1f,%1.1f,%1.1f,0)\n" % (meanXCenRot, meanYCenRot,
                                                     4.0*maxEllA, 4.0*maxEllB))
    fout.close()

    # 2) Creation of sum of normalised images
    dataList = []
    for band in bandList:
        norm_fits("rotated_%s.fits"%(band), "norm_tmp_%s.fits"%(band))
        normHDU = pyfits.open("norm_tmp_%s.fits"%(band))
        dataList.append(normHDU[0].data)
        normHDU.close()
        os.remove("norm_tmp_%s.fits"%(band))
    sumNormData = zeros(dataList[0].shape)
    for data in dataList:
        sumNormData += data
    outHDU = pyfits.PrimaryHDU(data=sumNormData)
    if exists("sum_norm.fits"):
        os.remove("sum_norm.fits")
    outHDU.writeto("sum_norm.fits")

    # 3) Show box to user
    ds9Proc = subprocess.Popen(["ds9", "sum_norm.fits",
                                "-regions", "crop_box.reg",
                                "-scale", "log"])
    ds9Proc.wait()

    # 4) Get data from box
    for line in open("crop_box.reg"):
        if "box" in line:
            params = line[4:-2].split(",")
            cropWidth = float(params[2]) / 2.0
            cropHeight = float(params[3]) / 2.0
    # 5) Finally crop images
    for band in bandList:
        inFile = "rotated_%s.fits" % (band)
        outFile = "cropped_%s.fits" % (band)
        crop_fits(inFile, outFile,
                  meanXCenRot, meanYCenRot,
                  cropWidth, cropHeight)
    os.remove("crop_box.reg")
    os.remove("sum_norm.fits")
    os.remove("field.cat")

    # Create mask. At first, individual mask for every passband.
    # Then, merge all masks into one, so images in all filters
    # will be processed with the same mask
    for band in bandList:
        # Creation of automatic mask based on SExtractor results
        auto_mask(band)
    for i, band in enumerate(bandList):
        # show automatical generated mask to the user
        # so he can fix it if nessesary
        print "Check the mask for %s band" % (band)
        ds9Proc = subprocess.Popen(["ds9", "cropped_%s.fits"%(band),
                                    "-regions", "mask_%s.reg"%(band),
                                    "-scale", "log"])
        ds9Proc.wait()
        if i <= (len(bandList)-2):
            shutil.copy("mask_%s.reg"%(band), "mask_%s.reg"%(bandList[i+1]))
        regions2mask("mask_%s.reg" % (band), "mask_%s.dat" % (band), "cropped_%s.fits" % (band))
        fits_mask("mask_%s.dat" % (band), "cropped_%s.fits" % (band), "mask_%s.fits" % (band))
    # merge all mask in one
    masterMaskData = None
    for band in bandList:
        bandMaskHDU = pyfits.open("mask_%s.fits" % (band))
        bandMaskData = bandMaskHDU[0].data
        if masterMaskData is None:
            masterMaskData = zeros_like(bandMaskData)
        masterMaskData = masterMaskData + bandMaskData
    masterMaskHDU = pyfits.PrimaryHDU(data=masterMaskData)
    if exists("mask.fits"):
        os.remove("mask.fits")
    masterMaskHDU.writeto("mask.fits")

    # Create normalised fits files
    if args.norm:
        for band in bandList:
            norm_fits("cropped_%s.fits" % (band), "norm_%s.fits" % (band))

    # Create combined fits files
    if args.combined:
        for band in bandList:
            make_combined_fits("cropped_%s.fits" % (band), "mask.fits", "combined_%s.fits" % (band))

    os.remove("segm.fits")
    os.remove("galaxy_only.fits")
    os.remove("field.cat")
    if not exists("./%s_results/" % (objName)):
        os.makedirs("./%s_results" % (objName))
    if not exists("./%s_results/intermediate" % (objName)):
        os.makedirs("./%s_results/intermediate" % (objName))
    move("mask.fits", "./%s_results/mask.fits" % (objName))
    move("pa_gauss.png", "./%s_results/intermediate/pa_gauss.png" % (objName))
    for band in bandList:
        move("mask_%s.dat" % (band), "./%s_results/intermediate/mask_%s.dat" % (objName, band))
        move("rotated_%s.fits" % (band), "./%s_results/intermediate/rotated_%s.fits" % (objName, band))
        move("cropped_%s.fits" % (band), "./%s_results/cropped_%s.fits" % (objName, band))
        move("mask_%s.reg" % (band), "./%s_results/intermediate/mask_%s.reg" % (objName, band))
        move("mask_%s.fits" % (band), "./%s_results/intermediate/mask_%s.fits" % (objName, band))
        move("back_clean_%s.fits" % (band), "./%s_results/intermediate/back_clean_%s.fits" % (objName, band))
        if args.norm:
            move("norm_%s.fits" % (band), "./%s_results/norm_%s.fits" % (objName, band))
        if args.combined:
            move("combined_%s.fits" % (band), "./%s_results/combined_%s.fits" % (objName, band))
    shutil.make_archive("./%s_results/intermediate" % (objName),
                        "bztar",
                        base_dir="./%s_results/intermediate" % (objName))
    shutil.rmtree("%s_results/intermediate" % (objName))
