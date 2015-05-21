#include <string.h>
#include <stdio.h>
#include "math.h"
#include "fitsio.h"
#include "unistd.h"

#define PI 3.14159265

int crotima(char fitsName[], char outName[], double xOrig, double yOrig,
	    double theta, int workHDU, double newCoordsOfOrigin[]);

typedef struct {
  double x;
  double y;
} Point;

Point rotPoint(Point inPoint, Point centPoint, double angle)
{
  Point outPoint;
  double x1, y1;
  x1 = cos(angle)*(inPoint.x-centPoint.x) - sin(angle)*(inPoint.y-centPoint.y) + centPoint.x;
  y1 = sin(angle)*(inPoint.x-centPoint.x) + cos(angle)*(inPoint.y-centPoint.y) + centPoint.y;
  outPoint.x = x1;
  outPoint.y = y1;
  return outPoint;
} 

int crotima(char fitsName[], char outName[], double xOrig, 
	    double yOrig, double theta, int workHDU,
	    double newCoordsOfOrigin[])
{
  fitsfile *fptr, *outima;     
  long fpixel[2] = {1, 1};
  int i, j;
  int status = 0, numHDUs, hdutype=0, bitpix=0;
  int nkeys = 0;
  long naxes[2] = {0, 0};
  long naxesOut[2] = {0, 0};
  long xSize, ySize, outDimX, outDimY, imsize, sizeXout, sizeYout;
  double angle=PI*theta/180.0;
  long minX, maxX, minY, maxY;
  double costh, sinth;
  Point p1, p2, p3, p4;
  Point p1r, p2r, p3r, p4r;
  Point orig;
  orig.x=xOrig;
  orig.y=yOrig;
  double *inDataArr;
  double *outDataArr;
  double xPr, yPr, x, y, iDest;
  int xInd, yInd;
  char record[FLEN_CARD];
  double keyValue;

  // Open file
  fits_open_file(&fptr, fitsName, READONLY, &status);

  // Let's check if work hdu is exists
  fits_get_num_hdus(fptr, &numHDUs, &status);
  if (workHDU > numHDUs) return 1;

  // Set work HDU
  fits_movabs_hdu(fptr, workHDU, &hdutype, &status);

  // Get information about size of image and type of data
  fits_get_img_size(fptr, 2,  naxes, &status);
  fits_get_img_type(fptr, &bitpix, &status);
  xSize = naxes[0];
  ySize = naxes[1];

  // Read data from input file
  imsize = xSize*ySize;
  inDataArr = (double *) malloc(sizeof(double) * imsize);
  fits_read_pix(fptr, TDOUBLE, fpixel, imsize, NULL, inDataArr, NULL, &status);

  // Compute size of output image
  p1.x=0; p1.y=0;
  p2.x=xSize; p2.y=0;
  p3.x=xSize; p3.y=ySize;
  p4.x=0; p4.y=ySize;
  p1r = rotPoint(p1, orig, angle);
  p2r = rotPoint(p2, orig, angle);
  p3r = rotPoint(p3, orig, angle);
  p4r = rotPoint(p4, orig, angle);
  minX = (long)fmin(fmin(p1r.x, p2r.x), fmin(p3r.x, p4r.x));
  maxX = (long)fmax(fmax(p1r.x, p2r.x), fmax(p3r.x, p4r.x));
  minY = (long)fmin(fmin(p1r.y, p2r.y), fmin(p3r.y, p4r.y));
  maxY = (long)fmax(fmax(p1r.y, p2r.y), fmax(p3r.y, p4r.y));
  sizeXout = (maxX - minX);
  sizeYout = (maxY - minY);
  naxesOut[0] = sizeXout;  
  naxesOut[1] = sizeYout;
  newCoordsOfOrigin[0] = orig.x-minX;
  newCoordsOfOrigin[1] = orig.y-minY;
  // Create out image
  // Check if file with given name is already exists (and delete it if so)
  if (access(outName, F_OK) != -1 ) {
    // file exists
    printf("Warning (crotima): fits file with name \'%s\' is already exists. Removing it...\n\n", outName);
    unlink(outName);
  }
  fits_create_file(&outima, outName, &status);

  // Create new image unit in out fits file
  fits_create_img(outima, bitpix, 2, naxesOut, &status);
  outDataArr = (double *) calloc(sizeof(double), (sizeXout*sizeYout));
  // Perform rotation
  costh = cos(-angle);
  sinth = sin(-angle);
  for (i=minX; i<maxX; i++) {
    for (j=minY; j<maxY; j++) {
      xPr = costh*(i-orig.x) - sinth*(j-orig.y) + orig.x;
      yPr = sinth*(i-orig.x) + costh*(j-orig.y) + orig.y;
      xInd = (int)(floor(xPr));
      yInd = (int)(floor(yPr));
      if ((xInd > 0) && (xInd < xSize-1) && (yInd > 0) && (yInd < ySize-1)){
	x = xPr - floor(xPr);
	y = yPr - floor(yPr);
	iDest = ((1-x)*(1-y) * inDataArr[xInd + yInd*xSize]
		 +x*(1-y) * inDataArr[xInd+1 + yInd*xSize]
		 +y*(1-x)*inDataArr[xInd + (yInd+1)*xSize]
		 +x*y*inDataArr[xInd+1 + (yInd+1)*xSize]);
	outDataArr[i-minX + (j-minY)*sizeXout] = iDest;
      }
    }
  }

  // Flush rotated data into file
  fits_write_pixnull(outima, TDOUBLE, fpixel, sizeXout*sizeYout, outDataArr, NULL, &status);
  
  // Now we need to copy important header keywords from old (unrotated file)
  // to the rotated one
  const char * keywords[3];
  keywords[0] = "GAIN";
  keywords[1] = "M0";
  keywords[2] = "READOUT";
  for (i=0; i<=2; i++){
    status = 0;
    fits_read_key(fptr, TDOUBLE, keywords[i], &keyValue, NULL, &status);
    if (status == 0){
      fits_update_key(outima, TDOUBLE, keywords[i], &keyValue, NULL, &status);
    }
  }
  status = 0;
  fits_close_file(outima, &status);
  status = 0;
  fits_close_file(fptr, &status);
  free(inDataArr);
  free(outDataArr);
  return(0);
}

int crotima_crop(char fitsName[], char outName[], double xOrig, double yOrig, double theta, int workHDU)
{
  fitsfile *fptr, *outima, *outCropIma;
  long fpixel[2] = {1, 1};
  int i, j;
  int status = 0, numHDUs, hdutype=0, bitpix=0;
  long naxes[2] = {0, 0};
  long naxesOut[2] = {0, 0};
  long xSize, ySize, outDimX, outDimY, imsize, sizeXout, sizeYout;
  double angle=PI*theta/180.0;
  long minX, maxX, minY, maxY;
  double costh, sinth;
  Point p1, p2, p3, p4;
  Point p1r, p2r, p3r, p4r;
  Point orig;
  orig.x=xOrig;
  orig.y=yOrig;
  double *inDataArr;
  double *outDataArr;
  double *outDataArrCrop;
  double sizeXoutCrop, sizeYoutCrop;
  double xPr, yPr, x, y, iDest;
  int xInd, yInd;
  Point cenNew;
  long fCropPixel[2];
  long lCropPixel[2];
  long inc[2] = {1, 1};
  int anynul;
  double *outCropArray;

  // Open file
  fits_open_file(&fptr, fitsName, READONLY, &status);

  // Let's check if work hdu is exists
  fits_get_num_hdus(fptr, &numHDUs, &status);
  if (workHDU > numHDUs) return 1;

  // Set work HDU
  fits_movabs_hdu(fptr, workHDU, &hdutype, &status);

  // Get information about size of image and type of data
  fits_get_img_size(fptr, 2,  naxes, &status);
  fits_get_img_type(fptr, &bitpix, &status);
  xSize = naxes[0];
  ySize = naxes[1];

  // Read data from input file
  imsize = xSize*ySize;
  inDataArr = (double *) malloc(sizeof(double) * imsize);
  fits_read_pix(fptr, TDOUBLE, fpixel, imsize, NULL, inDataArr, NULL, &status);

  // Compute size of output image
  p1.x=0; p1.y=0;
  p2.x=xSize; p2.y=0;
  p3.x=xSize; p3.y=ySize;
  p4.x=0; p4.y=ySize;
  p1r = rotPoint(p1, orig, angle);
  p2r = rotPoint(p2, orig, angle);
  p3r = rotPoint(p3, orig, angle);
  p4r = rotPoint(p4, orig, angle);
  minX = (long)fmin(fmin(p1r.x, p2r.x), fmin(p3r.x, p4r.x));
  maxX = (long)fmax(fmax(p1r.x, p2r.x), fmax(p3r.x, p4r.x));
  minY = (long)fmin(fmin(p1r.y, p2r.y), fmin(p3r.y, p4r.y));
  maxY = (long)fmax(fmax(p1r.y, p2r.y), fmax(p3r.y, p4r.y));
  sizeXout = (maxX - minX);
  sizeYout = (maxY - minY);

  // Create out image
  // Check if file with given name is already exists (and delete it if so)
  if (access(outName, F_OK) != -1 ) {
    // file exists
    printf("Warning (crotima): fits file with name \'%s\' is already exists. Removing it...\n\n", outName);
    unlink(outName);
  }

  fits_create_file(&outima, outName, &status);
  fits_create_img(outima, bitpix, 2, naxesOut, &status);
  outDataArr = (double *) calloc(sizeof(double), (sizeXout*sizeYout));
  // Perform rotation
  costh = cos(-angle);
  sinth = sin(-angle);
  for (i=minX; i<maxX; i++) {
    for (j=minY; j<maxY; j++) {
      xPr = costh*(i-orig.x) - sinth*(j-orig.y) + orig.x;
      yPr = sinth*(i-orig.x) + costh*(j-orig.y) + orig.y;
      xInd = (int)(floor(xPr));
      yInd = (int)(floor(yPr));
      if ((xInd > 0) && (xInd < xSize-1) && (yInd > 0) && (yInd < ySize-1)){
	x = xPr - floor(xPr);
	y = yPr - floor(yPr);
	iDest = ((1-x)*(1-y) * inDataArr[xInd + yInd*xSize]
		 +x*(1-y) * inDataArr[xInd+1 + yInd*xSize]
		 +y*(1-x)*inDataArr[xInd + (yInd+1)*xSize]
		 +x*y*inDataArr[xInd+1 + (yInd+1)*xSize]);
	outDataArr[i-minX + (j-minY)*sizeXout] = iDest;
      }
    }
  }

  // Flush rotated data into file
  fits_write_pixnull(outima, TDOUBLE, fpixel, sizeXout*sizeYout, outDataArr, NULL, &status);

  // Crop image
  cenNew.x = orig.x - minX;
  cenNew.y = orig.y - minY;
  // Make cropped image
  sizeXoutCrop = cenNew.x * 2;
  sizeYoutCrop = cenNew.y * 2;
  fCropPixel[0] = 10;//orig.x-cenNew.x;
  fCropPixel[1] = 10;//orig.y-cenNew.y;
  lCropPixel[0] = 20;//orig.x+cenNew.x;
  lCropPixel[1] = 20;//orig.y+cenNew.y;
  double nullval = 0.0;
  outCropArray = (double *) calloc(sizeof(double), sizeXoutCrop*sizeYoutCrop);
  fits_read_subset(outima, TDOUBLE, fCropPixel, lCropPixel, inc,
		   &nullval, outCropArray, &anynul, &status);
  printf("%i\n", status);
  FILE * stream = fopen("error.txt", "w");
  fits_report_error(stream, status);

  // Create cropped out image
  // Check if file with given name is already exists (and delete it if so)
  char outCropName[9] = "crop.fits";

  if (access(outCropName, F_OK) != -1 ) {
    // file exists
    printf("Warning (crotima): fits file with name \'%s\' is already exists. Removing it...\n\n", outName);
    unlink(outCropName);
  }

  fits_create_file(&outCropIma, outCropName, &status);
  fits_create_img(outCropIma, bitpix, 2, naxesOut, &status);
  fits_write_pixnull(outCropIma, TDOUBLE, fpixel, sizeXoutCrop*sizeYoutCrop, outCropArray, NULL, &status);
  fits_close_file(outima, &status);
  fits_close_file(fptr, &status);
  fits_close_file(outCropIma, &status);
  free(inDataArr);
  free(outDataArr);
  free(outCropArray);
  return(0);
}


int stretchima(char fitsName[], char outName[], double scale, int workHDU)
{
  fitsfile *fptr, *outima;     
  long fpixel[2] = {1, 1};
  int i, j, xInd, yInd;
  int status = 0, numHDUs, hdutype=0, bitpix=0;
  long naxes[2] = {0, 0};
  long naxesOut[2] = {0, 0};
  long xSize, ySize, outDimX, outDimY, imsize, sizeXout, sizeYout;
  long minX, maxX, minY, maxY;
  double *inDataArr;
  double *outDataArr;
  double xPr, yPr, x1, x2, iDest;

  // Open file
  fits_open_file(&fptr, fitsName, READONLY, &status);

  // Let's check if work hdu is exists
  fits_get_num_hdus(fptr, &numHDUs, &status);
  if (workHDU > numHDUs) return 1;

  // Set work HDU
  fits_movabs_hdu(fptr, workHDU, &hdutype, &status);

  // Get information about size of image and type of data
  fits_get_img_size(fptr, 2,  naxes, &status);
  fits_get_img_type(fptr, &bitpix, &status);
  xSize = naxes[0];
  ySize = naxes[1];

  // Read data from input file
  imsize = xSize*ySize;
  inDataArr = (double *) malloc(sizeof(double) * imsize);
  fits_read_pix(fptr, TDOUBLE, fpixel, imsize, NULL, inDataArr, NULL, &status);

  // Compute parameters of the output image
  sizeXout = (long) ceil(xSize * scale);
  sizeYout = ySize; // y-axis does not change
  naxesOut[0] = sizeXout;  
  naxesOut[1] = sizeYout;

  // Create out image
  // Check if file with given name is already exists (and delete it if so)
  if (access(outName, F_OK) != -1 ) {
    // file exists
    printf("Warning (crotima): fits file with name \'%s\' is already exists. Removing it...\n\n", outName);
    unlink(outName);
  }
  fits_create_file(&outima, outName, &status);
  fits_create_img(outima, bitpix, 2, naxesOut, &status);
  outDataArr = (double *) malloc(sizeof(double) * (sizeXout*sizeYout));
  // Fill array with out data by zeros
  for (i=0; i <sizeXout*sizeYout; i++)
    outDataArr[i] = 0.0;
  // Perform stratching
  for (i=0; i<sizeXout; i++) {
    for (j=0; j<sizeYout; j++) {
      xPr = ((double) i)/scale;
      yPr = j;
      xInd = (int)(floor(xPr));
      yInd = (int)(floor(yPr));
      if ((xInd > 0) && (xInd < xSize-1) && (yInd > 0) && (yInd < ySize-1)){
	x1 = xInd+1-xPr;
	x2 = 1-x1;
	iDest = x1 * inDataArr[xInd+yInd*xSize] + x2 * inDataArr[xInd+1+yInd*xSize];
	outDataArr[i + j*sizeXout] = iDest;
      }
    }
  }

  // Flush stratched data into file
  fits_write_pixnull(outima, TDOUBLE, fpixel, sizeXout*sizeYout, outDataArr, NULL, &status);
  fits_close_file(outima, &status);
  fits_close_file(fptr, &status);
  free(inDataArr);
  free(outDataArr);
  return(0);
}


int deproject(char infile[], char outfile[], double xCen, double yCen, double incl, double posang, int workHDU)
{
  /* If 0<incl<1 than it is treated as an axis ratio, if incl>1
     tnan it is regarded as an angle in degrees*/

  char tmpFile[25]={"temp_fits_file_sdfvs.fits"}; 
  double stretchParam, inclRadians;

  // Check if temporary file exists
  if (access(tmpFile, F_OK) != -1 )
    unlink(tmpFile);

  // Rotate file first
  double * newCoordsOfOrigin;
  newCoordsOfOrigin = (double *) malloc(sizeof(double)*2);
  crotima(infile, tmpFile, xCen, yCen, posang, workHDU, newCoordsOfOrigin);
  // Stretch now
  if (incl < 1.0)
    stretchParam = 1.0/incl;
  else
    {
      inclRadians = incl*PI/180.0;
      stretchParam = 1/cos(inclRadians);
    }
  stretchima(tmpFile, outfile, stretchParam, workHDU);
  unlink(tmpFile);
  return(0);
}


int main(int argc, char *argv[])
{
  int returnCode=0;
  deproject(argv[1], "depr.fits", 323.0, 323.0, 0.7, 110, 1);
  return(returnCode);
}
