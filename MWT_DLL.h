/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_DLL
#define MWT_DLL

// NOTE: Always update me!!!
#define REVISION 975


#include "MWT_Library.h"

// The 1d array handle structure for LabView
typedef struct {
  int dimSizes[1];
  int elt[1];
} TD1;

// The 2d array handle structure for LabView
typedef struct {
	// The dimensions of the array
	// Index 0 == x coordinate (coarse, x values advance lines)
	// Index 1 == y coordinate (fine, y values are adjacent)
	int dimSizes[2];
	// The actual array, this is a 2d array enfolded in to a 1d array.
  short elt[1];
} TD2;

typedef TD2 **TD2Hdl;
typedef TD1 **TD1Hdl;

extern "C" {
int MWT_reportRevisionNumber() { return REVISION; }

int MWT_getNewHandle();
int MWT_killHandle(int handle);

int MWT_setDate(int handle,int year,int month,int day,int hour,int minute,int second);
int MWT_borrowDate(int handle,int donor_handle);

int MWT_setOutput(int handle,char *path,char *prefix,bool save_obj,bool save_ref,bool save_im);
int MWT_beginOutput(int handle);
int MWT_setAndBeginOutput(int handle,char *path,char *prefix,bool save_obj,bool save_ref,bool save_im);

int MWT_setImageInfo(int handle,int bits,int width,int height);

int MWT_setRectangle(int handle,int x1,int y1,int x2,int y2);
int MWT_setEllipse(int handle,int centerx,int centery,int radiusx,int radiusy);
int MWT_addRectangle(int handle,int left,int right,int top,int bottom);
int MWT_addEllipse(int handle,int centerx,int centery,int radiusx,int radiusy);
int MWT_cutRectangle(int handle,int left,int right,int top,int bottom);
int MWT_cutEllipse(int handle,int centerx,int centery,int radiusx,int radiusy);

int MWT_showROI(int handle,TD2Hdl im, int bin=1);  

int MWT_setDancerBorderSize(int handle,int border);
int MWT_setRefIntensityThreshold(int handle,int intensity_low,int intensity_high);
int MWT_addReferenceObjectLocation(int handle,int x,int y);
int MWT_removeLastReferenceObject(int handle);

int MWT_setOutputType( int handle, int type);

int MWT_setUpdateBandNumber(int handle,int num_bands);

int MWT_scanRefs(int handle, TD2Hdl im, int bin=1);
int MWT_showRefs(int handle, TD2Hdl im, int bin=1);

int MWT_enableSkeletonization(int handle, bool enable); 
int MWT_enableOutlining(int handle, bool enable);


int MWT_setObjectIntensityThresholds(int handle,int intensity_to_fill, int intensity_of_new);
int MWT_setObjectSizeThresholds(int handle, int small_size,int small_good_size, int large_good_size, int large_size);
int MWT_setObjectPersistenceThreshold(int handle,int frames);
int MWT_setAdaptationRate(int handle,int alpha);

int MWT_scan(int handle,TD2Hdl im, int bin=1);  
int MWT_showScan(int handle,TD2Hdl im, int bin=1); 

int MWT_setDivisionImageCorrectionAlgorithm( int handle );
int MWT_setSubtractionImageCorrectionAlgorithm( int handle );
int MWT_reportImageCorrectionAlgorithm(int handle);
float MWT_generateDivisionCorrection(int handle, TD1Hdl coords);
int MWT_loadImage(int handle,TD2Hdl im,float time, int bin=1); 
int MWT_prepareImagePieces( int handle, float time);
int MWT_getNextPieceCoords( int handle, TD1Hdl coords); 
int MWT_loadThisImagePiece( int handle, TD2Hdl im, int x, int y, int bin=1);
int MWT_showLoaded(int handle, TD2Hdl im, int bin=1);
int MWT_markEvent(int handle, int event_number);
int MWT_processImage(int handle);
int MWT_checkErrors(int handle);
int MWT_complete(int handle);
  
int MWT_resizeRescale( int handle, TD2Hdl src, TD2Hdl dst, int dest_width, int dest_height);
int MWT_resizeRescaleMemory( int handle, TD2Hdl dest, int dest_width, int dest_height);
int MWT_resizeRescaleFixed( int handle, TD2Hdl dst, int dest_width, int dest_height);
int MWT_showResults(int handle, TD2Hdl im, int bin=1);

int MWT_reportNumber(int handle);
int MWT_reportNumberPersisting(int handle);
float MWT_reportPersistence(int handle);
float MWT_reportSpeed(int handle);
float MWT_reportAngularSpeed(int handle);
float MWT_reportLength(int handle);
float MWT_reportRelativeLength(int handle);
float MWT_reportWidth(int handle);
float MWT_reportRelativeWidth(int handle);
float MWT_reportAspect(int handle);
float MWT_reportRelativeAspect(int handle);
float MWT_reportEndWiggle(int handle);
float MWT_reportObjectPixelCount(int handle);
}
#endif
