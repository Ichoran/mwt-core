/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include "MWT_DLL.h"



/****************************************************************
              The One and Only Global Variable
****************************************************************/

TrackerLibrary the_library;


/****************************************************************
             Utility methods Not Attached to a Handle
****************************************************************/

void MWT_widen_image_8_to_16(char* in, int inStride, short* out, int outStride, int nx, int ny) {
  Image::copy8to16(in, inStride, out, outStride, nx, ny);
}


/****************************************************************
      Procedural Wrappers to Global TrackerLibrary Variable
****************************************************************/

int MWT_getNewHandle()
{
  return the_library.getNewHandle();
}

int MWT_killHandle(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.killHandle(handle);
}

int MWT_setDate(int handle,int year,int month,int day,int hour,int minute,int second)
{
  if( handle < 1 ) return -3;

  return the_library.setDate(handle,year,month,day,hour,minute,second);
}

int MWT_borrowDate(int handle,int donor_handle)
{
  if( handle < 1 || donor_handle < 1 ) return -3;

  return the_library.borrowDate(handle,donor_handle);
}

int MWT_setOutput(int handle,char *path,char *prefix,bool save_obj,bool save_ref,bool save_im)
{
  if( handle < 1 || path == NULL || prefix == NULL ) return -3;
  
  return the_library.setOutput(handle,path,prefix,save_obj,save_ref,save_im);
}

int MWT_beginOutput(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.beginOutput(handle);
}

int MWT_setAndBeginOutput(int handle,char *path,char *prefix,bool save_obj,bool save_ref,bool save_im)
{
  if( handle < 1 || path == NULL || prefix == NULL) return -3;
  
  return the_library.setAndBeginOutput(handle,path,prefix,save_obj,save_ref,save_im);
}

int MWT_setRectangle(int handle,int x1,int y1,int x2,int y2)
{
  if( handle < 1 ) return -3;

  return the_library.setRectangle(handle,x1,x2,y1,y2);
}

int MWT_setEllipse(int handle,int centerx,int centery,int radiusx,int radiusy)
{
  if( handle < 1 ) return -3;
  
  return the_library.setEllipse(handle,centerx,centery,radiusx,radiusy);
}

int MWT_addRectangle(int handle,int x1,int y1,int x2,int y2)
{
  if( handle < 1 ) return -3;
  
  return the_library.addRectangle(handle,x1,x2,y1,y2);
}

int MWT_addEllipse(int handle,int centerx,int centery,int radiusx,int radiusy)
{
  if( handle < 1 ) return -3;
  
  return the_library.addEllipse(handle,centerx,centery,radiusx,radiusy);
}

int MWT_cutRectangle(int handle,int left,int right,int top,int bottom)
{
  if( handle < 1 ) return -3;
  
  return the_library.cutRectangle(handle,left,right,top,bottom);
}

int MWT_cutEllipse(int handle,int centerx,int centery,int radiusx,int radiusy)
{
  if( handle < 1 ) return -3;
  
  return the_library.cutEllipse(handle,centerx,centery,radiusx,radiusy);
}

int MWT_showROI(int handle, TD2Hdl im, int bin)
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;
  
  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.showROI(handle,img);
}

int MWT_setDancerBorderSize(int handle,int border)
{
  if( handle < 1 || border < 1 ) return -3;
  
  return the_library.setDancerBorderSize(handle,border);
}

int MWT_setImageInfo(int handle,int bits,int width,int height)
{
  if( handle < 1 || bits < 1 || width < 1 || height < 1 ) return -3;
  
  return the_library.setImageInfo(handle,bits,width,height);
}

int MWT_setRefIntensityThreshold(int handle,int intensity_low,int intensity_high)
{
  if( handle < 1 ) return -3;

  return the_library.setRefIntensityThreshold(handle,intensity_low,intensity_high);
}

int MWT_addReferenceObjectLocation(int handle,int x,int y)
{
  if( handle < 1 ) return -3;
  
  return the_library.addReferenceObjectLocation(handle,x,y);
}

int MWT_removeLastReferenceObject(int handle)
{
  if( handle < 1 ) return -3;

  return the_library.removeLastReferenceObject(handle);
}

int MWT_setOutputType(int handle, int type) 
{
	if( handle < 1 ) return -3;
	
	if( type ) 
	{
		return the_library.setCombineBlobs( handle, true );
	}
	else 
	{
		return the_library.setCombineBlobs( handle, false );
	}
}

int MWT_setUpdateBandNumber(int handle, int num_bands) {
  if( handle < 1 ) return -3;
  
  return the_library.setUpdateBandNumber(handle,num_bands);
}


int MWT_scanRefs(int handle, TD2Hdl im, int bin) 
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;

  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.scanRefs(handle,img);
}

int MWT_showRefs(int handle, TD2Hdl im, int bin)
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;
  
  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.showRefs(handle,img,true);
}

int MWT_enableSkeletonization( int handle, bool enable ) 
{
  if( handle < 1 ) return -3;
  
  return the_library.enableSkeletonization(handle,enable); 
}

int MWT_enableOutlining(int handle, bool enable)
{
  if( handle < 1 ) return -3;
  
  return the_library.enableOutlining(handle, enable);
}

int MWT_setObjectIntensityThresholds(int handle,int intensity_to_fill,int intensity_of_new)
{  
  if( handle < 1 ) return -3;
  
  /* NOTE 
              2* is to scale the incoming intensity in to the intensity scores used in the library
   */
  return the_library.setObjectIntensityThresholds(handle,2*intensity_to_fill,2*intensity_of_new);
}

int MWT_setObjectSizeThresholds(int handle, int small_size,int small_good_size, int large_good_size, int large_size)
{
  if( handle < 1 || small_size < 1 || small_good_size < 1 || large_good_size < 1 || large_size < 1 ) return -3;
  
  return the_library.setObjectSizeThresholds(handle,small_size,small_good_size,large_good_size,large_size);
}

int MWT_setObjectPersistenceThreshold(int handle,int frames)
{
  if( handle < 1 || frames < 1 ) return -3;
  
  return the_library.setObjectPersistenceThreshold(handle,frames);
}

int MWT_setAdaptationRate(int handle,int alpha)
{
  if( handle < 1 || alpha < 1 ) return -3;
  
  return the_library.setAdaptationRate(handle,alpha);
}

int MWT_scan(int handle, TD2Hdl im, int bin)
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;

	Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.scanObjects(handle,img);
}

int MWT_showScan(int handle, TD2Hdl im, int bin)
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;
  
  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.showObjects(handle,img);
}

float MWT_generateDivisionCorrection( int handle, TD1Hdl coords ) 
{
	if( handle < 1 || coords == NULL ) return -3;
	else if( (*coords)->dimSizes[0] < NUM_PTS_TO_SAMPLE*2 ) return -4;
	
	// return immediately if we're not doing division correction.
	if( the_library.reportImageCorrectionAlgorithm(handle) == 0 ) return -5;
		
	int *array = (int*)malloc(2*NUM_PTS_TO_SAMPLE*sizeof(int));
	if( array == NULL ) return -3;
	array[0] = -1; // set dummy value so later method knows to fill array.
	
	// If coord array has been filled with data, copy it. 
	if( (*coords)->elt[0] >= 0 ) {
		for( int i = 0; i < NUM_PTS_TO_SAMPLE * 2; i++ ) array[i] = (*coords)->elt[i];
	}
	
	float ret = the_library.generateBackgroundAverageIntensity( handle, array );
	
	if( (*coords)->elt[0] < 0 ) {
		for( int i = 0; i < 2*NUM_PTS_TO_SAMPLE; i++ ) 
		{
			(*coords)->elt[i] = array[i];
		}	
	}
	free(array);
	return ret;
}

int MWT_setDivisionImageCorrectionAlgorithm( int handle )
{
	if( handle < 1 ) return -3;
	
	return the_library.setDivisionImageCorrectionAlgorithm( handle );
}

int MWT_setSubtractionImageCorrectionAlgorithm( int handle )
{
	if( handle < 1 ) return -3;
	
	return the_library.setSubtractionImageCorrectionAlgorithm( handle );
}

int MWT_reportImageCorrectionAlgorithm(int handle) 
{
	if( handle < 1 ) return -3;
	
	return the_library.reportImageCorrectionAlgorithm(handle);
}

int MWT_loadImage(int handle,TD2Hdl im,float time, int bin)
{
  if( handle < 1 || im == NULL || time < 0) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;
  
  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
	return the_library.loadImageAsPieces(handle,img,time);
}

int MWT_prepareImagePieces( int handle, float time ) 
{
  if( handle < 1 || time < 0 ) return -3;  

  return the_library.prepareImagePieces(handle,time);
}

int MWT_getNextPieceCoords( int handle, TD1Hdl coords ) 
{
  if( handle < 1 || coords == NULL ) return -3;
  else if( (*coords)->dimSizes[0] < 4 ) return -4;
  
  int ret;
  Rectangle rect;
	ret = the_library.getNextPieceCoords(handle,rect);
	// NOTE: +1 to the far point corrects for LabViews desire to return data from near to far-1.
  (*coords)->elt[0] = rect.near.x;
	(*coords)->elt[1] = rect.near.y;
  (*coords)->elt[2] = rect.far.x +1;
  (*coords)->elt[3] = rect.far.y +1;
  return ret;
}

int MWT_loadThisImagePiece( int handle, TD2Hdl im, int x, int y, int bin ) {
	if( handle < 1 || im == NULL )	return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 )	return -4;
	Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
	img.bounds += Point(x,y);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
	
  return the_library.loadThisImagePiece( handle, img );
}

int MWT_markEvent(int handle, int event_number) 
{
  if( handle < 1  ) return -3;

  return the_library.markEvent(handle,event_number);
}

int MWT_processImage(int handle)
{
  if( handle < 1 ) return -3;
  
	return the_library.processImage(handle);
}

int MWT_showLoaded(int handle, TD2Hdl im, int bin) 
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;
  
  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.showLoaded(handle,img);
}

// Should be called when display is enabled and it is set to "Fixed."
int MWT_showResults(int handle,TD2Hdl im, int bin)
{
  if( handle < 1 || im == NULL ) return -3;
  else if( (*im)->dimSizes[0] < 1 || (*im)->dimSizes[1] < 1 ) return -4;

  Image img((*im)->elt,Point((*im)->dimSizes[0],(*im)->dimSizes[1]),false);
  img.bin = bin;
  img.depth = the_library.getBitDepth(handle);
  return the_library.showResults(handle,img);
}

int MWT_resizeRescale( int handle, TD2Hdl s, TD2Hdl dst, int dest_width, int dest_height) 
{
  if( handle < 1 || s == NULL || dst == NULL ) return -3;
  else if( (*s)->dimSizes[0] < 1 || (*s)->dimSizes[1] < 1 ) return -4;
  else if( (*dst)->dimSizes[0] < 1 || (*dst)->dimSizes[1] < 1 ) return -4;
  else if( dest_width < 1 || dest_height < 1 ) return -5;
  
  Image src((*s)->elt,Point((*s)->dimSizes[0],(*s)->dimSizes[1]),false);
  Image dest((*dst)->elt,Point((*dst)->dimSizes[0],(*dst)->dimSizes[1]),false);
  src.depth = dest.depth = the_library.getBitDepth(handle);
  Rectangle rect( Point(0,0) , Point(dest_width-1,dest_height-1));
  return  the_library.resizeRescale(handle,src,dest,rect);
}

int MWT_resizeRescaleMemory( int handle, TD2Hdl dst, int dest_width, int dest_height) 
{
  if( handle < 1 || dst == NULL ) return -3;
  else if( (*dst)->dimSizes[0] < 1 || (*dst)->dimSizes[1] < 1 ) return -4;
  else if( dest_width < 1 || dest_height < 1 ) return -5;

  Image dest((*dst)->elt,Point((*dst)->dimSizes[0],(*dst)->dimSizes[1]),false);
  dest.depth = the_library.getBitDepth(handle);
  Rectangle rect( Point(0,0) , Point(dest_width-1,dest_height-1));
  return  the_library.resizeRescaleMemory(handle,dest,rect);
}

int MWT_resizeRescaleFixed( int handle, TD2Hdl dst, int dest_width, int dest_height) 
{
  if( handle < 1 || dst == NULL ) return -3;
  else if( (*dst)->dimSizes[0] < 1 || (*dst)->dimSizes[1] < 1 ) return -4;
  else if( dest_width < 1 || dest_height < 1 ) return -5;
  
  Image dest((*dst)->elt,Point((*dst)->dimSizes[0],(*dst)->dimSizes[1]),false);  
  dest.depth = the_library.getBitDepth(handle);
  Rectangle rect( Point(0,0) , Point(dest_width-1,dest_height-1));
  return  the_library.resizeRescaleFixed(handle,dest,rect);
}


int MWT_checkErrors(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.checkErrors(handle);
}

int MWT_complete(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.complete(handle);
}

int MWT_reportNumber(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportNumber(handle);
}

int MWT_reportNumberPersisting(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportNumberPersisting(handle);
}

float MWT_reportPersistence(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportPersistence(handle);
}

float MWT_reportSpeed(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportSpeed(handle);
}

float MWT_reportAngularSpeed(int handle)
{
  if( handle < 1 ) return -3;
  
  float return_value = the_library.reportAngularSpeed(handle) * (180/3.14159265); /* convert radians to degrees */
  return return_value;
}

float MWT_reportLength(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportLength(handle);
}

float MWT_reportRelativeLength(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportRelativeLength(handle);
}

float MWT_reportWidth(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportWidth(handle);
}

float MWT_reportRelativeWidth(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportRelativeWidth(handle);
}

float MWT_reportAspect(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportAspect(handle);
}

float MWT_reportRelativeAspect(int handle)
{
  if( handle < 1 ) return -3;
  
  return the_library.reportRelativeAspect(handle);
}

float MWT_reportEndWiggle(int handle) 
{
  if( handle < 1 ) return -3;
  
  float return_value = the_library.reportEndWiggle(handle) * (180/3.14159265); /* convert radians to degrees */
  return return_value;
}

float MWT_reportObjectPixelCount(int handle)
{
	if( handle < 1 ) return -3;	
	return the_library.reportObjectPixelCount(handle);
}
