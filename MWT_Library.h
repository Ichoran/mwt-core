/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_LIBRARY
#define MWT_LIBRARY


#include <time.h>
#include <stdio.h>


#include "MWT_Geometry.h"
#include "MWT_Lists.h"
#include "MWT_Storage.h"
#include "MWT_Image.h"
#include "MWT_Align.h"
#include "MWT_Blob.h"
#include "MWT_Model.h"



/****************************************************************
             Helper class to store summary data
****************************************************************/

// Stores data that isn't maintained by Performance
class SummaryData
{
public:
  int frame_number;
  float frame_time;
  bool initialized;
  
  int n_dancers_tracked;
  int n_dancers_good;
  float dancer_persistence;
  float dancer_speed;
  float dancer_angularspeed;
  float dancer_length;
  float dancer_relativelength;
  float dancer_width;
  float dancer_relativewidth;
  float dancer_aspect;
  float dancer_relativeaspect;
  float dancer_endwiggle;
  float dancer_pixelcount;
  float dancer_jitter_x;
  float dancer_jitter_y;
  int dancer_shift_x;
  int dancer_shift_y;
	
  int n_speed_updates;
  float update_frequency;
  
  ManagedList<int> event_list;
  
  SummaryData() : frame_number(-1),event_list(0) { } // Don't call this, just keeps compiler happy
  SummaryData(int frame_n,float frame_t,Storage< Listable<int> >* listore)
    : frame_number(frame_n),frame_time(frame_t),initialized(false),
      n_dancers_tracked(0),n_dancers_good(0),dancer_persistence(0.0),
      dancer_speed(0.0),dancer_angularspeed(0.0),
      dancer_length(0.0),dancer_relativelength(0.0),
      dancer_width(0.0),dancer_relativewidth(0.0),
      dancer_aspect(0.0),dancer_relativeaspect(0.0),
      dancer_endwiggle(0.0),dancer_pixelcount(0.0),
      dancer_jitter_x(0.0),dancer_jitter_y(0.0),dancer_shift_x(0),dancer_shift_y(0),
      n_speed_updates(0),update_frequency(1.0),
      event_list(listore)
  { }
  ~SummaryData() { }
  
  // Output to file
  void fprint(FILE* f , ManagedList<BlobOriginFate>& mlbof, ManagedList<BlobOutputFate>& mlbf);
};



/****************************************************************
     Classes to wrap one (or more) Performances in a Library
****************************************************************/

// Stores data about one performance
class TrackerEntry
{
public:
  Performance performance;
  Listable<int> *handle;
  ManagedList<Profile> reference_edges;
  ManagedList<Point> reference_objects;
  ManagedList<Point> working;
  
  Profile *cached_profile;
  int cached_profile_index;
  char* tracker_string;
  char* path_string;
  char* prefix_string;
  bool save_objects;
  bool save_refs;
  bool save_images;
  struct tm *output_date;
  int image_bits;  // Will need to use this when wrapping images from LabView--they don't say how many bits the camera has!
  int image_width;
  int image_height;
  
  // Error checking--make sure the appropriate things get initialized
  bool image_info_known;
  bool output_info_known;
  bool border_size_known;
  bool reference_edges_final;
  bool reference_intensities_known;
  bool object_intensities_known;
  bool object_sizes_known;
  bool object_persistence_known;
  bool adaptation_rate_known;
  bool references_found;
  bool objects_found;
  bool output_started;
  bool image_loaded;
  bool statistics_ready;
  
  // Statistics
  float update_frequency;
  Storage< Listable<int> > eventstore;
  ManagedList<SummaryData> summary;
  
  TrackerEntry()  // Don't call this one--just here to keep compiler happy
    : handle(NULL),tracker_string(NULL),path_string(NULL),prefix_string(NULL),eventstore(0),summary(0)
  { }
  TrackerEntry(Listable<int> *new_handle,int buffer_size)  // Call this one using placement new 
    : performance(new_handle->data,buffer_size),handle(new_handle),reference_edges(16,true),reference_objects(16,false),working(16,false),
      cached_profile(NULL),cached_profile_index(-1),
      tracker_string(NULL),path_string(NULL),prefix_string(NULL),output_date(NULL),
      update_frequency(1.0),eventstore(buffer_size,false) , summary(buffer_size,true)
  {
    image_info_known = false;
    output_info_known = false;
    border_size_known = false;
    reference_edges_final = false;
    reference_intensities_known = false;
    object_intensities_known = false;
    object_sizes_known = false;
    object_persistence_known = false;
    adaptation_rate_known = false;
    references_found = false;
    objects_found = false;
    output_started = false;
    image_loaded = false;
    statistics_ready = false;
  }
  ~TrackerEntry()
  {
    if (tracker_string!=NULL) { delete[] tracker_string; tracker_string=NULL; }
    if (path_string!=NULL) { delete[] path_string; path_string=NULL; }
    if (prefix_string!=NULL) { delete[] prefix_string; prefix_string=NULL; }
    if (output_date!=NULL) { delete output_date; output_date=NULL; }
  }
  void setTrackerName(const char *s) { if (tracker_string!=NULL) { delete[] tracker_string; } tracker_string = new char[strlen(s)+1]; strcpy(tracker_string, s); }
  void setPath(const char *s) { if (path_string!=NULL) { delete[] path_string; } path_string = new char[strlen(s)+1]; strcpy(path_string,s); }
  void setPrefix(const char *s) { if (prefix_string!=NULL) { delete[] prefix_string; } prefix_string = new char[strlen(s)+1]; strcpy(prefix_string,s); }
  void setDate(int year,int month,int day,int hour,int minute,int second)
  {
    if (output_date!=NULL) delete output_date;
    output_date = new struct tm;
    output_date->tm_year = year;
    output_date->tm_mon = month-1;
    output_date->tm_mday = day;
    output_date->tm_hour = hour;
    output_date->tm_min = minute;
    output_date->tm_sec = second;
  }
};


// Stores data about all performances with a C-ish handle-based interface
class TrackerLibrary
{
private:
  void swap(int& a,int& b) { int i=a; a=b; b=i; }  // Used to fix inputs that are entered backwards
public:
  static const int MAX_TRACKER_HANDLES = 999;
  static const int DEFAULT_DEFAULT_SIZE = 256; 
  
  TrackerEntry** all_trackers;
  int default_buffer_size;
  ManagedList<int> active_handles;
  Storage< Listable<Strip> > lsstore;
  
  TrackerLibrary();
  ~TrackerLibrary();
  
  // Creating new performances
  int getNewHandle();
  int killHandle(int victim);
  
  // Preparing output
	int setCombineBlobs( int handle, bool type );
  int setDate(int handle,int year,int month,int day,int hour,int minute,int second);
  int borrowDate(int handle,int donor_handle);
  int setAllDatesToMine(int handle);
  int setOutput(int handle,const char *path,const char *prefix,bool save_obj,bool save_ref,bool save_im);
  int beginOutput(int handle);
  int setAndBeginOutput(int handle,const char *path,const char *prefix,bool save_obj,bool save_ref,bool save_im)
  {
    int i = setOutput(handle,path,prefix,save_obj,save_ref,save_im);
    if (i==0) return beginOutput(handle);
    else return i;
  }
  
  // General image setup
  int setImageInfo(int handle,int bits,int width,int height);
  int getBitDepth(int handle);
  int setDancerBorderSize(int handle,int border);
  int setUpdateBandNumber(int handle,int n_bands);
  
	
  // Preparing and getting feedback on the region of interest
  int setRectangle(int handle,int left,int right,int top,int bottom);
  int setEllipse(int handle,int centerx,int centery,int radiusx,int radiusy);
  int addRectangle(int handle,int left,int right,int top,int bottom);
  int addEllipse(int handle,int centerx,int centery,int radiusx,int radiusy);
  int cutRectangle(int handle,int left,int right,int top,int bottom);
  int cutEllipse(int handle,int centerx,int centery,int radiusx,int radiusy);
  int showROI(int handle,Image& im);
  // 8-bit image version
  int showROI8(int handle,Image8& im);
  
  // Utility image processing functions
  int resizeRescale(int handle,Image& source,Image& dest,Rectangle dest_selection);
  int resizeRescaleMemory(int handle,Image& dest,Rectangle dest_selection);
  int resizeRescaleFixed(int handle,Image& dest,Rectangle dest_selection);
  // 8-bit image versions
  int resizeRescale8(int handle,Image8& source,Image8& dest,Rectangle dest_selection);
  int resizeRescaleMemory8(int handle,Image8& dest,Rectangle dest_selection);
  int resizeRescaleFixed8(int handle,Image8& dest,Rectangle dest_selection);

  // Finding edges for jitter alignment
private:
  int checkInImageBounds(int handle, const Rectangle &fromImage, const Rectangle &search, const Point &size);
public:
  float addBestXEdgesInImage(int handle, Image& im, const Rectangle &search, const Point &size);
  float addBestYEdgesInImage(int handle, Image& im, const Rectangle &search, const Point &size);
  float addBestXEdgesInImage8(int handle, Image8& im, const Rectangle &search, const Point &size);
  float addBestYEdgesInImage8(int handle, Image8& im, const Rectangle &search, const Point &size);
  int loadedProfiles(int handle);
private:
  Profile* getProfileIfValid(int handle, int profile);
public:
  int profileDirection(int handle, int profile);
  int profileCorner_x0(int handle, int profile);
  int profileCorner_x1(int handle, int profile);
  int profileCorner_y0(int handle, int profile);
  int profileCorner_y1(int handle, int profile);
  int profileEdgeCount(int handle, int profile);
  float profileEdgeLocation(int handle, int profile, int edge);
  float globalProfileEdgeLocation(int handle, int profile, int edge);
  float profileEdgeStrength(int handle, int profile, int edge);
  int profileEdgePolarity(int handle, int profile, int edge);
  int removeOneProfile(int handle, int profile);
  int removeAllProfiles(int handle);
  
  // Preparing and getting feedback on reference objects
  int setRefIntensityThreshold(int handle,int intensity_low,int intensity_high);
  int addReferenceObjectLocation(int handle,int x,int y);
  int removeLastReferenceObject(int handle);
  int scanRefs(int handle,Image& im);
  int showRefs(int handle,Image& im,bool show_as_white=true);
  // 8-bit image versions
  int scanRefs8(int handle,Image8& im);
  int showRefs8(int handle,Image8& im,bool show_as_white=true);
  
  // Preparing and getting feedback on moving objects
  int setObjectIntensityThresholds(int handle,int intensity_to_fill,int intensity_of_new);
  int setObjectSizeThresholds(int handle,int small_size,int small_good_size,int large_good_size,int large_size);
  int setObjectPersistenceThreshold(int handle,int frames);
  int setAdaptationRate(int handle,int alpha);
  int scanObjects(int handle,Image& im);
  int showObjects(int handle,Image& im);
  // 8-bit image versions
  int scanObjects8(int handle,Image8& im);
  int showObjects8(int handle,Image8& im);
  
  // Controls for what statistical data to collect
  int setVelocityIntegrationTime(int handle,float interval);
  int enableSkeletonization(int handle,bool enable=true);
  int enableOutlining(int handle,bool enable=true);
  
  // Running the main loop
  int prepareImagePieces(int handle,float time);
  int getNextPieceCoords(int handle,Rectangle& coords);
  int loadThisImagePiece(int handle,Image& im); 
  int loadImage(int handle,Image& im,float time);  // Prepares, gets, and loads all pieces at once
  int loadImageAsPieces(int handle,Image& im,float time);  // Inefficient version of loadImage with an extra copy of image pieces (useful for testing)
  int showLoaded(int handle,Image& im);
  int markEvent(int handle,int event_number);  // Refers to currently loaded image
  int processImage(int handle);
  int showResults(int handle,Image& im);
  int checkErrors(int handle);
  int complete(int handle);  // Must call this to save summary information!
  // 8-bit image versions
  int loadThisImagePiece8(int handle,Image8& im); 
  int loadImage8(int handle,Image8& im,float time);  // Prepares, gets, and loads all pieces at once
  int loadImageAsPieces8(int handle,Image8& im,float time);  // Inefficient version of loadImage with an extra copy of image pieces (useful for testing)
  int showLoaded8(int handle,Image8& im);
  int showResults8(int handle,Image8& im);
  
  // Image correction algorithm
  int setDivisionImageCorrectionAlgorithm( int handle );
  int setSubtractionImageCorrectionAlgorithm( int handle );
  int reportImageCorrectionAlgorithm( int handle );
  double generateBackgroundAverageIntensity( int handle, int *array );
	
  // Summary statistics
  int reportNumber(int handle);  
  int reportNumberPersisting(int handle);  
  float reportPersistence(int handle);  
  float reportSpeed(int handle);
  float reportAngularSpeed(int handle);
  float reportLength(int handle);  
  float reportRelativeLength(int handle);
  float reportWidth(int handle);  // Reports the mean width in pixels
  float reportRelativeWidth(int handle);  // Reports the current width / mean width
  float reportAspect(int handle);  // Reports width/length
  float reportRelativeAspect(int handle);  // Reports current (width/length) / mean (width/length)
  float reportEndWiggle(int handle);  // Angle of turned head or tail
  float reportObjectPixelCount(int handle); // Reports mean pixel count of found objects.
};



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int mwt_test_library(int bin);

#endif

