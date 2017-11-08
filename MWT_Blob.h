/* Copyright 2007-2010 Rex A. Kerr, Nicholas A. Swierczek, and HHMI Janelia
 * Copyright 2016 Rex A. Kerr and Calico Life Sciences
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_BLOB
#define MWT_BLOB

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#include "MWT_Geometry.h"
#include "MWT_Lists.h"
#include "MWT_Storage.h"
#include "MWT_Image.h"
#include "MWT_Align.h"

#define NUM_PTS_TO_SAMPLE 100

// (We need this later.)
class Dancer;
class Performance;



/****************************************************************
                 Utility Classes for Blobs
****************************************************************/
// A little class that tells us where blobs were output
class BlobOutputFate
{
public:
  int frame;
  int blob_ID;
  int file_ID;
  long byte_offset;
  BlobOutputFate() { byte_offset = -1; }
  BlobOutputFate( int f, int idb, int idf, unsigned long bo) : frame(f), blob_ID(idb), file_ID(idf), byte_offset(bo) { }
};

// A little class that tells us what happened to blobs
class BlobOriginFate
{
public:
  int frame;
  int ID_origin;
  int ID_fate;
  BlobOriginFate() { }
  BlobOriginFate(int f,int ido,int idf) : frame(f),ID_origin(ido),ID_fate(idf) { }
};


// A little error handling class
class IOErrorHandler
{
public:
  int ID;
  char *what;
  
  IOErrorHandler() : ID(-1),what(NULL) { }
  IOErrorHandler(int i,const char *msg) : ID(i) { int L=strlen(msg); what = new char[L+1]; strcpy(what,msg); }
  IOErrorHandler(int i,char *msg,bool dupe) : ID(i) { if (dupe) { int L=strlen(msg); what = new char[L+1]; strcpy(what,msg); } else what=msg; }
  ~IOErrorHandler() { if (what!=NULL) { delete[] what; what=NULL; } }
  IOErrorHandler& operator=(IOErrorHandler& ioeh)
  {
    ID = ioeh.ID; 
    if (what) delete[] what; 
    what = ioeh.what; 
    ioeh.what = NULL; 
    return *this;
  }
};


// A little formatted filename class
class FilenameComponent
{
public:
  static const int MAX_WIDTH = 16; // Must be at least the number of digits used by a signed integer (i.e. 11, for negative two billion)
  enum Identity { empty , frame_n , dancer_n , perf_n , text };
  Identity identity;
  int value;
  int width;
  char* format;
  FilenameComponent() : identity(empty),value(0),width(0),format(NULL) { }
  FilenameComponent(const FilenameComponent& fc);
  ~FilenameComponent() { if (format!=NULL) { delete[] format; format=NULL; } }
  void showNumber(int val,int wid);
  void showFrame(int val,int wid) { identity=frame_n; showNumber(val,wid); }
  void showDancer(int val,int wid) { identity=dancer_n; showNumber(val,wid); }
  void showPerf(int val,int wid) { identity=perf_n; showNumber(val,wid); }
  void showNameDate(const char* trackerName, struct tm& date);
  void showText(const char *t);
  void numberToText();
  int toString(char *s,int n);
  static void compact(ManagedList<FilenameComponent>& fc);
  static int safeLength(ManagedList<FilenameComponent>& fc);
  static void toString(ManagedList<FilenameComponent>& fc,char *s,int n);
  static void loadValues(ManagedList<FilenameComponent>& fc,int frame_id,int dancer_id,int perf_id);
};


// A little class to accumulate statistics
class Datum
{
public:
  int number;
  float sum;
  float sumsq;
  Datum() : number(0),sum(0.0),sumsq(0.0) {}
  float mean() { return (number==0) ? 0.0 : sum/number; }
  float std() { return (number==0) ? 0.0 : (number==1) ? sum/number : sqrt( sum*sum/(number*number) - sumsq/number ); }
  float sem() { return (number==0) ? 0.0 : (number==1) ? sum/number : sqrt( (sum*sum/(number*number) - (sumsq/number)) / (number-1) ); }
  void add(float value) { number++; sum += value; sumsq += value*value; }
};


// A little class to compactly store contours
class PackedContour
{
public:
  short x_start;
  short y_start;
  volatile int size;  // Length of contour including start point
  unsigned char *bits;
  PackedContour() : x_start(0),y_start(0),size(0),bits(NULL) { }
  PackedContour(Contour& c);
  // Warning: something VERY sketchy is going on here; destructor gets called twice!
  // `volatile int size` provides a workaround, but why is it called twice???
  ~PackedContour() { if (bits != NULL && size>0) { delete[] bits; bits=NULL; size = 0; } }
  Contour* unpack(Contour* c);
  void printData(FILE* f,bool add_newline);
};




/****************************************************************
    Blobs, Moving Blobs, and Regions with Many Moving Blobs
****************************************************************/

// A single object that stands out from the background
class Blob
{
public:
  static const int SKELETON_SIZE = 11;  // Number of points in a skeleton; 11 points means object is divided into 10 segments

  int frame;
  double time;
  int pixel_count;  // Will be lost from stats when history is cleared, so keep it here
  Dancer* dancer;   // We will use member variables from dancer in order to avoid duplicating them each frame
  Listable<FloodData>* stats; // We'll only have one, but we make it listable so we can just grab it from the floodfill list
  FPoint jitter;    // Jitter in the image
  Image *im;
  Listable<Contour>* skeleton; // Listable for easier disposal
  Listable<Contour>* outline; // Listable for easier disposal
  Listable<PackedContour>* packedline;  // Listable for easier disposal
  
  
  Blob(int F=0,double T=0.0) : frame(F),time(T),pixel_count(0),dancer(NULL),stats(NULL),jitter(0,0),im(NULL),skeleton(NULL),outline(NULL),packedline(NULL) { }
  ~Blob();
  
  inline FloodData& data() { return stats->data; }
  inline Mask& mask() { return stats->data.stencil; }
  inline Contour& edge() { return outline->data; }
  inline Contour& spine() { return skeleton->data; }
  
  bool operator<(const Blob& b) const
  {
    if (stats==NULL) return true;
    else if (b.stats==NULL) return false;
    else return stats->data.stencil.bounds < b.stats->data.stencil.bounds;
  }
  
  //Image processing stuff
  void init(int f,double t,Dancer &d) { frame=f; time=t; dancer=&d; }
  void tossStats();
  Rectangle findNextROI(Image* bg,int border) const;
  void clip(const Blob& last,Image* fg,Image* bg,int border);
  void clip8(const Blob& last, Image8* fg, Image* bg, int border);
  int find(Mask* old_mask,Mask* exclusion_mask);
  void adoptCandidate();
  void flush();
  
  // Stats
  void extraStats(Range fill_I , FPoint previous_direction);
  void packOutline();
  
  // Output stuff
  bool print(FILE* f,const char* prefix);
};


// Follows one blob over time
class Dancer
{
public:
  static constexpr float MIN_OVERLAP_RATIO = 0.5;

  Performance* performance;  // Will use member variables relating to the whole arena, e.g. for reporting errors

  // Parameters specifying how to find the blob
  bool blob_is_dark;
  DualRange fill_I;
  DualRange fill_size;
  int border;
  bool find_skel;
  bool find_outline;
  
  // Data retention/saving
  int n_keep_full;
  int n_long_enough;
  char* data_fname;
  ManagedList<FilenameComponent>* img_namer;
  
  // Data management
  Storage< Listable<FloodData> > floodstore;
  Storage< Stackable<Strip> > ssstore;
  Storage< Listable<Strip> > lsstore;
  Storage< Listable<Point> > lpstore;
  Storage< Listable<Contour> > contourstore;
  Storage< Listable<PackedContour> > packedstore;
  
  // The actual data--movement of blob over time
  ManagedList<Blob> movie;
  Listable<Blob> *clear_til;
  int n_cleared;
  
  // Used to find the next blob
  ManagedList<FloodData> candidates;
  
  // Summary information
  int ID;
  Range frames;
  FRange times;
  float last_speed;
  float last_speed_time;
  float last_angularspeed;
  Datum accumulated_length;
  Datum accumulated_width;
  Datum accumulated_aspect;
  bool validated;
  
  Dancer() :  // If this constructor gets called somehow, we'll be using placement new to overwrite it anyway
    floodstore(0),ssstore(0),lsstore(0),lpstore(0),contourstore(0),packedstore(0),
    movie( (Storage< Listable<Blob> >*)NULL ),
    candidates( (Storage< Listable<FloodData> >*)NULL)
  { }
  Dancer(int n) :
    performance(NULL) ,
    blob_is_dark(true) ,
    fill_I( Range(1,Image::DEFAULT_GRAY-1) , Range(1,Image::DEFAULT_GRAY-1) ) ,
    fill_size( Range(1,3840*2400) , Range(1,3840*2400) ) ,
    border(2) , find_skel(false) , find_outline(false) , 
    n_keep_full(0) , n_long_enough(2) ,
    data_fname(NULL) , img_namer(NULL) ,
    floodstore(n,true) , ssstore(n,false) , lsstore(n,false) , lpstore(n,false) , contourstore(n,true) , packedstore(n,true) ,
    movie(n,true) , clear_til(NULL) , n_cleared(0) ,
    candidates(&floodstore) ,
    ID(0) , frames(0,0) , times(0,0) , 
    last_speed(0.0) , last_speed_time(0.0) , last_angularspeed(0.0) ,
    accumulated_length() , accumulated_width() , accumulated_aspect() ,
    validated(false)
  {	}
  Dancer(Performance *p,int n,bool bid,const DualRange& fI,const DualRange& fs,int b,int nkf,int nle,int id,int f,double t) :
    performance(p) ,
    blob_is_dark(bid) ,
    fill_I(fI) , fill_size(fs) ,
    border(b) , find_skel(false) , find_outline(false) , 
    n_keep_full(nkf) , n_long_enough( (nle<2)?2:nle ) ,
    data_fname(NULL) , img_namer(NULL) ,
    floodstore(n,true) , ssstore(n,false) , lsstore(n,false) , lpstore(11*n,false) , contourstore(n,true) , packedstore(n,true) ,
    movie(n,true) , clear_til(NULL) , n_cleared(0) ,
    candidates(&floodstore) ,
    ID(id) , frames(f,f) , times(t,t) , 
    last_speed(0.0) , last_speed_time(t) , last_angularspeed(0.0) ,
    accumulated_length() , accumulated_width() , accumulated_aspect() ,
    validated(false)
	{ }
  ~Dancer();
  
  bool operator<(const Dancer& d) const
  {
    if (movie.size==0) return true;
    else if (d.movie.size==0) return false;
    else return (movie.list_tail->data < d.movie.list_tail->data);  // Compare most recent position of dancers--can't use t() because it breaks const
  }
    
  // Image processing
  inline void resetMovie() { movie.flush(); clear_til=NULL; n_cleared=0; }
  void tossCandidates();
  void setFirst(Image* im,FloodData *fd,int frame,double time);
  Blob* makeFirst(int frame,double time);
  void readyAnother(Image *fg,Image *bg,int frame,double time);
  void readyAnother8(Image8 *fg, Image *bg, int frame, double time);
  bool findAnother(bool best_guess, Mask* exclusion_mask, FPoint jitter);
  void validate();
  void invalidate();
  
  // Cleanup and output
  void complain(const char* topic,const char* message);
  void tossFilenames();
  void tidyHistoryLeaving(int n_keep);
  void tidyHistory() { tidyHistoryLeaving(n_keep_full); }
  bool printData(FILE* f,const char* prefix);
  
  // Statistics reporting
  float recentSpeed(float min_interval = 0.0);  // Also calculates angular speed, so ask for it too!
  float recentAngularSpeed(float min_interval = 0.0) { recentSpeed(min_interval); return last_angularspeed; }
};


// Follows a bunch of moving blobs (dancers) and resolves conflicts
class Performance
{
public:
  static const int DEFAULT_BG_DEPTH = 15;
  static const int DEFAULT_ADAPT_RATE = 5;
  static const int DEFAULT_BUFFER_SIZE = 128;
  static const int DEFAULT_SCAN_BANDS = 16;
  static const int MAX_BLOBS_PER_FILE = 1000;
  
  enum ImageLoadState { no_state,check_sitter,load_sitter,check_dancer,load_dancer,check_band,load_band,all_loaded };

  // Flag for how output is handled.
  bool combine_blobs;
  int last_blobs_fid;
  int output_count;
  
  // Data for finding blobs--will load into each dance
  bool blob_is_dark;
  DualRange fill_I;
  DualRange fill_size;
  DualRange ref_I;
  int border;
  bool find_dancer_skel;
  bool find_dancer_edge;
  FPoint jitter;
  Point ijitter;
  Profile edge_profile;
  
  // Default data retention/saving policies
  int n_keep_full;
  int n_long_enough;

  FILE* output_file;
  char* base_directory;
  char* prefix_name;
  char* path_date_connector;

  ManagedList<FilenameComponent>* blobs_fname;
  ManagedList<FilenameComponent>* dance_fname;
  ManagedList<FilenameComponent>* sit_fname;
  ManagedList<FilenameComponent>* img_fname;
  ManagedList<BlobOutputFate> output_fates;

  
  // Image data
  int bg_depth;        // Depth in bits of background image
  int adapt_rate;      // Background is updated with a decay rate of 2^(-adapt_rate)
  Image* foreground;   // Main image (local copy)
  Image* background;   // Reference image to subtract from main image
  Mask* full_area;     // Place to do background subtraction and find objects
  Mask* danger_zone;   // Objects that intersect this should be excluded (inverted full_area + border)
  int n_scan_bands;    // Divide up image into this many strips for background updating
  Image* band;         // A strip for background updating and worm detection
  Mask* band_area;     // The mask covering the strip that excludes existing worms
  
  // Tracking data
  ManagedList<Dancer> sitters;          // Reference object(s)
  ManagedList<Dancer> dancers;          // Moving object(s)
  ManagedList<FloodData> candidates;    // New objects detected in image
  Storage< Listable<Strip> > lsstore;   // Storage for image masks
  Storage< Stackable<Strip> > ssstore;  // Storage for flood fills
  
  bool use_division; // Flag for type of image correction to use. TRUE for division, FALSE for subtraction.
	
  // Summary information
  struct tm* date;
  int ID;
  int next_dancer_ID;
  int current_frame;
  double current_time;
  int dance_buf_size;
  ImageLoadState load_state;
  int expected_n_dancers;       // Used to pad extended filename with zeros (per dancer count)
  int expected_n_sitters;       // Same (per sitter)
  int expected_n_performances;  // Used to pad root filename with zeros (performance count)
  int expected_n_frames;        // Used to pad with zeros (image filenames only, frame count)
  ManagedList<BlobOriginFate> fates; 
  ManagedList<IOErrorHandler> errors;
  
  Performance(int id,int n) :
    // Output flags
		combine_blobs(false), 
    last_blobs_fid(0),
    output_count(0),

    // Data for finding blobs
    blob_is_dark(true) ,
    fill_I( Range(1,Image::DEFAULT_GRAY-1) , Range(1,Image::DEFAULT_GRAY-1) ),
    fill_size( Range(1,3840*2400) , Range(1,3840*2400) ),
    ref_I( Range(1,Image::DEFAULT_GRAY/2) , Range(1,Image::DEFAULT_GRAY/2) ),
    border(2),
    find_dancer_skel(false),
    find_dancer_edge(false),
    jitter(0, 0),
    ijitter(0, 0),
    edge_profile(Rectangle(0, 0, 0, 0), Profile::OverX),

    // Data retention/saving
    n_keep_full(2),
    n_long_enough(2),
    output_file(NULL),
    base_directory(NULL),
    prefix_name(NULL),
    path_date_connector(NULL),
    blobs_fname(NULL),
    dance_fname(NULL),
    sit_fname(NULL),
    img_fname(NULL),
    output_fates(n,false),

    // Image data 
    bg_depth(DEFAULT_BG_DEPTH),
    adapt_rate(DEFAULT_ADAPT_RATE),
    foreground(NULL),
    background(NULL),
    full_area(NULL),
    danger_zone(NULL),
    n_scan_bands(DEFAULT_SCAN_BANDS),
    band(NULL),
    band_area(NULL),

    // Tracking data
    sitters(n,true),
    dancers(n,true),
    candidates(n,true),
    lsstore(DEFAULT_BUFFER_SIZE*n,false),
    ssstore(DEFAULT_BUFFER_SIZE,false),
    use_division(false),

    // Summary information
    date(NULL),
    ID(id),
    next_dancer_ID(1),
    current_frame(0),
    current_time(0.0),
    dance_buf_size(DEFAULT_BUFFER_SIZE),
    load_state(no_state) ,
    expected_n_dancers(99999),
    expected_n_sitters(9),
    expected_n_performances(1),
    expected_n_frames(999999),
    fates(n,false), 
    errors(16,true)
  {
    path_date_connector = strdup("/");
  }
  Performance(): Performance(0, 2) {}
  ~Performance();

  // Setting parameters
  void setCombineBlobs( bool flag )
  {
      combine_blobs = flag;
  }

  FILE* getFileHandle();

  void addOutputFate(int frame, int idb, long byte_offset); 
  
  void unSetROI()
  { 
    if (full_area!=NULL) { delete full_area; full_area=NULL; }
    if (danger_zone!=NULL) { delete danger_zone; danger_zone=NULL; }
  }
  void setROI(const Rectangle& r);
  void setROI(const Ellipse& e);
  void setROI(Mask& m);
  void setNotInROI(const Image *im);
  
  // Image scanning
  int scanImage(Image* im);  // Lower level, called by latter functions
  void adoptScan(Image* im); // Also lower level
  int initialScan(Image *fg,double time);
  int initialRefs(Image *fg,ManagedList<Point>& locations,double time);
  // 8 bit versions
  int initialScan8(Image8 *fg, double time);
  int initialRefs8(Image8 *fg, ManagedList<Point>& locations, double time);
  
private:
  void calculateJitterDeltaSort(float delta, float *xedge, int &nxe, float *yedge, int &nye);
  void calculateJitterGivenDeltaSorted(float *xedge, int nxe, float *yedge, int nye);
public:
  void calculateJitter(Image* fg, ManagedList<Profile> &edges);
  void calculateJitter8(Image8* fg, ManagedList<Profile> &edges);
  int anticipateNext(double time);
  bool findNextItemBounds(Rectangle& im_bound);
  void loadNextSingleItem(Image* fg);
  void loadNextSingleItem8(Image8* fg);
  void readyNext(Image *fg, ManagedList<Profile> &edges, double time);
  void readyNext8(Image8 *fg, ManagedList<Profile> &edges, double time);
  int findNext();
  
	// Image correction 
	double generateBackgroundAverageIntensity( int *array );	
	
  // Output (mostly done during scanning--this is just for final cleanup)
  void enableOutput(Dancer& d,bool sitting=false);
  bool prepareOutput(const char* trackerName, const char* path,const char *prefix,bool save_dance,bool save_sit,bool save_img,struct tm* date_to_use);
  bool logErrors(char* err_fname);
  bool finishOutput();
  void imprint(Image* im,short borderI,int borderW,short maskI,int maskW,short dancerI,bool show_dancer,
    short sitterI,bool show_sitter,short dcenterI,int dcenterR);
  void imprint8(Image8* im, uint8_t borderI, int borderW, uint8_t maskI, int maskW, uint8_t dancerI, bool show_dancer,
    uint8_t sitterI, bool show_sitter, uint8_t dcenterI, int dcenterR);
};



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_blob_misc();
int test_mwt_blob_blob();
int test_mwt_blob_dancer();
int test_mwt_blob_performance();
int test_mwt_blob();

#endif

