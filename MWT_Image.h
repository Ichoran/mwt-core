/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */

#ifndef MWT_IMAGE
#define MWT_IMAGE

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "MWT_Geometry.h"
#include "MWT_Lists.h"
#include "MWT_Storage.h"


/****************************************************************
                      Pieces of an Image
****************************************************************/

// Holds a strip that marks a single portion of one line of an image.
class Strip
{
public:
  short x;
  short y0,y1;
  short id;
  
  Strip() {}
  Strip(short X,short Y0,short Y1) : x(X),y0(Y0),y1(Y1) { }
  Strip(short X,short Y0,short Y1,short ID) : x(X),y0(Y0),y1(Y1),id(ID) { }
  
  inline void set(short X,short Y0,short Y1) { x=X; y0=Y0; y1=Y1; }
  inline void set(short X,short Y0,short Y1,short ID) { x=X; y0=Y0; y1=Y1; id=ID; }
  
  inline short length() const { return 1+y1-y0; }
  
  // Tell if two strips touch each other
  inline bool overlaps(const Strip& s) const { return (x==s.x && !(y0>s.y1 || y1<s.y0)); }
  inline bool beyond(const Strip& s) const { return (x>s.x || (x==s.x && y0>s.y1)); }
  inline bool before(const Strip& s) const { return (x<s.x || (x==s.x && y1<s.y0)); }
  // Same thing, except with an additional offset
  inline bool overlaps(const Strip& s,const Point& p) { return (x==s.x+p.x && !(y0>s.y1+p.y || y1<s.y0+p.y)); }
  inline bool beyond(const Strip& s,const Point& p) { return (x>s.x+p.x || (x==s.x+p.x && y0>s.y1+p.y)); }
  inline bool before(const Strip& s,const Point& p) { return (x<s.x+p.x || (x==s.x+p.x && y1<s.y0+p.y)); }
  // Same thing, except comparing a strip to points only
  inline bool contains(const Point& p) const { return (x==p.x && y0<=p.y && y1>=p.y); }
  inline bool beyond(const Point& p) const { return (x>p.x || y0>p.y); }
  inline bool before(const Point& p) const { return (x<p.x || y1<p.y); }

  bool operator==(const Strip& s) const { return (x==s.x && y0==s.y0 && y1==s.y1); }
  bool operator!=(const Strip& s) const { return (x!=s.x || y0!=s.y0 || y1!=s.y1); }
  bool operator<(const Strip& s) const { return (x<s.x || (x==s.x && y0<s.y0)); }
  inline Strip& operator+=(const Point& p) { x+=p.x; y0+=p.y; y1+=p.y; return *this; }
  inline Strip& operator-=(const Point& p) { x-=p.x; y0-=p.y; y1-=p.y; return *this; }
  inline Strip operator+(const Point& p) const { return Strip(x+p.x,y0+p.y,y1+p.y,id); }
  inline Strip operator-(const Point& p) const { return Strip(x-p.x,y0-p.y,y1-p.y,id); }
};
inline Strip operator+(const Point& p,const Strip& s) { return s+p; }



/****************************************************************
                         An Image Mask
****************************************************************/

/* NOTE **
        * Masks are sorted sets of strips.  Strips aren't sorted when they
        * are added, so be sure to sort() after all strips have been
        * added!  Also, make sure that new strips don't overlap existing
        * ones--if they might, put the new strip into a mask and use +=.
** NOTE **
        * There are no operators that return a temporary copy.  If you need
        * a new mask that is formed from two others, use imitate() on the
        * first and then the +=, -=, etc. set of operators.
** NOTE **
        * The internal lines iterator is used by the overlap/contains
        * functions.  Make your own Listable<Strip> if you want one
        * that will be preserved across calls.
*********/

// Strips of pixels designating the extent of an object or image mask
class Mask
{
public:
  ManagedList<Strip> lines;
  int pixel_count;
  Rectangle bounds;
  
  // Create empty masks or those with a geometrical shape
  Mask(int nlines) : lines(nlines,false),pixel_count(0) { }
  Mask(Storage< Listable<Strip> >* store) : lines(store),pixel_count(0) { }
  Mask(const Rectangle& rect , Storage< Listable<Strip> >* store);
  Mask(const Ellipse& ellip , Storage< Listable<Strip> >* store);
  Mask(const Mask& existing , ManagedList<Strip>::CreationType action) // Careful with this one!  Shares storage with existing.
    : lines(existing.lines,action)
  {
    if (action==ManagedList<Strip>::SubordinateCopy)
    {
      pixel_count = existing.pixel_count;
      bounds = existing.bounds;
    }
    else pixel_count = 0;
  }
  ~Mask() { }
  
  // Shallow transfer of lines
  Mask& operator=(Mask& m) { lines=m.lines; pixel_count=m.pixel_count; bounds=m.bounds; return *this; }
  // Deep copy of lines
  Mask& imitate(const Mask& m) { lines.imitate(m.lines); pixel_count=m.pixel_count; bounds=m.bounds; return *this; }
  // Throw away all lines
  inline void flush() { lines.flush(); pixel_count = 0; }
  
  // Methods to make mask adopt geometrical shapes
  void rectangularize(const Rectangle& rect);
  void ellipticize(const Ellipse& ellip);
  
  // Shortcut methods for moving through and removing the strips in the mask
  inline void start() { lines.start(); }
  inline Strip& i() { return lines.i(); }
  inline Strip& ii(Listable<Strip>* ls) { return ls->data; }  // Syntactic sugar for using a separate iterator
  inline bool advance() { return lines.advance(); }
  inline bool retreat() { return lines.retreat(); }
  inline void end() { lines.end(); }
  inline void Backspace() { pixel_count -= i().length(); lines.Backspace(); }
  inline void Delete() { pixel_count -= i().length(); lines.Delete(); }
  
  // Be sure to sort after adding strips if you don't add them in order
  // add adds to end of list, push adds to front
  inline void addStrip(const Strip& s)
  {
    lines.Append(s);
    pixel_count += s.length();
  }
  inline void addStrip(short x,short y0,short y1)
  {
    Strip& s = *( lines.Append() );
    s.set(x,y0,y1);
    pixel_count += s.length();
  }
  inline void pushStrip(const Strip& s)
  {
    lines.Push(s);
    pixel_count += s.length();
  }
  inline void pushStrip(short x,short y0,short y1)
  {
    Strip& s = *( lines.Push() );
    s.set(x,y0,y1);
    pixel_count += s.length();
  }
  inline void sort() { lines.mergeSort(); }
  
  // Handy utility classes
  void findBounds(); // Finds bounding box
  void tidy();       // Removes any overlapping strips and corrects pixel count (must be sorted)
  
  // Everything after this point assumes that the strips are sorted and nonoverlapping and will keep them that way
  
  // Check whether a point or sorted list of points are in a mask
  bool contains(const Point& p);
  Stackable<Point>* contains(Stackable<Point>*& sp);
  Listable<Point>* contains(Listable<Point>*& sp);
  int contains(ManagedList<Point>& to_test,ManagedList<Point>& inside);  // Everything starts in to_test, those remaining are out
  
  // See if there is any overlap at all
  bool overlaps(Mask& m,Point offset);
  bool overlaps(Mask& m) { return overlaps(m,Point(0,0)); }
  
  // Return (fraction we are overlapped by them , fraction they are overlapped by us)
  FPoint overlapRatio(Mask& m);
  
  // Various trimming functions; only use after strips are sorted (these keep the strips sorted)
  void cropTo(const Rectangle& rect);      // Make the mask fit inside the rectangle
  int countOverlap(Mask& m);               // Count how many pixels overlap in the two masks
  Mask& operator+=(Mask& m);               // Union of the two masks
  Mask& operator-=(Mask& m);               // Remove the second mask
  Mask& operator-=(const Rectangle &rect); // Remove any of the mask in a rectangle
  Mask& operator*=(Mask& m);               // Intersection
  Mask& operator*=(const Rectangle &rect) { cropTo(rect); return *this; }  // Intersection with a rectangle--same as cropTo
  
  // Translate the mask by a vector (sorted strips not strictly necessary)
  Mask& operator+=(const Point& p)
  {
    lines.start();
    while (lines.advance()) lines.i() += p;
    return *this;
  }
  Mask& operator-=(const Point& p)
  {
    lines.start();
    while (lines.advance()) lines.i() -= p;
    return *this;
  }
  
  // Morphological operations normally applied to images
  void extractEdge(Mask& edges);
  void invert(const Rectangle& universe);
  void dilate(int n);
  void dilate(Mask& structuringElt);
  void erode(int n);
  void erode(Mask& structuringElt);
  void close(int n) { dilate(n); erode(n); }
  void close(Mask& structuringElt) { dilate(structuringElt); erode(structuringElt); }
  void open(int n) { erode(n); dilate(n); }
  void open(Mask& structuringElt) { erode(structuringElt); dilate(structuringElt); }
  
  // Text output of mask (with newline)
  void println();
};



/****************************************************************
                      A Contour of Points
****************************************************************/


// A single contour of points (assumed to be closed)
class Contour
{
protected:
  class Scorer  // Used for creating approximate contours
  {
  public:
    Listable<Point> *pt;
    float score;
    Scorer() { }
    bool operator<(const Scorer& s) { return (score < s.score); }
  };
public:
  ManagedList<Point> boundary;
  Rectangle bounds;
  
  Contour(int npix=0) : boundary(npix,false) { }
  Contour(Storage< Listable< Point> >* store) : boundary(store) { }
  Contour(const Rectangle& rect) : boundary(2*(rect.width()+rect.height()) , false)
  { findContour(rect); }
  Contour(const Ellipse& ellip) : boundary(4*(ellip.axes.x+ellip.axes.y) , false)
  { findContour(ellip); }
  Contour(Mask& m) : boundary(2*m.lines.size , false)
  { findContour(m); }
  ~Contour() { }
  
  // Iterator access
  inline Listable<Point>* getIterator() { return boundary.current; }
  inline void setIterator(Listable<Point>* lp) { boundary.current = lp; }
  // Iterators that make contour look linear
  inline void start() { boundary.start(); }
  inline bool advance() { return boundary.advance(); }
  inline bool retreat() { return boundary.retreat(); }
  inline void end() { boundary.end(); }
  // Iterators that make contour look circular
  inline bool first()
  { 
    if (boundary.size<=0) return false;
    else { boundary.current = boundary.list_head; return true; }
  }
  inline bool next()
  {
    if (boundary.advance()) return true;
    else { boundary.current = boundary.list_head; return false; }
  }
  inline bool prev()
  {
    if (boundary.retreat()) return true;
    else { boundary.current = boundary.list_tail; return false; }
  }
  inline bool last()
  {
    if (boundary.size<=0) return false;
    else { boundary.current = boundary.list_tail; return true; }
  }
  
  // Accessors
  inline Point& h() { return boundary.h(); }
  inline Point& t() { return boundary.t(); }
  inline Point& i() { return boundary.i(); }
  
  // Functions that insert single pixels
  inline void Push(const Point& p) { boundary.Push(p); }
  inline void Append(const Point& p) { boundary.Append(p); }
  inline void Insert(const Point& p) { boundary.Insert(p); }
  inline void Tuck(const Point& p) { boundary.Tuck(p); }
  inline void Backspace() { boundary.Backspace(); }
  inline void Delete() { boundary.Delete(); }
  
  //  Functions that build/operate on the whole contour
  inline int size() const { return boundary.size; }
  inline void flush() { boundary.flush(); }
  inline Contour& imitate(const Contour& c) { if (c.size()>0) { boundary.imitate(c.boundary); } bounds=c.bounds; return *this; }
  void findContour(const Rectangle& r,bool full=true);  // Can be either full or empty representation initially (default full)
  void findContour(const Ellipse& e);
  void findContour(Mask& m);
  void fillContour();
  void emptyContour();
  
  // Convert contour into various approximations of the full contour; assumes a single closed contour
  void approximateSpine(FPoint centroid,FPoint head_dir,int n_backbone_points);  // Rough linear spine down the center of the object
  void approximateContour(float max_err);        // Cubic fit to sets of 4 points will be within max_err pixels of original
  
  // Statistics
  void findBounds();
};



/****************************************************************
             Image Statistics and Data Structures
****************************************************************/

class Image;  // We'll be defining this later but point to it now


// Ranges of numbers
class Range : public Point
{
public:
  Range() { }
  Range(int L,int H) : Point(L,H) { }
  inline int length() { return 1+y-x; }
  inline int& lo() { return x; }
  inline int& hi() { return y; }
  inline bool excludes(int I) const { return (I<x || I>y); }
};

class FRange : public FPoint
{
public:
  FRange() { }
  FRange(float L,float H) : FPoint(L,H) { }
  inline float length() { return y-x; }
  inline float& lo() { return x; }
  inline float& hi() { return y; }
  inline bool excludes(float f) const { return (f<x || f>y); }
};

class DualRange
{
public:
  Range start;
  Range keep;
  DualRange() { }
  DualRange(Range R) : start(R) , keep(R) { }
  DualRange(Range S,Range K) : start(S),keep(K) { }
  DualRange(int l,int h) : start(l,h),keep(l,h) { }
  DualRange(int lS,int hS,int lK,int hK) : start(lS,hS),keep(lK,hK) { }
};

class DualFRange
{
public:
  FRange start;
  FRange keep;
  DualFRange() { }
  DualFRange(FRange S,FRange K) : start(S),keep(K) { }
  DualFRange(float lS,float hS,float lK,float hK) : start(lS,hS),keep(lK,hK) { }
};


// Useful stuff to know while we are filling
class FloodInfo
{
public:
  Image *im;          // Image we're trying to fill
  Range T;       // Range of valid values for fill
  Rectangle bounds;   // Spatial extent of fill
  long long cx,cy;    // Sum of intensity-weighted x and y coordinates
  long long sumI;     // Total intensity in object
  int pixel_count;    // Keep track if we don't have a FloodData to do it
  bool local_storage; // Set if we allocated our own storage
  Storage< Stackable<Strip> >* store;   // Place to get our filled strips from
  
  FloodInfo(Image* image , Storage< Stackable<Strip> >* storage);
  ~FloodInfo() { if (local_storage) delete store; }
  
  void floodPixel(int x,int y,short I);
  void freshen();
};


// Useful stuff to know about something we've filled
class FloodData
{
public:
  float score;             // Used/reused for various things
  FPoint centroid;         // Center of intensity
  FPoint major,minor;      // Least-squares axes (normally not filled in)
  Mask stencil;            // Actual places where pixels are
  float long_axis;         // Length of object along major axis
  float short_axis;        // Length of object along minor axis
  
  FloodData(Storage< Listable<Strip> >* sls=NULL) : stencil(sls) { }
  ~FloodData() { }
  
  bool operator<(const FloodData& fd) const { return score < fd.score; }
  void findCentroid(const Image& im,Range threshold);  // Find centroid (if not already done by flood algorithm)
  void principalAxes(const Image& im,Range threshold);  // Fill in major & minor axis vectors + lengths
  void duplicate(const FloodData& fd);
};


// Standard tag for TIFF files
class TiffIFD
{
public:
  short ifd_tag;
  short ifd_type;
  int ifd_number;
  int ifd_value;
  TiffIFD(short tag,short typ,int num,int val) : ifd_tag(tag),ifd_type(typ),ifd_number(num),ifd_value(val) {}
};




/****************************************************************
                        The Image Itself
****************************************************************/


/* HOWTO **
**** How Image deals with memory ****
  Images can have either one or two arrays of pixels; these are arrays of
  shorts either allocated with new (reclaim with delete[]), or allocated
  elsewhere.  If one image is assigned to another, the one on the left will
  inherit the memory.  If an image is managing its own memory, it will
  delete[] the array(s) when the destructor is called.
** HOWTO */
 
// Holds image data (as shorts), optionally with binning

class Image8; // Forward declaration of 8-bit images (needed so we can use them)

class Image
{
public:
  static const short DEFAULT_BIT_DEPTH = 10;
  static const short DEFAULT_GRAY = 1 << (DEFAULT_BIT_DEPTH-1);   // Ten bits
  
  enum ScaleType { Subsample , LinearFit };
  
  short *pixels;
  
  Rectangle bounds;
  Point size;
  
  short bin;
  short depth;
  bool owns_pixels;
  
  bool divide_bg; 

  // Creating and disposing of images
  Image() : pixels(NULL),bounds(),size(0,0),bin(0),depth(DEFAULT_BIT_DEPTH),owns_pixels(false),divide_bg(false) { }
  Image(short *raw_image,Point image_size, bool algo)
    : pixels(raw_image),bounds(Point(0,0),image_size-1),size(image_size),
      bin(0),depth(DEFAULT_BIT_DEPTH),owns_pixels(false),divide_bg(algo) { }
  Image(Point image_size, bool algo) : size(image_size),bin(0),depth(DEFAULT_BIT_DEPTH),owns_pixels(true),divide_bg(algo)
  {
    if (size.x<1) size.x=1;
    if (size.y<1) size.y=1;
    bounds = Rectangle(Point(0,0),size-1);
    pixels = new short[ bounds.area() ];
  }
  Image(Rectangle image_bounds, bool algo) : bounds(image_bounds),bin(0),depth(DEFAULT_BIT_DEPTH),owns_pixels(true),divide_bg(algo)
  {
    if (bounds.near.x > bounds.far.x) bounds.far.x = bounds.near.x;
    if (bounds.near.y > bounds.far.y) bounds.far.y = bounds.near.y;
    size = bounds.size();
    pixels = new short[ bounds.area() ];
  }
  Image(const Image& existing,const Rectangle &region, bool algo);
  ~Image()
  {
    if (owns_pixels && pixels!=NULL) { delete[] pixels; pixels=NULL; }
    owns_pixels = false;
  }
  
  // Stuff for binned images
  Rectangle getBounds() const { if (bin<=1) return bounds; else return Rectangle( (bounds.near+(bin-1))/bin , (bounds.far-(bin-1))/bin ); }
  int getWidth() const { if (bin<=1) return bounds.width(); else return 1 + (bounds.far.x-(bin-1))/bin - (bounds.near.x+(bin-1))/bin; }
  int getHeight() const { if (bin<=1) return bounds.height(); else return 1 + (bounds.far.y-(bin-1))/bin - (bounds.near.y+(bin-1))/bin; }
  inline int getGray() const { return 1<<(depth-1); }
  
  // Assigning images--pixels are assigned by reference, use clone() to duplicate, copy to transfer pixels
  Image& operator=(Image &im)
  {
    pixels = im.pixels;
    bounds = im.bounds;
    size = im.size;
    bin = im.bin;
    owns_pixels = im.owns_pixels;
    im.owns_pixels = false;
    return *this;
  }
  Image* clone()
  {
    Image* im = new Image(*this,getBounds(),divide_bg);
    im->bounds = bounds;
    return im;
  }

  // Creates an identical copy of the image with its own memory etc but everything else the same
  Image* identical() const;
  
  
  // Access to the underlying data
  inline short& raw(int x,int y) { return pixels[y+size.y*x]; }
  inline short& raw(Point p) { return pixels[p.y+size.y*p.x]; }
  inline short peek(int x,int y) const { return pixels[y+size.y*x]; }
  inline short peek(Point p) const { return pixels[p.y+size.y*p.x]; }
  // Access to the underlying data with offset
  inline short& rare(int x,int y) { return pixels[(y-bounds.near.y) + size.y*(x-bounds.near.x)]; }
  inline short& rare(Point p) { return pixels[(p.y-bounds.near.y) + size.y*(p.x-bounds.near.x)]; }
  inline short view(int x,int y) const { return pixels[(y-bounds.near.y) + size.y*(x-bounds.near.x)]; }
  inline short view(Point p) const { return peek(p-bounds.near); }
  // "Vector" access to underlying data (packed in ints)
  inline int& rareI(int x,int y) { return *((int*)(pixels + ((y-bounds.near.y) + size.y*(x-bounds.near.x)))); }
  inline int viewI(int x,int y) const { return *((int*)(pixels + ((y-bounds.near.y) + size.y*(x-bounds.near.x)))); }
  // Access to data in global coordinates, with binning as needed
  inline short get(int x,int y) const
  {
    if (bin<=1) return peek(x-bounds.near.x,y-bounds.near.y);
    else
    {
      int I=0;
      x = bin*x-bounds.near.x;
      y = bin*y-bounds.near.y;
      for (int i=0;i<bin;i++) for (int j=0;j<bin;j++) I += peek(x+i,y+j);
      return I/(bin*bin);
    }
  }
  inline short get(Point p) const
  {
    if (bin<=1) return peek(p-bounds.near);
    else
    {
      int I=0;
      Point q;
      p *= bin;
      p -= bounds.near;
      for (q.x=0;q.x<bin;q.x++) for (q.y=0;q.y<bin;q.y++) I += peek(p+q);
      return I/(bin*bin);
    }
  }
  void set(int x,int y,short I)
  {
    if (bin <= 1) raw(x-bounds.near.x,y-bounds.near.y) = I;
    else
    {
      x = bin*x-bounds.near.x;
      y = bin*y-bounds.near.y;
      for (int i=0;i<bin;i++) for (int j=0;j<bin;j++) raw(x+i,y+j) = I;
    }
  }
  void set(Point p,short I)
  {
    if (bin <= 1) raw(p-bounds.near) = I;
    else
    {
      Point q;
      p *= bin;
      p -= bounds.near;
      for (q.x=0;q.x<bin;q.x++) for (q.y=0;q.y<bin;q.y++) raw(p+q) = I;
    }
  }
  void set(Rectangle r,short I);
  void set(Mask &m,short I);
  void set(Contour& c,short I)
  {
    Point p;
    Rectangle b = getBounds();
    c.start();
    while (c.advance()) { if (b.contains(c.i())) set( c.i() , I ); }
  }
  void set(Point p1,Point p2,int width,short I); // Sets a line of specified radius (0=single pixel) from p1 up to but not including p2
  void set(Contour &c,int width,short I);        // Draws a contour as lines--if contour is not made of adjacent pixels, it's still connected
  
  // Operations to shift and scale image by a constant
  inline Image& operator=(short I) { for (int i=0;i<size.x*size.y;i++) pixels[i]=I; return *this; }
  inline Image& operator+=(short I) { for (int i=0;i<size.x*size.y;i++) pixels[i]+=I; return *this; }
  inline Image& operator+=(Point p) { bounds += p; return *this; }
  inline Image& operator-=(short I) { for (int i=0;i<size.x*size.y;i++) pixels[i]-=I; return *this; }
  inline Image& operator-=(Point p) { bounds -= p; return *this; }
  inline Image& operator<<=(short I) { for (int i=0;i<size.x*size.y;i++) pixels[i]<<=I; return *this; }
  inline Image& operator>>=(short I) { for (int i=0;i<size.x*size.y;i++) pixels[i]>>=I; return *this; }
  
  // Copy bits of one image to another in various ways
  void mimic(const Image& source , Rectangle my_region , Rectangle source_region , ScaleType method = Subsample);
  void mimic(const Image& source , ScaleType method = Subsample ) { mimic(source,getBounds(),source.getBounds(),method); }
  void copy(Point where,const Image& source,Point size,bool fix_depth=false);
  void copy(const Image& source,Mask& m,bool fix_depth=false);
  void copy(const Image& source,bool fix_depth=false) { copy( source.bounds.near , source , source.size , fix_depth ); }
  void adapt(Point where,const Image& im,Point size,int rate); // Adapting image has deeper bit depth
  void adapt(const Image& im,Mask &m,int rate);
  void diffCopy(Point where,const Image& source,Point size,const Image& bg);  // Diffcopied image is effectively one bit deeper than source, bg is much deeper
  void diffCopy(const Image& source,Mask& m,const Image& bg);
  void diffAdaptCopy(Point where,const Image& source,Point size,Image& bg,int rate);  // Same as diffCopy (but bg gets adapted)
  void diffAdaptCopy(const Image& source,Mask& m,Image& bg,int rate);

  // Same as above but using 8 bit images as a source
  void mimic8(const Image8& source , Rectangle my_region , Rectangle source_region , ScaleType method = Subsample);
  void mimic8(const Image8& source , ScaleType method = Subsample );
  void copy8(Point where,const Image8& source,Point size,bool fix_depth=false);
  void copy8(const Image8& source,Mask& m,bool fix_depth=false);
  void copy8(const Image8& source,bool fix_depth=false);
  void adapt8(Point where,const Image8& im,Point size,int rate); // Adapting image has deeper bit depth
  void adapt8(const Image8& im,Mask &m,int rate);
  void diffCopy8(Point where,const Image8& source,Point size,const Image& bg);  // Diffcopied image is effectively one bit deeper than source, bg is much deeper
  void diffCopy8(const Image8& source,Mask& m,const Image& bg);
  void diffAdaptCopy8(Point where,const Image8& source,Point size,Image& bg,int rate);  // Same as diffCopy (but bg gets adapted)
  void diffAdaptCopy8(const Image8& source,Mask& m,Image& bg,int rate);

  // Projections
  void maxOverX(Rectangle target, short *output);
  void maxOverY(Rectangle target, short *output);
  void minOverX(Rectangle target, short *output);
  void minOverY(Rectangle target, short *output);
  void meanOverX(Rectangle target, float* output);
  void meanOverY(Rectangle target, float* output);
  void deviationOverX(Rectangle target, const float* means, float* output);
  void deviationOverY(Rectangle target, const float* means, float* output);

  // Flood fills--note that resulting strips are unsorted
  void floodLine(FloodInfo& info,FloodData* data,Stackable<Strip>*& head,Stackable<Strip>*& tail);
  int floodFind(Point pt,DualRange threshold,Storage< Stackable<Strip> > *store,FloodData *result);
  int floodRect(DualRange threshold,Storage< Stackable<Strip> > *stripstore,Storage< Listable<Strip> > *storemask,
                ManagedList<FloodData>& result,Rectangle rect);
  int floodMask(DualRange threshold,Storage< Stackable<Strip> > *stripstore,Storage< Listable<Strip> > *storemask,
                ManagedList<FloodData>& result,Mask &m);
  
  // Output and testing
  int makeTiffHeader(unsigned char *buffer);
  int writeTiff(FILE *f);
  int writeTiff(const char *fname);
  void println() const;
};


/****************************************************************
                      An Eight-Bit Image
****************************************************************/


/* HOWTO **
**** How Image8 deals with memory ****
  Image8 can have either one or two arrays of pixels; these are arrays of
  bytes either allocated with new (reclaim with delete[]), or allocated
  elsewhere.  If one image is assigned to another, the one on the left will
  inherit the memory.  If an image is managing its own memory, it will
  delete[] the array(s) when the destructor is called.

  Unlike an Image16, all images must have a bit depth of exactly 8
** HOWTO */
 
// Holds image data (as bytes), optionally with binning
class Image8
{
public:
  
  enum ScaleType { Subsample , LinearFit };
  
  uint8_t *pixels;
  
  Rectangle bounds;
  Point size;
  
  short bin;
  bool owns_pixels;
  
  bool divide_bg;

  // Creating and disposing of images
  Image8() : pixels(NULL),bounds(),size(0,0),bin(0),owns_pixels(false),divide_bg(false) { }
  Image8(uint8_t *raw_image,Point image_size, bool algo)
    : pixels(raw_image),bounds(Point(0,0),image_size-1),size(image_size),
      bin(0),owns_pixels(false),divide_bg(algo) { }
  Image8(Point image_size, bool algo) : size(image_size),bin(0),owns_pixels(true),divide_bg(algo)
  {
    if (size.x<1) size.x=1;
    if (size.y<1) size.y=1;
    bounds = Rectangle(Point(0,0),size-1);
    pixels = new uint8_t[ bounds.area() ];
  }
  Image8(Rectangle image_bounds, bool algo) : bounds(image_bounds),bin(0),owns_pixels(true),divide_bg(algo)
  {
    if (bounds.near.x > bounds.far.x) bounds.far.x = bounds.near.x;
    if (bounds.near.y > bounds.far.y) bounds.far.y = bounds.near.y;
    size = bounds.size();
    pixels = new uint8_t[ bounds.area() ];
  }
  Image8(const Image8& existing,const Rectangle &region, bool algo);
  ~Image8()
  {
    if (owns_pixels && pixels!=NULL) { delete[] pixels; pixels=NULL; }
    owns_pixels = false;
  }
  
  // Stuff for binned images
  Rectangle getBounds() const { if (bin<=1) return bounds; else return Rectangle( (bounds.near+(bin-1))/bin , (bounds.far-(bin-1))/bin ); }
  int getWidth() const { if (bin<=1) return bounds.width(); else return 1 + (bounds.far.x-(bin-1))/bin - (bounds.near.x+(bin-1))/bin; }
  int getHeight() const { if (bin<=1) return bounds.height(); else return 1 + (bounds.far.y-(bin-1))/bin - (bounds.near.y+(bin-1))/bin; }
  
  // Assigning images--pixels are assigned by reference, use clone() to duplicate, copy to transfer pixels
  Image8& operator=(Image8 &im) {
    pixels = im.pixels;
    bounds = im.bounds;
    size = im.size;
    bin = im.bin;
    owns_pixels = im.owns_pixels;
    im.owns_pixels = false;
    return *this;
  }
  Image8* clone() {
    Image8* im = new Image8(*this,getBounds(),divide_bg);
    im->bounds = bounds;
    return im;
  }

  // Create an identical 16-bit-storage-depth image
  Image* identical16() const;
  
  
  // Access to the underlying data
  inline uint8_t& raw(int x,int y) { return pixels[y+size.y*x]; }
  inline uint8_t& raw(Point p) { return pixels[p.y+size.y*p.x]; }
  inline uint8_t peek(int x,int y) const { return pixels[y+size.y*x]; }
  inline uint8_t peek(Point p) const { return pixels[p.y+size.y*p.x]; }
  // Access to the underlying data with offset
  inline uint8_t& rare(int x,int y) { return pixels[(y-bounds.near.y) + size.y*(x-bounds.near.x)]; }
  inline uint8_t& rare(Point p) { return pixels[(p.y-bounds.near.y) + size.y*(p.x-bounds.near.x)]; }
  inline uint8_t view(int x,int y) const { return pixels[(y-bounds.near.y) + size.y*(x-bounds.near.x)]; }
  inline uint8_t view(Point p) const { return peek(p-bounds.near); }
  // "Vector" access to underlying data (packed in ints)
  inline unsigned short& rareS(int x,int y) { return *((unsigned short*)(pixels + ((y-bounds.near.y) + size.y*(x-bounds.near.x)))); }
  inline unsigned short viewS(int x,int y) const { return *((unsigned short*)(pixels + ((y-bounds.near.y) + size.y*(x-bounds.near.x)))); }
  // Access to data in global coordinates, with binning as needed
  inline uint8_t get(int x,int y) const {
    if (bin<=1) return peek(x-bounds.near.x,y-bounds.near.y);
    else {
      int I=0;
      x = bin*x-bounds.near.x;
      y = bin*y-bounds.near.y;
      for (int i=0;i<bin;i++) for (int j=0;j<bin;j++) I += peek(x+i,y+j);
      return I/(bin*bin);
    }
  }
  inline uint8_t get(Point p) const {
    if (bin<=1) return peek(p-bounds.near);
    else {
      int I=0;
      Point q;
      p *= bin;
      p -= bounds.near;
      for (q.x=0;q.x<bin;q.x++) for (q.y=0;q.y<bin;q.y++) I += peek(p+q);
      return I/(bin*bin);
    }
  }
  void set(int x,int y,uint8_t I) {
    if (bin <= 1) raw(x-bounds.near.x,y-bounds.near.y) = I;
    else {
      x = bin*x-bounds.near.x;
      y = bin*y-bounds.near.y;
      for (int i=0;i<bin;i++) for (int j=0;j<bin;j++) raw(x+i,y+j) = I;
    }
  }
  void set(Point p,uint8_t I) {
    if (bin <= 1) raw(p-bounds.near) = I;
    else {
      Point q;
      p *= bin;
      p -= bounds.near;
      for (q.x=0;q.x<bin;q.x++) for (q.y=0;q.y<bin;q.y++) raw(p+q) = I;
    }
  }
  void set(Rectangle r, uint8_t I);
  void set(Mask &m, uint8_t I);
  void set(Contour& c, uint8_t I) {
    Point p;
    Rectangle b = getBounds();
    c.start();
    while (c.advance()) { if (b.contains(c.i())) set( c.i() , I ); }
  }
  void set(Point p1, Point p2, int width, uint8_t I); // Sets a line of specified radius (0=single pixel) from p1 up to but not including p2
  void set(Contour &c, int width, uint8_t I);         // Draws a contour as lines--if contour is not made of adjacent pixels, it's still connected
  
  // Operations to shift and scale image by a constant
  inline Image8& operator=(uint8_t I) { for (int i=0;i<size.x*size.y;i++) pixels[i]=I; return *this; }
  inline Image8& operator+=(uint8_t I) { for (int i=0;i<size.x*size.y;i++) pixels[i]+=I; return *this; }
  inline Image8& operator+=(Point p) { bounds += p; return *this; }
  inline Image8& operator-=(uint8_t I) { for (int i=0;i<size.x*size.y;i++) pixels[i]-=I; return *this; }
  inline Image8& operator-=(Point p) { bounds -= p; return *this; }
  
  // Copy bits of one image to another in various ways
  void mimic(const Image8& source , Rectangle my_region , Rectangle source_region , ScaleType method = Subsample);
  void mimic(const Image8& source , ScaleType method = Subsample ) { mimic(source,getBounds(),source.getBounds(),method); }
  void copy(Point where,const Image8& source,Point size);
  void copy(const Image8& source,Mask& m);
  void copy(const Image8& source) { copy( source.bounds.near , source , source.size ); }
  
  // Similar, except use 16 bit image as a source
  void mimic16(const Image& source , Rectangle my_region , Rectangle source_region , ScaleType method = Subsample);
  void mimic16(const Image& source , ScaleType method = Subsample ) { mimic16(source,getBounds(),source.getBounds(),method); }
  void copy16(Point where,const Image& source,Point size,bool fix_depth=false);
  void copy16(const Image& source,Mask& m,bool fix_depth=false);
  void copy16(const Image& source,bool fix_depth=false) { copy16( source.bounds.near , source , source.size , fix_depth ); }
      
  // Projections
  void maxOverX(Rectangle target, uint8_t *output);
  void maxOverY(Rectangle target, uint8_t *output);
  void minOverX(Rectangle target, uint8_t *output);
  void minOverY(Rectangle target, uint8_t *output);
  void meanOverX(Rectangle target, float* output);
  void meanOverY(Rectangle target, float* output);
  void deviationOverX(Rectangle target, const float* means, float* output);
  void deviationOverY(Rectangle target, const float* means, float* output);

  // Output and testing
  int makeTiffHeader(unsigned char *buffer);
  int writeTiff(FILE *f);
  int writeTiff(const char *fname);
  void println() const;
};


/****************************************************************
                         Image Alignment
****************************************************************/
class ShiftWeight {
public:
  float shift;
  float weight;

  ShiftWeight() : shift(0), weight(0) {}
  ShiftWeight(float s, float w) : shift(s), weight(w) {}
};

class Feature1D {
public:
  float position;
  float strength;
  float diameter;
  int sign;

  Feature1D() : position(0), strength(0), diameter(0), sign(0) {}
  Feature1D(float p, float s, float d, int sgn): position(p), strength(s), diameter(d), sign(sgn) {}

  ShiftWeight operator-(Feature1D& that) {
    return ShiftWeight(position - that.position, (sign * that.sign) * fminf(strength, that.strength));
  }
};

#
class Profile {
public:
  enum Collapse { OverX, OverY };

  float *reference;
  float center;
  float jitter;
  int hist[64];
  Feature1D features[8];
  int n_features;
  Rectangle source;
  Collapse direction;
  float lateral_factor;
  int border;
  int n;

  Profile(const Rectangle src, Collapse dir);
  ~Profile() {
    if (reference != NULL) { 
      delete[] reference; 
      reference = NULL;
    }
  }

private:
  float find_mean(float *values, int count);
  ShiftWeight histify(float *values, int count);
  float find_center_via_hist(float *values, int count, ShiftWeight sw);
  void load_best_features(int m);
  void compute_new_values();
  float find_best_subpixel_shift(float *ref, float mean_r, float *probe, float mean_p, float scale, int n, int border);
  float find_best_shift(float *probe, float mean_p, float scale);
public
:
  /** Quality score for how useful of a profile we have for alignment */
  float quality();

  /** Adopt an already-computed profile */
  void adopt(float *profile, int length);

  /** Adopt a profile from the given image */
  void imprint(Image& frame);

  /** Adopt a profile from the given 8 bit image */
  void imprint8(Image8& frame);

  /** Find and adopt the best quality with in `search` distance plus or minus the current position.
    * The current position is updated to the best position.
    *
    * (Note: search is not exhaustive.)
    */
  void best_near(Image& frame, Point search);

  /** Find and adopt the best quality with in `search` distance plus or minus the current position.
    * The current position is updated to the best position.
    *
    * (Note: search is not exhaustive.)
    */
  void best_near8(Image8& frame, Point search);

  /** Given the profile in `probe`, find the offset between the stored profile and the probe profile. */
  float align(float *probe, float mean, float precision);

  /** Find the offset between the image in `frame` and the stored profile by computing and comparing the new profile. */
  float delta(Image& frame, float lateral);

  /** Find the offset between the 8 bit image in `frame` and the stored profile by computing and comparing the new profile. */
  float delta8(Image8& frame, float lateral);
};


/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_image_strip();
int test_mwt_image_mask();
int test_mwt_image_contour();
int test_mwt_image_image();
int test_mwt_image_image8();
int test_mwt_image();

#endif

