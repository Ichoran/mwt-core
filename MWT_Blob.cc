/* Copyright 2007-2010 Rex A. Kerr, Nicholas A. Swierczek, and HHMI Janelia
 * Copyright 2016 Rex A. Kerr and Calico Life Sciences
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <new>

#include "MWT_Blob.h"



/****************************************************************
                     Utility Class Methods
****************************************************************/

FilenameComponent::FilenameComponent(const FilenameComponent& fc) :
  identity(fc.identity),value(fc.value),width(fc.width)
{
  if (width==0 || fc.format==NULL) format = NULL;
  else
  {
    format = new char[width + 1];
    strncpy(format,fc.format,width);
    format[width] = 0;
  }
}

void FilenameComponent::showNumber(int val,int wid)
{
  value = val;
  width = wid;
  if (width<1) width=1;
  else if (width>MAX_WIDTH) width=MAX_WIDTH;
  if (format!=NULL) delete[] format;
  format = new char[8];
  sprintf(format,"%%0%dd",width);
}

void FilenameComponent::showNameDate(const char *trackerName, struct tm& date)
{
  if (format!=NULL) delete[] format; 
  // Make sure the buffer is long enough if you change the format string!
  int nameLen = (trackerName == NULL) ? 0 : (strlen(trackerName) + 1);
  format = new char[16 + nameLen];
  if (nameLen > 0) snprintf(format, nameLen, "%s_", trackerName);
  snprintf(format + nameLen, 16, "%04d%02d%02d_%02d%02d%02d",date.tm_year,date.tm_mon+1,date.tm_mday,date.tm_hour,date.tm_min,date.tm_sec);
  format[nameLen+15] = 0;
  width = strlen(format);  // Should always be nameLen+15, but just in case format string gets changed and counted wrong....
  identity=text;
}

void FilenameComponent::showText(const char *t)
{ 
  int n = strlen(t);
  if (format!=NULL) delete[] format;
  format=new char[n+1];
  memcpy(format,t,n);
  format[n] = 0;
  identity=text;
  width=n;
}

void FilenameComponent::numberToText()
{
  char *temp = new char[ MAX_WIDTH ];
  toString(temp,MAX_WIDTH);
  delete[] format;
  format = temp;
  width = strlen(format);
  identity = text;
}

int FilenameComponent::toString(char* s,int n)
{
  if (identity==empty) { s[0]=0; return 0; }
  else if (identity==text)
  {
    int i = strlen(format);
    if (i+1>n) i=n-1;
    memcpy(s,format,i);
    s[i] = 0;
    return i;
  }
  else return snprintf(s,n,format,value);
}

// Shorten the list of components into the most compact form possible by merging text fields
void FilenameComponent::compact(ManagedList<FilenameComponent>& mfc)
{
  Listable<FilenameComponent>* lfc;
  int n;
  char *c;
  mfc.start();
  while (mfc.advance())
  {
    if (mfc.i().identity!=text) continue;  // Can't merge someone else into a text field
    lfc = mfc.current;
    while (mfc.advance(lfc))
    {
      if (lfc->data.identity!=text) break;  // Can't merge a text field into us
      
      // String handling
      n = mfc.i().width + lfc->data.width;
      c = new char[ n + 1 ];
      memcpy(c,mfc.i().format,mfc.i().width);
      memcpy(c+mfc.i().width,lfc->data.format,lfc->data.width);
      c[n] = 0;
      delete[] mfc.i().format;
      mfc.i().format = c;
      mfc.i().width = n;
      
      mfc.Backspace(lfc); // Get rid of unnecessary field
    }
  }
}

// Figures out how big of a buffer you need to be completely safe when converting a bunch of filename bits into a string
int FilenameComponent::safeLength(ManagedList<FilenameComponent>& fc)
{
  int i=1;
  fc.start();
  while (fc.advance())
  {
    switch (fc.i().identity)
    {
      case empty: break;
      case text: i += fc.i().width; break;  // Easy--strings are as long as they are long.
      default:
        int j = (int)floor(log10(fabs((double)fc.i().value))+1);  // Number of digits in a numeric field
        if (fc.i().value<0) j++;                          // Need room for a minus sign
        if (j < fc.i().width) j = fc.i().width;           // We asked for padding, so give room for it
        i += j;
        break;
    }
  }
  return i;
}

void FilenameComponent::toString(ManagedList<FilenameComponent>& fc,char *s,int n)
{
  char *s_now = s;
  int n_left = n;
  int i;
  
  fc.start();
  while (fc.advance())
  {
    if (n_left <= 1) break;
    i = fc.i().toString(s_now,n_left);
    if (i>=0) { s_now += i; n_left -= i; } else break;
  }
  s[n-1] = 0;
}

void FilenameComponent::loadValues(ManagedList<FilenameComponent>& fc,int frame_id,int dancer_id,int perf_id)
{
  fc.start();
  while (fc.advance())
  {
    switch (fc.i().identity)
    {
      case frame_n: fc.i().value = frame_id; break;
      case dancer_n: fc.i().value = dancer_id; break;
      case perf_n: fc.i().value = perf_id; break;
      default: break;
    }
  }
}


// Loads an existing contour into a new packed contour
PackedContour::PackedContour(Contour& c) : x_start(0),y_start(0),size(0),bits(NULL)
{
  Point old;
  if (c.size()==0) return;
  
  c.start();
  if (c.advance())
  {
    old = c.i();
    x_start = old.x;
    y_start = old.y;
    size = 1;
    while (c.advance())
    {
      size += (old.x > c.i().x) ? old.x - c.i().x : c.i().x - old.x;
      size += (old.y > c.i().y) ? old.y - c.i().y : c.i().y - old.y;
      old = c.i();
    }
  }
  else return;
  
  int n = (size+3)/4 + 1;
  bits = new unsigned char[n];
  unsigned int accumulator = 0;
  int cycle_index = 0;
  int bits_index = 0;
  c.start();
  c.advance();
  old = c.i();
  while (c.advance())
  {
    while (old != c.i())
    {
      accumulator <<= 2;
      if (old.x > c.i().x) { old.x--; }
      else if (old.x < c.i().x) { old.x++; accumulator |= 0x1; }
      else if (old.y > c.i().y) { old.y--; accumulator |= 0x2; }
      else if (old.y < c.i().y) { old.y++; accumulator |= 0x3; }
      cycle_index++;
      if (cycle_index>=4)
      {
        if (bits_index >= n-1)
        {
          printf("Yikes!\n");
        }
        bits[bits_index++] = (unsigned char)accumulator;
        accumulator = 0;
        cycle_index = 0;
      }
    }
  }
  while (cycle_index != 0)  // Pack last byte with dummy bits
  {
    accumulator <<= 2;
    if ((cycle_index%1) == 1) accumulator |= 0x1;
    cycle_index++;
    if (cycle_index>=4)
    {
      if (bits_index >= n-1)
      {
        printf("Yikes!\n");
      }
      bits[bits_index++] = (unsigned char)accumulator;
      accumulator = 0;
      cycle_index = 0;
    }
  }
  while (bits_index < n) bits[bits_index++] = 0;
}

// Converts a packed contour into a regular one
Contour* PackedContour::unpack(Contour* c)
{
  if (c==NULL) c = new Contour( size+1 );
  Point p(x_start,y_start);
  c->Append( p );
  int bits_index = 0;
  int decumulator = 0;
  int dec_index = 0;
  int bitpair;
  for (int i=1 ; i<size ; i++)  // Start at 1 because (x_start,y_start) was 0
  {
    if (dec_index==0)
    {
      decumulator = bits[bits_index++];
      dec_index = 4;
    }
    dec_index--;
    bitpair = ( (decumulator >> (2*dec_index)) & 0x3 );
    if (bitpair&0x2)
    {
      if (bitpair&0x1) p.y++;
      else p.y--;
    }
    else
    {
      if (bitpair&0x1) p.x++;
      else p.x--;
    }
    c->Append(p);
  }
  return c;
}

// Writes to file a string containing the contour packed 6 bits at a time into ASCII, starting at '0'.
void PackedContour::printData(FILE* f, bool add_newline)
{
  int N = (size+2)/3;
  char output[ N + 1 ];
  int out_index = 0;
  int bits_index = 0;
  int accumulator = 0;
  int decumulator = 0;
  int acc_index = 0;
  int dec_index = 0;
  while (out_index < N)
  {
    if (dec_index==0)
    {
      decumulator = bits[bits_index++];
      dec_index = 4;
    }
    dec_index--;
    accumulator <<= 2;
    accumulator |= ( (decumulator >> (2*dec_index)) & 0x3 );
    acc_index++;
    if (acc_index>=3)
    {
      output[out_index++] = (char)('0' + accumulator);
      accumulator = 0;
      acc_index = 0;
    }
  }
  while (acc_index!=0)
  {
    accumulator <<= 2;
    if ((acc_index%1)==1) accumulator |= 0x1;
    acc_index++;
    if (acc_index>=3)
    {
      output[out_index++] = (char)('0' + accumulator);
      accumulator = 0;
      acc_index = 0;
    }
  }
  while (out_index<=N) output[out_index++] = 0;
  
  if (add_newline) fprintf(f,"%s\n",output);
  else fprintf(f,"%s",output);
}


/****************************************************************
                          Blob Methods
****************************************************************/

// Destructor refers to class Dancer, so keep code out of header
Blob::~Blob()
{
  if (im!=NULL) { delete im; im = NULL; }
  if (dancer!=NULL) tossStats();
  dancer = NULL;
}


// Get rid of any stats we might have
void Blob::tossStats()
{
  if (stats!=NULL)
  {
    dancer->floodstore.destroy(stats);
    stats = NULL;
  }
  if (skeleton!=NULL)
  {
    dancer->contourstore.destroy(skeleton);
    skeleton = NULL;
  }
  if (outline!=NULL)
  {
    dancer->contourstore.destroy(outline);
    outline = NULL;
  }
  if (packedline!=NULL)
  {
    dancer->packedstore.destroy(packedline);
    packedline = NULL;
  }
}


Rectangle Blob::findNextROI(Image* bg,int border) const
{
  if (stats==NULL) return Rectangle(0,0,0,0);
  Rectangle r = stats->data.stencil.bounds;
  r.expand(border);
  if (bg!=NULL) r *= bg->getBounds();
  return r;
}


// Clip out an expected blob given the last known location and a foreground image
// Background image is optional; can be subtracted.
void Blob::clip(const Blob& last,Image* fg,Image* bg,int border)
{
  if (im!=NULL) { delete im; im=NULL; }
  Rectangle r = last.findNextROI(bg,border);
  r *= fg->getBounds();
  if (r.area()<1) {
		return;
	}
  
  if (bg)
  {
    im = new Image(r,fg->divide_bg);
    im->depth = fg->depth+1;
    im->diffCopy(r.near,*fg,r.size(),*bg);
  }
  else im = new Image(*fg , r, fg->divide_bg);
}


// Clip out an expected blob given the last known location and a foreground image
// Background image is optional; can be subtracted.
void Blob::clip8(const Blob& last, Image8* fg, Image* bg, int border)
{
  if (im!=NULL) { delete im; im=NULL; }
  Rectangle r = last.findNextROI(bg,border);
  r *= fg->getBounds();
  if (r.area()<1) {
    return;
  }
  
  im = new Image(r, fg->divide_bg);
  if (bg)
  {
    im->depth = 9;
    im->diffCopy8(r.near, *fg, r.size(), *bg);
  }
  else im->copy8(r.near, *fg, r.size());
}


// Find blob(s) inside a clipped-out image and store potential matches inside dancer->candidates.
// This also uses the size threshold criteria of this blob's dancer and overlap with
// the previous blob.
int Blob::find(Mask* old_mask,Mask* exclusion_mask)
{
  if (im==NULL || dancer==NULL) {
		return 0;
	}
  
  ManagedList<FloodData>& found = dancer->candidates;
  found.flush();
  im->floodRect( dancer->fill_I , &(dancer->ssstore) , &(dancer->lsstore) , found , im->getBounds() );
  if (found.size==0) {
		return 0;
	}
  
  found.start();
  while (found.advance())
  {
		
    if (dancer->fill_size.keep.excludes(found.i().stencil.pixel_count)) {
			found.Backspace();  // Throw away mis-sized ones
		}
			
  }
  if (found.size==0) {
		return 0;
	}
  
  if (old_mask!=NULL)
  {
    found.start();
    while (found.advance())
    {
      if (!old_mask->overlaps(found.i().stencil,Point(0,0))) found.Backspace();  // Throw away nonoverlapping ones
      else if (exclusion_mask!=NULL && exclusion_mask->overlaps( found.i().stencil )) found.Backspace();  // Throw away ones that hit the border
    }
  }
  return found.size;
}


// Decide that the current candidate blob is us, and adopt its mask
void Blob::adoptCandidate()
{
  if (dancer==NULL || dancer->candidates.current==NULL) return;
  tossStats();
  
  // We can only do the following because the ManagedList uses the storage that our dancer provides!
  // Otherwise we'd have to make a copy of the data and then use Backspace() instead.
  stats = dancer->candidates.removeL();
  pixel_count = mask().pixel_count;
}


// Clean up image & outline if they're there
void Blob::flush()
{
  if (im!=NULL)
  {
    delete im;
    im = NULL;
  }
  if (outline!=NULL) packOutline();
}


// Compute any additional stats we need once we know it's worth it
void Blob::extraStats(Range fill_I , FPoint previous_direction)
{
  // Find major and minor axes
  if (im!=NULL) data().principalAxes( *im , fill_I );
  
  if ( data().major * previous_direction < 0 )
  {
    data().major = -data().major;
    data().minor = -data().minor;
  }
  
  if (mask().pixel_count==0 || dancer==NULL || !(dancer->find_skel || dancer->find_outline)) return;
  
  // Find contour around mask
  Contour* c;
  if (outline==NULL)
  {
    outline = dancer->contourstore.create(false);
    c = new( &(outline->data) ) Contour( &(dancer->lpstore) );
  }
  else
  {
    c = &(outline->data);
    c->flush();
  }
  c->findContour(mask());
  c->fillContour();
  
  if (dancer->find_outline)
  {
    if (dancer->find_skel)
    {
      if (skeleton==NULL)
      {
        skeleton = dancer->contourstore.create(false);
        new( &(skeleton->data) ) Contour( &(dancer->lpstore) );
      }
      else skeleton->data.flush();
      skeleton->data.imitate(outline->data);
    }
    
    float max_error = sqrt( mask().pixel_count )*0.04;
    if (max_error<1.0) max_error = 1.0;
    // Making an approximate contour is too expensive--just keep the whole thing but write it out in compact form.
    // (Also we switch to compact form when images are dumped.)
    // (Also, it may be broken now!  I've disabled unit tests!!!)
    // outline->data.approximateContour( max_error );
  }
  else // Swap outline into skeleton
  {
    Listable<Contour>* temp = outline;
    outline = skeleton;
    skeleton = temp;
  }
  
  if (dancer->find_skel)
  {
    skeleton->data.approximateSpine( data().centroid , data().major , SKELETON_SIZE );
  }
}

void Blob::packOutline()
{
  if (outline==NULL || dancer==NULL) return;
  packedline = dancer->packedstore.create(false);
  new( &(packedline->data) ) PackedContour( outline->data );
  dancer->contourstore.destroy(outline);
  outline = NULL;
}

  
// Output stats to one line of a text file, returns false if it fails
bool Blob::print(FILE* f,const char* prefix)
{
  int i;
  char tail_character;
  FloodData* d;
  
  if (stats==NULL) {
    i = -1;
  }
  else
  {
    tail_character = ((skeleton!=NULL && spine().size()==SKELETON_SIZE) || outline!=NULL || packedline!=NULL) ? ' ' : '\n';  // Are we going to print out shape(s)?
    if (prefix==NULL)
    {
      d = &(stats->data);
      i = fprintf(f,"%d %.3f  %.3f %.3f  %d  %.3f %.3f  %.3f  %.1f %.1f%c",
        frame,time,d->centroid.x,d->centroid.y,pixel_count,
        d->major.x,d->major.y,d->minor.length(),d->long_axis,d->short_axis,tail_character);
    }
    else
    {
      d = &(stats->data);
      i = fprintf(f,"%s %d %.3f  %.3f %.3f  %d  %.3f %.3f  %.3f  %.1f% .1f%c",
        prefix,frame,time,d->centroid.x,d->centroid.y,pixel_count,
        d->major.x,d->major.y,d->minor.length(),d->long_axis,d->short_axis,tail_character);
    }
    if (tail_character==' ' && skeleton!=NULL && spine().size()==SKELETON_SIZE) // Print skeleton
    {
      fprintf(f,"%%"); // Separator character for skeleton
      spine().start();
      while (spine().advance())
      {
        fprintf(f," %d %d",spine().i().x,spine().i().y);
      }
      if (outline!=NULL || packedline!=NULL) fprintf(f," ");
      else
      {
        fprintf(f,"\n");
        tail_character = '\n';
      }
    }
    if (tail_character==' ' && (outline!=NULL || packedline!=NULL)) // Print outline
    {
      if (packedline==NULL) packOutline();
      fprintf(f,"%%%% %d %d %d ",packedline->data.x_start,packedline->data.y_start,packedline->data.size);
      packedline->data.printData(f,true);
    }
  }
  return (i>=0);
}


/****************************************************************
                         Dancer Methods
****************************************************************/

Dancer::~Dancer()
{
  // First, figure out whether we need to write any files, write them if so, and leave error messages if something went wrong
    FILE *f = NULL;
    if( performance != NULL ) 
    {
      if( frames.lo() > 0 )
      {
        if( performance->combine_blobs ) {
          f = performance->getFileHandle();
					
          if( f ) {
            performance->addOutputFate(frames.hi(),ID,ftell(f));  
            fprintf(f,"%% %d\n",ID);
          }
        }
        else 
        {
          if (data_fname!=NULL )
          {
            f = fopen(data_fname,"w");
            if (!f) complain("Could not open file",data_fname);
          }          
        }
        if( f != NULL ) 
        {
          if (!f) complain("Could not open file",data_fname);
          else if (!printData(f,NULL)) complain("Incomplete write to file",data_fname);
          if( !performance->combine_blobs ) 
          {
            fclose(f);
          }
        }
      }
    }
  // Then get rid of everything that doesn't have its own destructor
  tidyHistoryLeaving(0);
  if (data_fname!=NULL) {
    delete[] data_fname;
  }
  if (img_namer!=NULL) {
    delete img_namer;
  }
  movie.flush();
}


// Set a newly found object to be the start of this dance
// FloodData must NOT use any of our storage!
void Dancer::setFirst(Image *im,FloodData* fd,int frame,double time)
{
  // Set up basic data
  frames = Range(frame,frame);
  times = FRange(time,time);
	
  // We'll need stats
  if (im) fd->principalAxes( (*im) , fill_I.keep );
  
  // Make a local copy of flood data
  candidates.flush();
  FloodData *new_fd = new( candidates.Append() ) FloodData(&lsstore);
  new_fd->duplicate(*fd);
  candidates.start(); candidates.advance();
  
  // Make the first blob
  resetMovie();
  Blob* b = new( movie.Append() ) Blob(frame,time);
  b->dancer = this;
  b->adoptCandidate();
  if (im) b->clip( *b , im , NULL , border );  // Clip out about self--works as long as we have stats
  
  // Make sure we know that we're all done
  validated = true;
}


// Create a blob to be a newly found object--need to fill in FloodData and Blob
Blob* Dancer::makeFirst(int frame,double time)
{
  frames = Range(frame,frame);
  times = FRange(time,time);
  candidates.flush();
  new( candidates.Append() ) FloodData(&lsstore);
  candidates.start(); candidates.advance();
  resetMovie();
  Blob *new_b = new( movie.Append() ) Blob(frame,time);
  new_b->dancer = this;
  new_b->adoptCandidate();
  validated = false;
  return new_b;
}


// Get ready to look for a new blob
void Dancer::readyAnother(Image* fg,Image* bg,int frame,double time)
{
  Blob& old_blob = movie.t();
  
  new( movie.Append() ) Blob(frame,time);
  
  movie.t().dancer = this;
  movie.t().clip( old_blob , fg , bg , border );
	
  frames.hi() = frame;
  times.hi() = time;
  validated = false;
}


// Get ready to look for a new blob
void Dancer::readyAnother8(Image8* fg, Image* bg, int frame, double time)
{
  Blob& old_blob = movie.t();
  
  new( movie.Append() ) Blob(frame,time);
  
  movie.t().dancer = this;
  movie.t().clip8( old_blob , fg , bg , border );
  
  frames.hi() = frame;
  times.hi() = time;
  validated = false;
}


// Actually look for a new blob
bool Dancer::findAnother(bool best_guess, Mask* exclusion_mask, FPoint jitter)
{
  if (validated) {
		return true;  // If we've already found one, we're done.
	}
  
  int i;
	if (movie.size==1) {
		i = movie.t().find(NULL,exclusion_mask);
	}
	
  else
  {
    movie.end();
    movie.retreat();
    i = movie.t().find(&movie.i().mask(),exclusion_mask);
  }
  
  if (i==0) {
		return false;
	}
  
  if (i==1)
  {
    candidates.start();
    candidates.advance();
    movie.t().adoptCandidate();
  }
  else
  {
    FPoint ratio;
    
    // Make "current" be the previous blob
    movie.end();
    movie.retreat();
    
    // Now compare each flooded object to the previous blob
    candidates.start();
    while (candidates.advance())
    {
      ratio = movie.i().mask().overlapRatio(candidates.i().stencil);
      candidates.i().score = sqrt(ratio.x*ratio.y);             // Geometric mean
      if (ratio.y < MIN_OVERLAP_RATIO) candidates.Backspace();  // New object doesn't overlap existing object enough
    }
    
    // Sort the list--highest overlap (end of list) is best match
    if (candidates.size>1) candidates.mergeSort();
    
    if (!best_guess) return false;
    else
    {
      candidates.end();
      movie.t().adoptCandidate();
    }
  }
  return true;
}


// Mark the current candidate as verified to be okay
// Also, compute any stats that we need to compute and do any needed cleanup/output
void Dancer::validate()
{
  if (validated) return;
  
  validated = true;
  
  // Find extra stats--need previous frame to keep direction vectors aligned
  if (movie.size<2) movie.t().extraStats( fill_I.keep , FPoint(0.0,0.0) );
  else
  {
    movie.end(); movie.retreat();
    if (movie.i().stats==NULL) movie.t().extraStats( fill_I.keep , FPoint(0.0,0.0) );
    else movie.t().extraStats( fill_I.keep , movie.i().data().major );
  }
  
  accumulated_length.add( movie.t().stats->data.long_axis );
  accumulated_width.add( movie.t().stats->data.short_axis );
  accumulated_aspect.add( movie.t().stats->data.short_axis / movie.t().stats->data.long_axis );
  
  if (frames.lo()==0)
  {
    movie.Behead();
    frames.lo() = movie.h().frame;
  }
  
  tidyHistory();
}


// Toss the current attempt, leaving a blank frame in the movie
void Dancer::invalidate()
{
  movie.Truncate();
  validated = true;  // Last one was validated
}


// Something to conveniently drop error messages into our performance
void Dancer::complain(const char* topic,const char* message)
{
  int tL = strlen(topic);
  int tM = strlen(message);
  char *tempstr = new char[ tL + tM + 2 ];
  strncpy(tempstr,topic,tL);
  tempstr[tL] = ' ';
  strncpy(tempstr+(tL+1),message,tM);
  tempstr[tL+tM+1] = 0;
  new( performance->errors.Append() ) IOErrorHandler(ID,tempstr,false);
  tempstr = NULL;
}


// Make sure we don't give any output--throw away output filenames
void Dancer::tossFilenames()
{
  if (data_fname!=NULL) delete[] data_fname;
  if (img_namer!=NULL) delete img_namer;
  data_fname = NULL;
  img_namer = NULL;
}


// Clean up our history.
void Dancer::tidyHistoryLeaving(int n_keep)
{
  if (movie.size==0 || movie.size<=n_keep+n_cleared) return;
  
  if (n_cleared==0)  // First time we're tidying
  {
    // Get started, point at the first frame in the movie
    movie.start(clear_til);
    if (!movie.advance(clear_til)) return;
    
    
    // Create directory--may result in empty directories if you set n_long_enough > n_keep
    // Why would you ever do that anyway?  You'll throw away images before you realize you wanted them!
    if (img_namer!=NULL)
    {
      FilenameComponent::loadValues(*img_namer , clear_til->data.frame , ID , performance->ID);
      int L = FilenameComponent::safeLength( *img_namer );
      char path_fname[ L ];
      FilenameComponent::toString(*img_namer , path_fname , L);
      char *c = strrchr(path_fname,'/');
      if (c==NULL) complain("Could not find path in",path_fname);
      else
      {
        *c = 0;
#ifdef WINDOWS
        if (mkdir(path_fname)!=0) complain("Could not create path for",path_fname);
#else
        if (mkdir(path_fname,01777)!=0) complain("Could not create path for",path_fname);
#endif
      }
    }
  }
  
  do
  {
    if (clear_til->data.im != NULL)  // We have an image
    {
      if (img_namer!=NULL && n_long_enough <= frames.length())  // Save it if we were supposed to
      {
        FilenameComponent::loadValues(*img_namer , clear_til->data.frame , ID , performance->ID);
        int L = FilenameComponent::safeLength( *img_namer );
        char extended_fname[ L ];
        FilenameComponent::toString(*img_namer , extended_fname , L);
        
        int i = clear_til->data.im->writeTiff(extended_fname);
        
        if (i!=0) complain("Could not write image file",extended_fname);
      }
      delete clear_til->data.im;
      clear_til->data.im = NULL;
      if (clear_til->data.outline!=NULL && n_keep!=0) clear_til->data.packOutline();
    }
    if (clear_til->data.packedline != NULL && n_keep == 0)
    {
      packedstore.destroy(clear_til->data.packedline);
      clear_til->data.packedline = NULL;
    }
    if (clear_til->data.outline != NULL && n_keep == 0)
    {
      contourstore.destroy(clear_til->data.outline);
      clear_til->data.outline = NULL;
    }
    if (clear_til->data.skeleton != NULL && n_keep == 0)
    {
       contourstore.destroy(clear_til->data.skeleton);
       clear_til->data.skeleton = NULL;
    }
    if (clear_til->data.stats != NULL)
    {
      if (clear_til->data.mask().pixel_count!=0) clear_til->data.mask().flush();
    }
    n_cleared++;
  } while (movie.advance(clear_til) && movie.size>n_keep+n_cleared);
}


// Print out the coordinates and other statistics for this dance in text
bool Dancer::printData(FILE* f,const char* prefix)
{
  if (movie.size <= 0) return true;
  
  bool ok = true;
  Listable<Blob>* lb;
  movie.start(lb);
  while (movie.advance(lb))
  {
    ok = lb->data.print(f,prefix);
    if (!ok) break;
  }
  return ok;
}


// Report velocity since the last time we were asked
float Dancer::recentSpeed(float min_interval)
{
  if (last_speed_time+min_interval >= times.hi() || movie.size<2) return last_speed;
  
  // Find endpoint of last speed measurement
  movie.end();
  while (movie.retreat())
  {
    if (movie.i().time <= last_speed_time) break;
  }
  if (movie.current==NULL) { movie.start(); movie.advance(); }
  
  if (movie.i().stats==NULL || movie.t().stats==NULL || movie.i().time >= movie.t().time)
  {
    last_speed_time = times.hi();
    return last_speed;
  }
  
  last_speed = (movie.t().data().centroid - movie.i().data().centroid).length();
  last_speed /= (movie.t().time - movie.i().time);
  last_angularspeed = (movie.t().data().major.cross( movie.i().data().major )).length();
  last_angularspeed /= movie.t().data().major.length()*movie.i().data().major.length();
  last_angularspeed = fabs( asin(last_angularspeed) );
  last_angularspeed /= (movie.t().time - movie.i().time);
  
  last_speed_time = movie.t().time;
  return last_speed;
}



/****************************************************************
                       Performance Methods
****************************************************************/

// Clean up everything we may have allocated
Performance::~Performance()
{
  if( combine_blobs ) {
    if( output_file != NULL ) {
      fclose(output_file);
      output_file = NULL;
    }
  }
  
  // Need to flush things that might require filenames or error handling or whatever
  candidates.flush();
  sitters.flush();
  dancers.flush();
  
  // Now toss everything we may explicitly have allocated
  if (date!=NULL) { delete date; date=NULL; }
  if (foreground!=NULL) { delete foreground; foreground=NULL; }
  if (background!=NULL) { delete background; background=NULL; }
  if (full_area!=NULL) { delete full_area; full_area=NULL; }
  if (danger_zone!=NULL) { delete danger_zone; danger_zone=NULL; }
  if (band!=NULL) { delete band; band=NULL; }
  if (band_area!=NULL) { delete band_area; band_area=NULL; }
  if (dance_fname!=NULL) { delete dance_fname; dance_fname=NULL; }
	if (blobs_fname!=NULL) { delete blobs_fname; blobs_fname=NULL; }
  if (sit_fname!=NULL) { delete sit_fname; sit_fname=NULL; }
  if (img_fname!=NULL) { delete img_fname; img_fname=NULL; }
  if (base_directory!=NULL) { delete[] base_directory; base_directory=NULL; }
  if (prefix_name!=NULL) { delete[] prefix_name; prefix_name=NULL; }
  if (path_date_connector!=NULL) { free(path_date_connector); }
}

FILE* Performance::getFileHandle()
{
  if( blobs_fname == NULL ) {
    output_file = NULL;
  }
  else 
  {
    if( output_count < MAX_BLOBS_PER_FILE ) {
      output_count++;
    }
    else {
      if( output_file != NULL ) 
      {
        fclose(output_file);
        output_file = NULL;        
      }
      output_count = 0;
    }    
    if( output_file == NULL )			
    {
      char *fname = NULL;
      int n;
      FilenameComponent::loadValues(*blobs_fname,current_frame,last_blobs_fid++,ID);
      n = FilenameComponent::safeLength(*blobs_fname);
      fname = new char[ n ];
      FilenameComponent::toString(*blobs_fname,fname,n);
      output_file = fopen(fname,"w");
			delete[] fname;
    }
  }
  return output_file;
}

void Performance::addOutputFate( int frame, int idb, long byte_offset ) 
{	
  output_fates.Append( BlobOutputFate(frame, idb, last_blobs_fid-1,byte_offset));
}


// Create the region of interest (if none is created, it will just be the whole field of view).
void Performance::setROI(const Rectangle& r)
{
  if( foreground != NULL ) { delete foreground; foreground = NULL; }
  if( background != NULL ) { delete background; background = NULL; }
  if (full_area!=NULL) delete full_area;
  if (danger_zone!=NULL) { delete danger_zone; danger_zone=NULL; }
  full_area = new Mask(r , &lsstore);
  full_area->findBounds();
}
void Performance::setROI(const Ellipse& e)
{
  if( foreground != NULL ) { delete foreground; foreground = NULL; }
  if( background != NULL ) { delete background; background = NULL; }
  if (full_area!=NULL) delete full_area;
  if (danger_zone!=NULL) { delete danger_zone; danger_zone=NULL; }
  full_area = new Mask(e, &lsstore);
  full_area->findBounds();
}
void Performance::setROI(Mask &m)
{
  if( foreground != NULL ) { delete foreground; foreground = NULL; }
  if( background != NULL ) { delete background; background = NULL; }
  if (full_area!=NULL) delete full_area;
  if (danger_zone!=NULL) { delete danger_zone; danger_zone=NULL; }
  full_area = new Mask(&lsstore);
  full_area->imitate(m);
  full_area->findBounds();
}

// Call this to set up a mask that detects when you're within the image but running too close to the edge of the ROI
void Performance::setNotInROI(const Image *im)
{
  if (im==NULL) return;
  if (danger_zone==NULL || !danger_zone->bounds.contains( im->getBounds() ))
  {
    if (danger_zone!=NULL) danger_zone->flush();
    else danger_zone = new Mask(&lsstore);
    danger_zone->imitate(*full_area);
    Rectangle r = im->getBounds();
    r.include( full_area->bounds );
    r.expand(2);
    danger_zone->invert(r);
    danger_zone->dilate(2);
    danger_zone->findBounds();
  }
}


// Scan an image, finding all blobs in it.  Return how many we've found.
// The image will usually be our foreground or band image.
int Performance::scanImage(Image *im)
{
  // Flood fill all possible candidates
  if (candidates.size>0) candidates.flush();
  im->floodMask( fill_I , &ssstore , &lsstore , candidates , (*full_area) );  // Be sure to recover FloodData items!
  if (candidates.size==0) return 0;

  // Throw away everything that's not the right size or overlaps the border
  candidates.start();
  while (candidates.advance())
  {
    fill_size.start.excludes( candidates.i().stencil.pixel_count);
    danger_zone->overlaps( candidates.i().stencil );
    if (fill_size.start.excludes( candidates.i().stencil.pixel_count )) candidates.Backspace();
    else if (danger_zone->overlaps( candidates.i().stencil )) candidates.Backspace();  
  }

  return candidates.size;
}


// Turn any scanned objects into dancers
void Performance::adoptScan(Image* im)
{
  if (candidates.size==0) return;
  
  Dancer* d;
  
  candidates.start();
  while (candidates.advance())
  {
    d = new ( dancers.Append() ) Dancer(this,dance_buf_size,blob_is_dark,fill_I,fill_size,border,n_keep_full,n_long_enough,
                                        next_dancer_ID++,current_frame,current_time);
    d->setFirst(im,&candidates.i(),current_frame,current_time);
    d->find_skel = find_dancer_skel;
    d->find_outline = find_dancer_edge;
  }
}
 

// Make an initial scan for blobs
int Performance::initialScan(Image *fg,double time)
{
  bool mask_was_set = false;
  // Start a new performance
  current_time = time;
  current_frame = 0;
  next_dancer_ID = 1;
  
  // If we found anything last time, better clean it up.
  if (dancers.size>0)  dancers.flush();
  if (full_area==NULL) full_area = new Mask( fg->getBounds() , &lsstore );
  else mask_was_set = true;
  if (background==NULL)  // First try ever--need to set up background image and exclusion zone
  {
    Rectangle r;
    if (foreground!=NULL)
    {
      delete foreground;
      foreground = NULL;
    }
    if (mask_was_set)
    {
      full_area->cropTo( fg->getBounds() );
      full_area->findBounds();
      r = full_area->bounds;
      r.expand(1);             // One-pixel safety margin around mask
      r *= fg->getBounds();    // But we still can't exceed the image bounds
    }
    else 
    {
      r = fg->getBounds();
    }
    background = new Image( (*fg) , r, fg->divide_bg );
    if (bg_depth > fg->depth) (*background) <<= (bg_depth - fg->depth);
    background->depth = bg_depth;
    foreground = new Image(fg->getBounds(), fg->divide_bg);
    foreground->depth = fg->depth+1;
    (*foreground) = 1<<fg->depth;  // Must all be gray, since we're subtracting ourselves
    return 0;
  }
  else
  {
    int n;
    
    /*
     * NOTE: This call corrects for any change that may have happened to the image size, perhaps due to a change in binning.
     *       There is a performance hit for this call, however, as it should only occur during initialization, it should be an allowable hit.
     *       A more intelligent implementation would be to move this call to the function that changes the libary's binning dimension. The problem there
     *       is that this function is based off an image which may or may not exist when the binning change function is called.
     */
    setNotInROI(foreground);
    foreground->diffAdaptCopy( (*fg) , (*full_area) , (*background) , adapt_rate );
    
    n = scanImage(foreground);
    if (n==0) return 0;
    adoptScan(foreground);
    return n;
  }
}
int Performance::initialScan8(Image8 *fg, double time)
{
  Image im(fg->getBounds(), false);  // Expensive, but initialization only
  im.depth = 8;
  im.copy8(*fg, true);
  return initialScan(&im, time);
}


// Find reference objects initially
// Removes any points that get filled--make a copy if you want to preserve point list!
// Returns the number of points successfully filled
int Performance::initialRefs(Image* fg,ManagedList<Point>& locations,double time)
{
  int i;
  int A = fg->getBounds().area();
  DualRange ref_size( Range(2,A) , Range(2,A) );
  FloodData* fd = NULL;
  Dancer* d;
  
  candidates.flush();
  if (sitters.size>0) sitters.flush();
  
  locations.start();
  while (locations.advance())
  {
    if (!fg->getBounds().contains(locations.i())) continue;
    if (fd==NULL) fd = new( candidates.Append() ) FloodData(&lsstore);
    i = fg->floodFind( locations.i() , ref_I , &ssstore , fd );
    if (!ref_size.start.excludes(i))
    {
      fd = NULL;  // If the size is OK, get a new one next time; otherwise reuse
      locations.Backspace();  // Turned point into flood, so throw it away
    }
  }
  if (fd!=NULL) candidates.Truncate();  // Last flood was bad, so throw away record
  
  if (candidates.size==0) return 0;
  
  i = 0;
  candidates.start();
  while (candidates.advance())
  {
    d = new ( sitters.Append() ) Dancer(this,dance_buf_size,true,ref_I,ref_size,border,2,2,++i,1,time);
    d->setFirst(fg , &candidates.i() , 1 , time );
  }
  
  return sitters.size;
}
int Performance::initialRefs8(Image8 *fg, ManagedList<Point>& locations, double time) {
  Image im(fg->getBounds(), false);
  im.depth = 8;
  im.copy8(*fg, true);
  return initialRefs(&im, locations, time);
}


void Performance::calculateJitterDeltaSort(float delta, float *xedge, int &nxe, float *yedge, int &nye) {
  float *eg;
  int nv;
  if (edge_profile.direction == Profile::OverY) {
    eg = xedge;
    nv = nxe;
    nxe += 1;
  }
  else {
    eg = yedge;
    nv = nye;
    nye += 1;
  }
  if (nv < 21) {
    // Insertion sort from small side
    while (nv > 0 && delta > eg[nv-1]) { eg[nv] = eg[nv-1]; nv--; }
    eg[nv] = delta;
  }
  else if (eg[0] > delta && delta < eg[20]) {
    int j;
    if (nv & 1) {
      // Half the time insertion sort from the right
      for (j = 20; j > 0 && delta > eg[j]; j--) eg[j] = eg[j-1];
    }
    else {
      // Half the time insertion sort from the left
      for (j = 0; j < 20 && delta < eg[j]; j++) eg[j] = eg[j+1];
    }
    eg[j] = delta;
  }
}

void Performance::calculateJitterGivenDeltaSorted(float *xedge, int nxe, float *yedge, int nye) {
  float jx = 0;
  float jy = 0;
  if (nxe > 0) jx = (nxe&1) ? xedge[nxe/2] : 0.5*(xedge[nxe/2] + xedge[nxe/2 - 1]);
  if (nye > 0) jy = (nye&1) ? yedge[nye/2] : 0.5*(yedge[nye/2] + yedge[nye/2 - 1]);
  jitter.set(jx, jy);
  if (fabsf(jx - ijitter.x) > 0.55) ijitter.x = (int)lrintf(jx);
  if (fabsf(jy - ijitter.y) > 0.55) ijitter.y = (int)lrintf(jy);
}

void Performance::calculateJitter(Image *fg, ManagedList<Profile> &edges) {
  if (edges.size == 0) {
    jitter.set(0, 0);
    ijitter.set(0, 0);
    return;
  }
  float edge[44]; memset(edge, 0, sizeof(float)*44);
  float *xedge = edge;
  float *yedge = edge + 22;
  int nxe = 0;
  int nye = 0;
  Listable<Profile> *lp = NULL;
  while (edges.advance(lp)) {
    edge_profile.imitate(lp->data);
    auto delta = edge_profile.delta(*fg, &(lp->data));
    if (isnan(delta)) continue;
    calculateJitterDeltaSort(delta, xedge, nxe, yedge, nye);
  }
  calculateJitterGivenDeltaSorted(xedge, nxe, yedge, nye);
}

void Performance::calculateJitter8(Image8 *fg, ManagedList<Profile> &edges) {
  if (edges.size == 0) {
    jitter.set(0, 0);
    ijitter.set(0, 0);
    return;
  }
  float edge[44]; memset(edge, 0, sizeof(float)*44);
  float *xedge = edge;
  float *yedge = edge + 22;
  int nxe = 0;
  int nye = 0;
  Listable<Profile> *lp = NULL;
  while (edges.advance(lp)) {
    edge_profile.imitate(lp->data);
    auto delta = edge_profile.delta8(*fg, &(lp->data));
    if (isnan(delta)) continue;
    calculateJitterDeltaSort(delta, xedge, nxe, yedge, nye);
  }
  calculateJitterGivenDeltaSorted(xedge, nxe, yedge, nye);
}

// Get ready to load in data from the next image (return how many items there will be)
int Performance::anticipateNext(double time)
{
  current_frame++;
  current_time = time;
  sitters.start();
  dancers.start();
  if (band_area!=NULL) band_area->flush();
  load_state = check_sitter;
  return sitters.size + dancers.size + 1;
}

// Find an item to load in the next image and report its bounds
bool Performance::findNextItemBounds(Rectangle& im_bound)
{
  switch (load_state)
  {
    case no_state:
      return false;
      break;

    case check_sitter:
      if (!sitters.advance())
      {
        load_state = check_dancer;
        return findNextItemBounds(im_bound);  // Recursive!
      }
      load_state = load_sitter;
      // Fall through!
    case load_sitter:
      im_bound = sitters.i().movie.t().findNextROI(NULL,sitters.i().border);
      break;

    case check_dancer:
      if (!dancers.advance())
      {
        load_state = check_band;
        return findNextItemBounds(im_bound);
        break;
      }
      load_state = load_dancer;
      // Fall through!
    case load_dancer:
      im_bound = dancers.i().movie.t().findNextROI(background,dancers.i().border);
      break;

    case check_band:
      // Need braces to enclose new local variables on stack
      {
        // Figure out how big the band should be
        Rectangle r = background->getBounds();

        int band_cols = (int) ceil( r.width()/(double)n_scan_bands );
        
        // Now find the position of these bands
        int i,j,x0,x1,y0,y1;
        j = n_scan_bands*2 + 1;  // Regular bands plus half-overlapping bands
        i = (current_frame*4) % j; // Go up in uneven pattern--4 is relatively prime to odd numbers, so this covers everyone
        x0 = (i-1)*band_cols/2 + r.near.x;
        x1 = x0 + band_cols - 1;
        if (x0 < r.near.x) x0 = r.near.x;
        else if (x0 > r.far.x) x0 = r.far.x;
        if (x1 > r.far.x) x1 = r.far.x;
        else if (x1 < r.near.x) x1 = r.near.x;
        
        // Create a mask to cover the area, minus where worms are
        if (band_area==NULL) band_area = new Mask(band_cols*2);
        else band_area->flush();
        full_area->start();
        while (full_area->advance())
        {
          if (full_area->i().x < x0) continue;
          if (full_area->i().x > x1) break;
          y0 = full_area->i().y0;
          y1 = full_area->i().y1;
          band_area->addStrip(full_area->i().x,y0,y1);
        }
        dancers.start();
        while(dancers.advance()) { *band_area -= dancers.i().movie.t().im->getBounds(); } 

        // Make the band the right size--ready to load image data
        if (band!=NULL)
        {
          if ( band->getBounds().width() != 1+x1-x0 || band->getBounds().height() != r.height() )
          {
            delete band;
            band = NULL;
          }
          else band->bounds = Rectangle( Point(x0,r.near.y) , Point(x1,r.far.y) );
        }
        if (band==NULL)
        {
          band = new Image( Rectangle( Point(x0,r.near.y) , Point(x1,r.far.y) ), use_division );
          band->bin = 0;
        }
      }
      
      load_state = load_band;
      
			// Can't fall through because of local variables

      im_bound = band->getBounds();
      break;

		case load_band:
			im_bound = band->getBounds();
      break;

    case all_loaded:
    default:
      return false;
      break;
  }
  return true;
}

// Load the bit of image for one item in the next image
void Performance::loadNextSingleItem(Image* fg)
{
  switch (load_state)
  {
    case no_state:
      break;
    case check_sitter:
      break;
    case check_dancer:
      break;
    case check_band:
      break;
    case load_sitter:
      sitters.i().readyAnother(fg,NULL,current_frame,current_time);
      load_state = check_sitter;
      break;
    case load_dancer: 
      dancers.i().readyAnother(fg,background,current_frame,current_time);
      load_state = check_dancer;
      break;
    case load_band:
      // We already prepared the image when we checked it, so we just need to load it
      band->depth = fg->depth+1;
      *band = 1 << (band->depth-1);  // Initialize to gray
      band->diffAdaptCopy(*fg,*band_area,*background,adapt_rate);
      load_state = all_loaded;
      break;
    case all_loaded:
    default:
      break;
  }
}

// Load the bit of image for one item in the next image
void Performance::loadNextSingleItem8(Image8* fg)
{
  switch (load_state)
  {
    case no_state:
      break;
    case check_sitter:
      break;
    case check_dancer:
      break;
    case check_band:
      break;
    case load_sitter:
      sitters.i().readyAnother8(fg, NULL, current_frame, current_time);
      load_state = check_sitter;
      break;
    case load_dancer: 
      dancers.i().readyAnother8(fg, background, current_frame, current_time);
      load_state = check_dancer;
      break;
    case load_band:
      // We already prepared the image when we checked it, so we just need to load it
      band->depth = 9;
      *band = 256;  // Initialize to gray
      band->diffAdaptCopy8(*fg, *band_area, *background, adapt_rate);
      load_state = all_loaded;
      break;
    case all_loaded:
    default:
      break;
  }
}
  

// Load data into appropriate images--use the one-item-at-a-time routines above
void Performance::readyNext(Image* fg, ManagedList<Profile> &edges, double time)
{
  calculateJitter(fg, edges);
  anticipateNext(time);
  
  fg->bounds -= ijitter;
  Rectangle single_bound;
  while (findNextItemBounds(single_bound)) loadNextSingleItem(fg);
  fg->bounds += ijitter;
}


// Load data into appropriate images--use the one-item-at-a-time routines above
void Performance::readyNext8(Image8* fg, ManagedList<Profile> &edges, double time)
{
  calculateJitter8(fg, edges);
  anticipateNext(time);
  
  fg->bounds -= ijitter;
  Rectangle single_bound;
  while (findNextItemBounds(single_bound)) loadNextSingleItem8(fg);
  fg->bounds += ijitter;
}


// Update positions of objects and find new objects; returns number of objects
int Performance::findNext()
// NOTE: Calls to ManagedList.Backspace() call the destructor contained object, which in the
//       case of Dancers, causes output to written to predefined files.
{
  bool found;
  Dancer* d;
  FPoint jit = jitter - ijitter;
  
  // First get reference objects
  sitters.start();
  while(sitters.advance())
  {
    found = sitters.i().findAnother(true, NULL, jit);
    if (found) sitters.i().validate();
    else sitters.i().invalidate();
    if (sit_fname!=NULL && current_frame<3) enableOutput(sitters.i(),true); // Only enable at beginning since we never lose these
  }
  
  // Then get moving objects
  dancers.start();
  while (dancers.advance())
  {
    //printf("I found a dancer!  %d spans %d - %d\n", dancers.i().ID, dancers.i().frames.lo(), dancers.i().frames.hi());
    found = dancers.i().findAnother(false, danger_zone, jit);
    if (!found)
    {
      if (dancers.i().candidates.size==0)  // Lost blob
      {
        if (dancers.i().frames.hi()>0)
        {
          fates.Append( BlobOriginFate(current_frame,dancers.i().ID,0) );
        }
        dancers.i().invalidate();  
        dancers.Backspace();

      }
      else  // Blob split, make new dancers out of them
      {
        dancers.i().candidates.start();
        while (dancers.i().candidates.advance())
        {
          d = new( dancers.Append() ) Dancer(this,dance_buf_size,blob_is_dark,fill_I,fill_size,border,n_keep_full,n_long_enough,
                                             next_dancer_ID++,current_frame,current_time);
          d->setFirst(dancers.i().movie.t().im,&dancers.i().candidates.i(),current_frame,current_time);
          d->find_skel = find_dancer_skel;
          d->find_outline = find_dancer_edge;
          fates.Append( BlobOriginFate(current_frame,dancers.i().ID,d->ID) );
        }

        dancers.i().invalidate();
        dancers.Backspace();       
      }
    }
  }
  // Then get new objects from strip
  if (candidates.size>0) candidates.flush();
  band->floodMask( fill_I , &ssstore , &lsstore , candidates , (*band_area) );
  candidates.start();
  while (candidates.advance())
  {
    if (!fill_size.start.excludes( candidates.i().stencil.pixel_count ))
    {
      d = new( dancers.Append() ) Dancer(this,dance_buf_size,blob_is_dark,fill_I,fill_size,border,n_keep_full,n_long_enough,
                                         next_dancer_ID++,0,current_time);
      d->setFirst(band,&candidates.i(),0,current_time);
      d->find_skel = find_dancer_skel;
      d->find_outline = find_dancer_edge;
    }
  }

  // Then handle any collisions
  Listable<Dancer>* ld;
  FloodData *fd,*fd2;
  Blob *b;
  Rectangle r;
  Image *im;
  if (dancers.size>1) dancers.mergeSort();
  dancers.start();
  while (dancers.advance())
  {
    fd = &dancers.i().movie.t().stats->data;
    ld = dancers.current;
    while (dancers.advance(ld))
    {
      // See if the two blobs overlap
      fd2 = &ld->data.movie.t().stats->data;
      if (fd->stencil.bounds.far.x < fd2->stencil.bounds.near.x) break;
      if (!fd->stencil.bounds.overlaps(fd2->stencil.bounds)) continue;
      if (!fd->stencil.overlaps(fd2->stencil,Point(0,0))) continue;
      
      // They do overlap--create merged dancer
      d = new( dancers.Tuck() ) Dancer(this,dance_buf_size,blob_is_dark,fill_I,fill_size,border,n_keep_full,n_long_enough,
                                         next_dancer_ID++,current_frame,current_time);
      d->find_skel = find_dancer_skel;     // Skeleton will probably be garbled because this is two dancers, but we'll do it anyway
      d->find_outline = find_dancer_edge;  // This should be OK even with two of them
					 
      // New dancer needs an image
      r = fd->stencil.bounds;
      r += fd2->stencil.bounds;
      im = new Image( r.size(), use_division );
      im->bounds += r.near;
      im->bin = ld->data.movie.t().im->bin;
      im->depth =  ld->data.movie.t().im->depth;
      *im = 1<<(im->depth-1);
      im->copy( fd->stencil.bounds.near , *dancers.i().movie.t().im , fd->stencil.bounds.size() );
      im->copy( fd2->stencil.bounds.near , *ld->data.movie.t().im , fd2->stencil.bounds.size() );
      
      // And it needs a blob that merges the two existing blobs
      b = d->makeFirst(current_frame,current_time);
      b->im = im; 
      b->stats->data.stencil.imitate( fd->stencil );
      b->stats->data.stencil += fd2->stencil;
      b->stats->data.stencil.findBounds();
      b->stats->data.findCentroid(*im,fill_I.keep);
      
      // And now we need to report what happened and throw away the existing dancers
      fates.Append( BlobOriginFate(current_frame,dancers.i().ID,d->ID) );
      fates.Append( BlobOriginFate(current_frame,ld->data.ID,d->ID) );
      dancers.i().movie.Truncate();
      ld->data.movie.Truncate();
      dancers.Destroy(ld);    // Remove 2nd dancer
      dancers.Backspace();   // Current is now our new blob (since we tucked it)
      fd = &dancers.i().movie.t().stats->data;  // Update flood data
      ld = dancers.current;                     // Start over again
    }
  }
  // Finally, validate our new list and set up output
  dancers.start();
  while (dancers.advance())
  {
    if (dancers.i().frames.lo()==0) fates.Append( BlobOriginFate(current_frame,0,dancers.i().ID) );
    dancers.i().validate();
    if (dance_fname!=NULL || img_fname!=NULL) enableOutput(dancers.i());
  }
  
  return dancers.size;
}


double Performance::generateBackgroundAverageIntensity( int *array ) 
{
	int lim = 0;
	long val = 0;
	int shift = ( bg_depth > foreground->depth ) ? bg_depth - foreground->depth : foreground->depth - bg_depth;

	if( array[0] < 0 ) // -1 being our dummy value, only fill the array once.
	{
		for( int a = 0, j = 0; a < NUM_PTS_TO_SAMPLE; a++, j += 2) 
		{			
			lim = 1 + (rand() % (full_area->lines.size-2));
			full_area->lines.start();
			
			for( int i = 0; i < lim; i++ ) full_area->lines.advance();
	
			array[j] = full_area->lines.i().x; 
			array[j+1] = full_area->lines.i().y0 + (rand() % (full_area->lines.i().y1 - full_area->lines.i().y0 - 1));
		}
	} 
	
	for( int i = 0; i < NUM_PTS_TO_SAMPLE*2; i+=2 ) {
		val += background->get(array[i],array[i+1]);
			//val += background->get(array[i],array[i+1]) >> (shift+1); // account for programatic discrepency.
	}
	
	return pow(2,-shift-1)*(((double)val)/NUM_PTS_TO_SAMPLE);
}


// Turn on output for a given dancer, if it's not already on.
void Performance::enableOutput(Dancer& d,bool sitting)
{
  int n = 0;
  if (!sitting && dance_fname!=NULL && d.data_fname==NULL)
  {
    if( !combine_blobs )
    {
      FilenameComponent::loadValues(*dance_fname,current_frame,d.ID,ID);
      n = FilenameComponent::safeLength(*dance_fname);
      d.data_fname = new char[ n ];
      FilenameComponent::toString(*dance_fname,d.data_fname,n);
    }
  }
  if (sitting && sit_fname!=NULL && d.data_fname==NULL)
  {
    FilenameComponent::loadValues(*sit_fname,current_frame,d.ID,ID);
    n = FilenameComponent::safeLength(*sit_fname);
    d.data_fname = new char[ n ];
    FilenameComponent::toString(*sit_fname,d.data_fname,n);
  }
  if (img_fname!=NULL && d.img_namer==NULL)
  {
    d.img_namer = new ManagedList<FilenameComponent>(16,true);
    FilenameComponent::loadValues(*img_fname,current_frame,d.ID,ID);
    img_fname->start();
    while (img_fname->advance())
    {
      d.img_namer->Append( img_fname->i() );
      if (d.img_namer->t().identity == FilenameComponent::perf_n) d.img_namer->t().numberToText();
      else if (d.img_namer->t().identity == FilenameComponent::dancer_n) d.img_namer->t().numberToText();
    }
    FilenameComponent::compact(*d.img_namer);
  }
}


// Create output strings and paths for this performance
bool Performance::prepareOutput(const char *trackerName, const char *path,const char *prefix,bool save_dance,bool save_sit,bool save_img,struct tm* date_to_use)
{
  FilenameComponent* fc;
  if (path==NULL) path = ".";
  if (prefix==NULL) prefix = "out";
  if (strchr(prefix,'/')!=NULL) prefix = "out";
  if (strchr(prefix,'\\')!=NULL) prefix = "out";
  if (access(path,W_OK) != 0) return false;
  if( blobs_fname!=NULL) { delete blobs_fname; blobs_fname=NULL; }
  if (dance_fname!=NULL) { delete dance_fname; dance_fname=NULL; }
  if (sit_fname!=NULL) { delete sit_fname; sit_fname=NULL; }
  if (img_fname!=NULL) { delete img_fname; img_fname=NULL; }
  if (date==NULL)
  {
    date = new struct tm();
    if (date_to_use==NULL)
    {
      time_t now = time(NULL);
#ifdef WINDOWS
      date = localtime(&now);
#else
      localtime_r( &now ,  date );
#endif      
    }
    else *date = *date_to_use;
  }
  
  // Generate directory and try to create it if it's not already there
  ManagedList<FilenameComponent> dir_name(4,true);
  fc = new( dir_name.Append() ) FilenameComponent(); fc->showText(path);
  fc = new( dir_name.Append() ) FilenameComponent(); fc->showText(path_date_connector);
  fc = new( dir_name.Append() ) FilenameComponent(); fc->showNameDate(trackerName, *date);
  fc = new( dir_name.Append() ) FilenameComponent(); fc->showText("/");
  int dirnamelen = FilenameComponent::safeLength(dir_name);
  if (base_directory!=NULL) delete[] base_directory;
  base_directory = new char[dirnamelen];
  FilenameComponent::toString(dir_name,base_directory,dirnamelen);
  if (access(base_directory,W_OK) != 0)
  {
#ifdef WINDOWS
    if (mkdir(base_directory) != 0) return false;  
#else
    if (mkdir(base_directory,01777) != 0) return false;
#endif      
  }
  
  // Save prefix just in case it's useful later
  int prefixlen = strlen(prefix);
  if (prefix_name!=NULL) delete[] prefix_name;
  prefix_name = new char[ prefixlen+1 ];
  strncpy(prefix_name,prefix,prefixlen);
  prefix_name[prefixlen] = 0;
  
  // Generate filename patterns for dancers, sitters, and images
  if (save_dance)
  {
    if( combine_blobs ) 
		{
      blobs_fname = new ManagedList<FilenameComponent>(12,true);
      fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showText(base_directory);
      fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showText(prefix);
      fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showText("_");
      if (expected_n_performances>1)
      {
        fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showPerf(ID,1+(int)floor(log10(expected_n_performances)));
        fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showText("_");
      }
      fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showDancer(0,1+(int)floor(log10(expected_n_dancers)));
      fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showText("k");
      fc = new( blobs_fname->Append() ) FilenameComponent(); fc->showText(".blobs");
      FilenameComponent::compact(*blobs_fname);
		}

      dance_fname = new ManagedList<FilenameComponent>(12,true);
      fc = new( dance_fname->Append() ) FilenameComponent(); fc->showText(base_directory);
      fc = new( dance_fname->Append() ) FilenameComponent(); fc->showText(prefix);
      fc = new( dance_fname->Append() ) FilenameComponent(); fc->showText("_");
      if (expected_n_performances>1)
      {
        fc = new( dance_fname->Append() ) FilenameComponent(); fc->showPerf(ID,1+(int)floor(log10(expected_n_performances)));
        fc = new( dance_fname->Append() ) FilenameComponent(); fc->showText("_");
      }
      fc = new( dance_fname->Append() ) FilenameComponent(); fc->showDancer(0,1+(int)floor(log10(expected_n_dancers)));
      fc = new( dance_fname->Append() ) FilenameComponent(); fc->showText(".blob");
      FilenameComponent::compact(*dance_fname);
  }
  if (save_sit)
  {
    sit_fname = new ManagedList<FilenameComponent>(12,true);
    fc = new( sit_fname->Append() ) FilenameComponent(); fc->showText(base_directory);
    fc = new( sit_fname->Append() ) FilenameComponent(); fc->showText(prefix);
    fc = new( sit_fname->Append() ) FilenameComponent(); fc->showText("_");
    if (expected_n_performances>1)
    {
      fc = new( sit_fname->Append() ) FilenameComponent(); fc->showPerf(ID,1+(int)floor(log10(expected_n_performances)));
      fc = new( sit_fname->Append() ) FilenameComponent(); fc->showText("_");
    }
    fc = new( sit_fname->Append() ) FilenameComponent(); fc->showDancer(0,1+(int)floor(log10(expected_n_sitters)));
    fc = new( sit_fname->Append() ) FilenameComponent(); fc->showText(".ref");
    FilenameComponent::compact(*sit_fname);
  }
  if (save_img)
  {
    img_fname = new ManagedList<FilenameComponent>(12,true);
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showText(base_directory);
    if (expected_n_performances>1)
    {
      fc = new( img_fname->Append() ) FilenameComponent(); fc->showPerf(ID,1+(int)floor(log10(expected_n_performances)));
      fc = new( img_fname->Append() ) FilenameComponent(); fc->showText("_");
    }
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showDancer(0,1+(int)floor(log10(expected_n_dancers)));
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showText("/");
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showText(prefix);
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showText("_");
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showFrame(0,1+(int)floor(log10(expected_n_frames)));
    fc = new( img_fname->Append() ) FilenameComponent(); fc->showText(".tif");
    FilenameComponent::compact(*img_fname);
  }
  return true;
}


// Write out any dangling output--for now, just delete all the dancers and the destructor will do the rest.
bool Performance::finishOutput()
{
  sitters.flush();
  dancers.flush();
  return true;
}


// Write some visual cues on an image that shows the current stuff being tracked
void Performance::imprint(
  Image* im, short borderI, int borderW, short maskI, int maskW, short dancerI, bool show_dancer,
  short sitterI, bool show_sitter, short dcenterI, int dcenterR)
{
  // WARNING--you need to MANUALLY keep this in sync with the 8 bit version!
  if (borderW > 0 && background!=NULL) {
    Rectangle r,R;
    R = background->getBounds();
    r = R; r.far.y = R.near.y; r.expand(borderW-1); im->set(r,borderI);
    r = R; r.far.x = R.near.x; r.expand(borderW-1); im->set(r,borderI);
    r = R; r.near.y = R.far.y; r.expand(borderW-1); im->set(r,borderI);
    r = R; r.near.x = R.far.x; r.expand(borderW-1); im->set(r,borderI);
  }
  if (maskW > 0 && full_area!=NULL) {
    Mask edgemask( *full_area , ManagedList<Strip>::SubordinateList );
    full_area->extractEdge(edgemask);
    edgemask.dilate(maskW);
    im->set(edgemask,maskI);
  }
  if (show_dancer) {
    Mask* m;
    dancers.start();
    while (dancers.advance()) {
      m = &dancers.i().movie.t().mask();
      im->set(*m,dancerI);
    }
  }
  if (show_sitter) {
    Mask* m;
    sitters.start();
    while (sitters.advance()) {     
      m = &sitters.i().movie.t().mask();
      im->set(*m,sitterI);
    }
  }
  if (dcenterR > 0) {
    Point p;
    Rectangle r;
    dancers.start();
    while (dancers.advance()) {
      p = dancers.i().movie.t().data().centroid.toPoint();
      r.near = p - (dcenterR-1);
      r.far = p + (dcenterR-1);
      im->set(r,dcenterI);
      
      if (dancers.i().find_skel && dancers.i().movie.t().skeleton != NULL) {
        Contour& c = dancers.i().movie.t().spine();
        r = im->getBounds();
        c.start();
        while (c.advance()) { if (r.contains( p + c.i() )) im->set(p+c.i() , dcenterI); }
      }
      if (dancers.i().find_outline && dancers.i().movie.t().outline != NULL) {
        Contour& c = dancers.i().movie.t().edge();
        r = im->getBounds();
        c.start();
        while (c.advance()) { if (r.contains( c.i() )) im->set(c.i(),dcenterI); }
      }
    }
  }
}

// Write some visual cues on an image that shows the current stuff being tracked
void Performance::imprint8(
  Image8* im, uint8_t borderI, int borderW, uint8_t maskI, int maskW, uint8_t dancerI, bool show_dancer,
  uint8_t sitterI, bool show_sitter, uint8_t dcenterI, int dcenterR)
{
  // WARNING--you need to MANUALLY keep this in sync with the 16 bit version!
  if (borderW > 0 && background!=NULL) {
    Rectangle r,R;
    R = background->getBounds();
    r = R; r.far.y = R.near.y; r.expand(borderW-1); im->set(r,borderI);
    r = R; r.far.x = R.near.x; r.expand(borderW-1); im->set(r,borderI);
    r = R; r.near.y = R.far.y; r.expand(borderW-1); im->set(r,borderI);
    r = R; r.near.x = R.far.x; r.expand(borderW-1); im->set(r,borderI);
  }
  if (maskW > 0 && full_area!=NULL) {
    Mask edgemask( *full_area , ManagedList<Strip>::SubordinateList );
    full_area->extractEdge(edgemask);
    edgemask.dilate(maskW);
    im->set(edgemask, maskI);
  }
  if (show_dancer) {
    Mask* m;
    dancers.start();
    while (dancers.advance()) {
      m = &dancers.i().movie.t().mask();
      im->set(*m, dancerI);
    }
  }
  if (show_sitter) {
    Mask* m;
    sitters.start();
    while (sitters.advance()) {     
      m = &sitters.i().movie.t().mask();
      im->set(*m, sitterI);
    }
  }
  if (dcenterR > 0) {
    Point p;
    Rectangle r;
    dancers.start();
    while (dancers.advance()) {
      p = dancers.i().movie.t().data().centroid.toPoint();
      r.near = p - (dcenterR-1);
      r.far = p + (dcenterR-1);
      im->set(r,dcenterI);
      
      if (dancers.i().find_skel && dancers.i().movie.t().skeleton != NULL) {
        Contour& c = dancers.i().movie.t().spine();
        r = im->getBounds();
        c.start();
        while (c.advance()) { if (r.contains( p + c.i() )) im->set(p+c.i() , dcenterI); }
      }
      if (dancers.i().find_outline && dancers.i().movie.t().outline != NULL) {
        Contour& c = dancers.i().movie.t().edge();
        r = im->getBounds();
        c.start();
        while (c.advance()) { if (r.contains( c.i() )) im->set(c.i(),dcenterI); }
      }
    }
  }
}


/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_blob_misc()
{
  IOErrorHandler ioeh1;
  IOErrorHandler* ioeh2 = new IOErrorHandler(1,"Funny!");
  ioeh1 = *ioeh2;
  delete ioeh2;
  if (strcmp(ioeh1.what,"Funny!")!=0) return 1;
  
  struct tm date;
  ManagedList<FilenameComponent> mfc(16,true);
  
  date.tm_sec=14;
  date.tm_min=5;
  date.tm_hour=13;
  date.tm_mday=12;
  date.tm_mon=11;
  date.tm_year=2007;
  
  mfc.Append( FilenameComponent() ); mfc.t().showText("Path/");
  mfc.Append( FilenameComponent() ); mfc.t().showNameDate(NULL, date);
  mfc.Append( FilenameComponent() ); mfc.t().showText("/");
  mfc.Append( FilenameComponent() ); mfc.t().showText("perf");
  mfc.Append( FilenameComponent() ); mfc.t().showPerf(0,1);
  mfc.Append( FilenameComponent() ); mfc.t().showText("/dance");
  mfc.Append( FilenameComponent() ); mfc.t().showDancer(0,5);
  mfc.Append( FilenameComponent() ); mfc.t().showText("_frame");
  mfc.Append( FilenameComponent() ); mfc.t().showFrame(0,5);
  mfc.Append( FilenameComponent() ); mfc.t().showText(".txt");
  FilenameComponent::compact(mfc);
  
  int L = FilenameComponent::safeLength(mfc);
  if (FilenameComponent::safeLength(mfc)!=53) return 2;
  
  char result[L];
  FilenameComponent::loadValues(mfc,717,2144,3);
  FilenameComponent::toString(mfc,result,L);
  
  if (!strcmp(result,"Path/20071112_130514/perf3/dance02144_frame00717.txt")) return 3;
  
  Contour c_one(16);
  Contour c_two(16);
  c_one.Append( Point(0,0) );
  c_one.Append( Point(0,1) );
  c_one.Append( Point(0,2) );
  c_one.Append( Point(1,2) );
  c_one.Append( Point(1,3) );
  c_one.Append( Point(2,3) );
  c_one.Append( Point(3,3) );
  c_one.Append( Point(3,2) );
  c_one.Append( Point(2,2) );
  c_one.Append( Point(2,1) );
  c_one.Append( Point(1,1) );
  c_one.Append( Point(1,0) );
  PackedContour pc(c_one);
  pc.unpack(&c_two);
  c_one.start();
  c_two.start();
  while (c_one.advance())
  {
    if (!c_two.advance()) return 4;
    if (c_one.i() != c_two.i()) return 5;
  }
  if (c_two.advance()) return 6;
  
  return 0;
}

int test_mwt_blob_blob() // Doesn't test find()--need Dancer for that
{
  Blob *b;
  
  Storage< Stackable<Strip> > sss(32);
  Storage< Listable<Strip> > sls(32);
  ManagedList<FloodData> mfd(8,true);
  ManagedList<Blob> mb(8,true);
 
  int i;
  Image* im = new Image( Point(25,25), false );
  im->depth = 12;
  *im = im->getGray();
  im->set( Rectangle( Point(10,5) , Point(14,19) ) , im->getGray()/2 );
  int g3_4 = (im->getGray()*3)/4;
  i = im->floodRect( DualRange(1,g3_4) , &sss , &sls , mfd , im->getBounds() );
  if (i!=1) return 1;
  
  b = new( mb.Append() ) Blob(1,0.02);
  b->dancer = NULL;
  b->im = im;
  b->stats = mfd.pop();
  b->pixel_count = b->mask().pixel_count;
  b = new( mb.Append() ) Blob(2,0.04);
  b->dancer = NULL;
  b->clip( mb.h() , mb.h().im , NULL , 2 );
  
  if (b->im->getBounds() != Rectangle( Point(8,3) , Point(16,21) )) return 2;
  if (b->im->get(11,11) != -im->getGray()/2 || b->im->get(9,4) != im->getGray()) return 3;
  if (mb.h().data().centroid != FPoint(12.0,12.0)) return 4;
  mb.h().extraStats( Range(1,g3_4) , FPoint(0.0,0.0) );  // Only finds reference object stats; others require dancer to be non-NULL
  if ( fabs(mb.h().data().major * FPoint(1.0,0.0)) > 1e-4 ) return 5;  // Major axis is Y axis, should be orthogonal to X
  
  FILE* f = fopen("test_blob.log","w");
  mb.h().print(f,"head blob : ");
  mb.h().print(f,NULL);
  fclose(f);
  
  return 0;
}

int test_mwt_blob_dancer() // Doesn't test output--need Performance for that
{
  Storage< Stackable<Strip> > sss(32);
  Storage< Listable<Strip> > sls(32);
  ManagedList<FloodData> mfd(4,true);
  ManagedList<Dancer> mld(2,true);
  Dancer* d = new( mld.Append() ) Dancer(NULL,32,true,DualRange(1,768),DualRange(1,999),3,5,4,1,1,0.02);
  
  Rectangle r1( 15,20 , 15,35 );
  Rectangle r2( 16,21 , 14,33 );
  Image im( Point(64,64), false );
  im.depth = 11;
  
  int G = im.getGray();
  int i;
  bool tf;
  
  im = G;
  im.set( r1 , G/2 );
  i = im.floodRect( DualRange(1,G-1) , &sss , &sls , mfd , im.getBounds() );
  if (i!=1) return 1;
  
  d->setFirst(&im , &mfd.h() , 1 , 0.02);
  if (d->movie.h().mask().pixel_count != r1.area()) return 1;
  if (d->movie.h().data().centroid != 0.5*(FPoint(r1.far)+FPoint(r1.near))) return 2;
  if ( fabs(d->movie.h().data().major * FPoint(1.0,0.0)) > 1e-4) return 3;
  
  im.set( r1 , G );
  im.set( r2 , G/2 );
  d->readyAnother(&im , NULL , 2 , 0.04);
  tf = d->findAnother(true, NULL, FPoint(0, 0));
  if (tf==false) return 4;
  d->validate();
  if (d->movie.t().mask().pixel_count != r2.area()) return 5;
  if (d->movie.t().data().centroid != 0.5*(FPoint(r2.far)+FPoint(r2.near))) return 6;
  
  return 0;
}

int test_mwt_blob_performance()
{
  ManagedList<Performance> mlp(2,true);
  Performance *p = new( mlp.Append() ) Performance(1,32);

  ManagedList<Point> refs(2);
  ManagedList<Point> saved_refs(2);

  struct tm date;
  date.tm_sec=14;
  date.tm_min=5;
  date.tm_hour=13;
  date.tm_mday=12;
  date.tm_mon=11;
  date.tm_year=2007;
  
  bool bad,good;
  int i,N;
  Rectangle r;
  Rectangle r1( 15,20 , 15,35 );
  Rectangle r2( 30,35 , 25,45 );
  Image im( Point(64,64), false );
  im.depth = 10;
  int W = (im.getGray()*2);      // White, or gray in fixed image
  int B = (W*3)/4;  // Background color
  saved_refs.Append( (r1.near+r1.far)/2 );
  refs.Append( saved_refs.h() );
  
  good = p->prepareOutput(NULL, ".","test",true,true,true,&date);
  if (!good) return 1;
  char fname[1024];
  FilenameComponent::toString(*p->dance_fname,fname,1024);
  if (strcmp(fname,"./20071212_130514/test_00000.blob")) return 2;
  FilenameComponent::toString(*p->sit_fname,fname,1024);
  if (strcmp(fname,"./20071212_130514/test_0.ref")) return 3;
  FilenameComponent::toString(*p->img_fname,fname,1024);
  if (strcmp(fname,"./20071212_130514/00000/test_000000.tif")) return 4;
  delete p->img_fname;
  p->img_fname = NULL;
  
  im = B;
  im.set(r1,B/2);
  im.set(r2,B/2);
  p->fill_I = DualRange(1,W-B/8);
  p->fill_size = DualRange( r2.area()/2 , r2.area()*2 );
  p->ref_I = DualRange(1,(3*B)/4);
  p->border = 5;
  p->adapt_rate = 5;
  p->setROI(Ellipse(Point(48,48),40));
  int n_blob = p->initialScan(&im , 0.00);
  int n_refs = p->initialRefs(&im , refs , 0.00);
  
  if (n_refs!=1) return 5;
  if (n_blob!=0) return 6;

  p->find_dancer_skel = true;
  p->find_dancer_edge = true;
  for (N=1;N<128;N++)
  {
    im.set(r1,B/2);
    im.set(r2,B);
    i = N%30;
    if (i>=15) i -= 30;
    if (i>=0) r2 += Point(1,1);
    else r2 -= Point(1,1);
    im.set(r2,B/2);
    if (refs.size==0) refs.Append( saved_refs.h() );  // initialRefs destroys list of points, so make and use a copy
    
    n_blob = p->initialScan(&im , 0.00);
    n_refs = p->initialRefs(&im , refs , 0.00);
  }
  
  if (n_refs!=1 || p->sitters.h().movie.h().mask().pixel_count!=r1.area()) return 7;
  if (n_blob!=1 || p->dancers.h().movie.h().mask().pixel_count!=r2.area()) return 8;

  ManagedList<Profile> no_edges(2, true);
  
  bad = false;
  for (N=128;N<256;N++)
  {
    im.set(r1,B/2);
    im.set(r2,B);
    i = N%30;
    if (i>=15) i -= 30;
    if (i>=0) r2 += Point(1,1);
    else r2 -= Point(1,1);
    im.set(r2,B/2);
    if (refs.size==0) refs.Append( saved_refs.h() );
    
    p->readyNext(&im , no_edges, 0.02*(N-127));
    i = p->findNext();
    if (i==0) bad=true;
    if ( (p->dancers.h().movie.t().data().centroid - 0.5*(FPoint(r2.far)+FPoint(r2.near))).length() > 2 ) bad = true;
  }
  if (bad) return 9;
  
  Image km(im.getBounds(), false);
  km.depth = 8;
  km.copy(im, true);
  Image8 jm(im.getBounds(), false);
  jm.copy16(km);
  for (int x = 0; x < im.size.x; x++) for (int y = 0; y < im.size.y; y++) if (jm.raw(x,y) != km.raw(x,y)) {
    printf("Copied: at %d %d we find %d != %d (originally %d)\n", x, y, jm.raw(x,y), km.raw(x,y), im.raw(x,y));
    return 10;
  }

  p->imprint( &im, 0, 3, W-1, 2, W-1, true, 0, true, 0, 2);

  p->imprint8(&jm, 0, 3, 255, 2, 255, true, 0, true, 0, 2);
  p->imprint( &km, 0, 3, 255, 2, 255, true, 0, true, 0, 2);
  for (int x = 0; x < im.size.x; x++) for (int y = 0; y < im.size.y; y++) if (jm.raw(x,y) != km.raw(x,y)) {
    printf("Imprinted: at %d %d we find %d != %d\n", x, y, jm.raw(x,y), km.raw(x,y));
    return 11;
  }

  im <<= 16-im.depth;
  i = im.writeTiff("performance_imprint.tiff");
  if (i!=0) return 12;

  i = jm.writeTiff("performance_imprint8.tiff");
  if (i != 0) return 13;
  
  good = p->finishOutput();
  if (!good) return 14;
  
  if (p->errors.size!=0)
  {
    p->errors.start();
    while (p->errors.advance())
    {
      printf("Error: #%d says %s\n",p->errors.i().ID,p->errors.i().what);
    }
    return 15;
  }
  
  return 0;
}

int test_mwt_blob()
{
  return test_mwt_blob_misc() + 10*test_mwt_blob_blob() + 1000*test_mwt_blob_dancer() + 100000*test_mwt_blob_performance();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_blob();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Blob test result is %d\n",i);
  return i>0;
}
#endif
