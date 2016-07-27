/* Copyright 2007-2010 Rex A. Kerr, Nicholas A. Swierczek, and HHMI Janelia
 * Copyright 2016 Rex A. Kerr and Calico Life Sciences
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include <math.h>
#include <stdio.h>

#include <new>

#include "MWT_Library.h"



/****************************************************************
                   Useful Utility Methods
****************************************************************/

// Write one line of summary data to disk
void SummaryData::fprint(FILE* f , ManagedList<BlobOriginFate>& mlbof , ManagedList<BlobOutputFate>& mlbf)
{
  // Advance to current frame in blob origin/fate list
  while (mlbof.current!=NULL && mlbof.i().frame<frame_number)
  {
    if (!mlbof.advance()) mlbof.current = NULL;
  }
  while( mlbf.current!=NULL && mlbf.i().frame<frame_number)
  {
    if(!mlbf.advance()) mlbf.current = NULL;
  }
  char last_character='\n';
  if (event_list.size>0 || (mlbof.current!=NULL && mlbof.i().frame==frame_number) || (mlbf.current!=NULL && mlbf.i().frame==frame_number) ) last_character = ' ';
  
  fprintf(f,"%d %.3f  %d %d %.1f  %.2f %.3f  %.1f %.3f  %.1f %.3f  %.3f %.3f  %.3f %.3f%c",
    frame_number,frame_time,
    n_dancers_tracked,n_dancers_good,dancer_persistence,
    dancer_speed,dancer_angularspeed,
    dancer_length,dancer_relativelength,
    dancer_width,dancer_relativewidth,
    dancer_aspect,dancer_relativeaspect,
    dancer_endwiggle,dancer_pixelcount,
    last_character
  );
  
  if (event_list.size > 0) // More to print
  {
    fprintf(f,"%%"); // Marker for event list
    event_list.start();
    while (event_list.advance()) { fprintf(f," 0x%X",event_list.i()); }
  }
  if (mlbof.current!=NULL && mlbof.i().frame==frame_number) // Yet more to print
  {
    fprintf(f," %%%%"); // Marker for worm origin/fate list
    while (mlbof.current!=NULL && mlbof.i().frame==frame_number)
    {
      fprintf(f," %d %d",mlbof.i().ID_origin,mlbof.i().ID_fate);
      if (!mlbof.advance()) mlbof.current=NULL;
    }
  }
  if( mlbf.current!=NULL && mlbf.i().frame==frame_number) // Even more to print
  {
    fprintf(f," %%%%%%"); // marker for worm output fate
    while( mlbf.current!=NULL && mlbf.i().frame==frame_number)
    {
      fprintf(f," %d %d.%ld",mlbf.i().blob_ID,mlbf.i().file_ID,mlbf.i().byte_offset);
      if( !mlbf.advance()) mlbf.current=NULL;
    }
  }
  if (last_character!='\n') fprintf(f,"\n");
}



/****************************************************************
                   TrackerLibrary Methods
****************************************************************/

// Constructor just creates an array of TrackerEntry pointers and sets them to NULL
TrackerLibrary::TrackerLibrary()
  : all_trackers(NULL),default_buffer_size(DEFAULT_DEFAULT_SIZE),active_handles(MAX_TRACKER_HANDLES+1,false),lsstore(1024,false)
{
  all_trackers = new TrackerEntry*[MAX_TRACKER_HANDLES+1];
  for (int i=0 ; i<=MAX_TRACKER_HANDLES ; i++) all_trackers[i] = NULL;
}

// Destructor deletes any handles that are still active
TrackerLibrary::~TrackerLibrary()
{
  if (all_trackers==NULL) return;
  
  if (active_handles.size>0)
  {
    active_handles.start();
    while (active_handles.advance())
    {
      if (all_trackers[active_handles.i()] == NULL) continue;
      delete all_trackers[active_handles.i()];
      all_trackers[active_handles.i()] = NULL;
      active_handles.Backspace();
    }
  }
  
  delete[] all_trackers;
  all_trackers=NULL;
}

// Find the smallest free handle, allocate memory and such for it, and return the handle, or -1 if none are free.
int TrackerLibrary::getNewHandle()
{
  int new_handle;
  if (active_handles.size==MAX_TRACKER_HANDLES) return -1;
  if (active_handles.size==0 || active_handles.h()!=1)
  {
    new_handle=1;
    active_handles.Push(new_handle);
    active_handles.start();
    active_handles.advance();
  }
  else
  {
    new_handle = active_handles.i();
    while (active_handles.advance())
    {
      if (active_handles.i()>new_handle+1) break;
      new_handle = active_handles.i();
    }
    new_handle++;
    if (new_handle<active_handles.i())
    {
      active_handles.Tuck(new_handle);
      active_handles.retreat();
    }
    else if (new_handle <= MAX_TRACKER_HANDLES)
    {
      active_handles.Append(new_handle);
      active_handles.advance();
    }
    else return -1;
  }
  
  if (all_trackers[new_handle] != NULL) delete all_trackers[new_handle];
  all_trackers[new_handle] = new TrackerEntry(active_handles.current , default_buffer_size);
  
  return new_handle;
}

// Delete a handle.  Return the number of the handle if it exists, or -1 if not.
int TrackerLibrary::killHandle(int victim)
{
  if (victim<1 || victim>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[victim]==NULL) return -1;
  active_handles.Destroy( all_trackers[victim]->handle );
  if (active_handles.size>0) { active_handles.start(); active_handles.advance(); }
  delete all_trackers[victim];
  all_trackers[victim] = NULL;
  return victim;
}

// Set the output type. This should only be done at initialization.
int TrackerLibrary::setCombineBlobs( int handle, bool type ) 
{
  if( handle < 1 || handle > MAX_TRACKER_HANDLES ) return -1;
  if( all_trackers[handle] == NULL ) return -1;
  all_trackers[handle]->performance.setCombineBlobs( type );
  return handle;
}


// Set the date in pieces.  Normally, don't do this--just let it pick its own date, or have one or all others adopt the date of one
int TrackerLibrary::setDate(int handle,int year,int month,int day,int hour,int minute,int second)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  all_trackers[handle]->setDate(year,month,day,hour,minute,second);
  return handle;  
}

// Set the date from another handle
int TrackerLibrary::borrowDate(int handle,int donor_handle)
{
  TrackerEntry *te;
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  if (donor_handle<1 || donor_handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[donor_handle]==NULL) return -1;
  
  te = all_trackers[donor_handle];
  if (te->output_date==NULL)
  {
    if (all_trackers[handle]->output_date!=NULL) delete all_trackers[handle]->output_date;
    all_trackers[handle]->output_date = NULL;
  }
  else
  {
    all_trackers[handle]->setDate(te->output_date->tm_year,te->output_date->tm_mon+1,te->output_date->tm_mday,
                                  te->output_date->tm_hour,te->output_date->tm_min,te->output_date->tm_sec);
  }
  
  return handle;  
}

// Make all handles have the same date as this one
int TrackerLibrary::setAllDatesToMine(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  Listable<int> *iterator;
  active_handles.start(iterator);
  while (active_handles.advance(iterator))
  {
    if (handle == iterator->data) continue;
    else
    {
      int i = borrowDate(handle,iterator->data);
      if (i != handle) return i;
    }
  }
  
  return handle;
}

// Set filename stuff, plus say what we want to output
int TrackerLibrary::setOutput(int handle,const char *path,const char *prefix,bool save_obj,bool save_ref,bool save_im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  te->setPath(path);
  te->setPrefix(prefix);
  te->save_objects = ( save_obj == 1 ) ? true : false;
  te->save_refs = ( save_ref == 1 ) ? true : false;
  te->save_images = (save_im == 1 ) ? true : false;
  te->output_info_known = true;
  
  return handle;  
}

// Actually create output directories and get ready to go.
int TrackerLibrary::beginOutput(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->output_info_known) return 0;
  bool good = te->performance.prepareOutput(te->path_string,te->prefix_string,te->save_objects,te->save_refs,te->save_images,te->output_date);
  if (!good) return 0;
  te->output_started = true;

  return handle;  
}


// Specify some important attributes of the image we'll be dealing with
int TrackerLibrary::setImageInfo(int handle,int bits,int width,int height)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];  
  if (bits<1) bits=1;
  if (bits>14) bits=14;
  te->image_bits = bits;
  if (width<1) width=1;
  te->image_width = width;
  if (height<1) height=1;
  te->image_height = height;
  te->image_info_known = true;
  
  return handle;  
}

// Query the bit depth that the tracker is expecting
// Does NOT return handle; returns bit depth that has been set, or -1 or 0 if error
int TrackerLibrary::getBitDepth(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  return te->image_bits;
}


// Sets the region of interest to a rectangle
int TrackerLibrary::setRectangle(int handle,int left,int right,int top,int bottom)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (right<left) swap(left,right);
  if (bottom<top) swap(top,bottom);
  if (left<0) left = 0;
  if (right>=te->image_width) right = te->image_width-1;
  if (top<0) top = 1;
  if (bottom>=te->image_height) bottom = te->image_height-1;

  te->performance.setROI( Rectangle( Point(left,top) , Point(right,bottom) ) );
  
  return handle;  
}

// Sets the region of interest to an ellipse
int TrackerLibrary::setEllipse(int handle,int centerx,int centery,int radiusx,int radiusy)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (radiusx<1) radiusx = 1;
  if (radiusy<1) radiusy = 1;
  if (radiusx>te->image_width+te->image_height) radiusx = te->image_width+te->image_height;
  if (radiusy>te->image_height+te->image_width) radiusy = te->image_height+te->image_width;
  if (centerx<0) centerx = 0;
  if (centery<0) centery = 0;
  if (centerx>=te->image_width) centerx = te->image_width-1;
  if (centery>=te->image_height) centery = te->image_height-1;

  te->performance.setROI( Ellipse( Point(centerx,centery) , Point(radiusx,radiusy) ) );
  
  return handle;  
}

// Add a rectangle to the existing region of interest
int TrackerLibrary::addRectangle(int handle,int left,int right,int top,int bottom)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;

  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (right<left) swap(left,right);
  if (bottom<top) swap(top,bottom);
  if (left<0) left = 0;
  if (right>=te->image_width) right = te->image_width-1;
  if (top<0) top = 1;
  if (bottom>=te->image_height) bottom = te->image_height-1;
  
  Performance& p = te->performance;
  if (p.full_area==NULL) p.setROI( Rectangle( Point(left,top) , Point(right,bottom) ) );
  else
  {
    Mask m( Rectangle( Point(left,top) , Point(right,bottom) ) , &lsstore );
    *p.full_area += m;
    p.full_area->findBounds();
  }

  return handle;  
}

// Add an ellipse to the region of interest
int TrackerLibrary::addEllipse(int handle,int centerx,int centery,int radiusx,int radiusy)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;

  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (radiusx<1) radiusx = 1;
  if (radiusy<1) radiusy = 1;
  if (radiusx>te->image_width+te->image_height) radiusx = te->image_width+te->image_height;
  if (radiusy>te->image_height+te->image_width) radiusy = te->image_height+te->image_width;
  if (centerx<0) centerx = 0;
  if (centery<0) centery = 0;
  if (centerx>=te->image_width) centerx = te->image_width-1;
  if (centery>=te->image_height) centery = te->image_height-1;

  Performance& p = te->performance;
  if (p.full_area==NULL) p.setROI( Ellipse( Point(centerx,centery) , Point(radiusx,radiusy) ) );
  else
  {
    Mask m( Ellipse( Point(centerx,centery) , Point(radiusx,radiusy) ) , &lsstore );
    *p.full_area += m;
    p.full_area->findBounds();
  }

  return handle;  
}

// Remove a rectangle from the region of interest.
int TrackerLibrary::cutRectangle(int handle,int left,int right,int top,int bottom)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;

  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (right<left) swap(left,right);
  if (bottom<top) swap(top,bottom);
  if (left<0) left = 0;
  if (right>=te->image_width) right = te->image_width-1;
  if (top<0) top = 1;
  if (bottom>=te->image_height) bottom = te->image_height-1;
  
  Performance& p = te->performance;
  
  if (p.full_area!=NULL) p.setROI( Rectangle(Point(0,0) , Point(te->image_width-1,te->image_height-1) ) );
  *p.full_area -= Rectangle( Point(left,top) , Point(right,bottom) );
  p.full_area->findBounds();

  return handle;  
}

// Remove an ellipse from the region of interest.
int TrackerLibrary::cutEllipse(int handle,int centerx,int centery,int radiusx,int radiusy)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;

  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (radiusx<1) radiusx = 1;
  if (radiusy<1) radiusy = 1;
  if (radiusx>te->image_width+te->image_height) radiusx = te->image_width+te->image_height;
  if (radiusy>te->image_height+te->image_width) radiusy = te->image_height+te->image_width;
  if (centerx<0) centerx = 0;
  if (centery<0) centery = 0;
  if (centerx>=te->image_width) centerx = te->image_width-1;
  if (centery>=te->image_height) centery = te->image_height-1;

  Performance& p = te->performance;

  if (p.full_area!=NULL) p.setROI( Rectangle(Point(0,0) , Point(te->image_width-1,te->image_height-1) ) );

  Mask m( Ellipse( Point(centerx,centery) , Point(radiusx,radiusy) ) , &lsstore );
  *p.full_area -= m;
  p.full_area->findBounds();

  return handle;  
}

// Draw the region onto an image
int TrackerLibrary::showROI(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
	short white = (1<<te->image_bits)-1;
	short color = te->performance.blob_is_dark ? 1 : white;  // Reference object will be oppositely colored from object
	im.divide_bg = te->performance.use_division;
  te->performance.imprint(&im , 0 , 3 , color , 2 , color , false , color , false , color , 0);
    
  return handle;  
}


// Copy and resize an image--doesn't actually need a handle
int TrackerLibrary::resizeRescale(int handle,Image& source,Image& dest,Rectangle dest_selection)
{
  dest.mimic(source,dest_selection,source.getBounds());
  return handle;
}

// Copy and resize the memory image--does need a handle
int TrackerLibrary::resizeRescaleMemory(int handle,Image& dest,Rectangle dest_selection)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  Performance& p = all_trackers[handle]->performance;
  if (p.background==NULL) return 0;
  dest.mimic( *p.background , dest_selection , p.background->getBounds() );
  return handle;
}

// Copy and resize the fixed image--needs a handle, and only works after a scanObjects
int TrackerLibrary::resizeRescaleFixed(int handle,Image& dest,Rectangle dest_selection)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  Performance& p = all_trackers[handle]->performance;
  if (p.foreground==NULL) return 0;
  dest.mimic( *p.foreground , dest_selection , p.foreground->getBounds() );
  return handle;
}


// Set the border around each tracked object so we can follow it moving
int TrackerLibrary::setDancerBorderSize(int handle,int border)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  if (border<1) return 0;
  
  TrackerEntry* te = all_trackers[handle];
  te->border_size_known = true;
  
  te->performance.border = border;
  
  return handle;  
}

// During load/process cycle, update image a strip at a time; this sets how many strips fit across the image (default is 16)
int TrackerLibrary::setUpdateBandNumber(int handle,int n_bands)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];

  if (n_bands<2) n_bands = 2;
  te->performance.n_scan_bands = n_bands;
  
  return handle;
}


// Set the intensity above or below which we will fill a reference object
int TrackerLibrary::setRefIntensityThreshold(int handle,int intensity_low,int intensity_high)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  if (intensity_low < 1) intensity_low = 1;
  if (intensity_high > (1<<te->image_bits)-1) intensity_high = (1<<te->image_bits)-1;
  if (intensity_low > intensity_high) intensity_low = intensity_high;
  te->reference_intensities_known = true;
  te->performance.ref_I = DualRange( Range(intensity_low,intensity_high) );
  
  return handle;  
}

// Add a new reference object; must have defined a region of interest
int TrackerLibrary::addReferenceObjectLocation(int handle,int x,int y)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  if (x<0 || y<0 || x>=te->image_width || y>=te->image_height) return 0;

  te->reference_objects.Append( Point(x,y) );
 
  return handle;
}

// Remove the last one (I guess we didn't like it)
int TrackerLibrary::removeLastReferenceObject(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (te->reference_objects.size==0) return 0;
  te->reference_objects.Truncate();
  
  return handle;  
}

// Try to find reference objects and return number found (-2 if initialization isn't done)
int TrackerLibrary::scanRefs(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->reference_intensities_known) return -2;
	im.divide_bg = te->performance.use_division;
  Performance& p = te->performance;
  te->working.imitate( te->reference_objects );
  p.initialRefs(&im , te->working , 0.0);
  te->references_found = true;
  
  return te->reference_objects.size - te->working.size;
}

// Show what reference objects were found
int TrackerLibrary::showRefs(int handle,Image& im,bool show_as_white)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->references_found) return 0;
	im.divide_bg = te->performance.use_division;
  short color = 1;
  if (show_as_white) color = (1<<te->image_bits)-1;
  te->performance.imprint(&im , 0 , 3 , (1<<te->image_bits)-1 , 2 , 0 , false , color , true , 0 , 0);
  
  return handle;  
}


// Set intensity thresholds for moving objects
int TrackerLibrary::setObjectIntensityThresholds(int handle,int intensity_to_fill,int intensity_of_new)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known) return 0;
  
  short white = 2*(1<<te->image_bits)-1;
  short gray = 1<<te->image_bits;
  if (intensity_to_fill == gray || intensity_of_new==gray) return 0;
  else if (intensity_to_fill > gray) // Bright blob
  {
    if (intensity_of_new < gray) return 0;
    if (intensity_of_new > white) intensity_of_new = white;
    if (intensity_of_new < intensity_to_fill) intensity_to_fill = intensity_of_new;
    te->performance.blob_is_dark = false;
    te->performance.fill_I = DualRange( Range(intensity_of_new,white) , Range(intensity_to_fill,white) );

  }
  else // Standard dark blob
  {
    if (intensity_of_new > gray) return 0;
    if (intensity_of_new < 1) intensity_of_new = 1;
    if (intensity_of_new > intensity_to_fill) intensity_to_fill = intensity_of_new;
    te->performance.blob_is_dark = true;
    te->performance.fill_I = DualRange( Range(1,intensity_of_new) , Range(1,intensity_to_fill) );
  }
  te->object_intensities_known = true;
  
  return handle;  
}

// Set how big the objects are allowed to be
int TrackerLibrary::setObjectSizeThresholds(int handle,int small_size,int small_good_size,int large_good_size,int large_size)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (small_size<1) small_size = 1;
  if (small_good_size<small_size) small_good_size = small_size;
  if (large_good_size<small_good_size) large_good_size = small_good_size;
  if (large_size<large_good_size) large_size = large_good_size;
  te->performance.fill_size = DualRange( Range(small_good_size,large_good_size) , Range(small_size,large_size) );
  te->object_sizes_known = true;
  
  return handle;  
}

// Set how long we need to track the objects in order to write them out
int TrackerLibrary::setObjectPersistenceThreshold(int handle,int frames)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  
  if (frames<2) frames=2;
  te->performance.n_long_enough = frames;
  
  if (frames<5) frames=5; // Guess we might need 5 for analysis purposes--revisit this later.
  te->performance.n_keep_full = frames;
  
  te->object_persistence_known = true;
  
  return handle;
}  


// Set rate of 2^-alpha for how much to update the background image from the current image
int TrackerLibrary::setAdaptationRate(int handle,int alpha)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];

  if (alpha<0) alpha = 0;
  if (alpha>12) alpha = 12;
  te->performance.adapt_rate = alpha;
  
  te->adaptation_rate_known = true;
  
  return handle;
}


// Finds moving objects in the image without tracking them
// Does NOT return the handle; returns the number of objects found
int TrackerLibrary::scanObjects(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->border_size_known || !te->object_intensities_known || 
      !te->object_sizes_known || !te->adaptation_rate_known)
  {
    return -2;
  }

  im.divide_bg = te->performance.use_division;
  int n = te->performance.initialScan(&im,0.0);
  te->objects_found = true;
  return n;
}

// Finds moving objects and writes stuff back to the image to show what happened
int TrackerLibrary::showObjects(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->objects_found) return 0;

  im.divide_bg = te->performance.use_division;
  short color = 0;
  short other = 0;
  if (te->performance.blob_is_dark) color = (1<<te->image_bits)-1;
  else other = (1<<te->image_bits)-1;
  te->performance.imprint(&im , 0 , 3 , (1<<te->image_bits)-1 , 2 , color , true , 0 , false , other , 2);  
  
  return handle;
}


// How long we let an object move before we calculate its velocity (default = 1.0s)
int TrackerLibrary::setVelocityIntegrationTime(int handle,float interval)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  if (interval<0.01) interval=0.01;
  
  TrackerEntry* te = all_trackers[handle];
  
  te->update_frequency = interval;

  return handle;
}

// Find skeletons for dancers (off by default, this turns it on with the default argument of enable=true)
int TrackerLibrary::enableSkeletonization(int handle,bool enable)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  Performance& p = all_trackers[handle]->performance;
  
  p.find_dancer_skel = enable;

  return handle;
}

// Find outlines for dancers (off by default, this turns it on with the default argument of enable=true)
int TrackerLibrary::enableOutlining(int handle,bool enable)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  Performance& p = all_trackers[handle]->performance;
  
  p.find_dancer_edge = enable;

  return handle;
}


// Prepare to load an image for processing a piece at a time; returns how many pieces will be needed
int TrackerLibrary::prepareImagePieces(int handle,float time)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->border_size_known || !te->object_intensities_known || 
      !te->object_sizes_known || !te->adaptation_rate_known || !te->object_persistence_known ||
      !te->reference_intensities_known)
  {
    return 0;
  }
	
  return te->performance.anticipateNext(time);
}

// Find the size of the next piece of the image to be loaded; returns the handle if more pieces remain, 0 otherwise
int TrackerLibrary::getNextPieceCoords(int handle,Rectangle& coords)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->border_size_known || !te->object_intensities_known || 
      !te->object_sizes_known || !te->adaptation_rate_known || !te->object_persistence_known ||
      !te->reference_intensities_known)
  {
    return 0;
  }
  
  if (te->performance.findNextItemBounds(coords)) return handle;
  else return 0;
}

// Load the next piece of the image, returns handle on success
int TrackerLibrary::loadThisImagePiece(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->border_size_known || !te->object_intensities_known || 
      !te->object_sizes_known || !te->adaptation_rate_known || !te->object_persistence_known ||
      !te->reference_intensities_known)
  {
    return 0;
  }
	
  Performance& p = te->performance;
	im.divide_bg = p.use_division;
  p.loadNextSingleItem(&im);
  if (p.load_state == Performance::all_loaded)
  {
    te->image_loaded = true;
    te->statistics_ready = false;
    SummaryData* sd = new( te->summary.Append() ) SummaryData( p.current_frame , p.current_time , &te->eventstore );
    sd->update_frequency = te->update_frequency;
  }
  return handle;
}


// Load an image for processing
int TrackerLibrary::loadImage(int handle,Image& im,float time)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->border_size_known || !te->object_intensities_known || 
      !te->object_sizes_known || !te->adaptation_rate_known || !te->object_persistence_known ||
      !te->reference_intensities_known || !te->objects_found)
  {
    return 0;
  }
  
  im.divide_bg = te->performance.use_division;	
  Performance& p = te->performance;
  p.readyNext(&im,time);
  te->image_loaded = true;
  te->statistics_ready = false;
  SummaryData* sd = new( te->summary.Append() ) SummaryData( p.current_frame , time , &te->eventstore );
  sd->update_frequency = te->update_frequency;
  
  return handle;
}


// Load an image for processing--the slow way, by copying out pieces
int TrackerLibrary::loadImageAsPieces(int handle,Image& im,float time)
{
  int howmany = prepareImagePieces(handle,time);
  if (howmany < 1) return howmany;
  int i;
  Rectangle piece_coords;
  for ( ; howmany>0 ; howmany--) {
    i = getNextPieceCoords(handle,piece_coords);
    if (i!=handle) return i;
    
    Image *temp = new Image(im,piece_coords,false);
    i = loadThisImagePiece(handle,*temp);
    if (i!=handle) return i;
    delete temp;
  }
  
  if (getNextPieceCoords(handle,piece_coords)!=0) return -2;  // Somehow we didn't get them all
  return handle;
}


// Show what we've loaded
int TrackerLibrary::showLoaded(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_loaded) return 0;
  im.divide_bg = te->performance.use_division;
  Performance& p = te->performance;
  
  if (p.band != NULL && p.band_area !=NULL)
  {
    im.copy(*p.band,*p.band_area,true);
  }
  if (p.dancers.size > 0)
  {
    p.dancers.start();
    while (p.dancers.advance())
    {
      if (p.dancers.i().movie.size==0) continue;
      if (p.dancers.i().movie.t().im==NULL) continue;
      im.copy( *(p.dancers.i().movie.t().im) , true );
    }
  }
  
  return handle; 
}

// Add an integer event number to our summary data to mark stimuli, etc.
int TrackerLibrary::markEvent(int handle,int event)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->output_started || te->summary.size==0) return 0;
  
  te->summary.t().event_list.Append(event);
  
  return handle; 
  
}

// Process the loaded image and calculate summary statistics
// Does NOT return the handle; returns the number of moving objects found
int TrackerLibrary::processImage(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_loaded || !te->output_started) return 0;
  
  Performance& p = te->performance;
  int n_dancers = p.findNext();
  te->image_loaded = false;
  
  // Need to calculate summary statistics
  SummaryData &sd = te->summary.t();
  
  sd.n_dancers_tracked = p.dancers.size;
  
  sd.n_dancers_good = 0;
  p.dancers.start();
  while (p.dancers.advance())
  {
    if (p.dancers.i().n_long_enough > p.dancers.i().frames.length()) continue;
    sd.n_dancers_good++;
  }
  
  Blob* b;
  Datum persist;
  Datum length;
  Datum rel_length;
  Datum width;
  Datum rel_width;
  Datum aspect;
  Datum rel_aspect;
	Datum pixelcount;
  p.dancers.start();
  while (p.dancers.advance())
  {
    b = &(p.dancers.i().movie.t());
    if (p.dancers.i().n_long_enough > p.dancers.i().frames.length()) continue;
    persist.add( p.dancers.i().times.length() );
    
    length.add( b->data().long_axis );
    if (p.dancers.i().accumulated_length.mean()==0) rel_length.add(0.0);
    else rel_length.add( b->data().long_axis / p.dancers.i().accumulated_length.mean() );
    
    width.add( b->data().short_axis );
    if (p.dancers.i().accumulated_width.mean()==0) rel_width.add(0.0);
    else rel_width.add( b->data().short_axis / p.dancers.i().accumulated_width.mean() );
    
    aspect.add( b->data().short_axis / b->data().long_axis );
    if (p.dancers.i().accumulated_aspect.mean()==0) rel_aspect.add(0.0);
    else rel_aspect.add( (b->data().short_axis / b->data().long_axis) / p.dancers.i().accumulated_aspect.mean() );
		
		pixelcount.add( b->data().stencil.pixel_count );		
  }
  sd.dancer_persistence = persist.mean();
  sd.dancer_length = length.mean();
  sd.dancer_relativelength = rel_length.mean();
  sd.dancer_width = width.mean();
  sd.dancer_relativewidth = rel_width.mean();
  sd.dancer_aspect = aspect.mean();
  sd.dancer_relativeaspect = rel_aspect.mean();
  sd.dancer_pixelcount = pixelcount.mean();
  if (sd.n_speed_updates * sd.update_frequency < sd.frame_time)
  {
    sd.n_speed_updates = (int)ceil( sd.frame_time / sd.update_frequency );
    Datum speed;
    Datum ang_speed;
    p.dancers.start();
    while (p.dancers.advance())
    {
      if (p.dancers.i().n_long_enough > p.dancers.i().frames.length()) continue;
      speed.add( p.dancers.i().recentSpeed( sd.update_frequency*0.5 ) );
      ang_speed.add( p.dancers.i().recentAngularSpeed( sd.update_frequency*0.5 ) );
    }
    sd.dancer_speed = speed.mean();
    sd.dancer_angularspeed = ang_speed.mean();
  }
  else
  {
    te->summary.end();
    if (te->summary.retreat())
    {
      sd.dancer_speed = te->summary.i().dancer_speed;
      sd.dancer_angularspeed = te->summary.i().dancer_angularspeed;
    }
  }
  
  if (p.find_dancer_skel)
  {
    // Figure out how far off the center line the nose or tail is pointing (report whichever deviates further)
    int i;
    float f;
    float nose_angle,tail_angle;
    FPoint nose,tail,noseless,tailless;  // Nose and tail are arbitrary; they refer to two ends, but may not correspond to organism's nose and tail
    Contour* c;
    Datum wiggle;
    p.dancers.start();
    while (p.dancers.advance())
    {
      if (p.dancers.i().n_long_enough > p.dancers.i().frames.length()) continue;
      if (p.dancers.i().movie.t().skeleton==NULL) continue;
      c = &(p.dancers.i().movie.t().spine());
      if (c->size() != Blob::SKELETON_SIZE) continue;
      
      c->start();
      for (i=0 ; (i+1)*5 < Blob::SKELETON_SIZE ; i++) c->advance();  // Pick off first 20% of dancer (nose)
      nose = c->h() - c->i();
      for ( ; (i+1)*3 < Blob::SKELETON_SIZE ; i++) c->advance();  // Pick off last 2/3 of dancer (all but nose)
      noseless = c->i() - c->t();
      f = nose.length();
      if (f==0) nose_angle = 0;
      else
      {
        nose /= f;
        f = noseless.length();
        if (f==0) nose_angle = 0;
        else
        {
          noseless /= f;
          nose_angle = nose * noseless;
        }
      }
      
      c->end();
      for (i=0 ; (i+1)*5 < Blob::SKELETON_SIZE ; i++) c->retreat();  // Pick off last 20% of dancer (tail)
      tail = c->t() - c->i();
      for ( ; (i+1)*3 < Blob::SKELETON_SIZE ; i++) c->retreat();  // Pick off first 2/3 of dancer (all but tail)
      tailless = c->i() - c->h();
      f = tail.length();
      if (f==0) tail_angle = 0;
      else
      {
        tail /= f;
        f = tailless.length();
        if (f==0) tail_angle = 0;
        else
        {
          tailless /= f;
          tail_angle = tail * tailless;
        }
      }
      
      if (nose_angle < tail_angle) wiggle.add( acos(nose_angle) );
      else wiggle.add( acos(tail_angle) );
    } 
    sd.dancer_endwiggle = wiggle.mean();
  }
  
  te->statistics_ready = true;

  return n_dancers;
}

// Write back to the image so we can see what happened
int TrackerLibrary::showResults(int handle,Image& im)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->image_info_known || !te->border_size_known || !te->object_intensities_known || 
      !te->object_sizes_known || !te->adaptation_rate_known || !te->object_persistence_known ||
      !te->reference_intensities_known)
  {
    return 0;
  }
  im.divide_bg = te->performance.use_division;
  Performance& p = te->performance;
  short white = (1<<te->image_bits)-1;
  short colorR = p.blob_is_dark ? 1 : white;  // Reference object will be oppositely colored from object
  short colorO = p.blob_is_dark ? white : 1;
  short colorC = p.blob_is_dark ? 1 : white;
  p.imprint(&im , 0 , 3 , white , 2 , colorO , true , colorR , true , colorC , 2);  
  return handle;
}

// Check for errors and clear the error buffer
// Does NOT return the handle; returns the number of errors (so we want 0!)
int TrackerLibrary::checkErrors(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  Performance& p = all_trackers[handle]->performance;
  int n = p.errors.size;
  p.errors.flush();  // Now that they know how many errors there are, we clear them.
  
  return n;
}

// Clean up everything after we're done, and write our summary information
int TrackerLibrary::complete(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->output_started) return 0;

  Performance& p = te->performance;
  p.finishOutput();
  te->output_started = false;
  
  FilenameComponent* fc;
  ManagedList<FilenameComponent> summary_fname(4,true);
  fc = new( summary_fname.Append() ) FilenameComponent(); fc->showText( p.base_directory );
  fc = new( summary_fname.Append() ) FilenameComponent(); fc->showText( p.prefix_name );
  fc = new( summary_fname.Append() ) FilenameComponent(); fc->showText( ".summary" );
  
  int summarylen = FilenameComponent::safeLength(summary_fname);
  char summary[ summarylen ];
  FilenameComponent::toString( summary_fname, summary , summarylen );
  
  FILE *f = fopen(summary,"w");
  if (!f) return 0;
  
  p.fates.start();
  p.fates.advance(); // Set current to first element of fates, if it exists
  p.output_fates.start();
  p.output_fates.advance();
  te->summary.start();
  while (te->summary.advance())
  {
    te->summary.i().fprint(f,p.fates,p.output_fates);  // This will move through fates list as needed
  }
  fclose(f);

  return handle;
}

int TrackerLibrary::setDivisionImageCorrectionAlgorithm( int handle )
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
	
	TrackerEntry* te = all_trackers[handle];
	te->performance.use_division = true;
	return handle; 		
}

int TrackerLibrary::setSubtractionImageCorrectionAlgorithm( int handle ) 
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
	
	TrackerEntry* te = all_trackers[handle];
	te->performance.use_division = false;
	return handle; 		
}

int TrackerLibrary::reportImageCorrectionAlgorithm( int handle ) 
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
	
	TrackerEntry* te = all_trackers[handle];
	
	return (te->performance.use_division) ? 1 : 0; 		
}


double TrackerLibrary::generateBackgroundAverageIntensity( int handle, int *array ) 
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
	
	TrackerEntry* te = all_trackers[handle];
	
	return te->performance.generateBackgroundAverageIntensity(array); 	
}

// Report the number of dancers being tracked
int TrackerLibrary::reportNumber(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().n_dancers_tracked;
}

// Report the number of dancers who have been around at least as long as persistence threshold
int TrackerLibrary::reportNumberPersisting(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().n_dancers_good;  
}

// Report the mean number of seconds the current dancers have been tracked
float TrackerLibrary::reportPersistence(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_persistence;
}

// Reports the mean speed (in pixels/second) since the last request
float TrackerLibrary::reportSpeed(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_speed;
}

// Reports the mean angular speed (in radians/second) since the last request
float TrackerLibrary::reportAngularSpeed(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_angularspeed;
}

// Reports the mean length in pixels
float TrackerLibrary::reportLength(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_length;
}

// Reports the mean (current length / mean length)
float TrackerLibrary::reportRelativeLength(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_relativelength;
}

// Reports the mean width in pixels
float TrackerLibrary::reportWidth(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_width;
}

// Reports the mean (current width / mean width)
float TrackerLibrary::reportRelativeWidth(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_relativewidth;
}

// Reports the aspect ratio (width / length)
float TrackerLibrary::reportAspect(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_aspect;
}

// Reports the mean (current aspect ratio / mean aspect ratio)
float TrackerLibrary::reportRelativeAspect(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_relativeaspect;
}

// Reports the mean angle between the direction of the head or tail and direction of body
float TrackerLibrary::reportEndWiggle(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_endwiggle;
}

// Reports the mean pixel count of found objects
float TrackerLibrary::reportObjectPixelCount(int handle)
{
  if (handle<1 || handle>MAX_TRACKER_HANDLES) return -1;
  if (all_trackers[handle]==NULL) return -1;
  
  TrackerEntry* te = all_trackers[handle];
  if (!te->statistics_ready) return -1;
  
  return te->summary.t().dancer_pixelcount;
}



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_library(int bin)
{
  TrackerLibrary a_library;
  
  // Unit tests only support up to binning of 3 (objects may be too small otherwise)
  if (bin<=0) bin=1;
  if (bin>3) bin=3;
  
  Image arena( Point(512,512), false );
  arena.depth = 10;
  int W = 1<<10;
  
  ModelCamera cam(arena);
#ifdef WINDOWS
  srand(3151);
#else
  srandom(3151);
#endif
  cam.scaled_noise = 1.0;
  cam.fixed_noise = 50;
  cam.scatterOffset(25);
  cam.rowsOffset(10);
  cam.scatterScale(0.05);
  cam.colsScale(0.05);

  ModelWorm worm1,worm2,worm3;
  worm1.setPose( 128+FPoint(60,40),FPoint(1,0),0 );
  worm1.setSize(40,2.5);
  worm1.setWiggle(4.0,30,0.7,-0.02);
  worm2.setPose( 128+FPoint(150,30),FPoint(-1,0),0 );
  worm2.setSize(45,3.0);
  worm2.setWiggle(4.5,32,0.68,0.02);
  worm3.setPose( 128+FPoint(125,230),FPoint(0,1),0);
  worm3.setSize(50,3.0);
  worm3.setWiggle(5.0,35,0.72,0.01);
  
  int i;
  int h1 = a_library.getNewHandle();
  int h2 = a_library.getNewHandle();
  a_library.setCombineBlobs(h1,true);
  a_library.setCombineBlobs(h2,true);
  if (h1!=1 || h2!=2) return 1;
  i = a_library.setImageInfo(h1,10,512/bin,512/bin); if (i!=h1) return 2;
  i = a_library.setImageInfo(h2,10,512/bin,512/bin); if (i!=h2) return 3;
  i = a_library.setDate(h1,2007,12,26,10,50,33); if (i!=h1) return 4;
  i = a_library.borrowDate(h2,h1); if (i!=h2) return 5;
  i = a_library.setOutput(h1,".","test",true,false,false); if (i!=h1) return 6;
  i = a_library.setOutput(h2,".","tset",false,true,false); if (i!=h2) return 7;
  i = a_library.setRectangle(h1,5/bin,(256+200)/bin,5/bin,(256+180)/bin); if (i!=h1) return 8;
  i = a_library.addRectangle(h1,100/bin,(250+256)/bin,80/bin,(256+250)/bin); if (i!=h1) return 9;
  i = a_library.setEllipse(h2,30/bin,(220+256)/bin,25/bin,25/bin); if (i!=h2) return 10;
  
  arena = (3*W)/4;
  arena.bin = bin;
  i = a_library.showROI(h1,arena); if (i!=h1) return 11;
  i = a_library.showROI(h2,arena); if (i!=h2) return 12;
  arena <<= 6; arena.writeTiff("test_two_ROIs.tiff");
  arena.bin = 0;
  
  i = a_library.setDancerBorderSize(h1,10/bin); i = a_library.setDancerBorderSize(h2,10/bin);
  i = a_library.setRefIntensityThreshold(h1,1,W/2); if (i!=h1) return 13;
  i = a_library.setRefIntensityThreshold(h2,1,W/2); if (i!=h2) return 14;
  i = a_library.addReferenceObjectLocation(h2,35/bin,(256+220)/bin); if (i!=h2) return 15;
  i = a_library.setObjectIntensityThresholds(h1,W-W/20,W-W/8); if (i!=h1) return 16;
  i = a_library.setObjectIntensityThresholds(h2,W-W/20,W-W/8); if (i!=h2) return 17;
  i = a_library.setObjectSizeThresholds(h1,30/(bin*bin),50/(bin*bin),10000/(bin*bin),15000/(bin*bin)); if (i!=h1) return 18;
  i = a_library.setObjectSizeThresholds(h2,30/(bin*bin),50/(bin*bin),10000/(bin*bin),15000/(bin*bin)); if (i!=h2) return 19;
  i = a_library.setObjectPersistenceThreshold(h1,8); if (i!=h1) return 20;
  i = a_library.setObjectPersistenceThreshold(h2,8); if (i!=h2) return 21;
  i = a_library.setAdaptationRate(h1,4); if (i!=h1) return 22;
  i = a_library.setAdaptationRate(h2,4); if (i!=h2) return 23;
  if (bin==1) { i = a_library.enableOutlining(h1,true); if (i!=h1) return 24; }
  else { i = a_library.enableSkeletonization(h1,true); if (i!=h1) return 24; }
  
  char s[1024];
  Mask m(128);
  for (int j=0;j<4;j++)
  {
    arena = (3*W)/4;
    worm1.imprint(arena,m,0.5,2);
    worm2.imprint(arena,m,0.5,2);
    worm3.imprint(arena,m,0.5,2);
    arena.set( Rectangle( Point(33,256+215) , Point(36,256+222) ) , W/3 );
    
    arena.bin = bin;
    i = a_library.scanObjects(h1,arena);
    if (j>2 && (i<2 || i>3)) return 25;
    i = a_library.scanRefs(h1,arena); if (i!=0) return 26;
    
    i = a_library.scanObjects(h2,arena); if (i!=0) return 27;
    i = a_library.scanRefs(h2,arena); if (i!=1) return 28;
    
    i = a_library.showObjects(h1,arena); if (i!=h1) return 29;
    i = a_library.showRefs(h1,arena); if (i!=h1) return 30;
    i = a_library.showObjects(h2,arena); if (i!=h2) return 31;
    i = a_library.showRefs(h2,arena); if (i!=h2) return 32;
    arena<<=6;
    sprintf(s,"test_lib_track_%02d.tiff",j);
    arena.writeTiff(s);
    arena.bin = 0;
    
    worm1.wiggle(0.4);
    worm2.wiggle(0.4);
    worm3.wiggle(0.4);
  }
  
  // Make sure error-checking is working; shouldn't be able to process images without loading them!
  i = a_library.processImage(h1); if (i!=0) return 33;
  i = a_library.processImage(h2); if (i!=0) return 34;
  
  i = a_library.beginOutput(h1); if (i!=h1) return 35;
  i = a_library.beginOutput(h2); if (i!=h2) return 36;
  i = a_library.setUpdateBandNumber(h1,12); if (i!=h1) return 37;
  i = a_library.setVelocityIntegrationTime(h1,1.2); if (i!=h1) return 38;
  
#ifndef PERFORMANCE_TEST
  Image temp_im( Point(300,280), false );
  temp_im.depth = arena.depth;
  temp_im = 0;
#endif
  
  for (int j=0;j<96;j++)
  {
    arena = (W*3)/4;
    worm1.imprint(arena,m,0.4,2);
    worm2.imprint(arena,m,0.4,2);
    worm3.imprint(arena,m,0.4,2);
    arena.set( Rectangle( Point(33,256+215) , Point(36,256+222) ) , W/3 );
    
    arena.bin = bin;
    if ((j%1)==0)
    {
      i = a_library.loadImage(h1,arena,0.4*j); if (i!=h1) return 39;  // Wholistic image loading
    }
    else
    {
      i = a_library.prepareImagePieces(h1,0.4*j); printf("%d\n",i);
      Rectangle r;
      while (a_library.getNextPieceCoords(h1,r)==h1)
      {
        i = a_library.loadThisImagePiece(h1,arena); if (i!=h1) return 40;
      }
    }
    i = a_library.loadImage(h2,arena,0.4*j); if (i!=h2) return 41;
    i = a_library.showLoaded(h1,arena); if (i!=h1) return 42;
    i = a_library.showLoaded(h2,arena); if (i!=h2) return 43;
    
    if (j==46)
    {
      i = a_library.markEvent(h2,1); if (i!=h2) return 44;
      i = a_library.markEvent(h1,2); if (i!=h1) return 45;
    }
    
    i = a_library.processImage(h1); if (i<2 || i>3) return 46;
    i = a_library.processImage(h2); if (i != 0) return 47;
    i = a_library.showResults(h1,arena); if (i!=h1) return 48;
    i = a_library.showResults(h2,arena); if (i!=h2) return 49;
    arena.bin = 0;

#ifndef PERFORMANCE_TEST    
    if ( (j%2) == 0)
    {
      printf("frame %02d: %d %d %.2f  %.2f %.2f  %.1f %.2f  %.1f %.2f  %.3f %.2f  %.2f\n",
        j+4,a_library.reportNumber(h1),a_library.reportNumberPersisting(h1),a_library.reportPersistence(h1),
        a_library.reportSpeed(h1),a_library.reportAngularSpeed(h1),
        a_library.reportLength(h1),a_library.reportRelativeLength(h1),
        a_library.reportWidth(h1),a_library.reportRelativeWidth(h1),
        a_library.reportAspect(h1),a_library.reportRelativeAspect(h1),
        a_library.reportEndWiggle(h1)
      );
    }
    
    arena.bin = bin;
    arena<<=6;
    sprintf(s,"test_lib_track_%02d.tiff",4+j);
    arena.writeTiff(s);
    arena.bin = 0;
    
    if (j>7 && j<10)
    {
      a_library.resizeRescale(h1,arena,temp_im,Rectangle(Point(0,0),Point(280,270)));
      sprintf(s,"test_%02d.tiff",4+j);
      temp_im.writeTiff(s);
    }
#endif

    worm1.wiggle(0.2);
    worm2.wiggle(0.2);
    worm3.wiggle(0.2);
  }

#ifdef PERFORMANCE_TEST
  arena = (W*3)/4;
  worm1.imprint(arena,m,0.4,2);
  worm2.imprint(arena,m,0.4,2);
  worm3.imprint(arena,m,0.4,2);
  Image newrena(arena,arena.getBounds());  
  for (i=0;i<1024*8;i++)
  {
    arena.copy( Point(0,0) , newrena , arena.size );
    arena.bin = bin;
    a_library.loadImage(h1,arena,0.4*96 + 0.1*i);
    a_library.processImage(h1);
    arena.bin = 0;
  }
#endif
  
  i = a_library.complete(h1); if (i!=h1) return 50;
  i = a_library.complete(h2); if (i!=h2) return 51;
  
  return 0;
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_library(1);
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Library test result is %d\n",i);
  return i>0;
}
#endif

