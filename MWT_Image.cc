/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include <time.h>
#include <string.h>
#include <stdio.h>

#ifdef ENABLE_SIMD
#include <smmintrin.h>
#endif

#include "MWT_Geometry.h"
#include "MWT_Lists.h"
#include "MWT_Storage.h"
#include "MWT_Image.h"



/****************************************************************
                      Image Mask Methods
****************************************************************/

// Make a bunch of strips that fill out a rectangle
Mask::Mask(const Rectangle& rect , Storage< Listable<Strip> >* store) : lines(store)
{
  rectangularize(rect);
}

// Make a bunch of strips that fill out an ellipse
Mask::Mask(const Ellipse& ellip , Storage< Listable<Strip> >* store) : lines(store)
{
  ellipticize(ellip);
}


// Set mask to a rectangular shape
void Mask::rectangularize(const Rectangle& rect)
{
  flush();
  for (int x = rect.near.x ; x <= rect.far.x ; x++) lines.Append( Strip(x,rect.near.y,rect.far.y) );
  bounds = rect;
  pixel_count = rect.area();
}

// Set mask to an elliptical shape
void Mask::ellipticize(const Ellipse& ellip)
{
  flush();
  Point p = Point(-ellip.axes.x,0);
  
  // First half of the ellipse plus the center line
  for ( ; p.x <= 0 ; p.x++ )
  {
    while (ellip.contains(p+ellip.center)) p.y++;
    p.y--;
    lines.Append( Strip( ellip.center.x+p.x , ellip.center.y-p.y , ellip.center.y+p.y ) );
    pixel_count += 2*p.y + 1;
  }
  
  // Second half is mirror image of first (skip center line)
  lines.end();
  lines.retreat();
  while (lines.current != NULL)
  {
    lines.Append( lines.i() + Point(2*(ellip.center.x - lines.i().x) , 0) );  // Reflect about center
    pixel_count += lines.i().length();
    lines.retreat();
  }
  bounds = ellip.bounds();
}


// Find the bounding box of a mask
void Mask::findBounds()
{
  Listable<Strip>* ls;
  
  lines.start(ls);
  if (lines.advance(ls)) { bounds = Rectangle(ls->data.x,ls->data.x,ls->data.y0,ls->data.y1); }
  while (lines.advance(ls))
  {
    bounds.includeX(ls->data.x);
    bounds.includeYY(ls->data.y0,ls->data.y1);
  }
}

// Clean up a mask that might have overlaps (assume sorted)
void Mask::tidy()
{
  Listable<Strip>* ls;
  
  pixel_count = 0;
  start();
  while (advance())
  {
    ls = lines.current;
    while (lines.advance(ls)) // Merge any later strips that overlap us or which touch in y
    {
      if (ls->data.x > i().x || ls->data.y0 > i().y1+1) break;
      if (i().y1 < ls->data.y1) i().y1 = ls->data.y1;   // Sort implies start is earlier, fix end if needed
      lines.Backspace(ls);
    }
    pixel_count += i().length();
  }
}


//Various point-inclusion functions; all assume that mask is sorted (points too, if there's a list)

// See whether a point is in the mask (slow, don't do this if you have many points)
bool Mask::contains(const Point& p)
{
  start();
  while (advance())
  {
    if (i().beyond(p)) return false;
    if (!i().before(p)) return true;
  }
  return false;
}

// Return the points contained in a mask leaving behind those that are not
Stackable<Point>* Mask::contains(Stackable<Point>*& sp)
{
  Stackable<Point> *temp,*in_h,*in_t,*out_h,*out_t;
  in_h = in_t = NULL;
  out_h = out_t = NULL;
  start();
  if (!advance()) return NULL;
  while (sp != NULL)
  {
    temp = Stackable<Point>::pop(sp);
    while (i().before(temp->data))
    {
      if (!advance())
      {
        if (out_t!=NULL) { out_t->next=temp; temp->next=sp; sp = out_h; return in_h; }
        else { temp->next=sp; sp = temp; return in_h; }
      }
    }
    if (i().beyond(temp->data)) temp->appendTo(out_h,out_t);
    else temp->appendTo(in_h,in_t);
  }
  sp = out_h;
  return in_h;
}

// Returns the points contained in a mask leaving behind those that are not
Listable<Point>* Mask::contains(Listable<Point>*& lp)
{
  Listable<Point> *temp,*ip,*in_h,*in_t;
  in_h = in_t = NULL;
  start();
  if (!advance()) return NULL;
  for ( ip=lp ; ip!=NULL ; )
  {
    while (i().before(ip->data)) { if (!advance()) return in_h; }
    if (i().beyond(ip->data)) ip=ip->next;
    else { temp=ip; ip = temp->extractR(); temp->appendTo(in_h,in_t); }
  }
  return in_h;
}

// Moves points that are contained in a mask from the first list to second and returns the number
// Assumes that the lists are clever enough to allocate their own memory
int Mask::contains(ManagedList<Point>& to_test,ManagedList<Point>& inside)
{
  int N = 0;
  start();
  if (!advance()) return 0;
  to_test.start();
  while (to_test.advance())
  {
    while (i().before(to_test.i())) { if (!advance()) return N; }
    if (!i().beyond(to_test.i()))
    {
      inside.Append(to_test.i());
      to_test.Backspace();
      N++;
    }
  }
  return N;
}

// End of various point-inclusion functions


// Fraction of two masks that overlap each other
FPoint Mask::overlapRatio(Mask& m)
{
  int N = 0;
  bool i_have_more,they_have_more;
  
  start();
  m.start();
  i_have_more = advance();
  they_have_more = m.advance();
  while (i_have_more && they_have_more)
  {
    if (i().before(m.i())) i_have_more = advance();
    else if (m.i().before(i())) they_have_more = m.advance();
    else
    {
      N += 1 + ((m.i().y1<i().y1) ? i().y1 : m.i().y1) - ((m.i().y0<i().y0) ? i().y0 : m.i().y0);
      if (i().y1 > m.i().y1) they_have_more = m.advance();
      else i_have_more = advance();
    }
  }
  
  if (N==0) return FPoint(0,0);
  else return FPoint( (float)N/(float)pixel_count , (float)N/(float)m.pixel_count );
}


// See if two masks overlap each other at all (assume they're sorted)
bool Mask::overlaps(Mask& m,Point offset)
{
  bool i_have_more,they_have_more;
  
  start();
  m.start();
  i_have_more = advance();
  they_have_more = m.advance();
  while (i_have_more && they_have_more)
  {
    if (i().x < m.i().x+offset.x) i_have_more = advance();
    else if (i().x > m.i().x+offset.x) they_have_more = m.advance();
    else if (i().y1 < m.i().y0+offset.y) i_have_more = advance();
    else if (i().y0 > m.i().y1+offset.y) they_have_more = m.advance();
    else return true;
  }
  return false;
}

// Make sure the mask fits within a bounding rectangle
void Mask::cropTo(const Rectangle& rect)
{
  start();
  while (advance())
  {
    if (i().x<rect.near.x) Backspace();
    else if (i().x>rect.far.x) Backspace();
    else if (i().y1<rect.near.y || i().y0>rect.far.y) Backspace();
    else
    {
      if (i().y0 < rect.near.y) { pixel_count -= rect.near.y - i().y0; i().y0 = rect.near.y; }
      if (i().y1 > rect.far.y) { pixel_count -= i().y1 - rect.far.y; i().y1 = rect.far.y; }
    }
  }
}


// Count the number of pixels in common between two masks.
// Strips must be sorted in both masks!
int Mask::countOverlap(Mask &m)
{
  int n_overlap = 0;
  bool i_have_more,they_have_more;
  
  start();
  m.start();
  i_have_more = advance();
  they_have_more = m.advance();
  while (i_have_more && they_have_more)
  {
    if (!i().overlaps(m.i()))
    {
      if (i() < m.i()) i_have_more = advance();
      else they_have_more = m.advance();
    }
    else
    {
      if (i().y1 < m.i().y1) // We end first
      {
        n_overlap += 1 + i().y1 - ((i().y0<m.i().y0) ? m.i().y0 : i().y0);
        i_have_more = advance();
      }
      else // They end first
      {
        n_overlap += 1 + m.i().y1 - ((i().y0<m.i().y0) ? m.i().y0 : i().y0);
        they_have_more = advance();
      }
    }
  }
  return n_overlap;
}


// Adds one mask to another
Mask& Mask::operator+=(Mask &m)
{
  bool i_have_more,they_have_more;
  
  start();
  m.start();
  i_have_more = advance();
  they_have_more = m.advance();
  while (i_have_more && they_have_more)
  {
    if (i().before(m.i())) i_have_more = advance();
    else if (i().beyond(m.i()))
    {
      pixel_count += m.i().length();
      lines.Tuck( m.i() );
      they_have_more = m.advance();
    }
    else
    {
      pixel_count -= i().length();
      if (i().y0 > m.i().y0) i().y0 = m.i().y0;
      if (i().y1 < m.i().y1) i().y1 = m.i().y1;
      pixel_count += i().length();
      they_have_more = m.advance();
    }
  }
  while (they_have_more)
  {
    lines.Append( m.i() );
    pixel_count += m.i().length();
    they_have_more = m.advance();
  }
  return *this;
}


// Removes one mask from another
Mask& Mask::operator-=(Mask &m)
{
  bool i_have_more,they_have_more;
  
  start();
  m.start();
  i_have_more = advance();
  they_have_more = m.advance();
  while (i_have_more && they_have_more)
  {
    while (they_have_more && m.i() < i())
    {
      if (m.i().x == i().x && m.i().y1 >= i().y0) break;
      they_have_more = m.advance();
    }
    if (m.i().x > i().x) i_have_more = advance();
    else if (m.i().y0 > i().y1) i_have_more = advance();
    else
    {
      if (m.i().y0 <= i().y0)
      {
        if (m.i().y1 >= i().y1)
        {
          pixel_count -= i().length();
          lines.Delete();
          if (lines.current==NULL) i_have_more = advance();
        }
        else
        {
          pixel_count -= 1+m.i().y1-i().y0;
          i().y0 = 1+m.i().y1;
        }
      }
      else if (m.i().y1 >= i().y1)
      {
        pixel_count -= 1+i().y1-m.i().y0;
        i().y1 = m.i().y0-1;
      }
      else
      {
        pixel_count -= m.i().length();
        lines.Insert( Strip(i().x , m.i().y1+1 , i().y1) );
        i().y1 = m.i().y0-1;
      }
    }
  }
  return *this;
}


// Removes a rectangle from a mask
Mask& Mask::operator-=(const Rectangle &rect)
{
  bool i_have_more;
  int x;
  
  start();
  i_have_more = advance();
  if (!i_have_more) return *this;
  
  for ( x = (i().x<rect.near.x) ? rect.near.x : i().x ; x<=rect.far.x && i_have_more ; x++ )
  {
    while ( (i().x < x || (i().x==x && i().y1 < rect.near.y)) && i_have_more) i_have_more = advance();
    while (i_have_more && i().x==x && !(i().y0 > rect.far.y))
    {
      if (rect.near.y <= i().y0)
      {
        if (rect.far.y >= i().y1)
        {
          pixel_count -= i().length();
          lines.Backspace();
        }
        else
        {
          pixel_count -= 1 + rect.far.y - i().y0;
          i().y0 = rect.far.y+1;
        }
      }
      else
      {
        if (rect.far.y >= i().y1)
        {
          pixel_count -= 1 + i().y1 - rect.near.y;
          i().y1 = rect.near.y-1;
        }
        else
        {
          pixel_count -= 1 + rect.far.y - rect.near.y;
          lines.Insert( Strip(x,rect.far.y+1,i().y1) );
          i().y1 = rect.near.y-1;
        }
      }
      i_have_more = advance();
    }
  }
  return *this;
}


// Finds the intersection of two masks
Mask& Mask::operator*=(Mask &m)
{
  bool i_have_more,they_have_more;
  
  start();
  m.start();
  i_have_more = advance();
  they_have_more = m.advance();
  while (i_have_more && they_have_more)
  {
    if (i().before(m.i()))
    {
      Backspace();
      i_have_more = advance();
    }
    else if (i().beyond(m.i())) they_have_more = m.advance();
    else
    {
      if (m.i().y0 > i().y0)
      {
        pixel_count -= m.i().y0-i().y0;
        i().y0 = m.i().y0;
      }
      if (m.i().y1 <= i().y1)
      {
        lines.Insert( Strip(i().x,m.i().y1+1,i().y1) );
        i().y1 = m.i().y1;
        they_have_more = m.advance();
      }
      i_have_more = advance();
    }
  }
  while (i_have_more)
  {
    pixel_count -= i().length();
    lines.Backspace();
    i_have_more = advance();
  }
  return *this;
}


// Morphological image operations on masks

// Copy an edge to a new mask
void Mask::extractEdge(Mask& edges)
{
  int y0,y1;
  Mask coverL( edges , ManagedList<Strip>::SubordinateList );
  Mask coverR( edges , ManagedList<Strip>::SubordinateList );
  Listable<Strip>* ls;
  edges.flush();
  
  start();
  while (advance())
  {
    // Special case--always add endpoints, so add whole strip if it's that short
    if (i().length() <= 2)
    {
      edges.addStrip(i());
      continue;
    }
    
    // Gather substrips that bound us to the left
    coverL.flush();
    ls = lines.current;
    while (lines.retreat(ls))
    {
      if (ii(ls).x < i().x-1) break;  // Too far
      if (ii(ls).x == i().x-1)
      {
        if (ii(ls).y1 < i().y0) break;  // Too far
        if (ii(ls).y0 > i().y1) continue;  // Not far enough
        y0 = i().y0;
        y1 = i().y1;
        if (ii(ls).y0 > y0) y0 = ii(ls).y0;
        if (ii(ls).y1 < y1) y1 = ii(ls).y1;
        coverL.pushStrip( i().x , y0 , y1 );
      }
      else continue;  // On same line--not far enough
    }
      
    // If completely bare, we add everything and are done
    if (coverL.lines.size==0)
    {
      edges.addStrip(i());
      continue;
    }
      
    // Gather substrips that bound us to the right
    coverR.flush();
    ls = lines.current;
    while (lines.advance(ls))
    {
      if (ii(ls).x > i().x+1) break;  // Too far
      if (ii(ls).x == i().x+1)
      {
        if (ii(ls).y0 > i().y1) break;  // Too far
        if (ii(ls).y1 < i().y0) continue;  // Not far enough
        y0 = i().y0;
        y1 = i().y1;
        if (ii(ls).y0 > y0) y0 = ii(ls).y0;
        if (ii(ls).y1 < y1) y1 = ii(ls).y1;
        coverR.addStrip( i().x , y0 , y1 );
      }
      else continue;  // Same line--not far enough
    }
    
    // If completely bare, again add everything
    if (coverR.lines.size==0)
    {
      edges.addStrip(i());
      continue;
    }
    
    // Find places with at least one side free and add those as an edge
    // coverL is reused as "bounded on both sides", coverR as "free on at least one side"
    coverL *= coverR;
    if (coverL.lines.size==0) edges.addStrip(i());  // Nothing sandwiched, add whole strip
    else
    {
      coverR.flush();
      coverR.addStrip(i());
      coverR -= coverL;
      if (coverR.lines.size==0)  // Everything sandwiched, only add ends
      {
        coverR.addStrip(i().x,i().y0,i().y0);
        coverR.addStrip(i().x,i().y1,i().y1);
      }
      else  // Add ends if we need to
      {
        if (coverR.lines.h().y0 > i().y0) // Need to cover small end
        {
          if (coverR.lines.h().y0 == i().y0+1)
          {
            coverR.lines.h().y0 = i().y0;   // Extend to end
            coverR.pixel_count++;
          }
          else  coverR.pushStrip( i().x , i().y0 , i().y0 ); // Add end
        }
        if (coverR.lines.t().y1 < i().y1) // Need to cover large end
        {
          if (coverR.lines.t().y1 == i().y1-1)
          {
            coverR.lines.t().y1 = i().y1;  // Extend to end
            coverR.pixel_count++;
          }
          else coverR.addStrip( i().x,i().y1,i().y1 );
        }
      }
      coverR.start();
      while (coverR.advance()) edges.addStrip(coverR.i());
    }
  }
}

// Invert a mask to a negative mask inside some rectangle
void Mask::invert(const Rectangle& universe)
{
  Listable<Strip> *ls;
  Point p;
  
  if (lines.size==0)
  {
    rectangularize(universe);
    return;
  }
  
  // Pick out where strips are not and place them at the beginning of the list (in order)
  ls = lines.list_head; 
  p.x = universe.near.x;
  p.y = universe.near.y;
  start();
  while (advance())
  {
    if (p.x > i().x) continue;
    while (p.x < i().x && p.x<=universe.far.x)
    {
      if (p.y <= universe.far.y) lines.Tuck(ls , Strip(p.x,p.y,universe.far.y));
      p.y = universe.near.y;
      p.x++;
    }
    if (p.x>universe.far.x) break;
    
    if (p.y < i().y0) lines.Tuck(ls , Strip(p.x,p.y,i().y0-1));
    p.y = i().y1+1;
  }
  
  // Add any strip-free region in the universe the last strip in the mask
  while (p.x <= universe.far.x)
  {
    if (p.y <= universe.far.y) lines.Tuck(ls , Strip(p.x,p.y,universe.far.y));
    p.y = universe.near.y;
    p.x++;
  }
  
  // Remove all the old strips
  lines.end();
  while (lines.current != ls) lines.Backspace();
  lines.Backspace();
  
  // Count the pixels in the new strips
  pixel_count = 0;
  start();
  while (advance()) pixel_count += i().length();
}
  

// Dilate an object using a square structuring element
void Mask::dilate(int n)
{
  if (n==0) return;
  if (n<0)
  {
    erode(-n);
    return;
  }
  
  n=n-1;
  int j;
  start();
  while (advance())
  {
    for (j = i().x-n ; j <= i().x+n ; j++) pushStrip( j , i().y0-n , i().y1+n );
    Backspace();
  }
  sort();
  tidy();
}

// Dilate an object using another mask as a structuring element
void Mask::dilate(Mask& structuringElt)
{
  start();
  while (advance())
  {
    structuringElt.start();
    while (structuringElt.advance())
    {
      pushStrip( i().x + structuringElt.i().x , i().y0 + structuringElt.i().y0 , i().y1 + structuringElt.i().y1 );
    }
  }
  sort();
  tidy();
}

// Erode an object using a square structuring element (erode by inverting, dilating, inverting)
void Mask::erode(int n)
{
  if (n==0) return;
  if (n<0)
  {
    dilate(-n);
    return;
  }
  
  // Find a safe box within which to do our inversion
  findBounds();
  Rectangle r = bounds;
  r.expand(n);
  invert(r);
  dilate(n);
  invert(r);
}

// Erode an object using another mask as a structuring element (invert, dilate, invert)
void Mask::erode(Mask& structuringElt)
{
  // Not sure if bounds are set--be on the safe side (slows it down, could optimize?)
  findBounds();
  structuringElt.findBounds();
  Rectangle r;
  
  // Find safe bounds for edge of structuring element
  Point near_delta,far_delta;
  near_delta = Point(-1,-1);
  if (structuringElt.bounds.near.x < 0) near_delta.x += structuringElt.bounds.near.x;
  if (structuringElt.bounds.near.y < 0) near_delta.y += structuringElt.bounds.near.y;
  far_delta = Point(1,1);
  if (structuringElt.bounds.far.x > 0) far_delta.x += structuringElt.bounds.far.x;
  if (structuringElt.bounds.far.y > 0) far_delta.y += structuringElt.bounds.far.y;
  
  r = bounds;
  r.near += near_delta; r.near += near_delta;  // Need a double-width safety margin in case structuring elt performs a shift
  r.far += far_delta; r.far += far_delta;
  invert(r);
  
  dilate(structuringElt);
  
  r.near -= near_delta;  // Shrink margin so structuring element can't have moved junk in from outside universe
  r.far -= far_delta;
  invert(r);
}
    
// End of morphological operations on masks


// Provide some text output so we can see mask for debugging purposes
void Mask::println()
{
  start();
  while (advance())
  {
    printf("%d  %d %d\n",i().x,i().y0,i().y1);
  }
  printf("  N=%d\n",pixel_count);
}



/****************************************************************
                     Image Contour Methods
****************************************************************/

// Build a contour from a rectangle--either just the four corners or the whole border
void Contour::findContour(const Rectangle& r,bool full)
{
  flush();
  if (full==false)
  {
    boundary.Append( r.near );
    boundary.Append( Point(r.near.x,r.far.y) );
    boundary.Append( r.far );
    boundary.Append( Point(r.far.y,r.near.x) );
  }
  else
  {
    Point p,q;
    p = q = r.firstBoundaryPoint();
    do
    {
      boundary.Append( p );
      p = r.nextBoundaryPoint(p);
    } while (p != q);
  }
}

// Build a contour from an ellipse.  Pixels will be adjacent; use emptyContour to sparsen if desired.
void Contour::findContour(const Ellipse& e)
{
  flush();
  Point p,q;
  p = q = e.firstBoundaryPoint();
  do
  {
    Append(p);
    p = e.nextBoundaryPoint(p);
  } while (p != q);
}

// Build a contour from a mask; assume that the mask is a single object.  Only the outer contour is grabbed.
// Contour is not formed from adjacent pixels.  Use fillContour to fill in the gaps.
// The contour is built starting at the left-most corner and traveling counterclockwise.  Thus, the mask
// must always be on the left-hand side of someone walking along the contour.
void Contour::findContour(Mask& m)
{
  int x_last; // The previous x coordinate of the contour
  int y;      // The current y coordinate of the contour (current x is always the x of the current strip)
  int i;      // Temporary variable used to test against ends of strips
  int dir;    // Direction along strips to travel to find more contour (1=travel in +y direction, -1 = -y)
  Listable<Strip>* ls;
  Listable<Strip>* next_strip;
  
  flush();
  if (m.pixel_count==0) return;
  
  m.start();
  m.advance();
  dir = 1;
  x_last = m.i().x;
  y = m.i().y0;
  do
  {
    // Add the current pixel unless it's already the last thing in the contour.
    // (Check size() to avoid null deref. on t().)
    if (size()==0 || t().x != m.i().x || t().y != y) Append( Point(m.i().x,y) );
    
    if (dir>0)
    {
      // Where should we look for the next strip?
      if (x_last<m.i().x) i = y+1;  // Next strip has to be forward a bit or it'd be part of the strip we were just in
      else i = y-1;                 // Next strip might be "behind" us.
      
      // Search for the next strip that might be to our right--in +y direction, this is to the -x direction
      next_strip = NULL;
      ls = m.lines.current;
      while (m.lines.retreat(ls))
      {
        if (ls->data.x < m.i().x-1) break;  // Too far, not even touching, irrelevant
        if (ls->data.x == m.i().x) continue;  // Same x strip as we're on, not far enough
        if (ls->data.y1<i || ls->data.y0 > m.i().y1+1) continue;  // Doesn't touch necessary part of existing strip
        next_strip = ls;
      }
      if (next_strip==NULL)  // Nobody is on our right--double back!
      {
        y = m.i().y1;
        dir = -1;
      }
      else  // Someone is on our right!  Let's jump.
      {
        if (next_strip->data.y0<y)  // Don't travel down our strip, just go up the other one.
        {
          y--;
          dir = -1;
        }
        else
        {
          if (next_strip->data.y0 > y+1) Append( Point(m.i().x,next_strip->data.y0-1) );  // We got this far in our strip
          y = next_strip->data.y0;  // Now move to the new one
        }
      }
    }
    else // Should be exactly the same as the above with all the directions flipped (see comments above, but flip left/right, up/down).
    {
      if (x_last>m.i().x) i = y-1;
      else i = y+1;
      
      next_strip = NULL;
      ls = m.lines.current;
      while (m.lines.advance(ls))
      {
        if (ls->data.x > m.i().x+1) break;
        if (ls->data.x == m.i().x) continue;
        if (ls->data.y0>i || ls->data.y1 < m.i().y0-1) continue;
        next_strip = ls;
      }
      if (next_strip==NULL)
      {
        y = m.i().y0;
        dir = 1;
      }
      else
      {
        if (next_strip->data.y1>y)
        {
          y++;
          dir = 1;
        }
        else
        {
          if (next_strip->data.y1 < y-1) Append( Point(m.i().x,next_strip->data.y1+1) );
          y = next_strip->data.y1;
        }
      } 
    }
    x_last = m.i().x;
    if (next_strip!=NULL) m.lines.current = next_strip;
  } while (dir!=1 || m.i()!=m.lines.h());
}


// Add pixels to a contour until pixels are adjacent (that is, contour has no gaps)
void Contour::fillContour()
{
  Listable<Point>* old_lp;
  Listable<Point>* new_lp;
  
  Point offset;
  Point p,q;
  
  if (size()<2) return;
  
  first();
  old_lp = getIterator();
  next();
  new_lp = getIterator();
  do
  {
    offset = new_lp->data - old_lp->data;  // Location of new pixel relative to old one
    if (offset*offset > 2)  // Distance is greater that sqrt(2), so we need to fill in pixels along the way
    {
      setIterator( old_lp );
      if (offset.x==0 || offset.y==0 || offset.x*offset.x == offset.y*offset.y)  // Along diagonal or x or y axis
      {
        p.x = (offset.x<0) ? -1 : (offset.x>0) ? 1 : 0;  // Travel this direction in x to reach next contour point
        p.y = (offset.y<0) ? -1 : (offset.y>0) ? 1 : 0;  // And this in y
        do
        {
          Insert( i()+p );  // Drop pixels as we go
          next();
          q = new_lp->data - i();  // Update how far is left to go
        } while (q*q > 2);
      }
      else // Need to find best set of pixels to match line
      {
        Point px;  // If we step in just the x-direction, go this way (x coord is +-1, other is 0)
        Point pxy; // If we step along a diagonal, this is the one to take
        Point py;  // If we step just in y, this is the right direction.
        FPoint o_hat,pf;
        float fx,fxy,fy;
        pxy.x = (offset.x<0) ? -1 : (offset.x>0) ? 1 : 0;
        pxy.y = (offset.y<0) ? -1 : (offset.y>0) ? 1 : 0;
        px.x = pxy.x;
        px.y = 0;
        py.x = 0;
        py.y = pxy.y;
        o_hat = FPoint(offset).unit();  // This vector points in the direction we wanted to go.
        
        // Now we walk along, trying to make sure that we keep pointing in the o_hat (i.e. original) direction
        do
        {
          pf = FPoint(i()+pxy).unit();
          fxy = pf*o_hat;  // Angle between original direction and new direction if we make an XY step
          pf = FPoint(i()+px).unit();
          fx = pf*o_hat;   // Angle between original and new if X step only
          if (fx<fxy) p = px;   // If X step is better, take it (if X is better than XY, then Y will be worse still)
          else
          {
            pf = FPoint(i()+py).unit();
            fy = pf*o_hat;  // Angle between original and new if Y step only
            if (fy<fxy) p = py;  // If Y step is better, take it
            else p = pxy;        // XY is better than either
          }
          Insert(i()+p);  // Drop in new pixel
          next();
          q = new_lp->data - i();  // Update how far is left to go
        } while (q*q > 2);
      }
    }
    
    // Jump to the next pair of points which might not be filled in
    setIterator( new_lp );
    old_lp = new_lp;
    next();
    new_lp = getIterator();
  } while (old_lp != boundary.list_head);
}

// Remove pixels until we leave only endpoints of lines
void Contour::emptyContour()
{
  Listable<Point>* old_lp; // Starting point
  Listable<Point>* the_lp; // Intermediate point--see if we can throw this away
  Listable<Point>* new_lp; // Ending point--advance this as far as we can
  
  Point offset;
  Point p,q;
  int i;
  bool advanced_old;
  
  if (size()<3) return;
  
  first();
  old_lp = getIterator();
  next();
  the_lp = getIterator();
  next();
  new_lp = getIterator();
  
  do
  {
    p = new_lp->data - the_lp->data;  // Direction from intermediate to end
    q = the_lp->data - old_lp->data;  // Direction from start to intermediate
    i = p*q;  // i = |p|*|q|*cos(angle between vectors)
    if (i*i == (p*p)*(q*q) && i>0)  // Same direction if cos(angle)==1, cheaper to square vectors than to find their length so check square
    {
      setIterator(the_lp);  // Don't need this guy
      Delete();
      setIterator(new_lp);  // Keep going from this one
      advanced_old = false;
    }
    else // Changed direction, have to leave intermediate pixel in place
    {
      old_lp = the_lp;
      advanced_old = true;
    }
    the_lp = new_lp;
    next();
    new_lp = getIterator();
  } while (!(advanced_old && old_lp == boundary.list_head));
}


// Generate a pretty crude spine down the center of a contour (not bad for worm-shaped things)
void Contour::approximateSpine(FPoint centroid,FPoint head_dir,int n_backbone_points)
{
  if (size()<3) return;
  
  FPoint p,q,r;
  int j,k;
  int N = size();
  float x[N];
  float y[N];
  float angle[N];
 
  // Pack points into local array (and convert to floats, and make centroid be 0,0)
  p = centroid;
  k = 0;
  start();
  while (advance())
  {
    x[k] = i().x - p.x;
    y[k] = i().y - p.y;
    k++;
  }
  
  // Find angles between points spaced a reasonable distance apart around the object
  int spacing = (N/(2*(n_backbone_points-1)));
  if (spacing<5) spacing = 5;
  for (k=0;k<N;k++)
  {
    p = FPoint(x[k],y[k]);
    j = (k+spacing); if (j>=N) j -= N;  // j %= N
    q = FPoint(x[j],y[j]);
    j = (k-spacing); if (j<0) j += N;  // j %= N
    r = FPoint(x[j],y[j]);
    q = (q-p).unit();  // Unit vector pointing towards later point on contour
    r = (r-p).unit();  // Unit vector pointing towards earlier point on contour
    
    // Want a measure that will be maximum when two vectors are close on an interior angle and minimum when they're close on exterior angle
    // (Assume contour was generated in a clockwise fashion.)
    // Could use signed angle (> pi means angle is convex viewed from inside), but atan2 is expensive.
    angle[k] = 1.0 + q*r;  // Dot product = cosine of angle, in range from 0 (antiparallel) to 2 (parallel)
    if (q.X(r) > 0) angle[k] = -angle[k];  // Flip if cross product is positive (i.e. exterior angle is smaller than interior angle)
  }
  
  // Find sharpest interior angle: probably a head or tail or other interesting feature
  float first_angle = angle[0];
  int first_index = 0;
  for (k=1;k<N;k++)
  {
    if (angle[k] > first_angle) { first_index=k; first_angle=angle[k]; }
  }
  
  // Decrease the weight of nearby interior angles, and look for another sharp interior angle (if we found tail, this would be head)
  angle[first_index] -= 4;
  for (k=1;k<N/2;k++)
  {
    j = (first_index+k); if (j>=N) j -= N;  // If is probably cheaper than mod N
    angle[j] *= (float)k/(float)((N+1)/2);
    j = (first_index-k); if (j<0) j += N;  // If is probably cheaper than mod N
    angle[j] *= (float)k/(float)((N+1)/2);
  }
  float second_angle = angle[0];
  int second_index = 0;
  for (k=1;k<N;k++)
  {
    if (angle[k] > second_angle) { second_index=k; second_angle=angle[k]; }
  }
  
  // Split worm into n_backbone_points-1 roughly equal segments by marking off border points 1/(n_backbone_points-1)
  // of the way from one sharp end to another
  int N_up = second_index-first_index;
  if (N_up>=N) N_up -= N;
  else if (N_up<0) N_up += N;
  int N_down = N - N_up;
  int up_index[n_backbone_points-2];
  int down_index[n_backbone_points-2];
  for (k=0;k<n_backbone_points-2;k++)
  {
    up_index[k] = first_index + (N_up*(k+1))/(n_backbone_points-1);
    if (up_index[k] >= N) up_index[k] -= N;
    down_index[k] = first_index - (N_down*(k+1))/(n_backbone_points-1);
    if (down_index[k] < 0) down_index[k] += N;
  }
  
  // Now create n_backbone_points skeleton points that are the two sharp ends and the average of the contour points we found on each side
  // Build forwards or backwards as needed so that the first element is in the "front" as determined by the direction of the head_dir vector
  bool head_first = ( head_dir * FPoint( x[first_index],y[first_index] ) >= head_dir * FPoint( x[second_index],y[second_index] ) );
  flush();
  p = FPoint( x[first_index] , y[first_index] );
  if (head_first) Append( p.toPoint() );
  else Push( p.toPoint() );
  for (k=0;k<n_backbone_points-2;k++)
  {
    p = 0.5 * (FPoint( x[ up_index[k] ] , y[ up_index[k] ] ) + FPoint( x[ down_index[k] ] , y[ down_index[k] ] ));
    if (head_first) Append( p.toPoint() );
    else Push( p.toPoint() );
  }
  p = FPoint( x[second_index] , y[second_index] );
  if (head_first) Append( p.toPoint() );
  else Push( p.toPoint() );  
}

// Generate a decent contour approximation via cubic segments (not spline, but center segment of 4 points)
void Contour::approximateContour(float max_err)
{
  if (max_err<=0) return;
  
  // Set up graph of x and y values of contour, plus indices to cubic fit points ("control" points, though it's not really a spline)
  // Also, wickedly repurpose the contour points so that x = index of point!
  int j,k,n,N;
  int x_coords[ size() ];
  int y_coords[ size() ];
  Listable<Scorer> *score_head,*score_tail;
#ifdef CLANG_WORKAROUND
  unsigned char clang_workaround_buffer[ size()*sizeof(Listable<Scorer>) ];
  Listable<Scorer> *scores = (Listable<Scorer>*) clang_workaround_buffer;
  for (int cwbi = 0; cwbi < size(); cwbi++) new (scores + cwbi) Listable<Scorer>();
#else
  Listable<Scorer> scores[ size() ];
#endif
  
  N = 0;
  start();
  while (advance())
  {
    x_coords[ N ] = i().x;
    y_coords[ N ] = i().y;
    i().x = N;
    i().y = N;
    N++;
  }
  
  // Values that we care about for the fit
  // We consider a cubic polynomial going through four points, and use it as a fit to a curve only between the central two points.
  // Thus, we are fitting an area between "L" and "R" points, and LL and RR are the points outside that.
  int xLL , xL , xR , xRR;      // X values
  int yLL , yL , yR , yRR;      // Y values
  int iLL , iL , iR , iRR;      // Index values
  int d_L2R;        // Distance between iL and iR
  float f_L2R;      // Inverse of distance between iL and iR
  int dL;           // Distance to left index value
  float t;          // Parameterized time (t=0 => at left index value, t=1 => at right (-1,+2 would be LL and RR))
  
  // Loop through, finding points that we can remove and removing them
  float value,error;
  bool keep_going = true;
  while (keep_going)
  {
    keep_going = false;
    
    if (size() <= 4) break;
    
    // Load up our scoring list
    n = 0;
    start();
    score_head = score_tail = NULL;
    while (advance())
    {
      i().y = i().x;  // We'll change i().y to mark a point that should be deleted
      scores[n].pushOnto(score_head,score_tail);
      scores[n].data.pt = getIterator();
      n++;
    }
    
    // Load the scores into each element of the list
    for (k=0;k<n;k++)
    {
      // Find points for fit (kinda messy because we have to worry about wrapping)
      iL = k-1;
      if (iL<0) { iL=n-1; iLL=n-2; }
      else { iLL=k-2; if (iLL<0) iLL=n-1; }
      iR = k+1;
      if (iR>=n) { iR=0; iRR=1; }
      else { iRR=k+2; if (iRR>=n) iRR=0; }
      
      // Got the raw indices, now load them using our reference tables
      iLL = scores[iLL].data.pt->data.x;
      iL = scores[iL].data.pt->data.x;
      iR = scores[iR].data.pt->data.x;
      iRR = scores[iRR].data.pt->data.x;
      xLL = x_coords[iLL]; yLL = y_coords[iLL];
      xL = x_coords[iL]; yL = y_coords[iL];
      xR = x_coords[iR]; yR = y_coords[iR];
      xRR = x_coords[iRR]; yRR = y_coords[iRR];
      d_L2R = iR-iL; if (d_L2R<0) d_L2R += N;
      f_L2R = 1.0/d_L2R;
      
      // Find maximum error in fit
      // Recall that a polynomial that goes through x0,x1,x2,x3 at times t0,t1,t2,t3 is given by
      // x0*(t-t1)*(t-t2)*(t-t3)/((t0-t1)*(t0-t2)*(t0-t3)) + (same thing for x1,x2,x3)
      // (just try plugging in t=t0 and note that everything cancels leaving x0*1).
      // We're using regularly spaced values of t, so we can precompute the (t0-t1)*... stuff. 
      error = 0;
      for (dL=1 ; dL<d_L2R ; dL++)
      {
        t = dL * f_L2R;  // Parameterize fitting region with a time between 0 and 1
        j = iL+dL;
        if (j>=N) j=0;
        
        // Make a prediction for the x coordinate
        value = ( -xLL*(t-2.0f) + xRR*(1.0f+t) )*t*(t-1.0f)*0.16666667f + ( xL*(t-1.0f) - xR*t )*(1.0f+t)*(t-2.0f)*0.5f;
        value -= x_coords[j];  // Subtract what we should have gotten--this is the error
        if (value<0) value = -value;  // Errors are always positive.
        if (value>error) error = value;  // We want to control maximum error, so check that
        
        // Same thing for the y coordinate
        value = ( -yLL*(t-2.0f) + yRR*(1.0f+t) )*t*(t-1.0f)*0.16666667f + ( yL*(t-1.0f) - yR*t )*(1.0f+t)*(t-2.0f)*0.5f;
        value -= y_coords[j];
        if (value<0) value = -value;
        if (value>error) error = value;
      }
      scores[k].data.score = error;
      if (error <= max_err) keep_going = true;
    }
    
    if (!keep_going) break;  // Didn't find anything to remove
    keep_going = false;
    
    // Sort scores and mark those to remove by dropping a -1 in the y coordinate
    // (only drop those that don't conflict with neighbors)
    Listable<Scorer>::mergeSort(score_head,score_tail);
    for ( ; score_head!=NULL && score_head->data.score <= max_err ; score_head = score_head->next)
    {
      setIterator( score_head->data.pt );
      next();
      if (i().y != -1)  // Previous one isn't removed
      {
        next();
        if (i().y != -1)  // Previous two aren't (LL and L will exist)
        {
          setIterator( score_head->data.pt );
          prev();
          if (i().y != -1)  // Next one isn't removed
          {
            prev();
            if (i().y != -1)  // Next two aren't (R and RR will exist)
            {
              // We have the surrounding points to make the approximation--remove this one.
              setIterator( score_head->data.pt );
              i().y = -1;
              keep_going = true;  // We removed some stuff, so we need to recompute and try to remove more next pass
            }
          }
        }
      }
    }
    
    if (!keep_going) break;  // Nothing changed, done.
    keep_going = false;
    
    // Now actually remove the marked items
    start();
    while (advance())
    {
      if (i().y==-1)
      {
        Backspace();
        keep_going = true;
      }
    }
  }
  
  // Put the remaining contour points back as x,y values
  start();
  while (advance())
  {
    k = i().x;
    i().x = x_coords[k];
    i().y = y_coords[k];
  }
}


// Find a bounding box
void Contour::findBounds()
{
  start();
  if (!advance()) return;
  
  bounds.near = bounds.far = i();
  
  while (advance()) bounds.include( i() );
}
 


/****************************************************************
                   Image Statistics Methods
****************************************************************/

// Typical version of the constructor
FloodInfo::FloodInfo(Image* image , Storage< Stackable<Strip> >* storage) : im(image),cx(0),cy(0),sumI(0)
{
  // Make bounds local but binned
  bounds = im->getBounds();
  bounds.far -= bounds.near;
  bounds.near = Point(0,0);
  pixel_count = 0;
  
  // Create storage if need be
  if (storage==NULL)
  {
    local_storage = true;
    store = new Storage< Stackable<Strip> >(im->size.x+im->size.y);
  }
  else
  {
    local_storage = false;
    store = storage;
  }
}


// Floods a pixel (makes it negative and counts it)
void FloodInfo::floodPixel(int x,int y,short I)
{
  short J = T.hi()-I;
  if (J > I-T.lo()) J = I-T.lo();
  im->set(x,y,-I);
  cx += J*x;
  cy += J*y;
  sumI += J;
  pixel_count++;
}


// Reset to start using again
void FloodInfo::freshen()
{
  cx = 0;
  cy = 0;
  sumI = 0;
  pixel_count = 0;
}


// Find centroid
void FloodData::findCentroid(const Image& im,Range threshold)
{
  short I,J;
  int y;
  long long cI,cx,cy;
  cI = cx = cy = 0;
  
  stencil.start();
  while (stencil.advance())
  {
    for (y=stencil.i().y0 ; y<=stencil.i().y1 ; y++)
    {
      I = im.get(stencil.i().x,y);
      if (I<0) I=-I;
      J = threshold.hi()-I;
      I -= threshold.lo();
      if (J<I) I=J;
      cI += I;
      cx += stencil.i().x*I;
      cy += y*I;
    }
  }
  
  centroid = FPoint( ((double)cx)/((double)cI) , ((double)cy)/((double)cI) );
}


// Find major and minor axes of an object
void FloodData::principalAxes(const Image& im,Range threshold)
{
  short I,J;
  int y;
  double a,b,f,L1,L2;
  double Sxx,Sxy,Syy,Si;
  
  // Least squares fit (standard formula)
  Si = 0.0;
  Sxx = Sxy = Syy = 0.0;
  stencil.start();
  while (stencil.advance())
  {
    a = stencil.i().x - centroid.x;
    for (y=stencil.i().y0 ; y<=stencil.i().y1 ; y++)
    {
      I = im.get(stencil.i().x,y);
      if (I<0) I=-I;
      J = threshold.hi()-I;
      I -= threshold.lo();
      if (J<I) I=J;
      b = y - centroid.y;
      Si += I;
      Sxx += I*a*a;
      Sxy += I*a*b;
      Syy += I*b*b;
    }
  }
  Sxx /= Si;
  Sxy /= Si;
  Syy /= Si;
  f = sqrt( (Sxx-Syy)*(Sxx-Syy)+4*Sxy*Sxy )/2;
  a = 0.5*(Sxx+Syy);
  b = 0.5*(Sxx-Syy);
  L1 = a + f;
  L2 = a - f;
  // If major axis poorly defined for the shape, arbitrarily choose X direction.
  if (fabs(f-b) < 1e-6) major = FPoint( 1.0 , 0.0 );  
  else
  {
    major = FPoint(Sxy/(f-b) , 1.0);
    if (fabs(major.x) < 1e-6) major.x = 0.0;
    major /= major.length();
  }
  // Likewise for minor axis
  if (fabs(b+f) < 1e-7) minor = FPoint( 1.0 , 0.0 );
  else
  {    
    minor = FPoint(Sxy/(b+f) , -1.0);
    if (fabs(minor.x) < 1e-6) minor.x = 0.0;
    minor /= minor.length();
  }
  
  // Find length along the least squares axes
  FPoint p;
  double maj_max,maj_min,min_max,min_min;
  maj_max = maj_min = min_max = min_min = 0.0;
  stencil.start();
  while (stencil.advance())
  {
    p.x = stencil.i().x - centroid.x;
    p.y = stencil.i().y0 - centroid.y;
    
    a = p*major;
    b = p*minor;
    if (a > maj_max) maj_max = a;
    else if (a < maj_min) maj_min = a;
    if (b > min_max) min_max = b;
    else if (b < min_min) min_min = b;
  }
  
  // Set data: axes (of standard deviation length) and length of object along two axes
  major *= sqrt(L1);
  minor *= sqrt(L2);
  long_axis = maj_max - maj_min;
  short_axis = min_max - min_min;
}


// Make a copy of a flood data object (only really needed to keep storage straight)
void FloodData::duplicate(const FloodData& fd)
{
  score = fd.score;
  centroid = fd.centroid;
  major = fd.major;
  minor = fd.minor;                         
  stencil.imitate(fd.stencil);
  long_axis = fd.long_axis;
  short_axis = fd.short_axis;
}



/****************************************************************
                         Image Methods
****************************************************************/

// Debugging methods to defer 8-bit input to a simple copy to
// scratch space plus a call to the corresponding 16 bit method
//#define DEFER_MIMIC
//#define DEFER_COPY
//#define DEFER_ADAPT
//#define DEFER_DIFFC
//#define DEFER_DAC


// Make a new image that copies a rectangular region of another
// If the old image is binned, it gets a copy of the binned (not original) version
Image::Image(const Image& existing,const Rectangle& region, bool algo) : bin(0),depth(existing.depth),divide_bg(algo)
{
  // WARNING!  You must keep this up to date with the 8 bit implementation BY HAND!!
  bounds = existing.getBounds();
  bounds.cropTo(region);
  if (bounds.isEmpty())
  {
    pixels=NULL;
    size=Point(0,0);
    owns_pixels=false;
    return;
  }
  
  size = (bounds.far-bounds.near)+1;
  pixels = new short[ bounds.area() ];
  owns_pixels = true;
  
  copy(bounds.near , existing , size);
}

Image* Image::identical() const {
  Image* im = new Image();
  im->pixels = new short[ bounds.area() ];
  im->bounds = bounds;
  im->size = size;
  im->bin = bin;
  im->depth = depth;
  im->owns_pixels = true;
  im->divide_bg = divide_bg;
  memcpy(im->pixels, pixels, bounds.area()*sizeof(short));
  return im;
}


// Draw rectangles on an image
void Image::set(Rectangle r,short I)
{
  // WARNING!  You must keep this up to date with the 8 bit implementation BY HAND!!
  Point p;
  
  if (bin <= 1)
  {
    if (!r.overlaps(bounds)) return;
    r.cropTo(bounds);
    r -= bounds.near;
  }
  else
  {
    Rectangle s = getBounds();
    if (!r.overlaps(s)) return;
    r.cropTo(s);
    r -= s.near;
    r *= bin;
    r.far += bin-1;
  }
  for (p.x=r.near.x ; p.x<=r.far.x ; p.x++) for (p.y=r.near.y ; p.y<=r.far.y ; p.y++) raw(p) = I;
}

// Draw a mask onto an image
void Image::set(Mask& m,short I)
{
  // WARNING!  You must keep this up to date with the 8 bit implementation BY HAND!!
  Point p;
  Rectangle b = getBounds();
  int x0,x1,y0,y1;
  m.start();
  while (m.advance())
  {
    p.x = m.i().x;
    if (p.x < b.near.x) continue;
    if (p.x > b.far.x) break;
    if (m.i().y1 < b.near.y || m.i().y0 > b.far.y) continue;
    y0 = m.i().y0;
    if (y0 < b.near.y) y0 = b.near.y;
    y1 = m.i().y1;
    if (y1 > b.far.y) y1 = b.far.y;
    if (bin>1)
    {
      x0 = p.x*bin;
      x1 = x0 + bin - 1;
      y0 *= bin;
      y1 *= bin;
      y1 += bin-1;
      x0 -= bounds.near.x;
      x1 -= bounds.near.x;
      y0 -= bounds.near.y;
      y1 -= bounds.near.y;
      for (p.x=x0 ; p.x<=x1 ; p.x++) for (p.y=y0 ; p.y<=y1 ; p.y++) raw(p) = I;
    }
    else for (p.y = y0 ; p.y <= y1 ; p.y++) rare(p) = I;  
  }
}


// Draw lines on an image--not efficient for single pixel width, but works
void Image::set(Point p1,Point p2,int width,short I)
{
  // WARNING!  You must keep this up to date with the 8 bit implementation BY HAND!!
  Rectangle r;
  
  if (width<1) return;
  else width--;
  
  if (p1.x==p2.x || p1.y==p2.y) // Single pixel or single horizontal or vertical line; easy.
  {
    if (p1.x==p2.x) r.near.x = r.far.x = p1.x;
    else if (p1.x<p2.x) { r.near.x = p1.x; r.far.x = p2.x-1; }
    else { r.near.x = p2.x+1; r.far.x = p1.x; }
    if (p1.y==p2.y) r.near.y = r.far.y = p1.y;
    else if (p1.y<p2.y) { r.near.y = p1.y; r.far.y = p2.y-1; }
    else { r.near.y = p2.y+1; r.far.y = p1.y; }
    r.near -= width;
    r.far += width;
    set(r,I);
  }
  else  // Diagonal, hard--draw a little rectangle of pixels around each point
  {
    r.near = p1;
    r.far = p1;
    r.expand(width);
    set(r,I);
    
    Point p,px,pxy,py,q;
    FPoint o_hat = FPoint(p2-p1).unit();  // Direction we want to go, dropping pixels as we do it
    FPoint pf;
    float fx,fxy,fy;
    
    // Figure out where we should go if we move only in X, only in Y, or in both X and Y (diagonal) this step
    pxy.x = (o_hat.x<0) ? -1 : (o_hat.x>0) ? 1 : 0;
    pxy.y = (o_hat.y<0) ? -1 : (o_hat.y>0) ? 1 : 0;
    px.x = pxy.x;
    px.y = 0;
    py.x = 0;
    py.y = pxy.y;
    
    q = p2-p1;
    while (q*q > 2)  // Farther than a diagonal to the last pixel--need to drop more pixels along the way
    {
      // Figure out if we have a better bearing after moving in X or in XY
      pf = FPoint(p1+pxy).unit();
      fxy = pf*o_hat;
      pf = FPoint(p1+px).unit();
      fx = pf*o_hat;
      pf = FPoint(p1+py).unit();
      if (fx > fxy)  // X was better, move in x and drop line of pixels (rest of block is already filled)
      {
        r.near.x = r.far.x = p1.x+px.x*(width+1);
        r.near.y = p1.y - width;
        r.far.y = p1.y + width;
        set(r,I);
        p1 += px;
      }
      else  // XY was better
      {
        // See if Y is better yet
        pf = FPoint(p1+py).unit();
        fy = pf*o_hat;
        if (fy > fxy) // Yes, move in Y and drop line of pixels
        {
          r.near.y = r.far.y = p1.y+py.y*(width+1);
          r.near.x = p1.x - width;
          r.far.x = p1.x + width;
          set(r,I);
          p1 += py;
        }
        else // XY was the way to go after all
        {
          r.near.x = r.far.x = p1.x+pxy.x*(width+1);
          r.near.y = p1.y - width + pxy.y;
          r.far.y = p1.y + width + pxy.y;
          set(r,I);
          if (width>0) // Can't just drop a line, so drop a full block of pixels
          {
            r.near.y = r.far.y = p1.y+pxy.y*(width+1);
            r.near.x = p1.x - width + pxy.x;
            r.far.x = p1.x + width + pxy.x - 1;
            set(r,I);
          }
          p1 += pxy;
        }
      }
      q = p2-p1;
    }
  }
}

// Draw a contour as a bunch of line segments
void Image::set(Contour& c,int width,short I)
{
  // WARNING!  You must keep this up to date with the 8 bit implementation BY HAND!!
  Point p0,p1,p2;
  c.first();
  p0 = p1 = c.i();
  if (!c.next()) { set(p1,p1,width,I); return; }
  p2 = c.i();
  do
  {
    set(p1,p2,width,I);
    p1 = p2;
    c.next();
    p2 = c.i();
  } while (p1!=p0);
}


// Copy a source image with resizing/rescaling so both images show the same scene
void Image::mimic(const Image& source,Rectangle my_region,Rectangle source_region,ScaleType method)
{
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->8 bits)
  // Convert region to point directly to image pixels
  if (bin > 1) my_region *= bin;
  if (source.bin > 1) source_region *= source.bin;
  if (!bounds.overlaps(my_region)) return;
  if (!source.bounds.overlaps(source_region)) return;
  my_region *= bounds;
  source_region *= source.bounds;
  my_region -= bounds.near;
  source_region -= source.bounds.near;
  
  float dx = (float)source_region.width() / (float)my_region.width();
  float dy = (float)source_region.height() / (float)my_region.height();
  if (method == Subsample) {  
    int i,j;     // Iterate over this image (muliply by dx, dy to get float coords)
    int si,sj;   // Integer version of source image iterators
    
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      si = source_region.near.x + (int)((i*dx)+1e-4);        // Make sure roundoff error doesn't get us!
      for (j=my_region.near.y ; j<=my_region.far.y ; j++) {
        sj = source_region.near.x + (int)((j*dy)+1e-4);
        raw(i,j) = source.peek(si, sj);
      }
      if (depth > source.depth)      for (j = my_region.near.y; j <= my_region.far.y; j++) raw(i, j) <<= (depth - source.depth);
      else if (depth < source.depth) for (j = my_region.near.y; j <= my_region.far.y; j++) raw(i, j) >>= (source.depth - depth);
    }
  }
  else if (method == LinearFit)
  {
    int i,j;               // Iterate over this image
    float fi0,fi1,fj0,fj1; // Real values over source image
    int si,sj;             // Integer values over source image
    float xfrac,yfrac;     // Width of strip in each direction
    float fI;              // Accumulated intensity
    float dx_i,dy_i;
    float shift_mult;      // Multiplicitive factor to take bit shift into account
    
    if (depth < source.depth) shift_mult = 1.0/(float)(1 << (source.depth-depth));
    else if (depth > source.depth) shift_mult = (float)(1 << (depth-source.depth));
    else shift_mult = 1.0;
    
    dx_i = 1.0/dx;
    dy_i = 1.0/dy;
    for (i=my_region.near.x ; i<=my_region.far.x ; i++)
    {
      fi0 = i*dx + source_region.near.x;
      fi1 = fi0 + dx;
      for (j=my_region.near.y ; j<=my_region.far.y ; j++)
      {
        fj0 = j*dy + source_region.near.y;
        fj1 = fj0 + dy;
        
        si = (int)(fi0+1e-4);
        if (si+1 > fi1) xfrac = fi1-fi0;
        else xfrac = (si+1)-fi0;
        xfrac *= dx_i;
        
        fI = 0.0;
        while (1) // Loop over source x; break out from middle
        {
          sj = (int)(fj0+1e-4);
          if (sj+1 > fj1) yfrac = fj1-fj0;
          else yfrac = (sj+1)-fj0;
          yfrac *= dy_i;
          
          while (1) // Loop over source y; break out from middle
          {
            fI += xfrac*yfrac*source.peek(si,sj);
            
            sj++;
            if (sj > fj1-1e-4) break;
            
            if (sj+1 > fj1) yfrac = (fj1-sj)*dy_i;
            else yfrac = dy_i;
          }
            
          si++;
          if (si > fi1-1e-4) break;
          
          if (si+1 > fi1) xfrac = (fi1-si)*dx_i;
          else xfrac = dx_i;
        }
        
        raw(i,j) = (short)(fI*shift_mult+0.5);
      }
    }
  }
  // No else clause--just do nothing for invalid image rescaling method
}

// Copy a source image with resizing/rescaling so both images show the same scene
void Image::mimic8(const Image8& source, Rectangle my_region, Rectangle source_region, ScaleType method)
{
#ifdef DEFER_MIMIC
  Image* temp = source.identical16();
  mimic(*temp, my_region, source_region, method);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->16 bits)
  // Convert region to point directly to image pixels
  if (bin > 1) my_region *= bin;
  if (source.bin > 1) source_region *= source.bin;
  if (!bounds.overlaps(my_region)) return;
  if (!source.bounds.overlaps(source_region)) return;
  my_region *= bounds;
  source_region *= source.bounds;
  my_region -= bounds.near;
  source_region -= source.bounds.near;
  
  float dx = (float)source_region.width() / (float)my_region.width();
  float dy = (float)source_region.height() / (float)my_region.height();
  if (method == Subsample) {  
    int i,j;     // Iterate over this image (muliply by dx, dy to get float coords)
    int si,sj;   // Integer version of source image iterators
    
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      si = source_region.near.x + (int)((i*dx)+1e-4);        // Make sure roundoff error doesn't get us!
      for (j=my_region.near.y ; j<=my_region.far.y ; j++) {
        sj = source_region.near.x + (int)((j*dy)+1e-4);
        raw(i,j) = source.peek(si, sj);
      }
      if (depth > 8)      for (j = my_region.near.y; j <= my_region.far.y; j++) raw(i, j) <<= (depth - 8);
      else if (depth < 8) for (j = my_region.near.y; j <= my_region.far.y; j++) raw(i, j) >>= (8 - depth);
    }
  }
  else if (method == LinearFit) {
    int i,j;               // Iterate over this image
    float fi0,fi1,fj0,fj1; // Real values over source image
    int si,sj;             // Integer values over source image
    float xfrac,yfrac;     // Width of strip in each direction
    float fI;              // Accumulated intensity
    float dx_i,dy_i;
    float shift_mult;      // Multiplicitive factor to take bit shift into account
    
    if (depth < 8) shift_mult = 1.0/(float)(1 << (8-depth));
    else if (depth > 8) shift_mult = (float)(1 << (depth-8));
    else shift_mult = 1.0;
    
    dx_i = 1.0/dx;
    dy_i = 1.0/dy;
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      fi0 = i*dx + source_region.near.x;
      fi1 = fi0 + dx;
      for (j=my_region.near.y ; j<=my_region.far.y ; j++) {
        fj0 = j*dy + source_region.near.y;
        fj1 = fj0 + dy;
        
        si = (int)(fi0+1e-4);
        if (si+1 > fi1) xfrac = fi1-fi0;
        else xfrac = (si+1)-fi0;
        xfrac *= dx_i;
        
        fI = 0.0;
        while (1) {      // Loop over source x; break out from middle
          sj = (int)(fj0+1e-4);
          if (sj+1 > fj1) yfrac = fj1-fj0;
          else yfrac = (sj+1)-fj0;
          yfrac *= dy_i;
          
          while (1) {      // Loop over source y; break out from middle
            fI += xfrac*yfrac*source.peek(si,sj);
            
            sj++;
            if (sj > fj1-1e-4) break;
            
            if (sj+1 > fj1) yfrac = (fj1-sj)*dy_i;
            else yfrac = dy_i;
          }
            
          si++;
          if (si > fi1-1e-4) break;
          
          if (si+1 > fi1) xfrac = (fi1-si)*dx_i;
          else xfrac = dx_i;
        }
        
        raw(i,j) = (short)(fI*shift_mult+0.5);
      }
    }
  }
  // No else clause--just do nothing for invalid image rescaling method
#endif
}
void Image::mimic8(const Image8& source , ScaleType method) { 
  mimic8(source,getBounds(),source.getBounds(),method); 
}



// Copy a source image onto (a part of) the existing image
void Image::copy(Point where,const Image& source,Point size,bool fix_depth)
{ 
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->8 bits)
  int x,y;
  if (bin <= 1 && source.bin <= 1) // Fast shortcut--grab memory directly
  {
    Point swhere = where - source.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    short *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    if (!fix_depth || depth==source.depth)
    {
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y)
      {
        if (size.y>3) memcpy(p,q,2*size.y);
        else for (y=0;y<size.y;y++) p[y] = q[y];
      }
    }
    else if (depth > source.depth)
    {
      short shift = depth - source.depth;
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y) for (y=0;y<size.y;y++) p[y] = q[y]<<shift;
    }
    else // depth < source.depth
    {
      short shift = source.depth - depth;
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y) for (y=0;y<size.y;y++) p[y] = q[y]>>shift;      
    }
  }
  else
  {
    Point stop = where + size;
    if (!fix_depth || depth==source.depth)
    {
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y));
    }
    else if (depth > source.depth)
    {
      short shift = depth - source.depth;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y)<<shift;
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y)<<shift);
    }
    else // depth < source.depth
    {
      short shift = source.depth - depth;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y)>>shift;
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y)>>shift);
    }
  }
}

// Same thing except copy using a mask
void Image::copy(const Image& source,Mask& m,bool fix_depth)
{
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->8 bits)
  int y,y0,y1;
  
  if (depth==source.depth) fix_depth = false;
  
  if (bin<=1 && source.bin<=1)
  {
    Rectangle safe = bounds * source.bounds;
    m.start();
    while (m.advance())
    {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1)
      {
        if (!fix_depth) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y);
        else if (depth>source.depth) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y)<<(depth-source.depth);
        else /*depth<source.depth*/ for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y)>>(source.depth-depth);
      }
    }
  }
  else
  {
    Rectangle safe = getBounds() * source.getBounds();
    m.start();
    while (m.advance())
    {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1)
      {
        if (!fix_depth)
        {
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y);
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y) );
        }
        else if (depth > source.depth)
        {
          short shift = depth - source.depth;
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y)<<shift;
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y)<<shift );          
        }
        else // depth < source.depth
        {
          short shift = source.depth - depth;
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y)>>shift;
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y)>>shift );          
        }
      }
    }
  }
}


// Copy a source image onto (a part of) the existing image
void Image::copy8(Point where, const Image8& source, Point size, bool fix_depth)
{ 
#ifdef DEFER_COPY
  Image* temp = source.identical16();
  copy(where, *temp, size, fix_depth);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->16 bits)
  int x,y;
  if (depth == 8) fix_depth = false;
  if (bin <= 1 && source.bin <= 1) { // Fast shortcut--grab memory directly
    Point swhere = where - source.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    uint8_t *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    if (!fix_depth) {
      for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = q[y];
    }
    else if (depth > 8) {
      short shift = depth - 8;
      for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = ((short)q[y])<<shift;
    }
    else { // depth < source.depth
      short shift = 8 - depth;
      for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = ((short)q[y])>>shift;      
    }
  }
  else
  {
    Point stop = where + size;
    if (!fix_depth) {
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y));
    }
    else if (depth > 8) {
      short shift = depth - 8;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = ((short)source.get(x,y))<<shift;
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , ((short)source.get(x,y))<<shift);
    }
    else // depth < 8
    {
      short shift = 8 - depth;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = ((short)source.get(x,y))>>shift;
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , ((short)source.get(x,y))>>shift);
    }
  }
#endif
}
void Image::copy8(const Image8& source, bool fix_depth) { 
  copy8( source.bounds.near , source , source.size , fix_depth );
}


// Same thing except copy using a mask
void Image::copy8(const Image8& source, Mask& m, bool fix_depth)
{
#ifdef DEFER_COPY
  Image* temp = source.identical16();
  copy(*temp, m, fix_depth);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->16 bits)
  int y,y0,y1;
  
  if (depth==8) fix_depth = false;
  
  if (bin<=1 && source.bin<=1) {
    Rectangle safe = bounds * source.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (!fix_depth)   for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y);
        else if (depth>8) for (y=y0;y<=y1;y++) rare(m.i().x,y) = ((short)source.view(m.i().x,y))<<(depth-8);
        else /*depth<8*/  for (y=y0;y<=y1;y++) rare(m.i().x,y) = ((short)source.view(m.i().x,y))>>(8-depth);
      }
    }
  }
  else {
    Rectangle safe = getBounds() * source.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (!fix_depth) {
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y);
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y) );
        }
        else if (depth > 8) {
          short shift = depth - 8;
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = ((short)source.get(m.i().x,y))<<shift;
          else for (y=y0;y<=y1;y++) set(m.i().x , y , ((short)source.get(m.i().x,y))<<shift );          
        }
        else // depth < source.depth
        {
          short shift = 8 - depth;
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = ((short)source.get(m.i().x,y))>>shift;
          else for (y=y0;y<=y1;y++) set(m.i().x , y , ((short)source.get(m.i().x,y))>>shift );          
        }
      }
    }
  }
#endif
}


// Adapt an image so that it is closer to a source image; requires rate>=0 and rate >= our depth - source depth
void Image::adapt(Point where,const Image& im,Point size,int rate)
{
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int x,y;
  bool irlz = false;
  int imrate = rate - depth + im.depth;
  if (imrate<0) { irlz=true; imrate = -imrate; }
  if (bin <= 1 && im.bin <= 1) {
    Point swhere = where - im.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    short *q = im.pixels + (swhere.x*im.size.y + swhere.y);
    if (imrate==0) for (x=0;x<size.x;x++,p+=this->size.y,q+=im.size.y) for (y=0;y<size.x;y++) p[y] += q[y] - (p[y]>>rate);
    else if (irlz) for (x=0;x<size.x;x++,p+=this->size.y,q+=im.size.y) for (y=0;y<size.x;y++) p[y] += (q[y]<<imrate) - (p[y]>>rate);
    else for (x=0;x<size.x;x++,p+=this->size.y,q+=im.size.y) for (y=0;y<size.x;y++) p[y] += (q[y]>>imrate) - (p[y]>>rate);
  }
  else {
    Point stop = where + size;
    if (bin<=1) {
      if (imrate==0) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) += im.get(x,y) - (rare(x,y)>>rate);
      else if (irlz) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) += (im.get(x,y)<<imrate) - (rare(x,y)>>rate);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) += (im.get(x,y)>>imrate) - (rare(x,y)>>rate);
    }
    else {
      short I;
      if (imrate==0) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) { I = get(x,y); set(x,y, I + im.get(x,y) - (I>>rate)); }
      else if (irlz) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) { I = get(x,y); set(x,y, I + (im.get(x,y)<<imrate) - (I>>rate)); }
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) { I = get(x,y); set(x,y, I + (im.get(x,y)>>imrate) - (I>>rate)); }
    }
  }
}

// Same thing except adapt using a mask
void Image::adapt(const Image& im,Mask &m,int rate)
{
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int y,y0,y1;
  bool irlz = false;
  int imrate = rate - depth + im.depth;
  if (imrate<0) { irlz=true; imrate = -imrate; }
  if (bin<=1 && im.bin<=1) {
    Rectangle safe = bounds * im.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0>y1) continue;
      if (imrate==0) for (y=y0;y<=y1;y++) rare(m.i().x,y) += im.view(m.i().x,y) - (rare(m.i().x,y)>>rate);
      else if (irlz) for (y=y0;y<=y1;y++) rare(m.i().x,y) += (im.view(m.i().x,y)<<imrate) - (rare(m.i().x,y)>>rate);
      else for (y=y0;y<=y1;y++) rare(m.i().x,y) += (im.view(m.i().x,y)>>imrate) - (rare(m.i().x,y)>>rate);
    }
  }
  else {
    short I;
    Rectangle safe = getBounds() * im.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (bin<=1) {
      	  if (imrate==0) for (y=y0;y<=y1;y++) rare(m.i().x,y) += im.get(m.i().x,y) - (rare(m.i().x,y)>>rate);
      	  else if (irlz) for (y=y0;y<=y1;y++) rare(m.i().x,y) += (im.get(m.i().x,y)<<imrate) - (rare(m.i().x,y)>>rate);
      	  else for (y=y0;y<=y1;y++) rare(m.i().x,y) += (im.get(m.i().x,y)>>imrate) - (rare(m.i().x,y)>>rate);
      	}
        else {
      	  if (imrate==0) for (y=y0;y<=y1;y++) { I = get(m.i().x,y); set(m.i().x , y , I + im.get(m.i().x,y) - (I>>rate) ); }
      	  else if (irlz) for (y=y0;y<=y1;y++) { I = get(m.i().x,y); set(m.i().x , y , I + (im.get(m.i().x,y)<<imrate) - (I>>rate) ); }
      	  else for (y=y0;y<=y1;y++) { I = get(m.i().x,y); set(m.i().x , y , I + (im.get(m.i().x,y)>>imrate) - (I>>rate) ); }
      	}
      }
    }
  }
}


// Adapt an image so that it is closer to a source image; requires rate>=0 and rate >= our depth - 8
void Image::adapt8(Point where, const Image8& im, Point size, int rate)
{
#ifdef DEFER_ADAPT
  Image* temp = im.identical16();
  adapt(where, *temp, size, rate);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the 16-bit-source implementation BY HAND!!
  int x,y;
  bool irlz = false;
  int imrate = rate - depth + 8;
  if (imrate<0) { irlz=true; imrate = -imrate; }
  if (bin <= 1 && im.bin <= 1) {
    Point swhere = where - im.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    uint8_t *q = im.pixels + (swhere.x*im.size.y + swhere.y);
    if (imrate==0) for (x=0;x<size.x;x++,p+=this->size.y,q+=im.size.y) for (y=0;y<size.x;y++) p[y] += ((short)q[y]) - (p[y]>>rate);
    else if (irlz) for (x=0;x<size.x;x++,p+=this->size.y,q+=im.size.y) for (y=0;y<size.x;y++) p[y] += (((short)q[y])<<imrate) - (p[y]>>rate);
    else for (x=0;x<size.x;x++,p+=this->size.y,q+=im.size.y) for (y=0;y<size.x;y++) p[y] += (((short)q[y])>>imrate) - (p[y]>>rate);
  }
  else {
    Point stop = where + size;
    if (bin<=1) {
      if (imrate==0) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) += ((short)im.get(x,y)) - (rare(x,y)>>rate);
      else if (irlz) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) += (((short)im.get(x,y))<<imrate) - (rare(x,y)>>rate);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) += (((short)im.get(x,y))>>imrate) - (rare(x,y)>>rate);
    }
    else {
      short I;
      if (imrate==0) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) { I = get(x,y); set(x,y, I + ((short)im.get(x,y)) - (I>>rate)); }
      else if (irlz) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) { I = get(x,y); set(x,y, I + (((short)im.get(x,y))<<imrate) - (I>>rate)); }
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) { I = get(x,y); set(x,y, I + (((short)im.get(x,y))>>imrate) - (I>>rate)); }
    }
  }
#endif
}

// Same thing except adapt using a mask
void Image::adapt8(const Image8& im, Mask &m, int rate)
{
#ifdef DEFER_ADAPT
  Image* temp = im.identical16();
  adapt(*temp, m, rate);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the 16-bit-source implementation BY HAND!!
  int y,y0,y1;
  bool irlz = false;
  int imrate = rate - depth + 8;
  if (imrate<0) { irlz=true; imrate = -imrate; }
  if (bin<=1 && im.bin<=1) {
    Rectangle safe = bounds * im.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0>y1) continue;
      if (imrate==0) for (y=y0;y<=y1;y++) rare(m.i().x,y) += ((short)im.view(m.i().x,y)) - (rare(m.i().x,y)>>rate);
      else if (irlz) for (y=y0;y<=y1;y++) rare(m.i().x,y) += (((short)im.view(m.i().x,y))<<imrate) - (rare(m.i().x,y)>>rate);
      else for (y=y0;y<=y1;y++) rare(m.i().x,y) += (((short)im.view(m.i().x,y))>>imrate) - (rare(m.i().x,y)>>rate);
    }
  }
  else {
    short I;
    Rectangle safe = getBounds() * im.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (bin<=1) {
          if (imrate==0) for (y=y0;y<=y1;y++) rare(m.i().x,y) += ((short)im.get(m.i().x,y)) - (rare(m.i().x,y)>>rate);
          else if (irlz) for (y=y0;y<=y1;y++) rare(m.i().x,y) += (((short)im.get(m.i().x,y))<<imrate) - (rare(m.i().x,y)>>rate);
          else for (y=y0;y<=y1;y++) rare(m.i().x,y) += (((short)im.get(m.i().x,y))>>imrate) - (rare(m.i().x,y)>>rate);
        }
        else {
          if (imrate==0) for (y=y0;y<=y1;y++) { I = get(m.i().x,y); set(m.i().x , y , I + ((short)im.get(m.i().x,y)) - (I>>rate) ); }
          else if (irlz) for (y=y0;y<=y1;y++) { I = get(m.i().x,y); set(m.i().x , y , I + (((short)im.get(m.i().x,y))<<imrate) - (I>>rate) ); }
          else for (y=y0;y<=y1;y++) { I = get(m.i().x,y); set(m.i().x , y , I + (((short)im.get(m.i().x,y))>>imrate) - (I>>rate) ); }
        }
      }
    }
  }
#endif
}



// Copy an image with background subtraction and a gray offset; requires bg depth >= source depth
// Does not scale the bit depth for the current image!  It's assumed to be one more than the source.  Make sure there's room!
// Division method adds 1 to fg & bg values to avoid errors with zeros (0/whatever = 0, whatever/0 = error).
void Image::diffCopy(Point where, const Image& source, Point size, const Image& bg)
{ 
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int x,y;
  short gray = 1 << source.depth;
  short shift = (source.depth<bg.depth) ? bg.depth-source.depth : 0;
  if (bin <= 1 && source.bin <= 1) { // Fast shortcut--grab memory directly
    int shiftmask = 0xFFFF>>shift;
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);

    Point swhere = where-source.bounds.near;
    Point bgwhere = where - bg.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    short *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    short *g = bg.pixels + (bgwhere.x*bg.size.y + bgwhere.y);
    
    if (divide_bg) { // This is 2*fg/(fg+bg) . . . slower but better for large changes in background!
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g += bg.size.y) for (y=0;y<size.y;y++)
          p[y] = ((1+(unsigned int)q[y])<<(source.depth+1)) / ((unsigned int)q[y] + (unsigned int)(g[y]>>shift) + 2);
    }          
    else { // This is just fg-bg . . . fast!  (Even faster when we use an int to do two short operations at once.)
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g += bg.size.y) {
        for (y=0;y<size.y-1;y+=2) *(int*)(p+y) = (dualgray + *(int*)(q+y)) - (((*(int*)(g+y))>>shift)&dualshiftmask);
        for (;y<size.y;y++) p[y] = q[y] - (g[y]>>shift) + gray;
      }
    }
  }
  else {
    Point stop = where + size;
    if (divide_bg) {
      unsigned int I_fg;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) {
        I_fg = 1 + (unsigned int)source.get(x,y);
        rare(x,y) = (I_fg<<(source.depth+1)) / (1 + I_fg + (unsigned int)(bg.get(x,y)>>shift));
      }
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) {
        I_fg = 1 + (unsigned int)source.get(x,y);
        set(x,y , (I_fg<<(source.depth+1)) / (1 + I_fg + (unsigned int)(bg.get(x,y)>>shift)) );
      }
    }
    else {
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y) - (bg.get(x,y)>>shift) + gray;
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y) - (bg.get(x,y)>>shift) + gray);
    }
  }
}

// Same thing except copy using a mask
void Image::diffCopy(const Image& source,Mask& m,const Image& bg)
{
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int y,y0,y1;
  unsigned int I_fg;
  short gray = 1 << source.depth;
  short shift = (source.depth<bg.depth) ? bg.depth-source.depth : 0;
  if (bin<=1 && source.bin<=1) {
    int shiftmask = 0xFFFF>>shift;
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);

    Rectangle safe = bounds * source.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (divide_bg) for (y=y0;y<=y1;y++) {
          I_fg = 1+(unsigned int)source.view(m.i().x,y);
          rare(m.i().x,y) = (I_fg<<(source.depth+1)) / (1 + I_fg + (unsigned int)(bg.view(m.i().x,y)>>shift));
        }
        else {
          for (y=y0;y<y1;y+=2) rareI(m.i().x,y) = (dualgray + source.viewI(m.i().x,y)) - ((bg.viewI(m.i().x,y)>>shift)&dualshiftmask);
          for (;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y) - (bg.view(m.i().x,y)>>shift) + gray;
        }
      }
    }
  }
  else {
    Rectangle safe = getBounds() * source.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (divide_bg) {
          if (bin<=1) for (y=y0;y<=y1;y++) {
            I_fg = 1+(unsigned int)source.get(m.i().x,y);
            rare(m.i().x,y) = (I_fg<<(source.depth+1))/(1 + I_fg + (unsigned int)(bg.get(m.i().x,y)>>shift));
          }
          else for (y=y0;y<=y1;y++) {
            I_fg = 1+(unsigned int)source.get(m.i().x,y);
            set(m.i().x,y , (I_fg<<(source.depth+1)) / (1 + I_fg + (unsigned int)(bg.get(m.i().x,y)>>shift)) );
          }
        }
        else {
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y) - (bg.get(m.i().x,y)>>shift) + gray;
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y) - (bg.get(m.i().x,y)>>shift) + gray );
        }
      }
    }
  }
}



// Copy an image with background subtraction and a gray offset; requires bg depth >= source depth
// Does not scale the bit depth for the background image!  It's assumed to be one more than the source.  Make sure there's room!
// Division method adds 1 to fg & bg values to avoid errors with zeros (0/whatever = 0, whatever/0 = error).
void Image::diffCopy8(Point where, const Image8& source, Point size, const Image& bg)
{
#ifdef DEFER_DIFFC
  Image* temp = source.identical16();
  diffCopy(where, *temp, size, bg);
  delete temp;
#else 
  // WARNING!  You must keep this up to date with the 16-bit-source implementation BY HAND!!
  int x,y;
  short gray = 256;
  short shift = (8<bg.depth) ? bg.depth-8 : 0;
  if (bin <= 1 && source.bin <= 1) { // Fast shortcut--grab memory directly
    int shiftmask = 0xFFFF>>shift;
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);

    Point swhere = where-source.bounds.near;
    Point bgwhere = where - bg.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    uint8_t *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    short *g = bg.pixels + (bgwhere.x*bg.size.y + bgwhere.y);
    
    if (divide_bg) { // This is 2*fg/(fg+bg) . . . slower but better for large changes in background!
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g += bg.size.y) for (y=0;y<size.y;y++)
          p[y] = ((1+(unsigned int)q[y])<<9) / ((unsigned int)q[y] + (unsigned int)(g[y]>>shift) + 2);
    }          
    else { // This is just fg-bg . . . fast!  (Even faster when we use an int to do two short operations at once.)
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g += bg.size.y) {
        for (y=0;y<size.y-1;y+=2) {
          unsigned short qy = *(unsigned short*)(q+y);
          *(int*)(p+y) = (dualgray + (qy & 0xFF) + ((int)(qy & 0xFF00) << 8)) - (((*(int*)(g+y))>>shift)&dualshiftmask);
        }
        for (;y<size.y;y++) p[y] = (short)q[y] - (g[y]>>shift) + gray;
      }
    }
  }
  else {
    Point stop = where + size;
    if (divide_bg) {
      unsigned int I_fg;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) {
        I_fg = 1 + (unsigned int)source.get(x,y);
        rare(x,y) = (I_fg<<9) / (1 + I_fg + (unsigned int)(bg.get(x,y)>>shift));
      }
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) {
        I_fg = 1 + (unsigned int)source.get(x,y);
        set(x,y , (I_fg<<9) / (1 + I_fg + (unsigned int)(bg.get(x,y)>>shift)) );
      }
    }
    else {
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y) - (bg.get(x,y)>>shift) + gray;
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y) - (bg.get(x,y)>>shift) + gray);
    }
  }
#endif
}

// Same thing except copy using a mask
void Image::diffCopy8(const Image8& source, Mask& m, const Image& bg)
{
#ifdef DEFER_DIFFC
  Image* temp = source.identical16();
  diffCopy(*temp, m, bg);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the 16-bit-source implementation BY HAND!!
  int y,y0,y1;
  unsigned int I_fg;
  short gray = 256;
  short shift = (8<bg.depth) ? bg.depth-8 : 0;
  if (bin<=1 && source.bin<=1) {
    int shiftmask = 0xFFFF>>shift;
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);

    Rectangle safe = bounds * source.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (divide_bg) for (y=y0;y<=y1;y++) {
          I_fg = 1+(unsigned int)source.view(m.i().x,y);
          rare(m.i().x,y) = (I_fg<<9) / (1 + I_fg + (unsigned int)(bg.view(m.i().x,y)>>shift));
        }
        else {
          for (y=y0;y<y1;y+=2) {
            unsigned short qy = source.viewS(m.i().x, y);
            rareI(m.i().x,y) = (dualgray + (qy & 0xFF) + ((int)(qy&0xFF00) << 8)) - ((bg.viewI(m.i().x,y)>>shift)&dualshiftmask);
          }
          for (;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y) - (bg.view(m.i().x,y)>>shift) + gray;
        }
      }
    }
  }
  else {
    Rectangle safe = getBounds() * source.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (divide_bg) {
          if (bin<=1) for (y=y0;y<=y1;y++) {
            I_fg = 1+(unsigned int)source.get(m.i().x,y);
            rare(m.i().x,y) = (I_fg<<9)/(1 + I_fg + (unsigned int)(bg.get(m.i().x,y)>>shift));
          }
          else for (y=y0;y<=y1;y++) {
            I_fg = 1+(unsigned int)source.get(m.i().x,y);
            set(m.i().x,y , (I_fg<<9) / (1 + I_fg + (unsigned int)(bg.get(m.i().x,y)>>shift)) );
          }
        }
        else {
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y) - (bg.get(m.i().x,y)>>shift) + gray;
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y) - (bg.get(m.i().x,y)>>shift) + gray );
        }
      }
    }
  }
#endif
}


// Copy an image with background subtraction and gray offset and adapt the background
void Image::diffAdaptCopy(Point where,const Image& source,Point size,Image& bg,int rate)
{ 
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int x,y;
  short gray = 1 << source.depth;
  short shift = (source.depth<bg.depth) ? bg.depth-source.depth : 0;
  short srcrate = rate - bg.depth + source.depth;  // Amount to bit shift each time we adapt
  bool sclz = false;  // Was srcrate less than zero?
  if (srcrate<0) { sclz=true; srcrate = -srcrate; }  // Yes, better use left shift instead of "negative right shift". 
  if (bin <= 1 && source.bin <= 1 && bg.bin <= 1) { // Fast shortcut--grab memory directly
    // Variables for "vector" processing of two shorts at a time using ints
    int shiftmask = 0xFFFF>>shift;
    int ratemask = 0xFFFF>>rate;
    int srcratemask = (sclz) ? 0 : ((0xFFFF>>srcrate)&0xFFFF);
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);
    int dualratemask = ratemask | (ratemask<<16);
    int dualsrcratemask = srcratemask | (srcratemask<<16);
    
    // Other task-specific variables
    Point swhere = where-source.bounds.near;
    Point bgwhere = where - bg.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    short *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    short *g = bg.pixels + (bgwhere.x*bg.size.y + bgwhere.y);
    if (sclz) {  // Background is deep--need to shift foreground left even given adaptation rate
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g+= bg.size.y) {
        if (divide_bg) {
          for (y=0;y<size.y;y++) { 
            g[y] += (q[y]<<srcrate) - (g[y]>>rate);
            p[y] = ((1+(unsigned int)q[y])<<(source.depth+1)) / (2 + (unsigned int)q[y] + (unsigned int)(g[y]>>shift));
          }
        }
        else {
          for (y=0;y<size.y-1;y+=2) {  // Two at a time with ints
            *((int*)(g+y)) = *((int*)(g+y)) + ((*(int*)(q+y))<<srcrate) - (((*(int*)(g+y))>>rate)&dualratemask);
            *((int*)(p+y)) = (dualgray + *((int*)(q+y))) - (((*(int*)(g+y))>>shift)&dualshiftmask);
          }
          for (;y<size.y;y++) { // Finish up with shorts (leave as for loop in case we switch to 64 bit mode and use long longs above)
            g[y] += (q[y]<<srcrate) - (g[y]>>rate);
            p[y] = q[y] - (g[y]>>shift) + gray;
          }
        }
      }
    }
    else { // Need to shift foreground right before updating background.
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g+= bg.size.y) {
        if (divide_bg) {
          for (y=0;y<size.y;y++) {
            g[y] += (q[y]>>srcrate) - (g[y]>>rate);
            p[y] = ((1+(unsigned int)q[y])<<(source.depth+1)) / (2 + (unsigned int)q[y] + (unsigned int)(g[y]>>shift));
          }
        }
        else {
          for (y=0;y<size.y-1;y+=2) {  // Batch process with ints
            *((int*)(g+y)) = *((int*)(g+y)) + (((*(int*)(q+y))>>srcrate)&dualsrcratemask) - (((*(int*)(g+y))>>rate)&dualratemask);
            *((int*)(p+y)) = (dualgray + *((int*)(q+y))) - (((*(int*)(g+y))>>shift)&dualshiftmask);
          }
          for (;y<size.y;y++) { // Finish with shorts 
            g[y] += (q[y]>>srcrate) - (g[y]>>rate);
            p[y] = q[y] - (g[y]>>shift) + gray;
          }
        }
      }
    }
  }
  else {
    short I,J;
    Point stop = where + size;
    for (x=where.x;x<stop.x;x++) {
      for (y=where.y;y<stop.y;y++) {
      	I = bg.get(x,y);
      	J = source.get(x,y);
      	if (srcrate==0) I += J - (I>>rate);
      	else if (sclz) I += (J<<srcrate) - (I>>rate);
      	else I += (J>>srcrate) - (I>>rate);
      	bg.set(x,y,I);
        if (divide_bg) set(x,y , ((1+(unsigned int)J)<<(source.depth+1)) / (2 + (unsigned int)J + (unsigned int)(I>>shift)));
      	else set(x,y , J - (I>>shift) + gray);
      }
    }
  }
}

// Same thing except copy using a mask
void Image::diffAdaptCopy(const Image& source,Mask& m,Image& bg,int rate)
{
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int y,y0,y1;
  int x,I_fg;
  short gray = 1 << source.depth;
  short shift = (source.depth<bg.depth) ? bg.depth-source.depth : 0;
  short srcrate = rate - bg.depth + source.depth;
  bool sclz = false;
  if (srcrate<0) { sclz=true; srcrate = -srcrate; }
  
  
  if (bin<=1 && source.bin<=1 && bg.bin<=1) {
    // Variables for "vector" processing of two shorts at a time using ints
    int shiftmask = 0xFFFF>>shift;
    int ratemask = 0xFFFF>>rate;
    int srcratemask = (sclz) ? 0 : ((0xFFFF>>srcrate)&0xFFFF);
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);
    int dualratemask = ratemask | (ratemask<<16);
    int dualsrcratemask = srcratemask | (srcratemask<<16);  
  
    Rectangle safe = bounds * source.bounds * bg.bounds;
    
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0>y1) continue;
      x = m.i().x;
      if (sclz) {
        if (divide_bg) {
          for (y=y0;y<=y1;y++) {
            bg.rare(x,y) += (source.view(x,y)<<srcrate) - (bg.rare(m.i().x,y)>>rate);
            I_fg = 1+(unsigned int)source.view(x,y);
            rare(x,y) = (I_fg<<(source.depth+1)) / ( 1 + I_fg + (unsigned int)(bg.rare(x,y)>>shift) );
          }
        }
        else {
          for (y=y0;y<y1;y+=2) { // "Vector" mode
            bg.rareI(x,y) = (bg.rareI(x,y) + (source.viewI(x,y)<<srcrate)) - ((bg.rareI(x,y)>>rate)&dualratemask);
            rareI(x,y) = (dualgray + source.viewI(x,y)) - ((bg.rareI(x,y)>>shift)&dualshiftmask);
          }
          for (;y<=y1;y++) { // Individual mode to finish any danglers on the end
      	    bg.rare(x,y) += (source.view(x,y)<<srcrate) - (bg.rare(x,y)>>rate);
	          rare(x,y) = source.view(x,y) - (bg.rare(x,y)>>shift) + gray;
          }
      	}
      }
      else {
        if (divide_bg) {
          for (y=y0;y<=y1;y++) {
            bg.rare(m.i().x,y) += (source.view(m.i().x,y)>>srcrate) - (bg.rare(m.i().x,y)>>rate);
            I_fg = 1+(unsigned int)source.view(m.i().x,y);
            rare(m.i().x,y) = (I_fg<<(source.depth+1)) / ( 1 + I_fg + (unsigned int)(bg.rare(m.i().x,y)>>shift) );
          }
        }
        else {
          for (y=y0;y<y1;y+=2) { // "Vector" mode
            bg.rareI(x,y) = (bg.rareI(x,y) + ((source.viewI(x,y)>>srcrate)&dualsrcratemask)) - ((bg.rareI(x,y)>>rate)&dualratemask);
            rareI(x,y) = (dualgray + source.viewI(x,y)) - ((bg.rareI(x,y)>>shift)&dualshiftmask);
          }
          for (;y<=y1;y++) { // Finish danglers
            bg.rare(m.i().x,y) += (source.view(m.i().x,y)>>srcrate) - (bg.rare(m.i().x,y)>>rate);
            rare(m.i().x,y) = source.view(m.i().x,y) - (bg.rare(m.i().x,y)>>shift) + gray;
          }
        }
      }
    }
  }
  else {
    short I,J;
    Rectangle safe = getBounds() * source.getBounds() * bg.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0 > y1) continue;
      
      if (sclz) {
        for (y=y0;y<=y1;y++) {
          I = bg.get(m.i().x,y);
          J = source.get(m.i().x,y);
          I += (J<<srcrate) - (I>>rate);
          bg.set(m.i().x,y,I);
          if (divide_bg) set(m.i().x,y , ((1+(unsigned int)J)<<(source.depth+1)) / (2 + (unsigned int)J + (unsigned int)(I>>shift)));
          else set(m.i().x , y , J - (I>>shift) + gray );
        }
      }
      else {
        for (y=y0;y<=y1;y++) {
          I = bg.get(m.i().x,y);
          J = source.get(m.i().x,y);
          I += (J>>srcrate) - (I>>rate);
          bg.set(m.i().x,y,I);
          if (divide_bg) set(m.i().x,y , ((1+(unsigned int)J)<<(source.depth+1)) / (2 + (unsigned int)J + (unsigned int)(I>>shift)));
          else set(m.i().x , y , J - (I>>shift) + gray );
        }
      }
    }
  }
}


// Copy an image with background subtraction and gray offset and adapt the background
void Image::diffAdaptCopy8(Point where, const Image8& source, Point size, Image& bg, int rate)
{ 
#ifdef DEFER_DAC
  Image* temp = source.identical16();
  diffAdaptCopy(where, *temp, size, bg, rate);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int x,y;
  short gray = 1 << 8;
  short shift = (8<bg.depth) ? bg.depth-8 : 0;
  short srcrate = rate - bg.depth + 8;  // Amount to bit shift each time we adapt
  bool sclz = false;  // Was srcrate less than zero?
  if (srcrate<0) { sclz=true; srcrate = -srcrate; }  // Yes, better use left shift instead of "negative right shift". 
  if (bin <= 1 && source.bin <= 1 && bg.bin <= 1) { // Fast shortcut--grab memory directly
    // Variables for "vector" processing of two shorts at a time using ints
    int shiftmask = 0xFFFF>>shift;
    int ratemask = 0xFFFF>>rate;
    int srcratemask = (sclz) ? 0 : ((0xFFFF>>srcrate)&0xFFFF);
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);
    int dualratemask = ratemask | (ratemask<<16);
    int dualsrcratemask = srcratemask | (srcratemask<<16);
    
    // Other task-specific variables
    Point swhere = where-source.bounds.near;
    Point bgwhere = where - bg.bounds.near;
    where -= bounds.near;
    short *p = pixels + (where.x*this->size.y + where.y);
    uint8_t *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    short *g = bg.pixels + (bgwhere.x*bg.size.y + bgwhere.y);
    if (sclz) {  // Background is deep--need to shift foreground left even given adaptation rate
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g+= bg.size.y) {
        if (divide_bg) {
          for (y=0;y<size.y;y++) { 
            g[y] += (((short)q[y])<<srcrate) - (g[y]>>rate);
            p[y] = ((1+(unsigned int)q[y])<<(8+1)) / (2 + (unsigned int)q[y] + (unsigned int)(g[y]>>shift));
          }
        }
        else {
          for (y=0;y<size.y-1;y+=2) {  // Two at a time with ints
            uint16_t qy = *((uint16_t*)(q+y));
            int qq = (((int)(qy & 0xFF00)) << 8) | (int)(qy & 0xFF);
            *((int*)(g+y)) = *((int*)(g+y)) + (qq<<srcrate) - (((*(int*)(g+y))>>rate)&dualratemask);
            *((int*)(p+y)) = (dualgray + qq) - (((*(int*)(g+y))>>shift)&dualshiftmask);
          }
          for (;y<size.y;y++) { // Finish up with shorts (leave as for loop in case we switch to 64 bit mode and use long longs above)
            g[y] += (((short)q[y])<<srcrate) - (g[y]>>rate);
            p[y] = q[y] - (g[y]>>shift) + gray;
          }
        }
      }
    }
    else { // Need to shift foreground right before updating background.
      for (x=0 ; x<size.x ; x++ , q += source.size.y , p += this->size.y , g+= bg.size.y) {
        if (divide_bg) {
          for (y=0;y<size.y;y++) {
            g[y] += (((short)q[y])>>srcrate) - (g[y]>>rate);
            p[y] = ((1+(unsigned int)q[y])<<(8+1)) / (2 + (unsigned int)q[y] + (unsigned int)(g[y]>>shift));
          }
        }
        else {
          for (y=0;y<size.y-1;y+=2) {  // Batch process with ints
            uint16_t qy = *((uint16_t*)(q+y));
            int qq = (((int)(qy & 0xFF00)) << 8) | (int)(qy & 0xFF);
            *((int*)(g+y)) = *((int*)(g+y)) + ((qq>>srcrate)&dualsrcratemask) - (((*(int*)(g+y))>>rate)&dualratemask);
            *((int*)(p+y)) = (dualgray + qq) - (((*(int*)(g+y))>>shift)&dualshiftmask);
          }
          for (;y<size.y;y++) { // Finish with shorts 
            g[y] += (((short)q[y])>>srcrate) - (g[y]>>rate);
            p[y] = ((short)q[y]) - (g[y]>>shift) + gray;
          }
        }
      }
    }
  }
  else {
    short I,J;
    Point stop = where + size;
    for (x=where.x;x<stop.x;x++) {
      for (y=where.y;y<stop.y;y++) {
        I = bg.get(x,y);
        J = source.get(x,y);
        if (srcrate==0) I += J - (I>>rate);
        else if (sclz) I += (J<<srcrate) - (I>>rate);
        else I += (J>>srcrate) - (I>>rate);
        bg.set(x,y,I);
        if (divide_bg) set(x,y , ((1+(unsigned int)J)<<(8+1)) / (2 + (unsigned int)J + (unsigned int)(I>>shift)));
        else set(x,y , J - (I>>shift) + gray);
      }
    }
  }
#endif
}

// Same thing except copy using a mask
void Image::diffAdaptCopy8(const Image8& source, Mask& m, Image& bg, int rate)
{
#ifdef DEFER_DAC
  Image* temp = source.identical16();
  diffAdaptCopy(*temp, m, bg, rate);
  delete temp;
#else
  // WARNING!  You must keep this up to date with the 8-bit-source implementation BY HAND!!
  int y,y0,y1;
  int x,I_fg;
  short gray = 1 << 8;
  short shift = (8<bg.depth) ? bg.depth-8 : 0;
  short srcrate = rate - bg.depth + 8;
  bool sclz = false;
  if (srcrate<0) { sclz=true; srcrate = -srcrate; }
  
  if (bin<=1 && source.bin<=1 && bg.bin<=1) {
    // Variables for "vector" processing of two shorts at a time using ints
    int shiftmask = 0xFFFF>>shift;
    int ratemask = 0xFFFF>>rate;
    int srcratemask = (sclz) ? 0 : ((0xFFFF>>srcrate)&0xFFFF);
    int dualgray = ((int)gray) | (((int)gray)<<16);
    int dualshiftmask = shiftmask | (shiftmask<<16);
    int dualratemask = ratemask | (ratemask<<16);
    int dualsrcratemask = srcratemask | (srcratemask<<16);  
  
    Rectangle safe = bounds * source.bounds * bg.bounds;

    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0>y1) continue;
      x = m.i().x;
      if (sclz) {
        if (divide_bg) {
          for (y=y0;y<=y1;y++) {
            bg.rare(x,y) += (((short)source.view(x,y))<<srcrate) - (bg.rare(m.i().x,y)>>rate);
            I_fg = 1+(unsigned int)source.view(x,y);
            rare(x,y) = (I_fg<<(8+1)) / ( 1 + I_fg + (unsigned int)(bg.rare(x,y)>>shift) );
          }
        }
        else {
          for (y=y0;y<y1;y+=2) { // "Vector" mode
            uint16_t sy = source.viewS(x,y);
            int iyy = (((int)(sy & 0xFF00)) << 8) | (sy & 0xFF);
            bg.rareI(x,y) = (bg.rareI(x,y) + (iyy<<srcrate)) - ((bg.rareI(x,y)>>rate)&dualratemask);
            rareI(x,y) = (dualgray + iyy) - ((bg.rareI(x,y)>>shift)&dualshiftmask);
          }
          for (;y<=y1;y++) { // Individual mode to finish any danglers on the end
            bg.rare(x,y) += (((short)source.view(x,y))<<srcrate) - (bg.rare(x,y)>>rate);
            rare(x,y) = ((short)source.view(x,y)) - (bg.rare(x,y)>>shift) + gray;
          }
        }
      }
      else {
        if (divide_bg) {
          for (y=y0;y<=y1;y++) {
            bg.rare(m.i().x,y) += (((short)source.view(m.i().x,y))>>srcrate) - (bg.rare(m.i().x,y)>>rate);
            I_fg = 1+(unsigned int)source.view(m.i().x,y);
            rare(m.i().x,y) = (I_fg<<(8+1)) / ( 1 + I_fg + (unsigned int)(bg.rare(m.i().x,y)>>shift) );
          }
        }
        else {
          for (y=y0;y<y1;y+=2) { // "Vector" mode
            uint16_t sy = source.viewS(x,y);
            int iyy = (((int)(sy & 0xFF00)) << 8) | (sy & 0xFF);
            bg.rareI(x,y) = (bg.rareI(x,y) + ((iyy>>srcrate)&dualsrcratemask)) - ((bg.rareI(x,y)>>rate)&dualratemask);
            rareI(x,y) = (dualgray + iyy) - ((bg.rareI(x,y)>>shift)&dualshiftmask);
          }
          for (;y<=y1;y++) { // Finish danglers
            bg.rare(m.i().x,y) += (((short)source.view(m.i().x,y))>>srcrate) - (bg.rare(m.i().x,y)>>rate);
            rare(m.i().x,y) = ((short)source.view(m.i().x,y)) - (bg.rare(m.i().x,y)>>shift) + gray;
          }
        }
      }
    }
  }
  else {
    short I,J;
    Rectangle safe = getBounds() * source.getBounds() * bg.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0 > y1) continue;
      
      if (sclz) {
        for (y=y0;y<=y1;y++) {
          I = bg.get(m.i().x,y);
          J = source.get(m.i().x,y);
          I += (J<<srcrate) - (I>>rate);
          bg.set(m.i().x,y,I);
          if (divide_bg) set(m.i().x,y , ((1+(unsigned int)J)<<(8+1)) / (2 + (unsigned int)J + (unsigned int)(I>>shift)));
          else set(m.i().x , y , J - (I>>shift) + gray );
        }
      }
      else {
        for (y=y0;y<=y1;y++) {
          I = bg.get(m.i().x,y);
          J = source.get(m.i().x,y);
          I += (J>>srcrate) - (I>>rate);
          bg.set(m.i().x,y,I);
          if (divide_bg) set(m.i().x,y , ((1+(unsigned int)J)<<(8+1)) / (2 + (unsigned int)J + (unsigned int)(I>>shift)));
          else set(m.i().x , y , J - (I>>shift) + gray );
        }
      }
    }
  }
#endif
}


// Flooding--iterative method to add one line at a time
void Image::floodLine(FloodInfo &info,FloodData *data,Stackable<Strip>*& head,Stackable<Strip>*& tail)
{
  int y,y0,y1;
  short I;
  Rectangle stop = getBounds();
  
  Stackable<Strip> *ss;
  Strip s = head->data;
  
  for (y=s.y0 ; y<=s.y1 ; y++)
  {
    // Fill current pixel
    I = info.im->get(s.x,y);
    if (info.T.excludes(I)) continue;
    info.floodPixel(s.x,y,I);
    
    // Fill pixels in the -y direction
    for (y0=y-1 ; y0 >= stop.near.y ; y0--)
    {
      I = info.im->get(s.x,y0);
      if (info.T.excludes(I)) break;
      info.floodPixel(s.x,y0,I);
    }
    
    // Fill pixels in the +y direction
    for (y1=y+1 ; y1 <= stop.far.y ; y1++)
    {
      I = info.im->get(s.x,y1);
      if (info.T.excludes(I)) break;
      info.floodPixel(s.x,y1,I);
    }
    
    // Overshot by one in each direction, so retract
    y0++;
    y1--;
    
    // Now we've filled a strip, so save our work if we can
    if (data!=NULL) data->stencil.addStrip(s.x,y0,y1);
    
    // Try to fill adjacent lines to this one (breadth first)
    if (s.id>=0 && s.x<stop.far.x)  // We were going forward--continue!
    {
      ss = info.store->create();
      ss->data.set(s.x+1,y0,y1,1);
      ss->appendTo(head,tail);
      
      if (s.id>0 && s.x>stop.near.x) // Try going backwards but avoid places we've been
      {
        if (y0 < s.y0)
        {
          ss = info.store->create();
          ss->data.set(s.x-1,y0,s.y0-1,-1);
          ss->appendTo(head,tail);
        }
        if (y1 > s.y1)
        {
          ss = info.store->create();
          ss->data.set(s.x-1,s.y1+1,y1,-1);
          ss->appendTo(head,tail);
        }
      }
    }
    if (s.id<=0 && s.x>stop.near.x) // Going backwards--continue!
    {
      ss = info.store->create();
      ss->data.set(s.x-1,y0,y1,-1);
      ss->appendTo(head,tail);
      
      if (s.id<0 && s.x<stop.far.x) // Try going forwards but avoid places we've been
      {
        if (y0 < s.y0)
        {
          ss = info.store->create();
          ss->data.set(s.x+1,y0,s.y0-1,1);
          ss->appendTo(head,tail);
        }
        if (y1 > s.y1)
        {
          ss = info.store->create();
          ss->data.set(s.x+1,s.y1+1,y1,1);
          ss->appendTo(head,tail);
        }
      }
    }
    
    // We already got up to y1 and we know that y1+1 is bad; keep going from there
    y = y1+1;
  }
  
  // We took care of first item on the list, so dispose of it
  info.store->destroy( Stackable<Strip>::pop(head) );
  if (head==NULL) tail=NULL;
}


// Fill a region of an image (turning values there negative) in some range
// Returns the number of pixels filled (0 means none)
int Image::floodFind(Point pt,DualRange threshold,Storage< Stackable<Strip> > *store,FloodData *result)
{
  short I = get(pt);
  if (I<=0 || threshold.start.excludes(I))
  {
    if (result!=NULL) result->stencil.flush();
    return 0;
  }
  else
  {
    Point p = pt - getBounds().near;
    FloodInfo info(this,store);
    Stackable<Strip> *head,*tail;
    
    info.T = threshold.keep;
    if (info.T.lo()<1) info.T.lo()=1;
    if (info.T.hi()<info.T.lo()) info.T.hi()=info.T.lo();
    if (threshold.start.lo()<1) threshold.start.lo()=1;
    if (threshold.start.hi()<threshold.start.lo()) threshold.start.hi()=threshold.start.lo();
    
    head = info.store->create();
    head->next=NULL;
    head->data.set(p.x,p.y,p.y,0);
    tail=head;
    
    if (result!=NULL) result->stencil.flush();
    
    while (head!=NULL) floodLine(info,result,head,tail);
    
    if (result!=NULL)
    {
      result->centroid = FPoint( ((double)info.cx)/((double)info.sumI) , ((double)info.cy)/((double)info.sumI) );
      result->stencil.sort();
      result->stencil.findBounds();
    }
    return info.pixel_count;
  }
}


// Fill an entire image in some range
int Image::floodRect(DualRange threshold,Storage< Stackable<Strip> > *stripstore,Storage< Listable<Strip> > *storemask,
                     ManagedList<FloodData>& result,Rectangle rect)
{
  Point pt;
  FloodInfo info(this,stripstore);
  FloodData* data;
  Stackable<Strip> *head,*tail;
  short I;
  int N_filled = 0;
  
  rect *= bounds;
  info.T=threshold.keep;
  if (info.T.lo()<1) info.T.lo()=1;
  if (info.T.hi()<info.T.lo()) info.T.hi()=info.T.lo();
  if (threshold.start.lo()<1) threshold.start.lo()=1;
  if (threshold.start.hi()<threshold.start.lo()) threshold.start.hi()=threshold.start.lo();
  
  for (pt.x = rect.near.x ; pt.x <= rect.far.x ; pt.x++)
  {
    for (pt.y = rect.near.y ; pt.y <= rect.far.y ; pt.y++)
    {
      I = get(pt);
      if (threshold.start.excludes(I)) continue;
      
      head = info.store->create();
      head->next = NULL;
      head->data.set(pt.x,pt.y,pt.y,0);
      tail=head;
      
      data = new( result.Append() ) FloodData(storemask);
      
      info.freshen();
      floodLine(info,data,head,tail);
      if (data->stencil.lines.size>0) pt.y=data->stencil.lines.h().y1+1;  // Skip what we already filled
      while (head!=NULL) floodLine(info,data,head,tail);
      
      if (info.pixel_count==0) result.Backspace();
      else
      {
        data->centroid = FPoint( ((double)info.cx)/((double)info.sumI) , ((double)info.cy)/((double)info.sumI) );
        data->stencil.sort();
        data->stencil.findBounds();
        N_filled++;
      }
    }
  }
  return N_filled;
}


// Fill an image starting from only the masked pieces
int Image::floodMask(DualRange threshold,Storage< Stackable<Strip> > *stripstore,Storage< Listable<Strip> > *storemask,
                     ManagedList<FloodData>& result,Mask &m)
{
  Point pt;
  FloodInfo info(this,stripstore);
  FloodData* data;
  Stackable<Strip> *head,*tail;
  short I;
  int N_filled = 0;
  int y0,y1;
  
  info.T=threshold.keep;
  if (info.T.lo()<1) info.T.lo()=1;
  if (info.T.hi()<info.T.lo()) info.T.hi()=info.T.lo();
  if (threshold.start.lo()<1) threshold.start.lo()=1;
  if (threshold.start.hi()<threshold.start.lo()) threshold.start.hi()=threshold.start.lo();
  
  m.start();
  while (m.advance())
  {
    if (m.i().beyond( bounds.far )) break;
    if (m.i().before( bounds.near )) continue;
    pt.x = m.i().x;
    y0 = (m.i().y0 < bounds.near.y) ? bounds.near.y : m.i().y0;
    y1 = (m.i().y1 > bounds.far.y) ? bounds.far.y : m.i().y1;
    
    if (y1-y0 < 0) continue;
    
    for (pt.y = y0 ; pt.y <= y1 ; pt.y++)
    {
      I = get(pt);
      if (threshold.start.excludes(I)) continue;
      
      head = info.store->create();
      head->next = NULL;
      head->data.set(pt.x,pt.y,pt.y,0);
      tail=head;
      
      data = new( result.Append() ) FloodData(storemask);
      
      info.freshen();
      floodLine(info,data,head,tail);
      if (data->stencil.lines.size>0) pt.y=data->stencil.lines.h().y1+1;  // Skip what we already filled
      while (head!=NULL) floodLine(info,data,head,tail);
      
      if (info.pixel_count==0) result.Backspace();
      else
      {
        data->centroid = FPoint( ((double)info.cx)/((double)info.sumI) , ((double)info.cy)/((double)info.sumI) );
        data->stencil.sort();
        data->stencil.findBounds();        
        N_filled++;
      }
    }
  }
  return N_filled;
}


// Finds the maximum value for each Y over all X in the rectangle
void Image::maxOverX(Rectangle target, short *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y ; y += 4) {
      __m64 max4 = *((__m64*)&rare(safe.near.x, y));
      for (int x = safe.near.x+1; x <= safe.far.x; x++) {
        __m64 next4 = *((__m64*)&rare(x, y));
        max4 = _mm_max_pi16(max4, next4);
      }
      *((__m64*)(output + (y - safe.near.y))) = max4;
    }
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) output[y-safe.near.y] = get(safe.near.x, y);
    for (int x = safe.near.x+1; x <= safe.far.x; x++) {
      for (int y = safe.near.y; y <= safe.far.y; y++) {
        auto p = output[y - safe.near.y];
        auto q = get(x, y);
        if (q > p) output[y - safe.near.y] = q;
      }
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the maximum value for each X over all Y in the rectangle
void Image::maxOverY(Rectangle target, short *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m64 accum = *((__m64*)&rare(x, safe.near.y));
      for (int y = safe.near.y+4; y <= safe.far.y; y += 4) {
        accum = _mm_max_pi16(accum, *((__m64*)&rare(x, y)));
      }
      short one = ((short*)&accum)[0];
      short two = ((short*)&accum)[1];
      short three = ((short*)&accum)[2];
      short four = ((short*)&accum)[3];
      if (three < four) three = four;
      if (one < two) one = two;
      output[x - safe.near.x] = (three < one) ? one : three;
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      short o = get(x, safe.near.y);
      for (int y = safe.near.y+1; y <= safe.far.y; y++) {
        auto q = get(x, y);
        if (q > o) o = q;
      }
      output[x - safe.near.x] = o;
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the minimum value for each Y over all X in the rectangle
void Image::minOverX(Rectangle target, short *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y ; y += 4) {
      __m64 min4 = *((__m64*)&rare(safe.near.x, y));
      for (int x = safe.near.x+1; x <= safe.far.x; x++) {
        __m64 next4 = *((__m64*)&rare(x, y));
        min4 = _mm_min_pi16(min4, next4);
      }
      *((__m64*)(output + (y - safe.near.y))) = min4;
    }
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) output[y-safe.near.y] = get(safe.near.x, y);
    for (int x = safe.near.x+1; x <= safe.far.x; x++) {
      for (int y = safe.near.y; y <= safe.far.y; y++) {
        auto p = output[y - safe.near.y];
        auto q = get(x, y);
        if (q < p) output[y - safe.near.y] = q;
      }
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the minimum value for each X over all Y in the rectangle
void Image::minOverY(Rectangle target, short *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m64 accum = *((__m64*)&rare(x, safe.near.y));
      for (int y = safe.near.y+4; y <= safe.far.y; y += 4) {
        accum = _mm_min_pi16(accum, *((__m64*)&rare(x, y)));
      }
      short one = ((short*)&accum)[0];
      short two = ((short*)&accum)[1];
      short three = ((short*)&accum)[2];
      short four = ((short*)&accum)[3];
      if (three > four) three = four;
      if (one > two) one = two;
      output[x - safe.near.x] = (three > one) ? one : three;
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      short o = get(x, safe.near.y);
      for (int y = safe.near.y+1; y <= safe.far.y; y++) {
        auto q = get(x, y);
        if (q < o) o = q;
      }
      output[x - safe.near.x] = o;
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean value for each Y averaged over all of X in the rectangle
void Image::meanOverX(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y; y += 4) {
      __m128 sum4 = _mm_cvtpi16_ps(*(__m64*)(&rare(safe.near.x, y)));
      for (int x = safe.near.x+1; x <= safe.far.x; x++) {
        __m128 next4 = _mm_cvtpi16_ps(*(__m64*)(&rare(x, y)));
        sum4 = _mm_add_ps(sum4, next4);
      }
      sum4 = _mm_div_ps(sum4, _mm_set_ps1(safe.width() < 1 ? 1 : safe.width()));
      *(__m128*)(output + (y-safe.near.y)) = sum4;
    }  
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) {
      int sum = 0;
      for (int x = safe.near.x; x <= safe.far.x; x++) sum += get(x, y);
      output[y - safe.near.y] = (float)((double)sum / safe.width());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean value for each X averaged over all of Y in the rectangle
void Image::meanOverY(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m128 accum = _mm_cvtpi16_ps(*(__m64*)(&rare(x, safe.near.y)));
      for (int y = safe.near.y+4; y <= safe.far.y; y += 4) {
        __m128 data4 = _mm_cvtpi16_ps(*(__m64*)(&rare(x, y)));
        accum = _mm_add_ps(accum, data4);
      }
      float sum = ((float*)(&accum))[0] + ((float*)(&accum))[1] + ((float*)(&accum))[2] + ((float*)(&accum))[3];
      output[x - safe.near.x] = sum / safe.height();
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      int sum = 0;
      for (int y = safe.near.y; y <= safe.far.y; y++) sum += get(x, y);
      output[x - safe.near.x] = (float)((double)sum / safe.height());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean square value for each Y averaged over all of X in the rectangle
void Image::meanSqOverX(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y; y += 4) {
      __m128 sum4 = _mm_cvtpi16_ps(*(__m64*)(&rare(safe.near.x, y)));
      sum4 = _mm_mul_ps(sum4, sum4);
      for (int x = safe.near.x+1; x <= safe.far.x; x++) {
        __m128 next4 = _mm_cvtpi16_ps(*(__m64*)(&rare(x, y)));
        sum4 = _mm_add_ps(sum4, _mm_mul_ps(next4, next4));
      }
      sum4 = _mm_div_ps(sum4, _mm_set_ps1(safe.width() < 1 ? 1 : safe.width()));
      *(__m128*)(output + (y-safe.near.y)) = sum4;
    }  
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) {
      double sum = 0;
      for (int x = safe.near.x; x <= safe.far.x; x++) { int v = get(x, y); sum += v*v; }
      output[y - safe.near.y] = (float)(sum / safe.width());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean square value for each X averaged over all of Y in the rectangle
void Image::meanSqOverY(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m128 accum = _mm_cvtpi16_ps(*(__m64*)(&rare(x, safe.near.y)));
      accum = _mm_mul_ps(accum, accum);
      for (int y = safe.near.y+4; y <= safe.far.y; y += 4) {
        __m128 data4 = _mm_cvtpi16_ps(*(__m64*)(&rare(x, y)));
        accum = _mm_add_ps(accum, _mm_mul_ps(data4, data4));
      }
      float sum = ((float*)(&accum))[0] + ((float*)(&accum))[1] + ((float*)(&accum))[2] + ((float*)(&accum))[3];
      output[x - safe.near.x] = sum / safe.height();
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      double sum = 0;
      for (int y = safe.near.y; y <= safe.far.y; y++) { int v = get(x, y); sum += v*v; }
      output[x - safe.near.x] = (float)(sum / safe.height());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the standard deviation for each Y across all X in the rectangle, given the means
void Image::deviationOverX(Rectangle target, const float* means, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y; y += 4) {
      __m128 sumsq4 = _mm_set_ps1(0.0);
      __m128 meanys = *((__m128*)(means + (y-safe.near.y)));
      for (int x = safe.near.x; x <= safe.far.x; x++) {
        __m128 next4 = _mm_cvtpi16_ps(*(__m64*)(&rare(x, y)));
        __m128 delta = _mm_sub_ps(next4, meanys);
        sumsq4 = _mm_add_ps(sumsq4, _mm_mul_ps(delta, delta));
      }
      sumsq4 = _mm_div_ps(sumsq4, _mm_set_ps1(safe.width() < 2 ? 1 : (safe.width() - 1)));
      *(__m128*)(output + (y-safe.near.y)) = _mm_sqrt_ps(sumsq4);
    }
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) {
      double ssqe = 0;
      auto m = means[y - safe.near.y];
      for (int x = safe.near.x; x <= safe.far.x; x++) {
        auto q = (float)get(x, y) - m;
        ssqe += (double)(q*q);
      }
      auto w = safe.width() - 1; if (w < 1) w = 1;
      output[y - safe.near.y] = (float)(sqrt(ssqe / w));
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the standard deviation for each X across all Y in the rectangle, given the means
void Image::deviationOverY(Rectangle target, const float* means, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m128 sumsq4 = _mm_set_ps1(0.0);
      __m128 meanxs = _mm_set_ps1(*(means + (x - safe.near.x)));
      for (int y = safe.near.y; y <= safe.far.y; y += 4) {
        __m128 next4 = _mm_cvtpi16_ps(*(__m64*)(&rare(x, y)));
        __m128 delta = _mm_sub_ps(next4, meanxs);
        sumsq4 = _mm_add_ps(sumsq4, _mm_mul_ps(delta, delta));
      }
      output[x - safe.near.x] = ((float*)(&sumsq4))[0] + ((float*)(&sumsq4))[1] + ((float*)(&sumsq4))[2] + ((float*)(&sumsq4))[3];
    }
    __m128 n = _mm_set_ps1((float)(safe.height() - 1));
    int x = safe.near.x;
    for ( ; x <= safe.far.x-3; x += 4) {
      __m128 four = *((__m128*)(output + (x - safe.near.x)));
      *(__m128*)(output + (x - safe.near.x)) = _mm_sqrt_ps(_mm_div_ps(four, n));
    }
    for (; x <= safe.far.x; x++) {
      _mm_store_ss(output + (x-safe.near.x), _mm_sqrt_ss(_mm_load_ps1(output + (x - safe.near.x))));
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      int sum = 0;
      auto m = means[x - safe.near.x];
      for (int y = safe.near.y; y <= safe.far.y; y++) {
        auto q = (float)get(x, y) - m;
        sum += (double)(q*q);
      }
      auto h = safe.height() - 1; if (h < 1) h = 1;
      output[x - safe.near.x] = (float)(sqrt(sum / h));
    }
#ifdef ENABLE_SIMD
  }
#endif
}



// Creates an appropriate TIFF header for this image; returns header size
int Image::makeTiffHeader(unsigned char *buffer)
{
  int bits = (depth<=8) ? 8 : 16;
  int bytes = (bits>>3)*getHeight()*getWidth();
  
  short* sbuf = (short*)buffer;
  int* ibuf = (int*)buffer;
  TiffIFD* tbuf = (TiffIFD*)(buffer + 10);
  
  sbuf[0] = 0x0102;
  if (buffer[0]==2)  buffer[0] = buffer[1] = 'I';         // Little endian
  else buffer[0] = buffer[1] = 'M';                       // Big endian 
  sbuf[1] = 42;                                           // Tiff file identifier
  ibuf[1] = 8;                                            // Offset of start of IFD
  sbuf[4] = 10;                                           // Ten entries in IFD
  tbuf[0] = TiffIFD(256,4,1,getHeight());                 // Image width (binned)--our data is row-minor so it's flipped
  tbuf[1] = TiffIFD(257,4,1,getWidth());                  // Image height (binned)
  tbuf[2] = TiffIFD(258,4,1,bits);                        // Bit depth (8 or 16)
  tbuf[3] = TiffIFD(259,4,1,1);                           // No compression
  tbuf[4] = TiffIFD(262,4,1,1);                           // Zero is black
  tbuf[5] = TiffIFD(273,4,1,30+10*sizeof(TiffIFD));       // Offsets of image
  tbuf[6] = TiffIFD(278,4,1,getWidth());                  // All of image in one strip
  tbuf[7] = TiffIFD(279,4,1,bytes);                       // Size of image
  tbuf[8] = TiffIFD(282,5,1,14+10*sizeof(TiffIFD));       // X resolution (dummy value)
  tbuf[9] = TiffIFD(283,5,1,22+10*sizeof(TiffIFD));       // Y resolution (dummy value)
  ibuf = (int*)(buffer + (10 + 10*sizeof(TiffIFD)));
  ibuf[0] = 0;                                            // No more IFDs
  ibuf[1] = 72;                                           // X resolution is 72 pixels...
  ibuf[2] = 1;                                            // ...per one inch
  ibuf[3] = 72;                                           // Y resolution is 72 pixels...
  ibuf[4] = 1;                                            // ...per one inch
  return 30 + 10*sizeof(TiffIFD);
}


// Write an image to a tiff file given an open file
int Image::writeTiff(FILE *f)
{
  if (!f) return 1;

  int i,x,y;
  int pixbytes;  
  unsigned char header[(46 + 10*sizeof(TiffIFD)) << 2];  // Way more space than we need
  unsigned char* data;
  
  if (bin<=1 && depth>8)
  {
    data = (unsigned char*)pixels;
    pixbytes = bounds.area() * sizeof(short);
  }
  else
  {
    Rectangle r = getBounds();
    i = 0;
    if (depth<=8)
    {
      data = new unsigned char[r.area()];
      for (x=r.near.x;x<=r.far.x;x++)
      {
        for (y=r.near.y;y<=r.far.y;y++)
        {
          data[i++] = (unsigned char)get(x,y);
        }
      }
      pixbytes = i;
    }
    else
    {
      short *tempdata = new short[r.area()];
      for (x=r.near.x;x<=r.far.x;x++)
      {
        for (y=r.near.y;y<=r.far.y;y++)
        {
          tempdata[i++] = get(x,y);
        }
      }
      data = (unsigned char*)tempdata;
      pixbytes = i * sizeof(short);
    }
  }
  
  i = makeTiffHeader(header);
  x = fwrite(header,i,1,f);
  y = fwrite(data,1,pixbytes,f);
  
  if (!(bin<=1)) delete[] data;
  
  if (x!=1 || y != pixbytes) return 1;

  return 0;
}


// Write an image to a tiff file given the filename
int Image::writeTiff(const char *fname)
{
  int i;
  FILE* f = fopen(fname,"wb");
  if (!f) return 1;
  i = writeTiff(f);
  fclose(f);
  return i;
}


// Write an image in "raw16" format to an open file.
// Format is width (2 bytes), height (2 bytes), data X-major.
// Everything is packed in native order.
int Image::writeRaw16(FILE *f) {
  if (!f) return 1;
  short header[] = { (short)bounds.width(), (short)bounds.height() };
  if (fwrite(header, 2, 2, f) != 2) return 1;
  if (bin <= 1 && bounds.near.x == 0 && bounds.near.y == 0 && size.x == bounds.width() && size.y == bounds.height()) {
    if (fwrite(pixels, 2, size.x*size.y, f) != (unsigned int)(size.x*size.y)) return 1;
  }
  else {
    short *buffer = new short[size.x * size.y];
    int i = 0;
    for (int x = 0; x < size.x; x++) for (int y = 0; y < size.y; y++) buffer[i++] = get(x, y);
    bool problem = fwrite(buffer, 2, size.x*size.y, f) != (unsigned int)(size.x*size.y);
    delete[] buffer;
    if (problem) return 1;
  }
  return 0;
}


// Write an image to a raw16 file given the filename
int Image::writeRaw16(const char* fname) {
  int i;
  FILE *f = fopen(fname, "wb");
  if (!f) return 1;
  i = writeRaw16(f);
  fclose(f);
  return i;
}


// Read an image in "raw16" format to a new object allocated with `new`.
// Format as in `writeRaw16(FILE*)`.
// Returns `NULL` if anything goes wrong.
Image* Image::readRaw16(FILE *f) {
  short size[2];
  if (fread(size, 2, 2, f) != 2) return NULL;
  unsigned int width = size[0];
  unsigned int height = size[1];
  if (width < 0 || height < 0 || width > 32767 || height > 32767) return NULL;
  Image *img = new Image(Point(width, height), false);
  if (fread(img->pixels, 2, width*height, f) != width*height) {
    delete[] img;
    return NULL;
  }
  return img;
}


// Read an image in raw16 format given a filename; `NULL` indicates something went wrong.
Image* Image::readRaw16(const char* fname) {
  FILE *f = fopen(fname, "rb");
  if (!f) return NULL;
  Image* ans = readRaw16(f);
  fclose(f);
  return ans;
}


// Write a text version of this image
void Image::println() const
{
  printf("%d,%d %d,%d\n",bounds.near.x,bounds.near.y,bounds.far.x,bounds.far.y);
  Rectangle b = getBounds();
  printf("%d,%d %d,%d\n",b.near.x,b.near.y,b.far.x,b.far.y);
  for (int i=b.near.x;i<=b.far.x;i++)
  {
    for (int j=b.near.y;j<=b.far.y;j++)
    {
      printf("%02x",get(i,j)&0xFF /*(char)(get(i,j)+' ')*/);
    }
    printf("\n");
  }
  printf("\n");
}


/****************************************************************
                         Image8 Methods
****************************************************************/


// Make a new image8 that copies a rectangular region of another
// If the old image8 is binned, it gets a copy of the binned (not original) version
Image8::Image8(const Image8& existing, const Rectangle& region, bool algo) : bin(0),divide_bg(algo)
{
  // WARNING!  You must keep this up to date with the 16 bit implementation BY HAND!!
  bounds = existing.getBounds();
  bounds.cropTo(region);
  if (bounds.isEmpty()) {
    pixels=NULL;
    size=Point(0,0);
    owns_pixels=false;
    return;
  }
  
  size = (bounds.far-bounds.near)+1;
  pixels = new uint8_t[ bounds.area() ];
  owns_pixels = true;
  
  copy(bounds.near , existing , size);
}


Image* Image8::identical16() const {
  Image* im = new Image();
  im->pixels = new short[ bounds.area() ];
  im->bounds = bounds;
  im->size = size;
  im->bin = bin;
  im->depth = 8;
  im->owns_pixels = true;
  im->divide_bg = divide_bg;
  int N = bounds.area();
  short *p = im->pixels;
  uint8_t *q = pixels;
  while (N > 0) { *p = *q; p++; q++; N--; }
  return im;
}

// Draw rectangles on an image8
void Image8::set(Rectangle r, uint8_t I)
{
  // WARNING!  You must keep this up to date with the 16 bit implementation BY HAND!!
  Point p;
  
  if (bin <= 1) {
    if (!r.overlaps(bounds)) return;
    r.cropTo(bounds);
    r -= bounds.near;
  }
  else {
    Rectangle s = getBounds();
    if (!r.overlaps(s)) return;
    r.cropTo(s);
    r -= s.near;
    r *= bin;
    r.far += bin-1;
  }
  for (p.x=r.near.x ; p.x<=r.far.x ; p.x++) for (p.y=r.near.y ; p.y<=r.far.y ; p.y++) raw(p) = I;
}


// Draw a mask onto an Image8
void Image8::set(Mask& m, uint8_t I)
{
  // WARNING!  You must keep this up to date with the 16 bit implementation BY HAND!!
  Point p;
  Rectangle b = getBounds();
  int x0,x1,y0,y1;
  m.start();
  while (m.advance()) {
    p.x = m.i().x;
    if (p.x < b.near.x) continue;
    if (p.x > b.far.x) break;
    if (m.i().y1 < b.near.y || m.i().y0 > b.far.y) continue;
    y0 = m.i().y0;
    if (y0 < b.near.y) y0 = b.near.y;
    y1 = m.i().y1;
    if (y1 > b.far.y) y1 = b.far.y;
    if (bin>1) {
      x0 = p.x*bin;
      x1 = x0 + bin - 1;
      y0 *= bin;
      y1 *= bin;
      y1 += bin-1;
      x0 -= bounds.near.x;
      x1 -= bounds.near.x;
      y0 -= bounds.near.y;
      y1 -= bounds.near.y;
      for (p.x=x0 ; p.x<=x1 ; p.x++) for (p.y=y0 ; p.y<=y1 ; p.y++) raw(p) = I;
    }
    else for (p.y = y0 ; p.y <= y1 ; p.y++) rare(p) = I;  
  }
}


// Draw lines on an Image8--not efficient for single pixel width, but works
void Image8::set(Point p1, Point p2, int width, uint8_t I)
{
  // WARNING!  You must keep this up to date with the 16 bit implementation BY HAND!!
  Rectangle r;
  
  if (width<1) return; else width--;
  
  if (p1.x==p2.x || p1.y==p2.y) {   // Single pixel or single horizontal or vertical line; easy.
    if (p1.x==p2.x)       r.near.x =         r.far.x = p1.x;
    else if (p1.x<p2.x) { r.near.x = p1.x;   r.far.x = p2.x-1; }
    else                { r.near.x = p2.x+1; r.far.x = p1.x; }
    if (p1.y==p2.y)       r.near.y =         r.far.y = p1.y;
    else if (p1.y<p2.y) { r.near.y = p1.y;   r.far.y = p2.y-1; }
    else                { r.near.y = p2.y+1; r.far.y = p1.y; }
    r.near -= width;
    r.far += width;
    set(r,I);
  }
  else {   // Diagonal, hard--draw a little rectangle of pixels around each point
    r.near = p1;
    r.far = p1;
    r.expand(width);
    set(r,I);
    
    Point p,px,pxy,py,q;
    FPoint o_hat = FPoint(p2-p1).unit();  // Direction we want to go, dropping pixels as we do it
    FPoint pf;
    float fx,fxy,fy;
    
    // Figure out where we should go if we move only in X, only in Y, or in both X and Y (diagonal) this step
    pxy.x = (o_hat.x<0) ? -1 : (o_hat.x>0) ? 1 : 0;
    pxy.y = (o_hat.y<0) ? -1 : (o_hat.y>0) ? 1 : 0;
    px.x = pxy.x;
    px.y = 0;
    py.x = 0;
    py.y = pxy.y;
    
    q = p2-p1;
    while (q*q > 2) {   // Farther than a diagonal to the last pixel--need to drop more pixels along the way
      // Figure out if we have a better bearing after moving in X or in XY
      pf = FPoint(p1+pxy).unit();
      fxy = pf*o_hat;
      pf = FPoint(p1+px).unit();
      fx = pf*o_hat;
      pf = FPoint(p1+py).unit();
      if (fx > fxy) {  // X was better, move in x and drop line of pixels (rest of block is already filled)
        r.near.x = r.far.x = p1.x+px.x*(width+1);
        r.near.y = p1.y - width;
        r.far.y = p1.y + width;
        set(r,I);
        p1 += px;
      }
      else {   // XY was better
        // See if Y is better yet
        pf = FPoint(p1+py).unit();
        fy = pf*o_hat;
        if (fy > fxy) {  // Yes, move in Y and drop line of pixels
          r.near.y = r.far.y = p1.y+py.y*(width+1);
          r.near.x = p1.x - width;
          r.far.x = p1.x + width;
          set(r,I);
          p1 += py;
        }
        else {  // XY was the way to go after all
          r.near.x = r.far.x = p1.x+pxy.x*(width+1);
          r.near.y = p1.y - width + pxy.y;
          r.far.y = p1.y + width + pxy.y;
          set(r,I);
          if (width>0) {   // Can't just drop a line, so drop a full block of pixels
            r.near.y = r.far.y = p1.y+pxy.y*(width+1);
            r.near.x = p1.x - width + pxy.x;
            r.far.x = p1.x + width + pxy.x - 1;
            set(r,I);
          }
          p1 += pxy;
        }
      }
      q = p2-p1;
    }
  }
}


// Draw a contour as a bunch of line segments
void Image8::set(Contour& c, int width, uint8_t I)
{
  // WARNING!  You must keep this up to date with the 16 bit implementation BY HAND!!
  Point p0,p1,p2;
  c.first();
  p0 = p1 = c.i();
  if (!c.next()) { set(p1,p1,width,I); return; }
  p2 = c.i();
  do {
    set(p1,p2,width,I);
    p1 = p2;
    c.next();
    p2 = c.i();
  } while (p1!=p0);
}


// Copy a source image with resizing/rescaling so both images show the same scene
void Image8::mimic(const Image8& source, Rectangle my_region, Rectangle source_region, ScaleType method)
{
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (16->8, 8->16, 16->16 bits)
  // Convert region to point directly to image pixels
  if (bin > 1) my_region *= bin;
  if (source.bin > 1) source_region *= source.bin;
  if (!bounds.overlaps(my_region)) return;
  if (!source.bounds.overlaps(source_region)) return;
  my_region *= bounds;
  source_region *= source.bounds;
  my_region -= bounds.near;
  source_region -= source.bounds.near;
  
  float dx = (float)source_region.width() / (float)my_region.width();
  float dy = (float)source_region.height() / (float)my_region.height();
  if (method == Subsample) {  
    int i,j;     // Iterate over this image (muliply by dx, dy to get float coords)
    int si,sj;   // Integer version of source image iterators
    
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      si = source_region.near.x + (int)((i*dx)+1e-4);        // Make sure roundoff error doesn't get us!
      for (j=my_region.near.y ; j<=my_region.far.y ; j++) {
        sj = source_region.near.x + (int)((j*dy)+1e-4);
        raw(i,j) = source.peek(si, sj);
      }
    }
  }
  else if (method == LinearFit) {
    int i,j;               // Iterate over this image
    float fi0,fi1,fj0,fj1; // Real values over source image
    int si,sj;             // Integer values over source image
    float xfrac,yfrac;     // Width of strip in each direction
    float fI;              // Accumulated intensity
    float dx_i,dy_i;
    
    dx_i = 1.0/dx;
    dy_i = 1.0/dy;
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      fi0 = i*dx + source_region.near.x;
      fi1 = fi0 + dx;
      for (j=my_region.near.y ; j<=my_region.far.y ; j++) {
        fj0 = j*dy + source_region.near.y;
        fj1 = fj0 + dy;
        
        si = (int)(fi0+1e-4);
        if (si+1 > fi1) xfrac = fi1-fi0;
        else xfrac = (si+1)-fi0;
        xfrac *= dx_i;
        
        fI = 0.0;
        while (1) {      // Loop over source x; break out from middle
          sj = (int)(fj0+1e-4);
          if (sj+1 > fj1) yfrac = fj1-fj0;
          else yfrac = (sj+1)-fj0;
          yfrac *= dy_i;
          
          while (1) {      // Loop over source y; break out from middle
            fI += xfrac*yfrac*source.peek(si,sj);
            
            sj++;
            if (sj > fj1-1e-4) break;
            
            if (sj+1 > fj1) yfrac = (fj1-sj)*dy_i;
            else yfrac = dy_i;
          }
            
          si++;
          if (si > fi1-1e-4) break;
          
          if (si+1 > fi1) xfrac = (fi1-si)*dx_i;
          else xfrac = dx_i;
        }
        
        raw(i,j) = (uint8_t)(fI + 0.5);
      }
    }
  }
  // No else clause--just do nothing for invalid image rescaling method
}

// Copy a source image with resizing/rescaling so both images show the same scene
void Image8::mimic16(const Image& source, Rectangle my_region, Rectangle source_region, ScaleType method)
{
  // WARNING!  You must keep this up to date with the other three implementations BY HAND!! (8->8, 8->16, 16->16 bits)
  // Convert region to point directly to image pixels
  if (bin > 1) my_region *= bin;
  if (source.bin > 1) source_region *= source.bin;
  if (!bounds.overlaps(my_region)) return;
  if (!source.bounds.overlaps(source_region)) return;
  my_region *= bounds;
  source_region *= source.bounds;
  my_region -= bounds.near;
  source_region -= source.bounds.near;
  
  float dx = (float)source_region.width() / (float)my_region.width();
  float dy = (float)source_region.height() / (float)my_region.height();
  if (method == Subsample) {  
    int i,j;     // Iterate over this image (muliply by dx, dy to get float coords)
    int si;      // Integer version of source image iterators
    
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      si = source_region.near.x + (int)((i*dx)+1e-4);        // Make sure roundoff error doesn't get us!
      if (source.depth == 8) {
        for (j=my_region.near.y ; j<=my_region.far.y ; j++)
          raw(i,j) = (uint8_t)source.peek(si, source_region.near.x + (int)((j*dy)+1e-4));
      }
      else if (source.depth > 8) {
        for (j=my_region.near.y ; j<=my_region.far.y ; j++)
          raw(i,j) = (uint8_t)(source.peek(si, source_region.near.x + (int)((j*dy)+1e-4)) >> (source.depth - 8));
      }
      else { //source.depth < 8
        for (j=my_region.near.y ; j<=my_region.far.y ; j++)
          raw(i,j) = (uint8_t)(source.peek(si, source_region.near.x + (int)((j*dy)+1e-4)) << (8 - source.depth));
      }
    }
  }
  else if (method == LinearFit) {
    int i,j;               // Iterate over this image
    float fi0,fi1,fj0,fj1; // Real values over source image
    int si,sj;             // Integer values over source image
    float xfrac,yfrac;     // Width of strip in each direction
    float fI;              // Accumulated intensity
    float dx_i,dy_i;
    float shift_mult;      // Multiplicitive factor to take bit shift into account
    
    if (8 < source.depth) shift_mult = 1.0/(float)(1 << (source.depth-8));
    else if (8 > source.depth) shift_mult = (float)(1 << (8-source.depth));
    else shift_mult = 1.0;
    
    dx_i = 1.0/dx;
    dy_i = 1.0/dy;
    for (i=my_region.near.x ; i<=my_region.far.x ; i++) {
      fi0 = i*dx + source_region.near.x;
      fi1 = fi0 + dx;
      for (j=my_region.near.y ; j<=my_region.far.y ; j++) {
        fj0 = j*dy + source_region.near.y;
        fj1 = fj0 + dy;
        
        si = (int)(fi0+1e-4);
        if (si+1 > fi1) xfrac = fi1-fi0;
        else xfrac = (si+1)-fi0;
        xfrac *= dx_i;
        
        fI = 0.0;
        while (1) {      // Loop over source x; break out from middle
          sj = (int)(fj0+1e-4);
          if (sj+1 > fj1) yfrac = fj1-fj0;
          else yfrac = (sj+1)-fj0;
          yfrac *= dy_i;
          
          while (1) {      // Loop over source y; break out from middle
            fI += xfrac*yfrac*source.peek(si,sj);
            
            sj++;
            if (sj > fj1-1e-4) break;
            
            if (sj+1 > fj1) yfrac = (fj1-sj)*dy_i;
            else yfrac = dy_i;
          }
            
          si++;
          if (si > fi1-1e-4) break;
          
          if (si+1 > fi1) xfrac = (fi1-si)*dx_i;
          else xfrac = dx_i;
        }
        
        raw(i,j) = (uint8_t)(fI*shift_mult+0.5);
      }
    }
  }
  // No else clause--just do nothing for invalid image rescaling method
}


// Copy a source image onto (a part of) the existing image
void Image8::copy(Point where, const Image8& source, Point size)
{ 
  int x,y;
  if (bin <= 1 && source.bin <= 1) { // Fast shortcut--grab memory directly
    Point swhere = where - source.bounds.near;
    where -= bounds.near;
    uint8_t *p = pixels + (where.x*this->size.y + where.y);
    uint8_t *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = q[y];
  }
  else {
    Point stop = where + size;
    if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = source.get(x,y);
    else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , source.get(x,y));
  }
}

// Same thing except copy using a mask
void Image8::copy(const Image8& source, Mask& m)
{
  int y,y0,y1;
  
  if (bin<=1 && source.bin<=1) {
    Rectangle safe = bounds * source.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.view(m.i().x,y);
    }
  }
  else {
    Rectangle safe = getBounds() * source.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y);
        else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y) );
      }
    }
  }
}

// Copy a source image onto (a part of) the existing image
void Image8::copy16(Point where, const Image& source, Point size, bool fix_depth)
{ 
  int x,y;
  if (8 == source.depth) fix_depth = false;
  if (bin <= 1 && source.bin <= 1) { // Fast shortcut--grab memory directly
    Point swhere = where - source.bounds.near;
    where -= bounds.near;
    uint8_t *p = pixels + (where.x*this->size.y + where.y);
    short *q = source.pixels + (swhere.x*source.size.y + swhere.y);
    if (!fix_depth) {
      for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = (uint8_t)q[y];
    }
    else if (source.depth < 8) {
      short shift = 8 - source.depth;
      for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = (uint8_t)(q[y]<<shift);
    }
    else { // source.depth > 8
      short shift = source.depth - 8;
      for (x=0 ; x<size.x ; x++, q += source.size.y, p += this->size.y) for (y=0; y<size.y; y++) p[y] = (uint8_t)(q[y]>>shift);      
    }
  }
  else {
    Point stop = where + size;
    if (!fix_depth) {
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = (uint8_t)source.get(x,y);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , (uint8_t)source.get(x,y));
    }
    else if (8 > source.depth) {
      short shift = 8 - source.depth;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = (uint8_t)(source.get(x,y)<<shift);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , (uint8_t)(source.get(x,y)<<shift));
    }
    else // 8 < source.depth
    {
      short shift = source.depth - 8;
      if (bin<=1) for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) rare(x,y) = (uint8_t)(source.get(x,y)>>shift);
      else for (x=where.x;x<stop.x;x++) for (y=where.y;y<stop.y;y++) set(x,y , (uint8_t)(source.get(x,y)>>shift));
    }
  }
}

// Same thing except copy using a mask
void Image8::copy16(const Image& source, Mask& m, bool fix_depth)
{
  int y,y0,y1;
  
  if (8 == source.depth) fix_depth = false;
  
  if (bin<=1 && source.bin<=1) {
    Rectangle safe = bounds * source.bounds;
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (!fix_depth) for (y=y0;y<=y1;y++) rare(m.i().x,y) = (uint8_t)source.view(m.i().x,y);
        else if (8 > source.depth) for (y=y0;y<=y1;y++) rare(m.i().x,y) = (uint8_t)(source.view(m.i().x,y)<<(8 - source.depth));
        else /*8 < source.depth*/ for (y=y0;y<=y1;y++) rare(m.i().x,y) = (uint8_t)(source.view(m.i().x,y)>>(source.depth - 8));
      }
    }
  }
  else {
    Rectangle safe = getBounds() * source.getBounds();
    m.start();
    while (m.advance()) {
      if (m.i().beyond( safe.far )) break;
      if (m.i().before( safe.near )) continue;
      y0 = (m.i().y0 < safe.near.y) ? safe.near.y : m.i().y0;
      y1 = (m.i().y1 > safe.far.y) ? safe.far.y : m.i().y1;
      if (y0<=y1) {
        if (!fix_depth) {
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = source.get(m.i().x,y);
          else for (y=y0;y<=y1;y++) set(m.i().x , y , source.get(m.i().x,y) );
        }
        else if (8 > source.depth) {
          short shift = 8 - source.depth;
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = (uint8_t)(source.get(m.i().x,y)<<shift);
          else for (y=y0;y<=y1;y++) set(m.i().x , y , (uint8_t)(source.get(m.i().x,y)<<shift) );          
        }
        else // source.depth > 8
        {
          short shift = source.depth - 8;
          if (bin<=1) for (y=y0;y<=y1;y++) rare(m.i().x,y) = (uint8_t)(source.get(m.i().x,y)>>shift);
          else for (y=y0;y<=y1;y++) set(m.i().x , y , (uint8_t)(source.get(m.i().x,y)>>shift) );          
        }
      }
    }
  }
}



// Finds the maximum value for each Y over all X in the rectangle
void Image8::maxOverX(Rectangle target, uint8_t *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y ; y += 8) {
      if (y+4 >= safe.far.y) {
        __m64 max8;
        __m64 next8;
        uint32_t last = *((uint32_t*)&rare(safe.near.x, y));
        *((uint32_t*)&max8) = last;
        *(((uint32_t*)&max8)+1) = last;
        *(((uint32_t*)&next8)+1) = last;
        for (int x = safe.near.x+1; x <= safe.far.x; x++) {
          *((uint32_t*)&next8) = *((uint32_t*)&rare(x, y));
          max8 = _mm_max_pu8(max8, next8);
        }
        *((uint32_t*)output) = *((uint32_t*)&max8);
      }
      else
      {
        __m64 max8 = *((__m64*)&rare(safe.near.x, y));
        for (int x = safe.near.x+1; x <= safe.far.x; x++) {
          __m64 next8 = *((__m64*)&rare(x, y));
          max8 = _mm_max_pu8(max8, next8);
        }
        *((__m64*)(output + (y - safe.near.y))) = max8;
      }
    }
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) output[y-safe.near.y] = get(safe.near.x, y);
    for (int x = safe.near.x+1; x <= safe.far.x; x++) {
      for (int y = safe.near.y; y <= safe.far.y; y++) {
        auto p = output[y - safe.near.y];
        auto q = get(x, y);
        if (q > p) output[y - safe.near.y] = q;
      }
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the maximum value for each X over all Y in the rectangle
void Image8::maxOverY(Rectangle target, uint8_t *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m64 accum;
      int y = safe.near.y;
      if (y+7 <= safe.far.y) {
        accum = *((__m64*)&rare(x, safe.near.y));
        y += 8;
      }
      else {
        uint32_t temp = *((uint32_t*)&rare(x, safe.near.y));
        *((uint32_t*)&accum) = temp;
        *(((uint32_t*)&accum)+1) = temp;
        y += 4;
      }
      for (; y+7 <= safe.far.y; y += 8) {
        accum = _mm_max_pu8(accum, *((__m64*)&rare(x, y)));
      }
      if (y+3 <= safe.far.y) {
        __m64 last;
        uint32_t temp = *((uint32_t*)&rare(x, y));
        *((uint32_t*)&last) = temp;
        *(((uint32_t*)&last)+1) = temp;
        accum = _mm_max_pu8(accum, last);
        y += 4;
      }
      accum = _mm_max_pu8(accum, _mm_shuffle_pi16(accum, 0x4E));
      accum = _mm_max_pu8(accum, _mm_shuffle_pi16(accum, 0xB1));
      uint16_t pair = *((uint16_t*)&accum);
      uint8_t lo = (uint8_t)pair;
      uint8_t hi = (uint8_t)(pair >> 8);
      output[x - safe.near.x] = (lo < hi) ? hi : lo;
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      short o = get(x, safe.near.y);
      for (int y = safe.near.y+1; y <= safe.far.y; y++) {
        auto q = get(x, y);
        if (q > o) o = q;
      }
      output[x - safe.near.x] = o;
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the minimum value for each Y over all X in the rectangle
void Image8::minOverX(Rectangle target, uint8_t *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y ; y += 8) {
      if (y+4 >= safe.far.y) {
        __m64 min8;
        __m64 next8;
        uint32_t last = *((uint32_t*)&rare(safe.near.x, y));
        *((uint32_t*)&min8) = last;
        *(((uint32_t*)&min8)+1) = last;
        *(((uint32_t*)&next8)+1) = last;
        for (int x = safe.near.x+1; x <= safe.far.x; x++) {
          *((uint32_t*)&next8) = *((uint32_t*)&rare(x, y));
          min8 = _mm_min_pu8(min8, next8);
        }
        *((uint32_t*)output) = *((uint32_t*)&min8);
      }
      else
      {
        __m64 min8 = *((__m64*)&rare(safe.near.x, y));
        for (int x = safe.near.x+1; x <= safe.far.x; x++) {
          __m64 next8 = *((__m64*)&rare(x, y));
          min8 = _mm_min_pu8(min8, next8);
        }
        *((__m64*)(output + (y - safe.near.y))) = min8;
      }
    }
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) output[y-safe.near.y] = get(safe.near.x, y);
    for (int x = safe.near.x+1; x <= safe.far.x; x++) {
      for (int y = safe.near.y; y <= safe.far.y; y++) {
        auto p = output[y - safe.near.y];
        auto q = get(x, y);
        if (q < p) output[y - safe.near.y] = q;
      }
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the minimum value for each X over all Y in the rectangle
void Image8::minOverY(Rectangle target, uint8_t *output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m64 accum;
      int y = safe.near.y;
      if (y+7 <= safe.far.y) {
        accum = *((__m64*)&rare(x, safe.near.y));
        y += 8;
      }
      else {
        uint32_t temp = *((uint32_t*)&rare(x, safe.near.y));
        *((uint32_t*)&accum) = temp;
        *(((uint32_t*)&accum)+1) = temp;
        y += 4;
      }
      for (; y+7 <= safe.far.y; y += 8) {
        accum = _mm_min_pu8(accum, *((__m64*)&rare(x, y)));
      }
      if (y+3 <= safe.far.y) {
        __m64 last;
        uint32_t temp = *((uint32_t*)&rare(x, y));
        *((uint32_t*)&last) = temp;
        *(((uint32_t*)&last)+1) = temp;
        accum = _mm_min_pu8(accum, last);
        y += 4;
      }
      accum = _mm_min_pu8(accum, _mm_shuffle_pi16(accum, 0x4E));
      accum = _mm_min_pu8(accum, _mm_shuffle_pi16(accum, 0xB1));
      uint16_t pair = *((uint16_t*)&accum);
      uint8_t lo = (uint8_t)pair;
      uint8_t hi = (uint8_t)(pair >> 8);
      output[x - safe.near.x] = (hi < lo) ? hi : lo;
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      short o = get(x, safe.near.y);
      for (int y = safe.near.y+1; y <= safe.far.y; y++) {
        auto q = get(x, y);
        if (q < o) o = q;
      }
      output[x - safe.near.x] = o;
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean value for each Y averaged over all of X in the rectangle
void Image8::meanOverX(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y; y += 4) {
      __m128 sum4 = _mm_set_ps1(0.0);
      for (int x = safe.near.x; x <= safe.far.x; x++) {
        uint8_t *p = &rare(x, y);
        __m128 next4 = _mm_cvtpi16_ps(_mm_set_pi16((short)p[3], (short)p[2], (short)p[1], (short)p[0]));
        sum4 = _mm_add_ps(sum4, next4);
      }
      sum4 = _mm_div_ps(sum4, _mm_set_ps1(safe.width() < 1 ? 1 : safe.width()));
      *(__m128*)(output + (y-safe.near.y)) = sum4;
    }  
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) {
      int sum = 0;
      for (int x = safe.near.x; x <= safe.far.x; x++) sum += get(x, y);
      output[y - safe.near.y] = (float)((double)sum / safe.width());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean value for each X averaged over all of Y in the rectangle
void Image8::meanOverY(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m128 accum = _mm_set_ps1(0.0);
      for (int y = safe.near.y; y <= safe.far.y; y += 4) {
        uint8_t *p = &rare(x, y);
        __m128 data4 = _mm_cvtpi16_ps(_mm_set_pi16((short)p[3], (short)p[2], (short)p[1], (short)p[0]));
        accum = _mm_add_ps(accum, data4);
      }
      float sum = ((float*)(&accum))[0] + ((float*)(&accum))[1] + ((float*)(&accum))[2] + ((float*)(&accum))[3];
      output[x - safe.near.x] = sum / safe.height();
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      int sum = 0;
      for (int y = safe.near.y; y <= safe.far.y; y++) sum += get(x, y);
      output[x - safe.near.x] = (float)((double)sum / safe.height());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean square value for each Y averaged over all of X in the rectangle
void Image8::meanSqOverX(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y; y += 4) {
      __m128 sum4 = _mm_set_ps1(0.0);
      for (int x = safe.near.x; x <= safe.far.x; x++) {
        uint8_t *p = &rare(x, y);
        __m128 next4 = _mm_cvtpi16_ps(_mm_set_pi16((short)p[3], (short)p[2], (short)p[1], (short)p[0]));
        sum4 = _mm_add_ps(sum4, _mm_mul_ps(next4, next4));
      }
      sum4 = _mm_div_ps(sum4, _mm_set_ps1(safe.width() < 1 ? 1 : safe.width()));
      *(__m128*)(output + (y-safe.near.y)) = sum4;
    }  
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) {
      double sum = 0;
      for (int x = safe.near.x; x <= safe.far.x; x++) { int v = get(x, y); sum += v*v; }
      output[y - safe.near.y] = (float)(sum / safe.width());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the mean square value for each X averaged over all of Y in the rectangle
void Image8::meanSqOverY(Rectangle target, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m128 accum = _mm_set_ps1(0.0);
      for (int y = safe.near.y; y <= safe.far.y; y += 4) {
        uint8_t *p = &rare(x, y);
        __m128 data4 = _mm_cvtpi16_ps(_mm_set_pi16((short)p[3], (short)p[2], (short)p[1], (short)p[0]));
        accum = _mm_add_ps(accum, _mm_mul_ps(data4, data4));
      }
      float sum = ((float*)(&accum))[0] + ((float*)(&accum))[1] + ((float*)(&accum))[2] + ((float*)(&accum))[3];
      output[x - safe.near.x] = sum / safe.height();
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      double sum = 0;
      for (int y = safe.near.y; y <= safe.far.y; y++) { int v = get(x, y); sum += v*v; }
      output[x - safe.near.x] = (float)(sum / safe.height());
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the standard deviation for each Y across all X in the rectangle, given the means
void Image8::deviationOverX(Rectangle target, const float* means, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int y = safe.near.y; y <= safe.far.y; y += 4) {
      __m128 sumsq4 = _mm_set_ps1(0.0);
      __m128 meanys = *((__m128*)(means + (y-safe.near.y)));
      for (int x = safe.near.x; x <= safe.far.x; x++) {
        uint8_t *p = &rare(x, y);
        __m128 next4 = _mm_cvtpi16_ps(_mm_set_pi16((short)p[3], (short)p[2], (short)p[1], (short)p[0]));
        __m128 delta = _mm_sub_ps(next4, meanys);
        sumsq4 = _mm_add_ps(sumsq4, _mm_mul_ps(delta, delta));
      }
      sumsq4 = _mm_div_ps(sumsq4, _mm_set_ps1(safe.width() < 2 ? 1 : (safe.width() - 1)));
      *(__m128*)(output + (y-safe.near.y)) = _mm_sqrt_ps(sumsq4);
    }
  }
  else {
#endif
    for (int y = safe.near.y; y <= safe.far.y; y++) {
      double ssqe = 0;
      auto m = means[y - safe.near.y];
      for (int x = safe.near.x; x <= safe.far.x; x++) {
        auto q = (float)get(x, y) - m;
        ssqe += (double)(q*q);
      }
      auto w = safe.width() - 1; if (w < 1) w = 1;
      output[y - safe.near.y] = (float)(sqrt(ssqe / w));
    }
#ifdef ENABLE_SIMD
  }
#endif
}

// Finds the standard deviation for each X across all Y in the rectangle, given the means
void Image8::deviationOverY(Rectangle target, const float* means, float* output) {
  Rectangle safe = target * getBounds();
#ifdef ENABLE_SIMD
  if (safe.height() % 4 == 0 && bin <= 1) {
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      __m128 sumsq4 = _mm_set_ps1(0.0);
      __m128 meanxs = _mm_set_ps1(*(means + (x - safe.near.x)));
      for (int y = safe.near.y; y <= safe.far.y; y += 4) {
        uint8_t *p = &rare(x, y);
        __m128 next4 = _mm_cvtpi16_ps(_mm_set_pi16((short)p[3], (short)p[2], (short)p[1], (short)p[0]));
        __m128 delta = _mm_sub_ps(next4, meanxs);
        sumsq4 = _mm_add_ps(sumsq4, _mm_mul_ps(delta, delta));
      }
      output[x - safe.near.x] = ((float*)(&sumsq4))[0] + ((float*)(&sumsq4))[1] + ((float*)(&sumsq4))[2] + ((float*)(&sumsq4))[3];
    }
    __m128 n = _mm_set_ps1((float)(safe.height() - 1));
    int x = safe.near.x;
    for ( ; x <= safe.far.x-3; x += 4) {
      __m128 four = *((__m128*)(output + (x - safe.near.x)));
      *(__m128*)(output + (x - safe.near.x)) = _mm_sqrt_ps(_mm_div_ps(four, n));
    }
    for (; x <= safe.far.x; x++) {
      _mm_store_ss(output + (x-safe.near.x), _mm_sqrt_ss(_mm_load_ps1(output + (x - safe.near.x))));
    }
  }
  else {
#endif
    for (int x = safe.near.x; x <= safe.far.x; x++) {
      int sum = 0;
      auto m = means[x - safe.near.x];
      for (int y = safe.near.y; y <= safe.far.y; y++) {
        auto q = (float)get(x, y) - m;
        sum += (double)(q*q);
      }
      auto h = safe.height() - 1; if (h < 1) h = 1;
      output[x - safe.near.x] = (float)(sqrt(sum / h));
    }
#ifdef ENABLE_SIMD
  }
#endif
}



// Creates an appropriate TIFF header for this image; returns header size
int Image8::makeTiffHeader(unsigned char *buffer)
{
  int bytes = getHeight()*getWidth();
  
  short* sbuf = (short*)buffer;
  int* ibuf = (int*)buffer;
  TiffIFD* tbuf = (TiffIFD*)(buffer + 10);
  
  sbuf[0] = 0x0102;
  if (buffer[0]==2)  buffer[0] = buffer[1] = 'I';         // Little endian
  else buffer[0] = buffer[1] = 'M';                       // Big endian 
  sbuf[1] = 42;                                           // Tiff file identifier
  ibuf[1] = 8;                                            // Offset of start of IFD
  sbuf[4] = 10;                                           // Ten entries in IFD
  tbuf[0] = TiffIFD(256,4,1,getHeight());                 // Image width (binned)--our data is row-minor so it's flipped
  tbuf[1] = TiffIFD(257,4,1,getWidth());                  // Image height (binned)
  tbuf[2] = TiffIFD(258,4,1,8);                           // Bit depth
  tbuf[3] = TiffIFD(259,4,1,1);                           // No compression
  tbuf[4] = TiffIFD(262,4,1,1);                           // Zero is black
  tbuf[5] = TiffIFD(273,4,1,30+10*sizeof(TiffIFD));       // Offsets of image
  tbuf[6] = TiffIFD(278,4,1,getWidth());                  // All of image in one strip
  tbuf[7] = TiffIFD(279,4,1,bytes);                       // Size of image
  tbuf[8] = TiffIFD(282,5,1,14+10*sizeof(TiffIFD));       // X resolution (dummy value)
  tbuf[9] = TiffIFD(283,5,1,22+10*sizeof(TiffIFD));       // Y resolution (dummy value)
  ibuf = (int*)(buffer + (10 + 10*sizeof(TiffIFD)));
  ibuf[0] = 0;                                            // No more IFDs
  ibuf[1] = 72;                                           // X resolution is 72 pixels...
  ibuf[2] = 1;                                            // ...per one inch
  ibuf[3] = 72;                                           // Y resolution is 72 pixels...
  ibuf[4] = 1;                                            // ...per one inch
  return 30 + 10*sizeof(TiffIFD);
}


// Write an image to a tiff file given an open file
int Image8::writeTiff(FILE *f)
{
  if (!f) return 1;

  int i,x,y;
  int pixbytes;  
  unsigned char header[(46 + 10*sizeof(TiffIFD)) << 2];  // Way more space than we need
  uint8_t* data;
  
  Rectangle r = getBounds();
  i = 0;
  data = new unsigned char[r.area()];
  for (x=r.near.x;x<=r.far.x;x++) for (y=r.near.y;y<=r.far.y;y++) data[i++] = get(x, y);
  pixbytes = i;
  
  i = makeTiffHeader(header);
  x = fwrite(header,i,1,f);
  y = fwrite(data,1,pixbytes,f);
  
  delete[] data;
  
  if (x != 1 || y != pixbytes) return 1;

  return 0;
}


// Write an image to a tiff file given the filename
int Image8::writeTiff(const char *fname)
{
  int i;
  FILE* f = fopen(fname,"wb");
  if (!f) return 1;
  i = writeTiff(f);
  fclose(f);
  return i;
}



// Write an image in "raw8" format to an open file.
// Format is width (2 bytes), height (2 bytes), data X-major.
// Everything is packed in native order.
int Image8::writeRaw8(FILE *f) {
  if (!f) return 1;
  short header[] = { (short)bounds.width(), (short)bounds.height() };
  if (fwrite(header, 2, 2, f) != 2) return 1;
  if (bin <= 1 && bounds.near.x == 0 && bounds.near.y == 0 && size.x == bounds.width() && size.y == bounds.height()) {
    if (fwrite(pixels, 1, size.x*size.y, f) != (unsigned int)(size.x*size.y)) return 1;
  }
  else {
    uint8_t *buffer = new uint8_t[size.x * size.y];
    int i = 0;
    for (int x = 0; x < size.x; x++) for (int y = 0; y < size.y; y++) buffer[i++] = get(x, y);
    bool problem = fwrite(buffer, 1, size.x*size.y, f) != (unsigned int)(size.x*size.y);
    delete[] buffer;
    if (problem) return 1;
  }
  return 0;
}


// Write an image to a raw8 file given the filename
int Image8::writeRaw8(const char* fname) {
  int i;
  FILE *f = fopen(fname, "wb");
  if (!f) return 1;
  i = writeRaw8(f);
  fclose(f);
  return i;
}


// Read an image in "raw8" format to a new object allocated with `new`.
// Format as in `writeRaw8(FILE*)`.
// Returns `NULL` if anything goes wrong.
Image8* Image8::readRaw8(FILE *f) {
  short size[2];
  if (fread(size, 2, 2, f) != 2) return NULL;
  unsigned int width = size[0];
  unsigned int height = size[1];
  if (width < 0 || height < 0 || width > 32767 || height > 32767) return NULL;
  Image8* img = new Image8(Point(width, height), false);
  if (fread(img->pixels, 1, width*height, f) != width*height) {
    delete[] img;
    return NULL;
  }
  return img;
}


// Read an image in raw8 format given a filename; `NULL` indicates something went wrong.
Image8* Image8::readRaw8(const char* fname) {
  FILE *f = fopen(fname, "rb");
  if (!f) return NULL;
  Image8* ans = readRaw8(f);
  fclose(f);
  return ans;
}


// Print a text version of the image
void Image8::println() const
{
  printf("%d,%d %d,%d\n",bounds.near.x,bounds.near.y,bounds.far.x,bounds.far.y);
  Rectangle b = getBounds();
  printf("%d,%d %d,%d\n",b.near.x,b.near.y,b.far.x,b.far.y);
  for (int i=b.near.x;i<=b.far.x;i++) {
    for (int j=b.near.y;j<=b.far.y;j++) printf("%02x",get(i,j));
    printf("\n");
  }
  printf("\n");
}


/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_image_strip()
{
  Strip sa[10];
  sa[0] = Strip(1,-1,2);
  sa[1] = Strip(1,2,5,0);
  sa[2].set(1,-5,0);
  sa[3].set(5,5,6,0);
  
  
  if (sa[0].length() + sa[1].length() + sa[2].length() + sa[3].length() != 16) return 1;
  if (sa[1].overlaps(sa[2]) || !sa[0].overlaps(sa[1]) || !sa[0].overlaps(sa[2])) return 2;
  if (!sa[3].beyond(sa[1]) || !sa[1].beyond(sa[2])) return 3;
  if (sa[3].before(sa[1]) || sa[1].before(sa[0])) return 4;
  if (!sa[0].overlaps(sa[3],Point(-4,-4)) || !sa[3].before(sa[1],Point(5,5))) return 5;
  
  sa[4] = sa[0]+Point(0,3);
  sa[5] = sa[1]-Point(0,3);
  if (sa[4]==sa[0] || sa[4]!=sa[1] || sa[5]==sa[1] || sa[5]!=sa[0]) return 6;
  if (!(sa[2]<sa[0] && sa[0]<sa[1] && sa[1]<sa[3])) return 7;
  
  return 0;
}

int test_mwt_image_mask()
{
  Ellipse ellip( Point(1,2) , 3 );
  ManagedList<Point> mlp(16),mlp1(16),mlp2(16),mlp3(16);
  Storage< Listable<Strip> > stor(16);
  Mask m1(32);
  Mask m2(&stor);
  Mask m3(Rectangle(Point(0,0),Point(5,5)),&stor);
  Mask m4(ellip,&stor);
  Mask m5(16);
  Mask m6(7);
  
  mlp.Append( Point(-3,-3) );
  mlp.Append( Point(-1,-1) );
  mlp.Append( Point(-1,0) );
  mlp.Append( Point(1,1) );
  mlp.Append( Point(3,3) );
  mlp.Append( Point(5,5) );
  mlp.Append( Point(7,7) );
  
  m1.imitate(m3);
  m1 += m4;
  m2.imitate(m4);
  m2 -= m3;
  m5.imitate(m3);
  m5 *= m4;
  m6.ellipticize(ellip);
  m6 *= Rectangle( Point(0,0) , Point(5,5) );
  
  mlp1.imitate(mlp);
  m3.contains(mlp1,mlp2);
  if (mlp1.size!=4 || mlp2.size!=3) return 1;
  mlp1.imitate(mlp);
  mlp2.flush();
  m4.contains(mlp1,mlp2);
  if (mlp1.size!=4 || mlp2.size!=3) return 2;
  mlp1.imitate(mlp);
  mlp2.flush();
  m1.contains(mlp1,mlp2);
  if (mlp1.size!=3 || mlp2.size!=4) return 3;
  mlp1.imitate(mlp);
  mlp2.flush();
  m6.contains(mlp1,mlp2);
  if (mlp1.size!=5 || mlp2.size!=2) return 4;
  
  if (m5.overlapRatio(m6) != FPoint(1.0,1.0)) return 5;
  if (m2.overlaps(m5,Point(0,0)) || !m2.overlaps(m6,Point(-2,-2))) return 6;
  
  m6.addStrip(Strip(3,-5,-3));
  m6.addStrip(Strip(1,6,7));
  m6.addStrip(Strip(5,-1,1));
  m6.sort();
  m6.start();
  m6.advance();
  Strip s = m6.i();
  while (m6.advance())
  {
    if (!s.before(m6.i())) return 7;
    s = m6.i();
  }
  
  if (m3.pixel_count != 36 || m4.pixel_count != 29 || m1.pixel_count != 43 || m2.pixel_count != 7 || m5.pixel_count != 22 || m6.pixel_count != 30) return 8; 
  if (m6.countOverlap(m1) !=  24 || m5.countOverlap(m4) != m5.pixel_count) return 9;
  
  Rectangle rect( Point(3,4) , Point(8,7) );
  m1.rectangularize( Rectangle( Point(-1,0) , Point(1,0) ) );
  m2.rectangularize( Rectangle( Point(0,-1) , Point(0,1) ) );
  
  m3.rectangularize( rect );
  m3.dilate(2);
  m4.rectangularize( rect.dup().expand(1) );
  if (m4.overlapRatio(m3) != FPoint(1.0,1.0)) return 10;
  
  m3.rectangularize( rect );
  m3.dilate(m1);
  m4.rectangularize( rect.dup().includeX(2).includeX(9) );
  if (m4.overlapRatio(m3) != FPoint(1.0,1.0)) return 11;
  
  m3.rectangularize( rect );
  m3.dilate(m2);
  m4.rectangularize( rect.dup().includeYY(3,8) );
  if (m4.overlapRatio(m3) != FPoint(1.0,1.0)) return 12;
  
  m3.rectangularize( rect );
  m4.flush();
  m3.extractEdge(m4);
  m5.flush();
  m5.addStrip( Strip(3,4,7) );
  for (int i=4 ; i<=7 ; i++)
  {
    m5.addStrip( Strip(i,4,4) );
    m5.addStrip( Strip(i,7,7) );
  }
  m5.addStrip( Strip(8,4,7) );
  if (m5.overlapRatio(m4) != FPoint(1.0,1.0)) return 13;
  
  m3.rectangularize( rect );
  m4.rectangularize( Rectangle( Point(0,0) , Point(10,10) ) );
  m4 -= m3;
  m3.invert( Rectangle(Point(0,0) , Point(10,10)) );
  if (m4.overlapRatio(m3) != FPoint(1.0,1.0)) return 14;
  
  m3.rectangularize( rect );
  m3.erode(2);
  m4.rectangularize( rect.dup().expand(-1) );
  if (m4.overlapRatio(m3) != FPoint(1.0,1.0)) return 15;
  
  m3.rectangularize( rect );
  m3.erode(m1);
  m4.rectangularize( Rectangle( Point(rect.near.x+1,rect.near.y) , Point(rect.far.x-1,rect.far.y) ) );
  if (m4.overlapRatio(m3) != FPoint(1.0,1.0)) return 16;
  
  return 0;
}

int test_mwt_image_contour()
{
  Mask m(7);
  m.ellipticize( Ellipse( Point(1,2) , 3 ) );
  m *= Rectangle( Point(0,0) , Point(5,5) );
  m.addStrip( Strip(3,-5,-3) );
  m.addStrip( Strip(1,6,7) );
  m.addStrip( Strip(5,-1,1) );
  m.sort();
  
  m.start();
  m.advance(); m.advance(); m.advance();
  m.i().y0 = 7; m.i().y1 = 8;
  
  Contour c( m );
  if (c.size()!=15) return 1;
  c.fillContour();
  if (c.size()!=20) return 2;
  
  Contour cc( 20 );
  cc.imitate(c);
  
  if (cc.size()!=20) return 3;
  
  // These tests were not for things that are used in the MWT
  // No reason they need to stay working, so disabled tests
  // Exercise caution if you need them again!
  //
  // c.emptyContour();
  // if (c.size()!=12) return 4;
  //
  // cc.approximateContour(1.0);
  // if (cc.size()!=6) { printf("%d\n",cc.size()); return 5; }
  
  return 0;
}

int test_mwt_image_image()
{
  int x,y;
  Mask m(32);
  Storage< Stackable<Strip> > ssstor(32);
  Storage< Listable<Strip> > slstor(32);
  ManagedList<FloodData> mlfd(32,true);
  Image im1(Point(10,11),false),im2(Point(6,6),false);
  im1 = 0;
  im2 = 5;
  im2 += Point(2,2);
  im1.copy( Point(2,2),im2,Point(6,6) );
  for (x=0;x<10;x++) for (y=0;y<11;y++) { if ((y<2 || y>8 || x<2 || x>8) && (im1.get(x,y)==5)) return 1; }
  
  im1.depth = im2.depth = 6;
  Image im3(im2,im2.getBounds(),false);
  im3 = 55;
  im3.adapt(Point(3,3),im2,Point(4,4),3);
  if (im3.get(2,2)!=55 || im3.get(3,3)!=49 || im3.get(6,6)!=49 || im3.get(7,7)!=55) return 2;
  
  im1.diffCopy(Point(4,4),im3,Point(2,2),im2);  // im1 might actually be 7 bits deep now, but we can ignore that here
  im1.diffAdaptCopy(Point(2,2),im3,Point(2,2),im2,3);
  if (im2.get(2,2)!=11 || im2.get(4,4)!=5) return 3;
#ifndef DIVIDE_BG
  if (im1.get(1,1)!=0 || im1.get(2,2)!=44+64 || im1.get(3,3)!=38+64 || im1.get(4,4)!=44+64 || im1.get(6,6)!=5) return 4;
#endif
  
  m.addStrip( Strip(2,4,8) );
  m.addStrip( Strip(3,8,8) );
  m.addStrip( Strip(4,8,8) );
  m.addStrip( Strip(5,8,8) );
  m.addStrip( Strip(6,5,8) );
  m.findBounds();
  Rectangle rrr = im1.getBounds();
  for (int yy = rrr.near.y ; yy <= rrr.far.y ; yy++)
  {
    for (int xx = rrr.near.x ; xx <= rrr.far.x ; xx++) printf("%02x ",im1.get(xx,yy));
    printf("\n");
  }
  im1.diffCopy(im1,m,im1);  // bg subtraction of self == set pixels to 64 in mask
  im1.set(3,3 , 38);
  x = im1.floodRect( DualRange(40,180,40,180) , &ssstor , &slstor , mlfd , im1.getBounds() );
  printf("%d %d\n",x,mlfd.h().stencil.pixel_count);
  if (x != 1 || mlfd.h().stencil.pixel_count!=19) return 5;
  
  //mlfd.flush();
  im1 = 0;
  im2 = 5;
  Mask m2(Rectangle(Point(2,2),Point(7,7)),&slstor);
  im1.copy(im2,m2);
  for (x=0;x<10;x++) for (y=0;y<11;y++) { if ((y<2 || y>8 || x<2 || x>8) && (im1.get(x,y)==5)) return 6; }
  
  im3 = 55;
  Mask m3(Rectangle(Point(3,3),Point(6,6)),&slstor);
  im3.adapt(im2,m3,3);
  if (im3.get(2,2)!=55 || im3.get(3,3)!=49 || im3.get(6,6)!=49 || im3.get(7,7)!=55) return 7;
  
  Mask m4(Rectangle(Point(4,4),Point(5,5)),&slstor);
  Mask m5(Rectangle(Point(2,2),Point(3,3)),&slstor);
  Image im4( im1.size*2, false );
  Image im5( im2 , im2.getBounds(), false );
  Image im6( im3 , im3.getBounds(), false );
  
  im4.bin = 2;
  im4.depth = im1.depth;
  im4.copy(im1);
  
  im1.diffCopy(im3,m4,im2);
  im4.diffCopy(im6,m4,im5);
  for (x=0;x<im1.size.x;x++) for (y=0;y<im1.size.y;y++) if (im1.get(x,y) != im4.get(x,y)) return 8;
  for (x=im2.bounds.near.x;x<=im2.bounds.far.x;x++) for (y=im2.bounds.near.y;y<=im2.bounds.far.y;y++) if (im2.get(x,y) != im5.get(x,y)) return 9;
  for (x=im3.bounds.near.x;x<=im3.bounds.far.x;x++) for (y=im3.bounds.near.y;y<=im3.bounds.far.y;y++) if (im3.get(x,y) != im6.get(x,y)) return 10;
  
  im1.diffAdaptCopy(im3,m5,im2,3);
  im4.diffAdaptCopy(im6,m5,im5,3);
  //for (x=0;x<im1.size.x;x++) for (y=0;y<im1.size.y;y++) printf("%d %d -> %d %d\n",x,y,(int)im2.get(x,y),(int)im5.get(x,y));
  for (x=0;x<im1.size.x;x++) for (y=0;y<im1.size.y;y++) if (im1.get(x,y) != im4.get(x,y)) return 11;
  for (x=im2.bounds.near.x;x<=im2.bounds.far.x;x++) for (y=im2.bounds.near.y;y<=im2.bounds.far.y;y++) if (im2.get(x,y) != im5.get(x,y)) return 12;
  
#ifndef DIVIDE_BG
  if (im1.get(1,1)!=0 || im1.get(2,2)!=44+64 || im1.get(3,3)!=38+64 || im1.get(4,4)!=44+64 || im1.get(6,6)!=5) return 13;
#endif
  
  im1.diffCopy(im1,m,im1);
  im1.set(3,3 , 38);
  x = im1.floodMask( DualRange(40,180,40,180) , &ssstor , &slstor , mlfd , m3 );
  if (x != 1 || mlfd.t().stencil.pixel_count!=19 || mlfd.size!=2) return 14;
  mlfd.t().principalAxes(im1,Range(40,180));
  
#ifdef UNIT_TEST_OWNER
  im1.depth = 16;
  im1.writeTiff("test_image.tiff");
#endif
  
  im1.depth = 5;
  im1 = 20;
  im1.set( Rectangle(Point(0,0) , Point(3,4)) , 30 );
  im2.mimic(im1);
  im2 -= Point(2,2);
  for (x=0;x<im2.size.x;x++)
  {
    for (y=0;y<im2.size.y;y++)
    {
      if (y<3 && x<3)
      {
        if (im2.get(x,y) != 30<<1) return 15;
      }
      else
      {
        if (im2.get(x,y) != 20<<1) return 16;
      }
    }
  }
  
  im1 = 20;
  im1.set( Rectangle(Point(0,0) , Point(4,5)) , 30 );
  im2.mimic(im1 , Image::LinearFit);
  for (x=0;x<im2.size.x;x++)
  {
    for (y=0;y<im2.size.y;y++)
    {
      if (y<3 && x<3)
      {
        if (im2.get(x,y) != 30<<1) return 17;
      }
      else if (y==3 && x<3)
      {
        if (im2.get(x,y) != 45) return 18;
      }
      else
      {
        if (im2.get(x,y) != 20<<1) return 19;
      }
    }
  }
  
  im1 = 20000;
  im1.depth = 16;
  im1.set( Rectangle(Point(0,0) , Point(4,5)) , 30000 );
  im1.bin = 2;
  if (im1.getBounds().area() != 25) return 20;
  if (im1.get(0,0) != 30000) return 21;
  if (im1.get(4,4) != 20000) return 22;
#ifdef UNIT_TEST_OWNER
  im1.writeTiff("test_bin2.tiff");
#endif
  im1.bin = 3;
  if (im1.getBounds().area() != 9) return 23;
  if (im1.get(0,0) != 30000) return 24;
  if (im1.get(2,2) != 20000) return 25;
#ifdef UNIT_TEST_OWNER
  im1.writeTiff("test_bin3.tiff");
#endif
  
#ifdef BENCHMARK
  Image huge_image(Point(1536,2048));
  Image big_image(Point(600,800));
  big_image.depth = huge_image.depth - 2;
  
  huge_image = 64;
  
  clock_t t0 = clock();
  
  for (x=0;x<100;x++) big_image.mimic(huge_image);
  
  clock_t t1 = clock();
  
  for (x=0;x<100;x++) big_image.mimic(huge_image,Image::LinearFit);
  
  clock_t t2 = clock();
  
  printf("%f per second\n", 100.0 / (( (double)t1-(double)t0 ) / (double)CLOCKS_PER_SEC) );
  printf("%f per second\n", 100.0 / (( (double)t2-(double)t1 ) / (double)CLOCKS_PER_SEC) );
#endif
  
  return 0;
}

int test_mwt_image_image8() {
  int x,y;
  Mask m(32);
  Storage< Stackable<Strip> > ssstor(32);
  Storage< Listable<Strip> > slstor(32);
  ManagedList<FloodData> mlfd(32,true);
  Image8 im1(Point(10,11),false),im2(Point(6,6),false);
  im1 = 0;
  im2 = 5;
  im2 += Point(2,2);
  im1.copy( Point(2,2),im2,Point(6,6) );
  for (x=0;x<10;x++) for (y=0;y<11;y++) { if ((y<2 || y>8 || x<2 || x>8) && (im1.get(x,y)==5)) return 1; }
  
  Image jm3(im2.getBounds(),false);
  jm3.depth = 8;
  jm3 = 55;
  jm3.adapt8(Point(3,3),im2,Point(4,4),3);
  if (jm3.get(2,2)!=55 || jm3.get(3,3)!=49 || jm3.get(6,6)!=49 || jm3.get(7,7)!=55) return 2;
  
  Image jm1(im1.getBounds(), false);
  jm1.depth = 8;
  jm1.copy8(im1);
  Image jm2(im2.getBounds(), false);
  jm2.depth = 8;
  jm2.copy8(im2);
  Image8 im3(jm3.getBounds(), false);
  im3.copy16(jm3);
  Image km1(jm1, jm1.getBounds(), false);
  Image km2(jm2, jm2.getBounds(), false);
  for (x = 0; x < jm1.size.x; x++) for (y = 0; y < jm1.size.y; y++) if (jm1.raw(x,y) != im1.raw(x,y) || km1.raw(x,y) != im1.raw(x,y)) return 4;
  for (x = 0; x < jm2.size.x; x++) for (y = 0; y < jm2.size.y; y++) if (jm2.raw(x,y) != im2.raw(x,y) || km2.raw(x,y) != im2.raw(x,y)) return 5;
  for (x = 0; x < jm3.size.x; x++) for (y = 0; y < jm3.size.y; y++) if (jm3.raw(x,y) != im3.raw(x,y)) return 6;
  jm1.diffCopy8(Point(4,4),im3,Point(2,2),jm2); 
  km1.diffCopy(Point(4,4),jm3,Point(2,2),jm2);
  for (x = 0; x < jm1.size.x; x++) for (y = 0; y < jm1.size.y; y++) if (jm1.raw(x,y) != km1.raw(x,y)) return 7;
  jm1.diffAdaptCopy8(Point(2,2),im3,Point(2,2),jm2,3);
  km1.diffAdaptCopy(Point(2,2),jm3,Point(2,2),km2,3);
  for (x = 0; x < jm2.size.x; x++) for (y = 0; y < jm2.size.y; y++) if (jm2.raw(x,y) != km2.raw(x,y)) return 8;
  for (x = 0; x < jm1.size.x; x++) for (y = 0; y < jm1.size.y; y++) if (jm1.raw(x,y) != km1.raw(x,y)) return 9;
  
  m.addStrip( Strip(2,2,5) );
  m.addStrip( Strip(3,4,6) );
  m.addStrip( Strip(4,2,4) );
  m.findBounds();
  for (x = 0; x < jm2.size.x; x++) for (y = 0; y < jm2.size.y; y++) jm2.raw(x,y) = km2.raw(x,y) = 0x50;
  
  jm1.diffCopy8(im3, m, jm2);
  km1.diffCopy(jm3, m, km2);
  jm1.diffAdaptCopy8(im3, m, jm2, 3);
  km1.diffAdaptCopy(jm3, m, km2, 3);
  for (x = 0; x < jm2.size.x; x++) for (y = 0; y < jm2.size.y; y++) if (jm2.raw(x,y) != km2.raw(x,y)) return 10;
  for (x = 0; x < jm1.size.x; x++) for (y = 0; y < jm1.size.y; y++) if (jm1.raw(x,y) != km1.raw(x,y)) return 11;

  return 0;
}

bool test_same_bytes(uint8_t *a, uint8_t *b, int n) {
  while (n > 0) {
    if (*a++ != *b++) return false;
    n--;
  }
  return true;
}

void test_print_bytes(uint8_t *a, uint8_t *b, int n) {
  while (n > 0) {
    printf("%02X %02X\n", *a++, *b++);
    n --;
  }
}

bool test_same_shorts(short *a, short *b, int n) {
  while (n > 0) {
    if (*a++ != *b++) return false;
    n--;
  }
  return true;
}

void test_print_shorts(short *a, short *b, int n) {
  while (n > 0) {
    printf("%04X %04X\n", *a++, *b++);
    n --;
  }
}

void test_print_floats(float *a, float *b, int n) {
  while (n > 0) {
    printf("%8.4f %8.4f\n", *a++, *b++);
    n --;
  }
}

bool test_same_floats(float *a, float *b, int n) {
  while (n > 0) {
    float err = fabsf(*a++ - *b++);
    if (err > 1e-2) return false;
    n--;
  }
  return true;
}

int test_mwt_image_flattening() {
  int test_image_height = 8;
  int test_image_width = 12;
  short test_image[] = {
    0x00, 0x00, 0x00, 0x00, 0x80, 0x80, 0x80, 0xA0,
    0x00, 0x00, 0x00, 0x00, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x00, 0x80, 0x80, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x10, 0x80, 0x80, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x20, 0x80, 0x80, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x30, 0x80, 0x80, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x40, 0x80, 0x80, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x50, 0x80, 0x80, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x60, 0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x70, 0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x80, 0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF,
    0x10, 0x80, 0xFF, 0xFF, 0xFF, 0x80, 0xFF, 0xFF
  };
  short test_image_max_y[] = {0xA0, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  short test_image_max_x[] = {0x10, 0x80, 0xFF, 0xFF, 0xFF, 0x80, 0xFF, 0xFF};
  short test_image_min_y[] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10};
  short test_image_min_x[] = {0x00, 0x00, 0x00, 0x00, 0x80, 0x80, 0x80, 0xA0};
  float test_image_mean_y[] = {68, 95.75, 127.75, 129.75, 131.75, 149.625, 151.625, 153.625, 171.5, 173.5, 175.5, 193.375};
  float test_image_mean_x[] = {1.333333, 58.666668, 117.25, 149, 202.083328, 128, 244.416672, 247.083328};
  float test_image_meansq_y[] = {9344.0, 20352.25, 24448.25, 24480.25, 24576.25, 30816.375, 31040.375, 31328.375, 37760.5, 38176.5, 38656.5, 44768.625};
  float test_image_meansq_x[] = {21.3333, 5717.3333, 17706.75, 29867.0, 44757.9167, 16384.0, 60971.5833, 61739.583};
  float test_image_sd_y[] = {73.4458, 113.0572, 96.3813, 93.4738, 90.826, 98.147, 95.918, 93.9771, 97.6773, 96.061, 94.7553, 91.8056};
  float test_image_sd_x[] = {4.6188, 49.8239, 65.72, 91.4489, 65.3959, 0, 36.6617, 27.4241};

  short result_short_y[12];
  short result_short_x[8];
  float result_float_y[12];
  float result_float_x[8];
  float result2_float_y[12];
  float result2_float_x[8];

  Image im(test_image, Point(test_image_width, test_image_height), false);

  im.maxOverY(Rectangle(0, test_image_width, 0, test_image_height), result_short_y);
  if (!test_same_shorts(result_short_y, test_image_max_y, test_image_width)) {
    test_print_shorts(result_short_y, test_image_max_y, test_image_width);
    return 1;
  }

  im.maxOverX(Rectangle(0, test_image_width, 0, test_image_height), result_short_x);
  if (!test_same_shorts(result_short_x, test_image_max_x, test_image_height)) {
    test_print_shorts(result_short_x, test_image_max_x, test_image_height);
    return 2;
  }

  im.minOverY(Rectangle(0, test_image_width, 0, test_image_height), result_short_y);
  if (!test_same_shorts(result_short_y, test_image_min_y, test_image_width)) {
    test_print_shorts(result_short_y, test_image_min_y, test_image_width);
    return 3;
  }

  im.minOverX(Rectangle(0, test_image_width, 0, test_image_height), result_short_x);
  if (!test_same_shorts(result_short_x, test_image_min_x, test_image_height)) {
    test_print_shorts(result_short_x, test_image_min_x, test_image_height);
    return 4;
  }

  im.meanOverY(Rectangle(0, test_image_width, 0, test_image_height), result_float_y);
  if (!test_same_floats(result_float_y, test_image_mean_y, test_image_width)) {
    test_print_floats(result_float_y, test_image_mean_y, test_image_width);
    return 5;
  }

  im.meanOverX(Rectangle(0, test_image_width, 0, test_image_height), result_float_x);
  if (!test_same_floats(result_float_x, test_image_mean_x, test_image_height)) {
    test_print_floats(result_float_x, test_image_mean_x, test_image_height);
    return 6;
  }

  im.meanSqOverY(Rectangle(0, test_image_width, 0, test_image_height), result2_float_y);
  if (!test_same_floats(result2_float_y, test_image_meansq_y, test_image_width)) {
    test_print_floats(result2_float_y, test_image_meansq_y, test_image_width);
    return 7;
  }
  for (int i = 0; i < 12; i++) result2_float_y[i] = (float)(sqrt((8.0/7)*(result2_float_y[i] - result_float_y[i]*result_float_y[i])));

  im.meanSqOverX(Rectangle(0, test_image_width, 0, test_image_height), result2_float_x);
  if (!test_same_floats(result2_float_x, test_image_meansq_x, test_image_height)) {
    test_print_floats(result2_float_x, test_image_meansq_x, test_image_height);
    return 8;
  }
  for (int i = 0; i < 8; i++) result2_float_x[i] = (float)(sqrt((12.0/11)*(result2_float_x[i] - result_float_x[i]*result_float_x[i])));

  im.deviationOverY(Rectangle(0, test_image_width, 0, test_image_height), test_image_mean_y, result_float_y);
  if (!test_same_floats(result_float_y, test_image_sd_y, test_image_width)) {
    test_print_floats(result_float_y, test_image_sd_y, test_image_width);
    return 9;
  }
  if (!test_same_floats(result_float_y, result2_float_y, test_image_width)) {
    test_print_floats(result_float_y, result2_float_y, test_image_width);
    return 10;
  }

  im.deviationOverX(Rectangle(0, test_image_width, 0, test_image_height), test_image_mean_x, result_float_x);
  if (!test_same_floats(result_float_x, test_image_sd_x, test_image_height)) {
    test_print_floats(result_float_x, test_image_sd_x, test_image_height);
    return 11;
  }
  if (!test_same_floats(result_float_x, result2_float_x, test_image_height)) {
    test_print_floats(result_float_x, result2_float_x, test_image_height);
    return 12;
  }

  uint8_t test_image8[] = {
    0x00, 0x00, 0x00, 0x00, 0x80, 0x80, 0x80, 0xA0,
    0x00, 0x00, 0x00, 0x00, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x00, 0x80, 0x80, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x10, 0x80, 0x80, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x20, 0x80, 0x80, 0x80, 0x80, 0xFF, 0xFF,
    0x00, 0x30, 0x80, 0x80, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x40, 0x80, 0x80, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x50, 0x80, 0x80, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x60, 0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x70, 0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF,
    0x00, 0x80, 0x80, 0xFF, 0xFF, 0x80, 0xFF, 0xFF,
    0x10, 0x80, 0xFF, 0xFF, 0xFF, 0x80, 0xFF, 0xFF
  };
  uint8_t test_image8_max_y[] = {0xA0, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
  uint8_t test_image8_max_x[] = {0x10, 0x80, 0xFF, 0xFF, 0xFF, 0x80, 0xFF, 0xFF};
  uint8_t test_image8_min_y[] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10};
  uint8_t test_image8_min_x[] = {0x00, 0x00, 0x00, 0x00, 0x80, 0x80, 0x80, 0xA0};
  uint8_t result_byte_y[12];
  uint8_t result_byte_x[8];

  Image8 im8(test_image8, Point(test_image_width, test_image_height), false);

  im8.maxOverY(Rectangle(0, test_image_width, 0, test_image_height), result_byte_y);
  if (!test_same_bytes(result_byte_y, test_image8_max_y, test_image_width)) {
    test_print_bytes(result_byte_y, test_image8_max_y, test_image_width);
    return 13;
  }

  im8.maxOverX(Rectangle(0, test_image_width, 0, test_image_height), result_byte_x);
  if (!test_same_bytes(result_byte_x, test_image8_max_x, test_image_height)) {
    test_print_bytes(result_byte_x, test_image8_max_x, test_image_height);
    return 14;
  }

  im8.minOverY(Rectangle(0, test_image_width, 0, test_image_height), result_byte_y);
  if (!test_same_bytes(result_byte_y, test_image8_min_y, test_image_width)) {
    test_print_bytes(result_byte_y, test_image8_min_y, test_image_width);
    return 15;
  }

  im8.minOverX(Rectangle(0, test_image_width, 0, test_image_height), result_byte_x);
  if (!test_same_bytes(result_byte_x, test_image8_min_x, test_image_height)) {
    test_print_bytes(result_byte_x, test_image8_min_x, test_image_height);
    return 16;
  }

  im8.meanOverY(Rectangle(0, test_image_width, 0, test_image_height), result_float_y);
  if (!test_same_floats(result_float_y, test_image_mean_y, test_image_width)) {
    test_print_floats(result_float_y, test_image_mean_y, test_image_width);
    return 17;
  }

  im8.meanOverX(Rectangle(0, test_image_width, 0, test_image_height), result_float_x);
  if (!test_same_floats(result_float_x, test_image_mean_x, test_image_height)) {
    test_print_floats(result_float_x, test_image_mean_x, test_image_height);
    return 18;
  }

  im8.meanSqOverY(Rectangle(0, test_image_width, 0, test_image_height), result2_float_y);
  if (!test_same_floats(result2_float_y, test_image_meansq_y, test_image_width)) {
    test_print_floats(result2_float_y, test_image_meansq_y, test_image_width);
    return 19;
  }
  for (int i = 0; i < 12; i++) result2_float_y[i] = (float)(sqrt((8.0/7)*(result2_float_y[i] - result_float_y[i]*result_float_y[i])));

  im8.meanSqOverX(Rectangle(0, test_image_width, 0, test_image_height), result2_float_x);
  if (!test_same_floats(result2_float_x, test_image_meansq_x, test_image_height)) {
    test_print_floats(result2_float_x, test_image_meansq_x, test_image_height);
    return 20;
  }
  for (int i = 0; i < 8; i++) result2_float_x[i] = (float)(sqrt((12.0/11)*(result2_float_x[i] - result_float_x[i]*result_float_x[i])));

  im8.deviationOverY(Rectangle(0, test_image_width, 0, test_image_height), test_image_mean_y, result_float_y);
  if (!test_same_floats(result_float_y, test_image_sd_y, test_image_width)) {
    test_print_floats(result_float_y, test_image_sd_y, test_image_width);
    return 21;
  }
  if (!test_same_floats(result_float_y, result2_float_y, test_image_width)) {
    test_print_floats(result_float_y, result2_float_y, test_image_width);
    return 22;
  }

  im8.deviationOverX(Rectangle(0, test_image_width, 0, test_image_height), test_image_mean_x, result_float_x);
  if (!test_same_floats(result_float_x, test_image_sd_x, test_image_height)) {
    test_print_floats(result_float_x, test_image_sd_x, test_image_height);
    return 23;
  }
  if (!test_same_floats(result_float_x, result2_float_x, test_image_height)) {
    test_print_floats(result_float_x, result2_float_x, test_image_height);
    return 24;
  }

  return 0;
}

int test_mwt_image()
{
  return 
    test_mwt_image_strip() + 
    10*test_mwt_image_mask() + 
    1000*test_mwt_image_contour() + 
    10000*test_mwt_image_image() +
    1000000*test_mwt_image_image8() +
    100000000*test_mwt_image_flattening();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_image();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Image test result is %d\n",i);
  return i>0;
}
#endif

