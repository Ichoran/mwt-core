/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <new>

#include "MWT_Model.h"



/****************************************************************
                   ModelWorm Class Methods
****************************************************************/

// Set various worm parameters with some sanity checking
void ModelWorm::setPose(FPoint c,FPoint A,float faze)
{
  center = c;
  if (A.length()<1e-6) bearing = FPoint(1,0);
  else bearing = A.unit();
  phase = faze;
}
void ModelWorm::setSize(float len,float thick)
{
  length = (len<1.0) ? 1.0 : len;
  thickness = (thick<0.1) ? 0.1 : thick;
  if (thickness*2.0 > length) thickness = 0.5* length;
}
void ModelWorm::setWiggle(float ampl,float wavel,float fr,float flx)
{
  amplitude = (ampl<0.0) ? 0.0 : ampl;
  wavelen = (wavel<1.0) ? 1.0 : wavel;
  freq = (fr<0.0) ? 0.0 : fr;
  flex = flx;
}


// Test whether a point is inside the worm
// This might make some sense if you read what the values are
bool ModelWorm::isAt(FPoint where)
{
  FPoint p,q;
  FPoint local;
  
  p = where - center;
  local.x = p*bearing;
  local.y = p*bearing.ccw();
  if (fabs(local.x) >= 0.5*length) return false;
  local.y -= 0.5*flex*local.x*local.x;
  local.y -= amplitude*sin( phase + two_pi*local.x/wavelen );
  if (local.x > 0.3*length) { if ( fabs(local.y) < thickness*(0.5*length-local.x)/(0.2*length) ) return true; }
  else if (local.x < -0.2*length) { if ( fabs(local.y) < thickness*(0.5*length+local.x)/(0.3*length) ) return true; }
  else if (fabs(local.y) < thickness) return true;
  
  return false;
}


// Make a mask that covers the worm
// Just rectangular for now, can elaborate later
void ModelWorm::makeMask(Mask& m)
{
  Rectangle r = getBounds();
  m.flush();
  for (int x = r.near.x ; x<=r.far.x ; x++) m.addStrip(x,r.near.y,r.far.y);
  m.findBounds();
}


// A bounding rectangle about the worm (conservative)
Rectangle ModelWorm::getBounds()
{
  FPoint p = center;
  FPoint q = 0.5*bearing*length;
  if (q.x<0) q.x = -q.x;
  if (q.y<0) q.y = -q.y;
  float f = fabs(0.5*flex*q.x*q.x) + amplitude + thickness;      // Approximate bulk of worm (width of body + width of wiggle + width of circular curve)
  Rectangle r = Rectangle( (p-q).toPoint() , (p+q).toPoint() );  // Approximate endpoints of worm
  r.expand( 1+(int)ceil(f) );  // If we've got endpoints and expand by the bulk, we should be safe.
  return r;
}


// Move the worm forwards based on its frequency--assumes small changes (trig functions ~ linear or quadratic)
void ModelWorm::wiggle(float dt)
{
  phase += dt*(two_pi*freq);  // How much we should advance in phase
  
  // Find new position
  float f = dt*freq*wavelen;  // How far we should move
  FPoint p = f*bearing + 0.5*flex*f*f*bearing.ccw();  // Direction--we go in a circle because of "flex"
  p *= f/p.length();  // Right distance and direction
  center += p;
  
  // Find new heading--we're pointing in a slightly different direction because of flex
  p = f*flex*bearing.ccw(); 
  bearing += p;
  bearing /= bearing.length();
}


// Move the worm a longer time by dividing it up in to safely-short changes
void ModelWorm::longWiggle(float T)
{
  float f = (freq/two_pi)*wavelen*flex*T;
  
  if (f < 1e-2) wiggle(T);
  else
  {
    int n = (int)ceil(f/1e-2);
    for (int i=0;i<n;i++) wiggle(T/n);
  }
}


// Place a worm onto an image
void ModelWorm::imprint(Image& im,Mask &m,float darkness,int aa)
{
  short I;
  int y,n;
  Strip s;
  float fI;
  FPoint p,q;
  float delta = 1.0/(float)aa;
  float nmax = aa*aa;
  
  makeMask(m);
  m.cropTo(im.getBounds());
  m.start();
  while (m.advance())
  {
    s = m.i();
    for (y=s.y0 ; y<=s.y1 ; y++)
    {
      q = FPoint(s.x,y);
      n = 0;
      for (p.x = -0.5*(1-delta) ; p.x < 0.5 ; p.x += delta)
      {
        for (p.y = -0.5*(1-delta) ; p.y < 0.5 ; p.y += delta)
        {
          if (isAt(p+q)) n++;
        }
      }
      if (n>0)
      {
        fI = 1.0*(nmax-n)/nmax + darkness*n/nmax;
        I = im.get(s.x,y);
        I = (short)( I*fI );
        im.set(s.x,y,I);
      }
    }
  }
}



/****************************************************************
                   ModelWorm Class Methods
****************************************************************/


// Turn a perfect image into a camera-corrupted image.
void ModelCamera::imprint(Image& im)
{
  if (size != im.size) return;
  short I;
  float f;
  FPoint p;
  for (int i = 0;i<size.x*size.y;i++)
  {
    f = im.pixels[i];
    f += offsets[i];
    f *= scales[i];
    p = randomGaussians();
    f += p.x*fixed_noise + p.y*sqrt(f)*scaled_noise;
    I = (short)(f+0.5);
    //if (I<min) I=min;
    //else if (I>max) I=max;
    im.pixels[i] = I;
  }
}


// A bunch of functions to introduce various types of errors
void ModelCamera::scatterOffset(float deviation)
{
  for (int i=0;i<size.x+size.y;i++) offsets[i] += (short)(deviation*randomGaussian()+0.5);
}
void ModelCamera::rowsOffset(float deviation)
{
  short s;
  for (int i=0;i<size.y;i++)
  {
    s = (short)(deviation*randomGaussian()+0.5);
    for (int j=0;j<size.x;j++) offsets[j*size.y + i] += s;
  }
}
void ModelCamera::colsOffset(float deviation)
{
  short s;
  for (int i=0;i<size.x;i++)
  {
    s = (short)(deviation*randomGaussian() + 0.5);
    for (int j=0;j<size.y;j++) offsets[i*size.y + j] += s;
  }
}
void ModelCamera::scatterScale(float deviation)
{
  for (int i=0;i<size.x+size.y;i++) scales[i] *= randomNonNegativeShiftedGaussian(deviation);
}
void ModelCamera::rowsScale(float deviation)
{
  float f;
  for (int i=0;i<size.y;i++)
  {
    f = randomNonNegativeShiftedGaussian(deviation);
    for (int j=0;j<size.x;j++) scales[j*size.y + i] *= f;
  }
}
void ModelCamera::colsScale(float deviation)
{
  float f;
  for (int i=0;i<size.x;i++)
  {
    f = randomNonNegativeShiftedGaussian(deviation);
    for (int j=0;j<size.y;j++) scales[i*size.y + j] *= f;
  }
}


/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_model_worm()
{
  Image im( Point(64,64), false );
  im.depth = 11;
  Mask m(128);
  im = (im.getGray()*3)/2;
  
  ModelWorm worm;
  
  worm.center = FPoint(32,24);
  worm.bearing = FPoint(0.8,0.4); worm.bearing /= worm.bearing.length();
  worm.amplitude = 2;
  worm.thickness = 2;
  worm.length = 30;
  worm.wavelen = 20;
  worm.freq = 0.5;
  worm.phase = 0.2;
  worm.flex = -0.02;
  
  worm.imprint(im , m , 0.6 , 3);
  worm.wiggle(1.0);
  worm.imprint(im , m , 0.6 , 3);
  
  im <<= (15-im.depth);
  int i = im.writeTiff("worm_imprint.tiff");
  if (i!=0) return 1;
  
  return 0;
}

int test_mwt_model_camera()
{  
  Image im( Point(64,64), false );
  im.depth = 11;
  Mask m(128);
  im = (im.getGray()*3)/2;
  
  ModelWorm worm;
  
  worm.center = FPoint(32,24);
  worm.bearing = FPoint(0.8,0.4); worm.bearing /= worm.bearing.length();
  worm.amplitude = 2;
  worm.thickness = 2;
  worm.length = 30;
  worm.wavelen = 20;
  worm.freq = 0.5;
  worm.phase = 0.2;
  worm.flex = -0.02;
  
  worm.imprint(im , m , 0.6 , 3);
  worm.wiggle(1.0);
  worm.imprint(im , m , 0.6 , 3);

#ifdef WINDOWS
  srand(1491);
#else
  srandom(1491);  
#endif  
  ModelCamera mc(im);
  mc.scaled_noise = 1.0;
  mc.fixed_noise = 50;
  mc.scatterOffset(10);
  mc.rowsOffset(25);
  mc.scatterScale(0.05);
  mc.colsScale(0.1);
  mc.imprint(im);
  
  im <<= (15-im.depth);
  int i = im.writeTiff("worm_noisy.tiff");
  if (i!=0) return 1;
  
  return 0;
}

int test_mwt_model()
{
  return test_mwt_model_worm() + 100*test_mwt_model_camera();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_model();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Model test result is %d\n",i);
  return i>0;
}
#endif

