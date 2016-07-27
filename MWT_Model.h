/* Copyright 2007-2010 Rex A. Kerr, Nicholas A. Swierczek, and HHMI Janelia
 * Copyright 2016 Rex A. Kerr and Calico Life Sciences
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_MODEL
#define MWT_MODEL

#include <math.h>
#include "MWT_Geometry.h"
#include "MWT_Lists.h"
#include "MWT_Storage.h"
#include "MWT_Image.h"


/****************************************************************
         Classes to Model Noisy Images of Moving Worms
****************************************************************/

// A pretend worm that we can put onto an image
class ModelWorm
{
public:
  static constexpr float two_pi = 6.283185308;
  
  FPoint center;    // Center of the body
  FPoint bearing;   // Direction that the worm's head is pointing (unit vector)
  float amplitude;  // Amplitude of sinusoidal bend, in pixels
  float thickness;  // Width of worm in pixels (front 20% of animal tapers to a point, as does back 30%)
  float length;     // Length of worm in pixels (half in front of and half behind the center)
  float wavelen;    // Wavelength of worm in pixels
  float freq;       // Number of cycles of movement per (model) second
  float phase;      // Phase in the current cycle
  float flex;       // Inverse of radius of curvature, approximately; as you go x forward from center, worm shifts in y by 0.5*x*x*flex.
  
  ModelWorm() { }
  
  void setPose(FPoint c,FPoint A,float faze);
  void setSize(float len,float thick);
  void setWiggle(float ampl,float wavel,float fr,float flx);
  bool isAt(FPoint where);
  void makeMask(Mask& m);
  Rectangle getBounds();
  void wiggle(float dt);
  void longWiggle(float T);
  void imprint(Image& im,Mask &m,float darkness,int aa);
};


// Applies various camera artifacts to an image
class ModelCamera
{
public:
  bool is_saved;
  float saved_rn;
  float fixed_noise;
  float scaled_noise;
  
  short min,max;
  Rectangle bounds;
  Point size;
  short *offsets;
  float *scales;

  ModelCamera(const Image& im) : 
    is_saved(false),saved_rn(0),fixed_noise(0),scaled_noise(0),
    min(1),max(2*im.getGray()-1),bounds(im.getBounds())
  {
    int i;
    size = bounds.size();
    offsets = new short[ size.x*size.y ];
    scales = new float[ size.x*size.y ];
    for (i=0;i<size.x*size.y;i++)
    {
      offsets[i] = 0;
      scales[i] = 1.0;
    }
  }
  ~ModelCamera()
  {
    if (offsets!=NULL) { delete[] offsets; offsets=NULL; }
    if (scales!=NULL) { delete[] scales; scales=NULL; }
  }
  
  FPoint randomGaussians()  // Polar Box-Muller
  {
#ifdef WINDOWS
    FPoint f( rand() , rand() );
#else
    FPoint f( random() , random() );
#endif    
    f /= (float)RAND_MAX;
    float L = f.length2();
    while (L>=1.0 || L<=1e-14)
    {
#ifdef WINDOWS
      f.x = rand();
      f.y = rand();
#else
      f.x = random();
      f.y = random();
#endif      
      f /= (float)RAND_MAX;
      L = f.length2();
    }
    f *= sqrt( -2.0f * log(L) / L );
    return f;
  }
  float randomGaussian()
  {
    if (is_saved) { is_saved=false; return saved_rn; }
    else { FPoint f = randomGaussians(); saved_rn = f.y; is_saved=true; return f.x; }
  }
  float randomNonNegativeShiftedGaussian(float deviation)
  {
    float f = randomGaussian();
    f *= deviation;
    f += 1.0;
    return (f<0) ? 0.0 : f;
  }

  void imprint(Image& im);
  inline void resetProblems() { for (int i=0;i<size.x*size.y;i++) { offsets[i]=0; scales[i]=1.0; } }
  void scatterOffset(float deviation);
  void rowsOffset(float deviation);
  void colsOffset(float deviation);
  void scatterScale(float deviation);
  void rowsScale(float deviation);
  void colsScale(float deviation);
};



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_model_worm();
int test_mwt_model_camera();
int test_mwt_model();

#endif

