/* Copyright (C) 2017 Rex A. Kerr and Calico Life Sciences
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */

#ifndef MWT_ALIGN
#define MWT_ALIGN

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "MWT_Geometry.h"
#include "MWT_Lists.h"
#include "MWT_Storage.h"
#include "MWT_Image.h"



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

  /** Find and adopt the best quality within `search` distance plus or minus the current position.
    * The current position is updated to the best position.
    *
    * (Note: search is not exhaustive.)
    */
  void best_near(Image& frame, Point search);

  /** Find and adopt the best quality within `search` distance plus or minus the current position.
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


int test_mwt_align();


#endif