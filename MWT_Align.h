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


class Profile {
public:
  enum Collapse { OverX, OverY };

  float *reference;
  float *squareref;
  float *buffer;
  float center;
  float jitter;
  int* hist;
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
    if (squareref != NULL) {
      delete[] squareref;
      squareref = NULL;
    }
    if (buffer != NULL) { 
      delete[] buffer; 
      buffer = NULL;
    }
    if (hist != NULL) {
      delete[] hist;
      hist = NULL;
    }
  }

private:
  float find_mean(float *values, int count);
  ShiftWeight histify(float *values, int count);
  float find_center_via_hist(float *values, int count, ShiftWeight sw);
  void load_best_features(int m);
  void compute_new_values();

  float sort_and_size(Feature1D *data, int k, int *ix);  // Sorts by position, returns mean diameter

  void set_margins(int dist, Collapse dir, Rectangle &sub, Rectangle &add, int &absdist);
  Rectangle constrain_source(Rectangle frameBounds, Rectangle regionBounds);
  float best_tiled_inside(Image& frame, Rectangle search);
  float best_tiled_inside8(Image8& frame, Rectangle search);
  void constrain_bounds_nearby(Rectangle &bounds);
public:


  /** Quality score for how useful of a profile we have for alignment */
  float quality();

  /** Adopt an already-computed profile */
  void adopt(float *profile, int length);

  /** Adopt a profile from the given image */
  void imprint(Image& frame);

  /** Adopt a profile from the given 8 bit image */
  void imprint8(Image8& frame);

  /** Updates the profile to be shifted from the current position in the direction of collapse as indicated */
  void slide(Image& frame, int distance);

  /** Updates the profile to be shifted from the current position in the direction of collapse as indicated */
  void slide8(Image8& frame, int distance);

  /** Updates the profile to be shifted from the current position in the direction of variation as indicated */
  void scroll(Image& frame, int distance);

  /** Updates the profile to be shifted from the current position in the direction of variation as indicated */
  void scroll8(Image8& frame, int distance);

  /** Find and adopt the best quality within `search`, keeping the dimensions of the Profile.
    * The current position is updated to the best position.  Distances and widths are rounded to
    * a multiple of four.
    *
    * The Profile will be imprinted on the best position when this routine is complete.
    *
    * (Note: the search is not exhaustive.)
    */
  void best_inside(Image& frame, Rectangle search);

  /** Find and adopt the best quality within `search`, keeping the dimensions of the Profile.
    * The current position is updated to the best position.  Distances and widths are rounded to
    * a multiple of four.
    *
    * The Profile will be imprinted on the best position when this routine is complete. 
    *
    * (Note: the search is not exhaustive.)
    */
  void best_inside8(Image8& frame, Rectangle search);

  /** Given the profile in `probe`, find the offset between the stored profile and the probe profile. */
  float align(Profile* that);

  /** Find the offset between the image in `frame` and the stored profile in `*that` by computing and comparing the new profile. */
  inline float delta(Image& frame, Profile* that) { imprint(frame); return align(that); }

  /** Find the offset between the 8 bit image in `frame` and the stored profile in `*that` by computing and comparing the new profile. */
  inline float delta8(Image8& frame, Profile* that) { imprint8(frame); return align(that); }

  friend int test_mwt_align_best();
};


int test_mwt_align();


#endif