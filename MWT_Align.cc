/* Copyright (C) 2017 Rex A. Kerr and Calico Life Sciences
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
#include "MWT_Align.h"


/****************************************************************
                         Image Alignment
****************************************************************/

Profile::Profile(const Rectangle src, Collapse dir) : 
  reference(NULL), center(0), n_features(0), source(src), direction(dir), lateral_factor(0.0)
{
  n = (dir == OverX) ? source.height() : source.width();
  border = (n < 10) ? 1 : n/10;
  if (n > 0) reference = new float[2*n];
}

float Profile::find_mean(float *values, int count) {
  double m = 0.0;
  for (int i = 0; i < count; i++) m += values[i];
  if (count > 0) m /= count;
  return (float)m;
}

ShiftWeight Profile::histify(float *values, int count) {
  if (count <= 0) return ShiftWeight(0, 0);
  if (count == 1) return ShiftWeight(values[0], 0);
  memset(hist, 0, 64*sizeof(int));
  float lo = values[0];
  float hi = 0;
  float last = lo;
  double sumsq = 0;
  for (int i = 0; i < count; i++) {
    float x = values[i];
    if (x < lo) lo = x;
    else if (x > hi) hi = x;
    float dx = x - last;
    sumsq += dx*dx;
    last = x;
  }
  jitter = (float)sqrt(sumsq/(count - 1));
  float delta = (hi - lo)/64;
  float scale = 0.999/delta;
  float shift = lo-0.0005;
  for (int i = 0; i < count ; i++) {
    float x = (values[i] - shift)*scale;
    hist[(int)floor(x)]++;
  }
  return ShiftWeight(lo, delta);
}

float Profile::find_center_via_hist(float *values, int count, ShiftWeight sw) {
  int n05 = ceil(0.05*count);
  int n95 = n05;
  int i05 = 0; for (; i05 < 64 && n05 > 0; i05++) n05 -= hist[i05];
  int i95 = 63; for (; i95 >= 0 && n95 > 0; i95--) n95 -= hist[i95];
  if (i05 == i95) return find_mean(values, count);
  else {
    float x05 = sw.shift + i05*sw.weight;
    float x95 = sw.shift + i95*sw.weight;
    return 0.5*(x05 + x95);
  }
}

void Profile::load_best_features(int m) {
  n_features = 0;

  if (m > 8) m = 8;
  if (m < 1) m = 1;
  float x;

  // Find a value that is not close to the center
  // a is always the LEFT-most point
  int a = 0;
  while (a < n && (x = reference[a] - center, x > -jitter && x < jitter)) a += 1;
  if (a >= n) return;

  // Find a value that is not close to the center in the opposite direction
  // b is the MIDDLE or MIDDLE-RIGHT point of some feature
  int b = a+1;
  float u = reference[a] - center;
  while (b < n && (x = reference[b] - center, (u < 0) ? (x < jitter) : (x > -jitter))) {
    if (u < 0) { if (x < u) u = x; }
    else { if (x > u) u = x; }
    b += 1;
  }
  if (b >= n) return;

  // Back up to find the tightest crossing
  // bb is the MIDDLE-LEFT point of some feature, or unused if the feature has no width
  int bb = b-1;
  while (bb > a && (x = reference[bb] - center, (u < 0) ? (x > -jitter) : (x < jitter))) bb -= 1;

  // Now go forwards again to find a value that comes back on the other side again
  // c is always the RIGHT-most point
  int c = b+1;
  float v = reference[b] - center;
  while (c < n && (x = reference[c] - center, (v < 0) ? (x < jitter) : (x > -jitter))) {
    if (v < 0) { if (x < v) v = x; }
    else { if (x > v) v = x; }
    c += 1;
  }

  // Main loop where we walk through the data
  while (true) {
    // Place the previous completed feature into the list
    // Uses insertion sort; this is usually fastest for small lists like this
    Feature1D feature((b+bb)*0.5, fminf(fabsf(u), fabsf(v)), b - bb, (u<0) ? 1 : -1);
    if (n_features < m || features[n_features-1].strength < feature.strength) {
      int j = n_features;
      if (n_features < m) n_features += 1;
      else j -= 1;
      while (j > 0 && features[j-1].strength < feature.strength) {
        features[j] = features[j-1];
        j -= 1;
      }
      features[j] = feature;
    }

    // Now we need to advance to the next feature
    a = b;
    u = v;
    b = c;
    if (b >= n || ((reference[a] < center) == (reference[b] < center))) break;

    // Back up to find tightest crossing
    bb = b - 1;
    while (bb > a && (x = reference[bb] - center, (u < 0) ? (x > -jitter) : (x < jitter))) bb -= 1;

    // Forwards again to find another crossing
    c = b+1;
    v = reference[b] - center;
    while (c < n && (x = reference[c] - center, (v < 0) ? (x < jitter) : (x > -jitter))) {
      if (v < 0) { if (x < v) v = x; }
      else { if (x > v) v = x; }
      c += 1;
    }
  }

  // Now we need to refine the position of all these features with a straight line fit
  // It is NOT important that a straight line be a particularly great fit; it's just a decent way to average crossing positions
  for (int i = 0; i < n_features; i++) {
    bb = (int)rint(features[i].position - 0.5*features[i].diameter);
    b = (int)rint(bb + features[i].diameter);
    u = features[i].strength;
    auto inv_span = 1.0/features[i].diameter;
    auto slope = (reference[b] - reference[bb])*inv_span;
    float score = 0;

    // Extend backwards away from core to places where slope is still steep enough
    for (a = bb-1; a >= 0 && (x = reference[a] - center, score = (slope < 0) ? u-x : x+u, score > 0); a--) {
      auto a_slope = (reference[bb] - reference[a])*inv_span;
      if (slope < 0) { if (a_slope > slope) score -= u*(1 - a_slope/slope); }
      else           { if (a_slope < slope) score -= u*(1 - a_slope/slope); }
      if (score <= 0) break;
    }
    if (a < 0) a = 0;
    // Don't actually retreat from original value!  That's ridiculous!
    while (a < bb-1 && (x = reference[bb] - reference[a], (x < 0) != (slope < 0))) a++;

    // Extends forwards away from core to places where slope is still steep enough
    for (c = b+1; c < n && (x = reference[c] - center, score = (slope < 0) ? x+u : u-x, score > 0); c++) {
      auto c_slope = (reference[c] - reference[b])*inv_span;
      if (slope < 0) { if (c_slope > slope) score -= u*(1 - c_slope/slope); }
      else           { if (c_slope < slope) score -= u*(1 - c_slope/slope); }
      if (score <= 0) break;
    }
    if (c >= n) c = n-1;
    // Don't actually retreat from original value!  That's ridiculous!
    while (c > b+1 && (x = reference[c] - reference[b], (x < 0) != (slope < 0))) c--;

    // Asymmetry may pull the called position out of place, so make symmetric
    if (bb - a > c - b) a = bb - (c - b);
    else                c = b + (bb - a);

    // If we have long tails, broaden the original diameter (symmetrically)
    x = c - a - 2*(b - bb);
    if (x > 2) {
      auto delta = (int)rint(x*0.25);
      b += delta;
      bb -= delta;
    }

    // Finally get the straight line fit, weighted less heavily towards the ends to be less affected by variability in cutoff position
    // Nomenclature of straight line fit is x(t) = p + q*t, centered at features[i].position and center
    float t0 = features[i].position;
    double sum_w = 0.0;
    double sum_wt = 0.0;
    double sum_wx = 0.0;
    double sum_wtt = 0.0;
    double sum_wtx = 0.0;
    if (a < bb) {
      int inv_na = 1.0/(bb - a);
      for (int i = a; i < bb; i++) {
        x = reference[i] - center;
        auto t = i - t0;
        double w = inv_na*(bb - i - 0.5);
        sum_w += w;
        sum_wt += w*t;
        sum_wx += w*x;
        sum_wtt += w*t*t;
        sum_wtx += w*t*x;
      }
    }
    for (int i = bb; i <= b; i++) {
      x = reference[i] - center;
      auto t = i - t0;
      sum_w += 1;
      sum_wt += t;
      sum_wx += x;
      sum_wtt += t*t;
      sum_wtx += t*x;
    }
    if (c > b) {
      int inv_nc = 1.0/(c - b);
      for (int i = b+1; i <=c; i++) {
        x = reference[i] - center;
        auto t = i - t0;
        double w = inv_nc*(0.5 + c - i);
        sum_w += w;
        sum_wt += w*t;
        sum_wx += w*x;
        sum_wtt += w*t*t;
        sum_wtx += w*t*x;
      }
    }
    // Do NOT need to calculate q or p, just -p/q which is where x = 0!
    auto denom = sum_w*sum_wtx - sum_wt*sum_wx;
    float tx0 = (denom == 0) ? 0 : (float)((sum_wt*sum_wtx - sum_wx*sum_wtt)/denom);

    // Now a little basic sanity checking and we're good to go.
    if (fabsf(tx0)*2 < features[i].diameter+2) {
      features[i].position += tx0;
      u = features[i].position - bb;
      v = b - features[i].position;
      if (u < v) u = v;
      if (u < 1) u = 1;
      features[i].diameter = u;
    }
  }
}

void Profile::compute_new_values() {
  ShiftWeight sw = histify(reference, n);
  center = find_center_via_hist(reference, n, sw);
  load_best_features(6);
}

void Profile::adopt(float *profile, int length) {
  n = (direction == OverX) ? source.height() : source.width();
  if (n > length) n = length;
  memcpy(reference, profile, n*sizeof(float));
  compute_new_values();
}

void Profile::imprint(Image& frame) {
  Rectangle safe = source * frame.getBounds();
  if (direction == OverX) {
    n = safe.height();
    frame.meanOverX(safe, reference);
  }
  else {
    n = safe.width();
    frame.meanOverY(safe, reference);
  }
  compute_new_values();
}

void Profile::imprint8(Image8& frame) {
  Rectangle safe = source * frame.getBounds();
  if (direction == OverX) {
    n = safe.height();
    frame.meanOverX(safe, reference);
  }
  else {
    n = safe.width();
    frame.meanOverY(safe, reference);
  }
  compute_new_values();
}



int test_mwt_align_features() {
  float data[] = { 0.1, 0.1, 0.1, 0.1, 0.5, 1.1, 1.7, 2, 2, 2, 2, 2, 2, 2, 1.8, 0.5, 0, 0, 0, 0, 0, 0 };
  Profile prof(Rectangle(0, 21, 0, 21), Profile::OverX);
  prof.adopt(data, prof.source.width());
  if (prof.center < 0.8 || prof.center > 1.3) return 1;
  if (prof.n_features != 2) return 2;
  if (prof.features[0].strength < prof.features[1].strength) return 3;
  //if (prof.features[0].position < 14.4 || prof.features[0].position > 14.7) return 4;
  //if (prof.features[1].position < 4.8 || prof.features[1].position > 5.2) return 5;

  float data2[128]; for (int i = 0; i < 128; i++) data2[i] = (float)sin(0.1*i);
  Profile prof2(Rectangle(0, 127, 0, 15), Profile::OverY);
  prof2.adopt(data2, 128);
  if (prof2.n_features < 4) return 6;
  for (int i = 0; i < prof2.n_features; i++) {
    auto x = prof2.features[i];
    printf("#%d: %f %f %f  %d\n", i, x.position, x.strength, x.diameter, x.sign);
    auto npi = x.position/31.4159;
    auto err = fabs(npi - rint(npi));
    // Note--data set has slight positive bias so expect upward slopes to be slightly right-shifted and downward slopes slightly left-shifted
    if (err > 0.01) return 10+i;
  }

  for (int j = 1; j < 10; j++) {
    float dataj[] = { -1, -1, -0.5, 0, 0.5, 1, 1, 1, 1, 0.5, 0, -0.5, -1, -1 };
    auto dx = (j - 5)*0.05;
    auto zeroL = 3 - dx*2;
    auto zeroR = zeroL + 7;
    for (int i = 1; i <= 5; i++) {
      auto x = dataj[i] + dx;
      if (x < -1) x = -1;
      if (x > 1) x = 1;
      dataj[i] = x;
    }
    for (int i = 8; i <= 12; i++) {
      auto x = dataj[i] - dx;
      if (x < -1) x = -1;
      if (x > 1) x = 1;
      dataj[i] = x;
    }
    Profile profj(Rectangle(0, 7, 1, 14), Profile::OverX);
    profj.adopt(dataj, profj.source.height());
    if (profj.n_features != 2) return 20+j;
    printf("Zeros: %f %f\n", zeroL, zeroR);
    for (int i = 0; i < profj.n_features; i++) {
      auto x = profj.features[i];
      if (fabsf(x.position - 3) < 2) { if (fabsf(x.position - zeroL) > 0.2) return 30+j; }
      if (fabsf(x.position - 10) < 2) { if (fabsf(x.position - zeroR) > 0.2) return 40+j; }
    }
  }
  return 0;
}

int test_mwt_align()
{
  return test_mwt_align_features();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_align();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Align test result is %d\n",i);
  return i>0;
}
#endif

