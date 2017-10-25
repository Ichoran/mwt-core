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
  reference(NULL), squareref(NULL), buffer(NULL), center(0), hist(NULL), n_features(0), source(src), direction(dir), lateral_factor(0.0)
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
  if (hist == NULL) hist = new int[64*sizeof(int)];
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
    if (x < 0) x = 0;
    hist[0x3F & (int)floor(x)]++;
  }
  return ShiftWeight(lo, delta);
}

float Profile::find_center_via_hist(float *values, int count, ShiftWeight sw) {
  if (hist == NULL) hist = new int[64*sizeof(int)];   // No data but at least not a segfault
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
    if (squareref != NULL) frame.meanSqOverX(safe, squareref);
  }
  else {
    n = safe.width();
    frame.meanOverY(safe, reference);
    if (squareref != NULL) frame.meanSqOverY(safe, squareref);
  }
  compute_new_values();
}

void Profile::imprint8(Image8& frame) {
  Rectangle safe = source * frame.getBounds();
  if (direction == OverX) {
    n = safe.height();
    frame.meanOverX(safe, reference);
    if (squareref != NULL) frame.meanSqOverX(safe, squareref);
  }
  else {
    n = safe.width();
    frame.meanOverY(safe, reference);
    if (squareref != NULL) frame.meanSqOverY(safe, squareref);
  }
  compute_new_values();
}

float Profile::sort_and_size(Feature1D *data, int k, int *ix) {
  if (k < 1) return 0.0;
  float sumw = data[0].diameter;
  for (int i = 0; i < k; i++) ix[i] = i;
  for (int i = 1; i < k; i++) {
    sumw += data[i].diameter;
    for (int j = i; j > 0 && data[ix[j]].position < data[ix[j-1]].position; j--) {
      auto temp = ix[j];
      ix[j] = ix[j-1];
      ix[j-1] = temp;
    }
  }
  return sumw;
}

float Profile::align(Profile* that) {
  if (n_features < 1) return NAN;
  if (that->n_features < 1) return NAN;
  if (n_features == 1) {
    float delta = features[0].position - that->features[0].position;
    for (int i = 1; i < that->n_features; i++) {
      float x = features[0].position - that->features[i].position;
      if (fabsf(x) < fabsf(delta)) { delta = x; }
    }
    return delta;
  }
  else if (that->n_features == 1) {
    float delta = features[0].position - that->features[0].position;
    for (int i = 1; i < n_features; i++) {
      float x = features[i].position - that->features[0].position;
      if (fabsf(x) < fabsf(delta)) delta = x;
    }
    return delta;
  }
  else {
    int i;

    // Need more extensive search if we're not within the feature widths
    int thisix[8];
    int thatix[8];
    // Insertion sort to order features by position (using indexing arrays)
    float meanw = 
      sort_and_size(this->features, this->n_features, thisix) +
      sort_and_size(that->features, that->n_features, thatix);
    meanw /= this->n_features + that->n_features;

    // Precompute spacings; we want these to match reasonably well
    float thisdx[7]; for (i = 1; i < this->n_features; i++) thisdx[i-1] = this->features[thisix[i]].position - this->features[thisix[i-1]].position;
    float thatdx[7]; for (i = 1; i < that->n_features; i++) thatdx[i-1] = that->features[thatix[i]].position - that->features[thatix[i-1]].position;
    float delta[8];
    int thisn = n_features-1;
    int thatn = that->n_features-1;
    float meandelta;
    if (thisn == thatn) {
      // Maybe they're all the same!
      meandelta = 0.0;
      for (i = 0; i < thisn; i++) delta[i] = thisdx[i] - thatdx[i];
      for (i = 0; i < thisn; i++) meandelta += fabsf(delta[i]);
      meandelta /= thisn;
      bool concordant = true;
      for (i = 0; concordant && i < thisn; i++) concordant = fabsf(delta[i]) <= meanw;
      if (concordant) {
        float dpos = 0.0;
        for (i = 0; i < n_features; i++) dpos += features[thisix[i]].position - that->features[thatix[i]].position;
        dpos /= n_features;
        return dpos;
      }
    }
    else meandelta = fabsf(features[0].position) + fabsf(features[n_features-1].position) + fabsf(that->features[0].position) + 1;

    // If we got here, then they're not all the same and we need to try dropping mismatching features
    Profile *more = (n_features < that->n_features) ? that : this;
    Profile *less = (n_features < that->n_features) ? this : that;
    float *moredx = (n_features < that->n_features) ? thatdx : thisdx;
    float *lessdx = (n_features < that->n_features) ? thisdx : thatdx;
    int moren = (thisn < thatn) ? thatn : thisn;
    int lessn = (thisn < thatn) ? thisn : thatn;
    if (lessn+2 >= moren) {
      float bestdelta = meandelta;
      int besti = -1;
      int bestj = -1;
      bool two = moren == lessn+2;
      bool any_concordant = false;
      for (int skippi = 0; skippi <= moren; skippi++) {
        for (int skippj = 0; skippj <= moren; skippj++) {
          if (two && skippj <= skippi) continue;
          if (moren == lessn+1) skippj = moren+1;
          int il = (skippi == 0) ? ((two && skippj == 1) ? 2 : 1) : 0;
          int ir = il+1; if (ir==skippi) ir++; if (two && ir==skippj) ir++;
          int jl = (!two && skippj == 0) ? 1 : 0;
          int jr = jl+1; if (!two && jr==skippj) jr++;
          int ndelta = 0;
          meandelta = 0.0;
          bool concordant = true;
          while (concordant && ir <= moren && jr <= lessn) {
            ndelta += 1;
            float di = moredx[il++]; for (; il < ir; il++) di += moredx[il];
            float dj = lessdx[jl++]; for (; jl < jr; jl++) dj += lessdx[jl];
            float d = fabs(di-dj);
            concordant = concordant && d <= meanw;
            meandelta += d;
            il = ir; ir += 1; if (ir==skippi) ir++; if (two && ir==skippj) ir++;
            jl = jr; jr += 1; if (!two && jr==skippj) jr++;
          }
          if (ndelta > 0 && concordant) {
            meandelta /= ndelta;
            if (meandelta < bestdelta) {
              bestdelta = meandelta;
              besti = skippi;
              bestj = skippj;
              any_concordant = true;
            }
          }
        }
      }
      if (any_concordant) {
        float dpos = 0.0;
        i = (besti == 0) ? ((two && bestj == 1) ? 2 : 1) : 0;
        int j = (!two && bestj==0) ? 1 : 0;
        int np = 0;
        int *moreix = (n_features < that->n_features) ? thatix : thisix;
        int *lessix = (n_features < that->n_features) ? thisix : thatix;
        while (i < more->n_features && j < less->n_features) {
          dpos += more->features[moreix[i]].position - less->features[lessix[j]].position;
          np += 1;
          i += 1; if (i == besti) i++; if (two && i == bestj) i++;
          j += 1; if (!two && j == bestj) j++;
        }
        dpos /= np;
        return (n_features < that->n_features) ? -dpos : dpos;
      }
    }

    // If we reach this point, we couldn't match the pattern at all!
    // So we'll just try matching the closest feature and hoping for the best.
    float best = 2*(features[0].position - that->features[0].position);
    for (i = 0; i < n_features; i++) {
      for (int j = 0; j < that->n_features; j++) {
        float p = features[i].position - that->features[j].position;
        if (fabsf(p) < fabsf(best)) best = p;
      }
    }
    return best;
  }
}

float Profile::quality() {
  if (n_features < 1) return 0;
  double q = 0.0;
  float c = 0.5*(n-1);
  for (int i = 0; i < n_features; i++) {
    Feature1D& fi = features[i];
    float centrality = 1.0 - fabsf(c - fi.position)/c;
    float strength = fi.strength;
    if (squareref != NULL) {
      int j0 = (int)floor(fi.position - fi.diameter); if (j0 < 0)  j0 = 0;
      int j1 = (int) ceil(fi.position + fi.diameter); if (j1 >= n) j1 = n-1;
      float svar = 0.0;
      for (int j = j0; j <= j1; j++) svar += squareref[j] - reference[j]*reference[j];
      svar /= 1 + j1 - j0;
      if (svar > 0 && svar*6 > fi.strength) {
        if (svar > strength) fi.strength *= 0.25;
        else strength *= 0.25 + (fi.strength/svar - 1)*0.15;
      }
    }
    q += strength*((centrality > 0.5) ? 1 : (3*centrality - 0.5));
  }
  return q;
}


void Profile::set_margins(int dist, Collapse dir, Rectangle &sub, Rectangle &add, int &absdist) {
  if (dist > 0) {
    absdist = dist;
    if (direction == dir) {
      add = Rectangle(source.far.x + 1, source.far.x + dist, source.near.y, source.far.y);
      sub = Rectangle(source.near.x, source.near.x + dist - 1, source.near.y, source.far.y);
    }
    else {
      add = Rectangle(source.near.x, source.far.x, source.far.y + 1, source.far.y + dist);
      sub = Rectangle(source.near.x, source.far.x, source.near.y, source.near.y + dist - 1);
    }
  }
  else {
    absdist = -dist;
    if (direction == dir) {
      add = Rectangle(source.near.x + dist, source.near.x - 1, source.near.y, source.far.y);
      sub = Rectangle(source.far.x + dist + 1, source.far.x, source.near.y, source.far.y);
    }
    else {
      add = Rectangle(source.near.x, source.far.x, source.near.y + dist, source.near.y - 1);
      sub = Rectangle(source.near.x, source.far.x, source.far.y + dist + 1, source.far.y);
    }
  }
}

void Profile::slide(Image &frame, int distance) {
  if (n*0.4 <= distance || distance == 0) {
    source += (direction == OverX) ? Point(distance, 0) : Point(0, distance);
    imprint(frame);
    return;
  }
  if (buffer == NULL) buffer = new float[n*2];
  Rectangle add, sub;
  int absdist;
  set_margins(distance, OverX, sub, add, absdist);
  if (direction==OverY) printf("%d to %d; %d to %d\n", add.near.x, add.far.x, add.near.y, add.far.y);
  sub = sub * frame.getBounds();
  add = add * frame.getBounds();
  float mult = absdist / (float)((direction == OverX) ? source.width() : source.height());

  if (direction == OverX) frame.meanOverX(sub, buffer);   else frame.meanOverY(sub, buffer);
  if (direction == OverX) frame.meanOverX(add, buffer+n); else frame.meanOverY(add, buffer+n);
  for (int i = 0; i < n; i++) reference[i] += mult*(buffer[n+i] - buffer[i]);

  if (squareref != NULL) {
    if (direction == OverX) frame.meanSqOverX(sub, buffer);   else frame.meanSqOverY(sub, buffer);
    if (direction == OverX) frame.meanSqOverX(add, buffer+n); else frame.meanSqOverY(add, buffer+n);
    for (int i = 0; i < n; i++) squareref[i] += mult*(buffer[n+i] - buffer[i]);
  }

  source += (direction == OverX) ? Point(distance, 0) : Point(0, distance);
  compute_new_values();
}

void Profile::slide8(Image8 &frame, int distance) {
  if (n*0.4 <= distance || distance == 0) {
    source += (direction == OverX) ? Point(distance, 0) : Point(0, distance);
    imprint8(frame);
    return;
  }
  if (buffer == NULL) buffer = new float[n*2];
  Rectangle add, sub;
  int absdist;
  set_margins(distance, OverX, sub, add, absdist);
  sub = sub * frame.getBounds();
  add = add * frame.getBounds();
  float mult = absdist / (float)((direction == OverX) ? source.width() : source.height());

  if (direction == OverX) frame.meanOverX(sub, buffer);   else frame.meanOverY(sub, buffer);
  if (direction == OverX) frame.meanOverX(add, buffer+n); else frame.meanOverY(add, buffer+n);
  for (int i = 0; i < n; i++) reference[i] += mult*(buffer[n+i] - buffer[i]);

  if (squareref != NULL) {
    if (direction == OverX) frame.meanSqOverX(sub, buffer);   else frame.meanSqOverY(sub, buffer);
    if (direction == OverX) frame.meanSqOverX(add, buffer+n); else frame.meanSqOverY(add, buffer+n);
    for (int i = 0; i < n; i++) squareref[i] += mult*(buffer[n+i] - buffer[i]);
  }

  source += (direction == OverX) ? Point(distance, 0) : Point(0, distance);
  compute_new_values();
}

void Profile::scroll(Image &frame, int distance) {
  if (n <= distance || distance == 0) {
    source += (direction == OverY) ? Point(distance, 0) : Point(0, distance);
    imprint(frame);
    return;
  }
  if (buffer == NULL) buffer = new float[n*2];
  Rectangle add, sub;
  int absdist;
  set_margins(distance, OverY, sub, add, absdist);
  add = add * frame.getBounds();

  if (direction == OverX) frame.meanOverX(add, buffer); else frame.meanOverY(add, buffer);
  if (distance > 0) {
    memmove(reference, reference+absdist, (n-absdist)*sizeof(float));
    memcpy(reference+n-absdist, buffer, absdist*sizeof(float));
  }
  else {
    memmove(reference+absdist, reference, (n-absdist)*sizeof(float));
    memcpy(reference, buffer, absdist*sizeof(float));
  }

  if (squareref != NULL) {
    if (direction == OverX) frame.meanSqOverX(add, buffer); else frame.meanSqOverY(add, buffer);
    if (distance > 0) {
      memmove(squareref, squareref+absdist, (n-absdist)*sizeof(float));
      memcpy(squareref+n-absdist, buffer, absdist*sizeof(float));
    }
    else {
      memmove(squareref+absdist, squareref, (n-absdist)*sizeof(float));
      memcpy(squareref, buffer, absdist*sizeof(float));
    }
  }

  source += (direction == OverY) ? Point(distance, 0) : Point(0, distance);
  compute_new_values();
}

void Profile::scroll8(Image8 &frame, int distance) {
  if (n <= distance || distance == 0) {
    source += (direction == OverY) ? Point(distance, 0) : Point(0, distance);
    imprint8(frame);
    return;
  }
  if (buffer == NULL) buffer = new float[n*2];
  Rectangle add, sub;
  int absdist;
  set_margins(distance, OverY, sub, add, absdist);
  add = add * frame.getBounds();

  if (direction == OverX) frame.meanOverX(add, buffer); else frame.meanOverY(add, buffer);
  if (distance > 0) {
    memmove(reference, reference+absdist, (n-absdist)*sizeof(float));
    memcpy(reference+n-absdist, buffer, absdist*sizeof(float));
  }
  else {
    memmove(reference+absdist, reference, (n-absdist)*sizeof(float));
    memcpy(reference, buffer, absdist*sizeof(float));
  }

  if (squareref != NULL) {
    if (direction == OverX) frame.meanSqOverX(add, buffer); else frame.meanSqOverY(add, buffer);
    if (distance > 0) {
      memmove(squareref, squareref+absdist, (n-absdist)*sizeof(float));
      memcpy(squareref+n-absdist, buffer, absdist*sizeof(float));
    }
    else {
      memmove(squareref+absdist, squareref, (n-absdist)*sizeof(float));
      memcpy(squareref, buffer, absdist*sizeof(float));
    }
  }

  source += (direction == OverY) ? Point(distance, 0) : Point(0, distance);
  compute_new_values();
}


Rectangle Profile::constrain_source(Rectangle frameBounds, Rectangle regionBounds) {
  Rectangle actual = frameBounds * regionBounds;
  if (actual.width() < source.width()) {
    source.near.x = actual.near.x;
    source.far.x  = actual.far.x;
  }
  if (actual.height() < source.height()) {
    source.near.y = actual.near.y;
    source.far.y  = actual.far.y;
  }
  if (actual.near.x > source.near.x)    source += Point(actual.near.x - source.near.x, 0);
  else if (actual.far.x < source.far.x) source -= Point(source.far.x - actual.far.x, 0);
  if (actual.near.y > source.near.y)    source += Point(0, actual.near.y - source.near.y);
  else if (actual.far.y < source.far.y) source -= Point(0, source.far.y - actual.far.y);
  int m = (direction == OverX) ? source.height() : source.width();
  if (m != n) {
    if (n > m) {
      if (reference != NULL) { delete[] reference; reference = new float[m]; }
      if (squareref != NULL) { delete[] squareref; squareref = new float[m]; }
      if (buffer != NULL) { delete[] buffer; buffer = new float[2*m]; }
    }
    n = m;
  }
  return actual;
}

float Profile::best_tiled_inside(Image& frame, Rectangle bounds) {
  bool dealloc = false;
  if (squareref == NULL) {
    squareref = new float[n];
    dealloc = true
;  }
  auto x0 = bounds.near.x; if ((x0 & 0x3) != 0) x0 += 4 - (x0 & 0x3);
  auto y0 = bounds.near.y; if ((y0 & 0x3) != 0) y0 += 4 - (y0 & 0x3);
  source += Point(x0, y0) - source.near;
  bool no_x = source.far.x > bounds.far.x;
  bool no_y = source.far.y > bounds.far.y;
  if (no_x) source += Point(bounds.near.x - x0, 0);
  if (no_y) source += Point(0, bounds.near.y - y0);
  auto corner = source.near;
  auto best_p = source.near;
  float best_q = 0;
  auto halfdelta = (source.size() / 8); halfdelta *= 4;
  if (halfdelta.x <= 1) halfdelta.x = 2;
  if (halfdelta.y <= 1) halfdelta.y = 2;
  auto delta = halfdelta * 2;
  for (auto c = corner; source.nearTo(c), source.far.x <= bounds.far.x; c.x += delta.x) {
    for (c.y = corner.y; source.nearTo(c), source.far.y <= bounds.far.y; c.y += delta.y) {
      imprint(frame);
      auto q = quality();
      printf("  %2d,%2d  %10.3f %2d %5.2f %7.2f\n", source.near.x, source.near.y, q, n_features, (n_features > 0) ? features[0].position : 0.0, (n_features > 0) ? features[0].strength : 0.0);
      if (q > best_q) {
        best_p = source.near;
        best_q = q;
      }
    }
    auto cc = Point(c.x, corner.y);
    cc += halfdelta;
    if (cc.x + source.width() <= bounds.far.x && cc.y + source.height() <= bounds.far.y) {
      for (; source.nearTo(cc), source.far.y <= bounds.far.y; cc.y += delta.y) {
        imprint(frame);
        auto q = quality();
        printf("  %2d,%2d  %10.3f %2d %5.2f %7.2f\n", source.near.x, source.near.y, q, n_features, (n_features > 0) ? features[0].position : 0.0, (n_features > 0) ? features[0].strength : 0.0);
        if (q > best_q) {
          best_p = source.near;
          best_q = q;
        }
      }
    }
  }
  source.nearTo(best_p);
  imprint(frame);
  if (dealloc) {
    delete[](squareref);
    squareref = NULL;
  }
  return best_q;
}

float Profile::best_tiled_inside8(Image8& frame, Rectangle bounds) {
  bool dealloc = false;
  if (squareref == NULL) {
    squareref = new float[n];
    dealloc = true
;  }
  auto x0 = bounds.near.x; if ((x0 & 0x3) != 0) x0 += 4 - (x0 & 0x3);
  auto y0 = bounds.near.y; if ((y0 & 0x3) != 0) y0 += 4 - (y0 & 0x3);
  source += Point(x0, y0) - source.near;
  bool no_x = source.far.x > bounds.far.x;
  bool no_y = source.far.y > bounds.far.y;
  if (no_x) source += Point(bounds.near.x - x0, 0);
  if (no_y) source += Point(0, bounds.near.y - y0);
  auto corner = source.near;
  auto best_p = source.near;
  float best_q = 0;
  auto halfdelta = (source.size() / 8); halfdelta *= 4;
  if (halfdelta.x <= 1) halfdelta.x = 2;
  if (halfdelta.y <= 1) halfdelta.y = 2;
  auto delta = halfdelta * 2;
  for (auto c = corner; source.nearTo(c), source.far.x <= bounds.far.x; c.x += delta.x) {
    for (c.y = corner.y; source.nearTo(c), source.far.y <= bounds.far.y; c.y += delta.y) {
      imprint8(frame);
      auto q = quality();
      if (q > best_q) {
        best_p = source.near;
        best_q = q;
      }
    }
    auto cc = Point(c.x, corner.y);
    cc += halfdelta;
    if (cc.x + source.width() <= bounds.far.x && cc.y + source.height() <= bounds.far.y) {
      for (; source.nearTo(cc), source.far.y <= bounds.far.y; cc.y += delta.y) {
        imprint8(frame);
        auto q = quality();
        if (q > best_q) {
          best_p = source.near;
          best_q = q;
        }
      }
    }
  }
  source.nearTo(best_p);
  imprint8(frame);
  if (dealloc) {
    delete[](squareref);
    squareref = NULL;
  }
  return best_q;
}


void Profile::constrain_bounds_nearby(Rectangle &bounds) {
  auto w = source.width(); w -= w % 4;
  auto h = source.height(); h -= h % 4;
  if (bounds.near.x < source.near.x - w) bounds.near.x = source.near.x - w;
  if (bounds.near.y < source.near.y - h) bounds.near.y = source.near.y - h;
  if (bounds.far.x > source.far.x + w)  bounds.far.x = source.far.x + w;
  if (bounds.far.y > source.far.y + source.height()) bounds.far.y = source.far.y + h;
  int x, y;
  if ((x = bounds.near.x & 0x3) != 0 && bounds.near.x + x <= source.near.x) bounds.near.x += x;
  if ((y = bounds.near.y & 0x3) != 0 && bounds.near.y + y <= source.near.y) bounds.near.y += y;
  if ((x = bounds.far.x & 0x3) != 0 && bounds.far.x - x >= source.far.x) bounds.far.x -= x;
  if ((y = bounds.far.y & 0x3) != 0 && bounds.far.y - y >= source.far.y) bounds.far.y -= y;
}


/*
void Profile::best_near(Image &frame, Rectangle search) {
  Rectangle bounds = constrain_source(frame.getBounds(), search);
  if (
    bounds.width() >= source.width() &&
    bounds.height() >= source.height() &&
    (bounds.width()/source.width())*(bounds.height()/source.height()) > 3
  ) {
    best_tile_near(frame, bounds);
    constrain_bounds_nearby(bounds);
  }
  int xlo = bounds.near.x - source.near.x;
  int xhi = bounds.far.x - source.far.x;
  int ylo = bounds.near.y - source.near.y;
  int yhi = bounds.far.y - source.far.y;
  int xbest = 0;
  int ybest = 0;
  float qbest = 0.0;
  if (squareref != NULL) squareref = new float[n];
  if (direction == OverX) {

    imprint(frame);
  }
}
*/


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

int test_mwt_align_align() {
  Rectangle r(0, 21, 0, 21);
  auto dir = Profile::OverX;

  Profile one_a(r, dir);
  one_a.n_features = 1;
  one_a.features[0] = Feature1D(10.7, 1, 5, 1);
  Profile one_b(r, dir);
  one_b.n_features = 1;
  one_b.features[0] = Feature1D(11.3, 1, 5, 1);
  if (!(fabsf(one_a.align(&one_b) + 0.6) < 0.01)) return 1;

  Profile two_a(r, dir);
  two_a.n_features = 2;
  two_a.features[0] = Feature1D(10.7, 1, 5, 1);
  two_a.features[1] = Feature1D(7.4, 0.8, 5, -1);
  Profile two_b(r, dir);
  two_b.n_features = 2;
  two_b.features[0] = Feature1D(11.3, 1, 5, 1);
  two_b.features[1] = Feature1D(7.99, 0.8, 5, -1);
  if (!(fabsf(two_a.align(&one_a)) < 0.01)) return 2;
  if (!(fabsf(one_a.align(&two_a)) < 0.01)) return 3;
  if (!(fabsf(two_a.align(&one_b) - one_a.align(&one_b)) < 0.01)) return 4;
  if (!(fabsf(one_a.align(&two_b) - one_a.align(&one_b)) < 0.01)) return 5;

  Profile three_a(r, dir);
  three_a.n_features = 3;
  three_a.features[0] = Feature1D(10.7, 1, 5, 1);
  three_a.features[1] = Feature1D(11.6, 0.9, 5, -1);
  three_a.features[2] = Feature1D(7.4, 0.8, 5, -1);
  if (!(fabsf(three_a.align(&two_b) - one_a.align(&one_b)) < 0.01)) return 6;
  if (!(fabsf(two_b.align(&three_a) + one_a.align(&one_b)) < 0.01)) return 7;

  Profile three_b(r, dir);
  three_b.n_features = 3;
  three_b.features[0] = Feature1D(11.3, 1, 5, 1);
  three_b.features[1] = Feature1D(12.21, 0.9, 5, -1);
  three_b.features[2] = Feature1D(19.7, 0.8, 3, 1);
  if (!(fabsf(three_a.align(&three_b) - one_a.align(&one_b)) < 0.01)) return 8;
  if (!(fabsf(three_b.align(&three_a) - one_b.align(&one_a)) < 0.01)) return 9;

  Profile four_a(r, dir);
  four_a.n_features = 4;
  four_a.features[0] = Feature1D(3.8, 1.1, 5, 1);
  four_a.features[1] = Feature1D(10.7, 1, 5, 1);
  four_a.features[2] = Feature1D(11.6, 0.9, 5, -1);
  four_a.features[3] = Feature1D(7.4, 0.8, 5, -1);
  if (!(fabsf(four_a.align(&two_b) - one_a.align(&one_b)) < 0.01)) return 10;
  if (!(fabsf(two_b.align(&four_a) + one_a.align(&one_b)) < 0.01)) return 11;

  Profile five_a(r, dir);
  five_a.n_features = 5;
  five_a.features[0] = Feature1D(10.7, 1, 5, 1);
  five_a.features[1] = Feature1D(13.9, 0.9, 5, -1);
  five_a.features[2] = Feature1D(16.4, 0.8, 5, 1);
  five_a.features[3] = Feature1D(7.41, 0.7, 5, -1);
  five_a.features[4] = Feature1D(1.5, 0.6, 3, 1);
  if (!(fabsf(five_a.align(&two_b) - one_a.align(&one_b)) < 0.1)) return 10;
  if (!(fabsf(two_b.align(&five_a) + one_a.align(&one_b)) < 0.1)) return 11;

  Profile zero_a(r, dir);
  zero_a.n_features = 0;

  if (five_a.quality() >= four_a.quality()) return 12;
  if (three_a.quality() >= four_a.quality()) return 13;
  if (three_b.quality() >= three_a.quality()) return 14;
  if (two_a.quality() >= three_a.quality()) return 15;
  if (one_a.quality() >= two_a.quality()) return 16;
  if (zero_a.quality() >= one_a.quality()) return 17;
  if (one_b.quality() != one_a.quality()) return 18;

  return 0;
}

int test_mwt_align_imprint() {
  char data[16][21] = {
   /*01234567891123456789*/
    "       !AWYZYWSYZYWZ",
    "       !AXZXWYWXYWXZ",
    "       @NZWYZYWYZWYZ",
    "        @NZYWWYZWYWZ",
    "       !AWYZYWSYZYWZ",
    "       !AXZXWYWXYWXZ",
    "       @NZWYZYWYZWYZ",
    "        @NZYWWYZWYWZ",
    "       !AWYZYWSYZYWZ",
    "       !AXZXWYWXYWXZ",
    "       @NZWYZYWYZWYZ",
    "        @NZYWWYZWYWZ",
    "       !AWYZYWSYZYWZ",
    "       !AXZXWYWXYWXZ",
    "       @NZWYZYWYZWYZ",
    "        @NZYWWYZWYWZ"
  };
  short *data_i16 = new short[20*16*sizeof(short)];
  uint8_t *data_u8 = new uint8_t[20*16*sizeof(uint8_t)];
  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 20; j++) {
      data_i16[i*20+j] = (short)data[i][j] * 64;
      data_u8[i*20+j] = (uint8_t)data[i][j];
    }
  }

  Image full_i16(data_i16, Point(16, 20), false); full_i16.owns_pixels = true;
  Image8 full_u8(data_u8, Point(16, 20), false); full_u8.owns_pixels = true;

  Image test_a_i16(full_i16, Rectangle(2, 13, 0, 11), false);
  Image8 test_a_u8(full_u8, Rectangle(2, 13, 0, 11), false);
  Profile p_a_i16(Rectangle(2, 9, 0, 11), Profile::OverX);
  Profile p_a_u8(Rectangle(2, 9, 0, 11), Profile::OverX);
  p_a_i16.imprint(test_a_i16);
  auto delta_a1 = p_a_u8.delta8(test_a_u8, &p_a_i16);
  auto delta_a2 = p_a_i16.delta(test_a_i16, &p_a_u8);
  if (fabsf(delta_a1) > 0.01) return 1;
  if (fabsf(delta_a2) > 0.01) return 2;
  if (fabsf(delta_a1 + delta_a2) > 0.001) return 3;

  Image test_b_i16(full_i16, Rectangle(2, 13, 2, 13), false);
  Profile p_b_i16(Rectangle(2, 9, 0, 11), Profile::OverX);
  auto delta_ab = p_b_i16.delta(test_b_i16, &p_a_i16);
  if (fabsf(delta_ab + 2) > 0.01) { printf("%f\n", delta_ab); return 4; }

  Image8 test_c_u8(full_u8, Rectangle(4, 15, 0, 11), false);
  Profile p_c_u8a(Rectangle(2, 9, 0, 11), Profile::OverX);
  Profile p_c_u8b(Rectangle(3, 10, 0, 11), Profile::OverX);
  auto delta_aca = p_c_u8a.delta8(test_c_u8, &p_a_u8);
  auto delta_cab = p_c_u8b.delta8(test_c_u8, &p_c_u8a);
  if (fabsf(delta_aca) > 0.1) { printf("%f\n", delta_aca); return 5; }
  if (fabsf(delta_cab) > 0.1) { printf("%f\n", delta_cab); return 6; }

  return 0;
}

int test_mwt_align_slide() {
  char data[16][21] = {
   /*01234567891123456789*/
    "        PPPPPPPPPPPP",
    "        PPPPPPPPPPPP",
    "        PPPPPPPPPPPP",
    "        PPPPPPPPPPPP",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "        XXXXXXXXXXXX",
    "                    ",
    "                    ",
    "                    ",
    "                    "
  };
  short *data_i16 = new short[20*16*sizeof(short)];
  uint8_t *data_u8 = new uint8_t[20*16*sizeof(uint8_t)];
  for (int i = 0; i < 16; i++) {
    for (int j = 0; j < 20; j++) {
      data_i16[i*20+j] = (short)data[i][j];
      data_u8[i*20+j] = (uint8_t)data[i][j];
    }
  }
  Image full_i16(data_i16, Point(16, 20), false); full_i16.owns_pixels = true;
  Image8 full_u8(data_u8, Point(16, 20), false); full_u8.owns_pixels = true;

  char center[]  = "    XXXXXXXX";
  char up_2[]    = "    VVVVVVVV";
  char down_2[]  = "    JJJJJJJJ";
  char left_2[]  = "      XXXXXX";
  char right_2[] = "  XXXXXXXXXX";
  char u1_r2[]   = "  WWWWWWWWWW";

  Profile p16(Rectangle(4, 11, 4, 15), Profile::OverX);
  p16.imprint(full_i16);    for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - center[i])  > 0.15) return 1;
  p16.slide(full_i16, -2);  for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - up_2[i])    > 0.15) return 2;
  p16.slide(full_i16, 2);   for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - center[i])  > 0.15) return 3;
  p16.slide(full_i16, 2);   for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - down_2[i])  > 0.15) return 4;
  p16.slide(full_i16, -2);  for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - center[i])  > 0.15) return 5;
  p16.scroll(full_i16, -2); for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - left_2[i])  > 0.15) return 6;
  p16.scroll(full_i16, 2);  for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - center[i])  > 0.15) return 7;
  p16.scroll(full_i16, 2);  for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - right_2[i]) > 0.15) return 8;
  p16.slide(full_i16, -1);  for (int i=0; i<p16.n; i++) if (fabsf(p16.reference[i] - u1_r2[i])   > 0.15) return 9;

  Profile p8(Rectangle(4, 11, 4, 15), Profile::OverX);
  p8.imprint8(full_u8);    for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - center[i])  > 0.15) return 11;
  p8.slide8(full_u8, -2);  for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - up_2[i])    > 0.15) return 12;
  p8.slide8(full_u8, 2);   for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - center[i])  > 0.15) return 13;
  p8.slide8(full_u8, 2);   for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - down_2[i])  > 0.15) return 14;
  p8.slide8(full_u8, -2);  for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - center[i])  > 0.15) return 15;
  p8.scroll8(full_u8, -2); for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - left_2[i])  > 0.15) return 16;
  p8.scroll8(full_u8, 2);  for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - center[i])  > 0.15) return 17;
  p8.scroll8(full_u8, 2);  for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - right_2[i]) > 0.15) return 18;
  p8.slide8(full_u8, -1);  for (int i=0; i<p8.n; i++) if (fabsf(p8.reference[i] - u1_r2[i])   > 0.15) return 19;

  char kenter[] = "JJJJJJJJ  ";
  char vp_2[]   = "DDJJJJJJJJ";
  char town_2[] = "JJJJJJ    ";
  char rite_2[] = "XXXXXXXX  ";
  char meft_2[] = "<<<<<<<<  ";
  char v2_m2[]  = "88<<<<<<<<";

  Profile q16(Rectangle(4, 13, 6, 13), Profile::OverY);
  q16.imprint(full_i16);    for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - kenter[i])  > 0.15) return 21;
  q16.scroll(full_i16, -2); for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - vp_2[i])    > 0.15) return 22;
  q16.scroll(full_i16, 2);  for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - kenter[i])  > 0.15) return 23;
  q16.scroll(full_i16, 2);  for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - town_2[i])  > 0.15) return 24;
  q16.scroll(full_i16, -2); for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - kenter[i])  > 0.15) return 25;
  q16.slide(full_i16, 2);   for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - rite_2[i])  > 0.15) return 26;
  q16.slide(full_i16, -2);  for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - kenter[i])  > 0.15) return 27;
  q16.slide(full_i16, -2);  for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - meft_2[i])  > 0.15) return 28;
  q16.scroll(full_i16, -2); for (int i=0; i<q16.n; i++) if (fabsf(q16.reference[i] - v2_m2[i])   > 0.15) return 29;

  Profile q8(Rectangle(4, 13, 6, 13), Profile::OverY);
  q8.imprint8(full_u8);    for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - kenter[i])  > 0.15) return 31;
  q8.scroll8(full_u8, -2); for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - vp_2[i])    > 0.15) return 32;
  q8.scroll8(full_u8, 2);  for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - kenter[i])  > 0.15) return 33;
  q8.scroll8(full_u8, 2);  for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - town_2[i])  > 0.15) return 34;
  q8.scroll8(full_u8, -2); for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - kenter[i])  > 0.15) return 35;
  q8.slide8(full_u8, 2);   for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - rite_2[i])  > 0.15) return 36;
  q8.slide8(full_u8, -2);  for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - kenter[i])  > 0.15) return 37;
  q8.slide8(full_u8, -2);  for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - meft_2[i])  > 0.15) return 38;
  q8.scroll8(full_u8, -2); for (int i=0; i<q8.n; i++) if (fabsf(q8.reference[i] - v2_m2[i])   > 0.15) return 39;

  return 0;
}

int test_mwt_align_best() {
  char data[26][41] = {
   /*0123456789112345678921234567893123456789*/
    "     !AWYZYWSYZYWZzjabgeacbaedcfzjervvih",
    "      !AWYZYWSYZYWZzjabgeacbaedcfzjzvcih",
    "       !AWYZYWSYZYWZzjabgeacbaedcfzjszih",
    "        !AWYZYWSYZYWZzjabgeacgaedcfzjfih",
    "         !AWYZYWSYZYWZzjabeacbgaedcfzjih",
    "          !AXZXWYWXYWXZzjjabgacbgaedcfih",
    "           @NZWYZYWYZWYZzabgecbgaedcfsih",
    "             @NZYWWYZWYWZabgeacbgaedfxih",
    "             !AWYZYWSYZYWZabgecbgaedcfih",
    "             !AXZXWYWXYWXZabgeacbaedcfih",
    "             @NZWYZYWYZWYZabgeacgaedcfih",
    "              @NZYWWYZWYWZbgeacbgaedcfih",
    "             !AWYZYWSYZYWZageacbgaedcfih",
    "             !AXZXWYWXYWXZabeacbgaedcfih",
    "             @NZWYZYWYZWYZabeacbgaedcfih",
    "              @NZYWWYZWYWZabgeabgaedcfih",
    "             !AWYZYWSYZYWZabgacbgaedcfih",
    "               !AXZXWYWXYWabgeabgaedcfih",
    "                @NZWYZYWYZabgecbgaedcfih",
    "                  @NZYWWYZabgeabgaedcfih",
    "                   @NZYWWYZabacbgaedcfih",
    "                    @NZYWWYZagebgaedcfih",
    "                     @NZYWWYZageacedcfih",
    "                      @NZYWWYZageacbgaih",
    "                       @NZYWWYZageacgaih",
    "                        @NZYWWYZagebgaih"
  };
  short *data_i16 = new short[26*40*sizeof(short)];
  short *data_i16t = new short[26*40*sizeof(short)];
  uint8_t *data_u8 = new uint8_t[26*40*sizeof(uint8_t)];
  uint8_t *data_u8t = new uint8_t[26*40*sizeof(uint8_t)];
  for (int i = 0; i < 26; i++) {
    for (int j = 0; j < 40; j++) {
      data_i16[i*40+j] = (short)data[i][j] * 64;
      data_i16t[j*26+i] = data_i16[i*40+j];
      data_u8[i*40+j] = (uint8_t)data[i][j];
      data_u8t[j*26+i] = (uint8_t)data[i][j];
    }
  }

  Image full_i16(data_i16, Point(26, 40), false); full_i16.owns_pixels = true;
  Image8 full_u8(data_u8, Point(26, 40), false); full_u8.owns_pixels = true;
  Image full_i16t(data_i16t, Point(40, 26), false); full_i16t.owns_pixels = true;
  Image8 full_u8t(data_u8t, Point(40, 26), false); full_u8t.owns_pixels = true;

  Profile search(Rectangle(0, 7, 0, 7), Profile::OverX);
  Rectangle over = full_i16.bounds;
  search.best_tiled_inside(full_i16, over);
  printf("%d,%d %dx%d\n", search.source.near.x, search.source.near.y, search.source.width(), search.source.height());

  Profile searcht(Rectangle(0, 7, 0, 7), Profile::OverY);
  Rectangle overt = full_i16t.bounds;
  searcht.best_tiled_inside(full_i16t, overt);
  printf("%d,%d %dx%d\n", search.source.near.x, search.source.near.y, search.source.width(), search.source.height());

  Profile search8(Rectangle(0, 7, 0, 7), Profile::OverX);
  Rectangle over8 = full_u8.bounds;
  search8.best_tiled_inside8(full_u8, over8);
  printf("%d,%d %dx%d\n", search.source.near.x, search.source.near.y, search.source.width(), search.source.height());

  Profile searcht8(Rectangle(0, 7, 0, 7), Profile::OverY);
  Rectangle overt8 = full_u8t.bounds;
  searcht8.best_tiled_inside8(full_u8t, overt8);
  printf("%d,%d %dx%d\n", search.source.near.x, search.source.near.y, search.source.width(), search.source.height());

  Point correct_tile(12, 12);
  Point size(8, 8);

  if (search.source.near != correct_tile) return 1;
  if (searcht.source.near != correct_tile) return 2;
  if (search8.source.near != correct_tile) return 3;
  if (searcht8.source.near != correct_tile) return 4;

  if (search.source.size() != size) return 5;
  if (searcht.source.size() != size) return 6;
  if (search8.source.size() != size) return 7;
  if (searcht8.source.size() != size) return 8;

  return 0;
}

int test_mwt_align()
{
  return test_mwt_align_features() +
    100*test_mwt_align_align() +
    10000*test_mwt_align_imprint() +
    1000000*test_mwt_align_slide() +
    100000000*test_mwt_align_best();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_align();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Align test result is %d\n",i);
  return i>0;
}
#endif

