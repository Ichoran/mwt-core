/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "MWT_Geometry.h"

int test_mwt_geometry_point()
{
  Point p1(0,1);
  Point p2 = p1.dup();
  if (p1!=p2 || !(p1==p2) || p1<p2 || p2<p1) return 1;
  
  p1 += Point(5,2);  // Should be (5,3)
  p1 -= Point(4,4);  // (1,-1)
  p1 *= 4;           // (4,-4)
  p1 /= 2;           // (2,-2)
  p1 += 3;           // (5,1)
  p1 -= 2;           // (3,-1)
  
  if (p1.distSq(Point(2,2)) != 10) return 2;
  if (p1.distSqX(Point(5,4)) != 4) return 3;
  if (p1.distSqY(Point(3,2)) != 9) return 4;
  
  p2 = (1+p2+Point(1,2)).cross(p2*2-Point(4,2)/2);

  if (p2.x+p2.y != 10) return 5;
  if (p1*p2 != -2) return 6;
  if (p1.X(p2) != 26) return 7;
  if (p1.cw()*p2.ccw() != (-p1)*p2) return 8;
  
  if (p1<p2) return 9;
  
  return 0;
}

int test_mwt_geometry_fpoint()
{
  FPoint p1(0,1);
  FPoint p2 = p1.dup();
  if (p1!=p2 || !(p1==p2) || p1<p2 || p2<p1) return 1;
  
  p1 += FPoint(5,2);  // Should be (5,3)
  p1 -= FPoint(4,4);  // (1,-1)
  p1 *= 4;           // (4,-4)
  p1 /= 2;           // (2,-2)
  p1 += 3;           // (5,1)
  p1 -= 2;           // (3,-1)
  
  if (p1.distSq(FPoint(2,2)) != 10) return 2;
  if (p1.distSqX(FPoint(5,4)) != 4) return 3;
  if (p1.distSqY(FPoint(3,2)) != 9) return 4;
  
  p2 = (1+p2+FPoint(1,2)).cross(p2*2-FPoint(4,2)/2);

  if (p2.x+p2.y != 10) return 5;
  if (p1*p2 != -2) return 6;
  if (p1.X(p2) != 26) return 7;
  if (p1.cw()*p2.ccw() != (-p1)*p2) return 8;
  
  if (p1<p2) return 9;
  
  if ( fabs(p2.length()-8.2462) > 0.0001 ) return 10;
  if ( p1.unit().distSq( FPoint(-0.3217505544) ) > 0.0000001 ) return 11;
  
  return 0;
}

int test_mwt_geometry_rectangle()
{
  Rectangle r1(0,5,0,5) , r2(Point(3,3),Point(4,6));
  
  if (r2<r1) return 1;
  if (!r1.contains(Point(5,3)) || !r2.containsY(5) || r2.containsX(2)) return 2;
  if (!r2.overlaps(r1)) return 3;
  if (r1.isEmpty()) return 4;
  if (r2.size().distSq(Point(0,0)) != 20) return 5;
  if (r1.area() != 36) return 6;
  if (r1<r1.bounds() || r1.bounds()<r1) return 7;
  
  r1.includeX(-1); r1.includeY(6);
  r2.includeXX(2,4); r2.includeYY(2,4);
  r1.include(Point(6,-1));
  if (!r1.contains(Point(6,6))) return 8;
  if (!r2.contains(2,2)) return 9;
  
  Rectangle r3(r2),r4(r2);
  r3.include(r1);
  r4 += r1;
  if (r3 != r4) return 10;
  r3 += Point(2,3);
  r4 -= Point(2,3);
  if (r3==r4 || !r3.overlaps(r4)) return 11;
  if ( (r2*r1) != r2) return 12;
  if ( (r1/Point(2,3)).area() != 12) return 13;
  
  Point p = r1.firstBoundaryPoint();
  Point q = p;
  do
  {
    if (!r1.contains(q)) return 14;
    if (r1.contains(q+Point(1,1)) && r1.contains(q-Point(1,1))) return 15;
    q = r1.nextBoundaryPoint(q);
  } while (p != q);
  
  return 0;
}

int test_mwt_geometry_ellipse()
{
  Point p,q;
  Ellipse e1( Point(4,4) , 5 ) , e2( Point(3,3) , Point(4,6) );
  
  if (e2.area() != 73) return 1;
  if (e1.contains( Point(8,0) ) || !e1.contains(Point(7,1))) return 2;
  if (e2.bounds() != Rectangle( Point(-1,-3) , Point(7,9) )) return 3;
  
  p = e2.firstBoundaryPoint();
  q = p;
  do
  {
    if (!e2.contains(q)) return 4;
    if (e2.contains(q+1) && e2.contains(q-1) && e2.contains(q+Point(1,-1)) && e2.contains(q+Point(-1,1))) return 5;
    q = e2.nextBoundaryPoint(q);
  } while (q != p);
  
  return 0;
}

int test_mwt_geometry()
{
  return test_mwt_geometry_point() + test_mwt_geometry_fpoint()*100 + test_mwt_geometry_rectangle()*10000 + test_mwt_geometry_ellipse()*1000000;
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_geometry();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Geometry test result is %d\n",i);
  return i>0;
}
#endif

