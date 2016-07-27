/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_GEOMETRY
#define MWT_GEOMETRY


#include <math.h>


/****************************************************************
                         Points (in 2D)
****************************************************************/

// A point in integer coordinates in Cartesian coordinate space
class Point
{
public:
  int x;
  int y;
  
  Point() { }
  Point(int X,int Y) : x(X),y(Y) { }
  
  Point dup() const { return Point(*this); }
  
  inline bool operator==(const Point& p) const { return (x==p.x && y==p.y); }
  inline bool operator!=(const Point& p) const { return (x!=p.x || y!=p.y); }
  inline bool operator<(const Point& p) const { return (x<p.x || (x==p.x && y<p.y)); }
  
  inline Point& operator+=(const Point& p) { x+=p.x; y+=p.y; return *this; }
  inline Point& operator+=(int i) { x+=i; y+=i; return *this; }
  inline Point& operator-=(const Point& p) { x-=p.x; y-=p.y; return *this; }
  inline Point& operator-=(int i) { x-=i; y-=i; return *this; }
  inline Point& operator*=(int i) { x*=i; y*=i; return *this; }
  inline Point& operator/=(int i) { x/=i; y/=i; return *this; }
  
  inline Point operator+(const Point& p) const { return Point(x+p.x,y+p.y); }
  inline Point operator+(int i) const { return Point(x+i,y+i); }
  inline Point operator-() const { return Point(-x,-y); }
  inline Point operator-(const Point& p) const { return Point(x-p.x,y-p.y); }
  inline Point operator-(int i) const { return Point(x-i,y-i); }
  inline int operator*(const Point& p) const { return x*p.x + y*p.y; }
  inline Point operator*(int i) const { return Point(i*x,i*y); }
  inline Point operator/(int i) const { return Point(x/i,y/i); }
  inline int X(const Point& p) const { return x*p.y - y*p.x; }
  inline Point cross(const Point& p) const { return Point(x*p.y,-y*p.x); }
  
  inline Point ccw() const { return Point(-y,x); }
  inline Point cw() const { return Point(y,-x); }

  inline int distSqX(const Point& p) const { return (x-p.x)*(x-p.x); }
  inline int distSqY(const Point& p) const { return (y-p.y)*(y-p.y); }
  inline int distSq(const Point& p) const { return (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y); }
  inline int length2() const { return x*x + y*y; }
};
inline Point operator+(int i,const Point& p) { return p+i; }
inline Point operator-(int i,const Point& p) { return p-i; }
inline Point operator*(int i,const Point &p) { return p*i; }


// A point in floating-point coordinates in Cartesian coordinate space
class FPoint
{
public:
  float x;
  float y;
  
  FPoint() { }
  FPoint(float X,float Y) : x(X),y(Y) { }
  FPoint(float theta) : x(cos(theta)),y(sin(theta)) { }
  FPoint(const Point &p) : x(p.x),y(p.y) { }
  
  FPoint dup() const { return FPoint(*this); }
  Point toPoint() const { return Point( (int)floor(x+0.5) , (int)floor(y+0.5) ); }
  
  inline bool operator==(const FPoint& p) const { return (x==p.x && y==p.y); }
  inline bool operator!=(const FPoint& p) const { return (x!=p.x || y!=p.y); }
  inline bool operator<(const FPoint &p) const { return (x<p.x || (x==p.x && y<p.y)); }
  
  inline FPoint& operator+=(const FPoint &f) { x+=f.x; y+=f.y; return *this; }
  inline FPoint& operator+=(float f) { x+=f; y+=f; return *this; }
  inline FPoint& operator-=(const FPoint &f) { x-=f.x; y-=f.y; return *this; }
  inline FPoint& operator-=(float f) { x-=f; y-=f; return *this; }
  inline FPoint& operator*=(float f) { x*=f; y*=f; return *this; }
  inline FPoint& operator/=(float f) { x/=f; y/=f; return *this; }

  inline FPoint operator+(const FPoint &p) const { return FPoint(x+p.x,y+p.y); }
  inline FPoint operator+(float f) const { return FPoint(x+f,y+f); }
  inline FPoint operator-() const { return FPoint(-x,-y); }
  inline FPoint operator-(const FPoint &p) const { return FPoint(x-p.x,y-p.y); }
  inline FPoint operator-(float f) const { return FPoint(x-f,y-f); }
  inline float operator*(const FPoint &p) const { return x*p.x + y*p.y; }
  inline FPoint operator*(float i) const { return FPoint(i*x,i*y); }
  inline FPoint operator/(float i) const { return FPoint(x/i,y/i); }
  inline float X(const FPoint& p) const { return x*p.y-y*p.x; }
  inline FPoint cross(const FPoint &p) const { return FPoint(x*p.y,-y*p.x); }

  inline FPoint ccw() const { return FPoint(-y,x); }
  inline FPoint cw() const { return FPoint(y,-x); }

  inline float distSqX(const FPoint &p) const { return (x-p.x)*(x-p.x); }
  inline float distSqY(const FPoint &p) const { return (y-p.y)*(y-p.y); }
  inline float distSq(const FPoint &p) const { return (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y); }
  inline float length() const { return sqrt(x*x+y*y); }
  inline float length2() const { return x*x + y*y; }
  inline FPoint unit() const { return (*this)/length(); }
};
inline FPoint operator+(float f,const FPoint& p) { return p+f; }
inline FPoint operator-(float f,const FPoint& p) { return p-f; }
inline FPoint operator*(float f,const FPoint& p) { return p*f; }



/****************************************************************
           Geometric Shapes: Rectangle, Circle, etc.
****************************************************************/

/* NOTE **
        * Adding rectangles together produces a rectangle enclosing both (union)
        * Multiplying rectangles produces a rectangle enclosing the intersection
        * Adding or subtracting a point to a rectangle translates the rectangle
        * Multiplying or dividing a rectangle by a point scales the rectangle
** NOTE **
        * nextBoundaryPoint is done when it returns the same point as firstBoundaryPoint
** NOTE */

// A rectangle in integer coordinates
class Rectangle
{
public:
  Point near;
  Point far;
  
  Rectangle() { }
  Rectangle(Point A,Point B) : near(A),far(B) { }
  Rectangle(int L,int R,int B,int U) : near(L,B),far(R,U) { }
  
  Rectangle dup() const { return Rectangle(near,far); }
  
  inline bool operator<(const Rectangle& r) const { return (near.x<r.near.x || (near.x==r.near.x && near.y<r.near.y)); }
  inline bool operator==(const Rectangle& r) const { return near==r.near && far==r.far; }
  inline bool operator!=(const Rectangle& r) const { return near!=r.near || far!=r.far; }
  
  inline bool containsX(int x) const { return (near.x<=x && x<=far.x); }
  inline bool containsY(int y) const { return (near.y<=y && y<=far.y); }
  inline bool contains(const Point& p) const { return (near.x<=p.x && p.x<=far.x && near.y<=p.y && p.y<=far.y); }
  inline bool contains(int x,int y) const { return (near.x<=x && x<=far.x && near.y<=y && y<=far.y); }
  inline bool contains(const Rectangle& r) const
  {
    return (near.x<=r.near.x && r.far.x<=far.x && near.y<=r.near.y && r.far.y<=far.y);
  }
  inline bool overlaps(const Rectangle& r) const
  {
    return !(far.x<r.near.x || near.x>r.far.x || far.y<r.near.y || near.y>r.far.y);
  }

  inline bool isEmpty() const { return (near.x>far.x || near.y>far.y); }
  inline int width() const { return 1+far.x-near.x; }
  inline int height() const { return 1+far.y-near.y; }
  inline Point size() const { return Point(1,1)+far-near; }
  inline int area() const { return width()*height(); }
  inline Rectangle bounds() { return dup(); }

  inline Rectangle& includeX(int x) { if (x<near.x) { near.x=x; } if (far.x<x) { far.x=x; } return *this; }
  inline Rectangle& includeY(int y) { if (y<near.y) { near.y=y; } if (far.y<y) { far.y=y; } return *this; }
  inline Rectangle& includeXX(int x0,int x1) { if (x0<near.x) { near.x=x0; } if (x1>far.x) { far.x=x1; } return *this; }
  inline Rectangle& includeYY(int y0,int y1) { if (y0<near.y) { near.y=y0; } if (y1>far.y) { far.y=y1; } return *this; }
  inline Rectangle& include(const Point& p) { includeX(p.x); includeY(p.y); return *this; }
  inline Rectangle& include(const Rectangle& r) { includeXX(r.near.x,r.far.x); includeYY(r.near.y,r.far.y); return *this; }
  
  inline Rectangle& operator+=(const Rectangle& r) { return include(r); }
  inline Rectangle& operator+=(const Point& p) { near+=p; far+=p; return *this; }
  inline Rectangle& operator-=(const Point& p) { near-=p; far-=p; return *this; }
  Rectangle& operator*=(const Rectangle& r)
  {
    if (r.near.x>near.x) near.x = r.near.x;
    if (r.far.x<far.x) far.x = r.far.x;
    if (r.near.y>near.y) near.y = r.near.y;
    if (r.far.y<far.y) far.y = r.far.y;
    return *this;
  }
  inline Rectangle& operator*=(const Point& p) { near.x*=p.x; near.y*=p.y; far.x*=p.x; far.y*=p.y; return *this; }
  inline Rectangle& operator*=(int k) { near*=k; far*=k; return *this; }
  inline Rectangle& operator/=(const Point& p) { near.x/=p.x; near.y/=p.y; far.x/=p.x; far.y/=p.y; return *this; }
  inline Rectangle& operator/=(int k) { near/=k; far/=k; return *this; }
  inline Rectangle& cropTo(const Rectangle& r) { return (*this) *= r; }
  inline Rectangle& expand(int i) { near-=i; far+=i; return *this; }
  inline Rectangle& adoptOrigin(const Point& p) { near -= p; far -= p; return *this; }
  
  inline Rectangle operator+(const Point& p) const { return Rectangle(near+p,far+p); }
  inline Rectangle operator-(const Point& p) const { return Rectangle(near-p,far-p); }
  inline Rectangle operator*(const Rectangle& r) const { return dup() *= r; }
  inline Rectangle operator*(const Point& p) const { return dup() *= p; }
  Rectangle operator*(int k) const { return Rectangle(near*k,far*k); }
  inline Rectangle operator/(const Point& p) const { return dup() /= p; }
  Rectangle operator/(int k) const { return Rectangle(near/k,far/k); }

  inline Point firstBoundaryPoint() const { return near; }
  Point nextBoundaryPoint(Point p) const
  {
    if (p.x==near.x)
    {
      if (p.y==far.y) p.x++;
      else p.y++;
    }
    else if (p.x==far.x)
    {
      if (p.y==near.y) p.x--;
      else p.y--;
    }
    else if (p.y==far.y) p.x++;
    else p.x--;
    return p;
  }
};
inline Rectangle operator+(const Point& p,const Rectangle& r) { return r+p; }
inline Rectangle operator*(int k,const Rectangle& r) { return r*k; }

/* NOTE **
        * nextBoundaryPoint is done when it returns the same point as firstBoundaryPoint
** NOTE */

// An ellipse (circles are a special case)
class Ellipse
{
public:
  Point center;
  Point axes;
  
  Ellipse() { }
  Ellipse(Point C,int R) : center(C),axes(R,R) { }
  Ellipse(Point C,Point A) : center(C),axes(A) { if (axes.x*axes.y<0) { if (axes.x<0) axes.x=-axes.x; else axes.y=-axes.y; } }
  
  // A point p is in an ellipse centered at 0 with axes of length a,b if (p.x/a)^2 + (p.y/b)^2 <= 1.
  // To avoid divisions, multiply through by a*a*b*b: b^2*p.x^2 + a^2*p.y^2 <= a^2*b^2
  bool contains(const Point& p) const
  {
    Point q = p-center;
    Point r = q.cross(axes);
    long long ll_rx2,ll_ry2,ll_axesx2,ll_axesy2;  // Need long long to avoid overflow
    ll_rx2 = r.x; ll_rx2 *= ll_rx2;
    ll_ry2 = r.y; ll_ry2 *= ll_ry2;
    ll_axesx2 = axes.x; ll_axesx2 *= ll_axesx2;
    ll_axesy2 = axes.y; ll_axesy2 *= ll_axesy2;
    return (ll_rx2 + ll_ry2 <= ll_axesx2*ll_axesy2);
  }
  Rectangle bounds() const { return Rectangle(center-axes,center+axes); }
  
  Point firstBoundaryPoint() const { return Point(center.x+axes.x,center.y); }
  Point nextBoundaryPoint(Point p) const
  {
    if (p.y>=center.y)
    {
      if (p.x>center.x)
      {
        p.y++;
        if (!contains(p)) p.x--;
        if (!contains(p)) p.y--;
      }
      else
      {
        p.x--;
        if (!contains(p)) p.y--;
        while (!contains(p)) p.x++;
      }
    }
    else
    {
      if (p.x<center.x)
      {
        p.y--;
        if (!contains(p)) p.x++;
        if (!contains(p)) p.y++;
      }
      else
      {
        p.x++;
        if (!contains(p)) p.y++;
        if (!contains(p)) p.x--;
      }
    }
    return p;
  }
  
  int area() const
  {
    if (axes.x <= 0 || axes.y <=0) return 0;
    
    Point p = Point(axes.x,0);
    int A = 2*p.x + 1;
    while (p.y < axes.y)
    {
      p.y++;
      while (!contains(p+center)) p.x--;
      A += 4*p.x+2;
    }
    return A;
  }
};



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_geometry_point();
int test_mwt_geometry_fpoint();
int test_mwt_geometry_rectangle();
int test_mwt_geometry_ellipse();
int test_mwt_geometry();

#endif

