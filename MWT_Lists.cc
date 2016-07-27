/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#include <string.h>
#include <stdio.h>

#include "MWT_Lists.h"
#include "MWT_Geometry.h"

int test_less_than_function(const Point& p,const Point& q,void *extra)
{
  return (p.y < q.y || (p.y==q.y && p.x<q.x));
}

int test_mwt_lists_stackable()
{
  Stackable<Point> points[10];
  Stackable<Point> *sp=NULL;
  Stackable<Point> *sph=NULL,*spt=NULL;
  int i;
  
  for (i=0;i<10;i++) points[i].data = Point(i,-i);
  for (i=0;i<10;i++)
  {
    if (i&0x1) points[i].pushOnto(sp);
    else if (i&0x2) points[i].pushOnto(sph,spt);
    else points[i].appendTo(sph,spt);
  }
  
  if (sp->next->next->next->data != Point(3,-3)) return 1;
  if (spt->data != Point(8,-8) || sph->next->data != Point(2,-2)) return 2;
  
  Stackable<Point>::pushList(sph,spt,sp);
  if (sph->next->next->next->data != Point(3,-3) || spt->data != Point(8,-8)) return 3;
  
  for (i=0;i<3;i++) sp = Stackable<Point>::pop(sph,spt);
  if (sph->data != Point(3,-3)) return 4;
  
  sph = Stackable<Point>::mergeSort(sph);
  for (spt=NULL,sp=sph;sp!=NULL;spt=sp,sp=sp->next)
  {
    if (spt!=NULL && spt->data.x >= sp->data.x) return 5;
    if (sp->data.x == 5 || sp->data.x == 7 || sp->data.x == 9) return 6;
  }
  if (sph->next->data.x != 1) return 7;
  
  sph = Stackable<Point>::mergeSort(sph,test_less_than_function,NULL);
  for (spt=NULL,sp=sph;sp!=NULL;spt=sp,sp=sp->next)
  {
    if (spt!=NULL && spt->data.x <= sp->data.x) return 8;
  }
  if (sph->next->data.x != 6) return 9;
  
  return 0;
}

int test_mwt_lists_listable()
{
  Listable<Point> points[10];
  Listable<Point> *lp=NULL,*lq=NULL;
  Listable<Point> *lph=NULL,*lpt=NULL;
  int i;
  
  for (i=0;i<10;i++) points[i].data = Point(i,-i);
  for (i=0;i<10;i++)
  {
    if (i&0x1) points[i].pushOnto(lph,lpt);
    else if (i&0x2) points[i].addAfter(lph);
    else points[i].appendTo(lph,lpt);
  }
  for (i=0,lp=lph; i<2 ; i++,lp=lp->next) {}
  lp=lp->extractL(lph,lpt);
  lpt->prev->prev->prev->prev->extractR();
  
  if (lp->data != Point(7,-7)) return 1;
  if (lph->next->next->next->data != Point(3,-3)) return 2;
  if (lpt->prev->prev->next->prev->prev->data != Point(2,-2)) return 3;
  
  Listable<Point>::mergeSort(lph,lpt);
  for (lq=NULL,lp=lph;lp!=NULL;lq=lp,lp=lp->next)
  {
    if (lq!=NULL && lq->data.x >= lp->data.x) return 4;
    if (lp->data.x == 5 || lp->data.x == 1) return 5;
  }
  if (lph->next->data.x != 2) return 6;
  
  lp=lq=NULL;
  points[1].pushOnto(lp,lq);
  points[5].appendTo(lp,lq);
  Listable<Point>::joinL(lph,lpt,lp,lq);
  if (lph->data.x != 1 || lph->next->next->next->data.x != 2) return 7;
  
  Listable<Point>::mergeSort(lph,lpt,test_less_than_function,NULL);
  for (lq=NULL,lp=lph;lp!=NULL;lq=lp,lp=lp->next)
  {
    if (lq!=NULL && lq->data.x <= lp->data.x) return 8;
  }
  if (lph->next->data.x != 8) return 9;
  
  return 0;
}

int test_mwt_lists()
{
  return test_mwt_lists_stackable() + 100*test_mwt_lists_listable();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_lists();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Lists test result is %d\n",i);
  return i>0;
}
#endif

