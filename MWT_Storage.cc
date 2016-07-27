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

#include "MWT_Storage.h"

int test_mwt_storage_arrays()
{
  int i,j;
  
  Array1D<double>& data = *Array1D<double>::createArray(1024);
  for (i=0;i<1024;i++) data[i] = i*2.718281828459045;
  for (i=0;i<1024;i++) if ( fabs(*(data+i) - i*2.718281828459045) > 1e-12 ) return 1;
  Array1D<double>::destroyArray(data);
  
  Array2D<short>& img = *Array2D<short>::createArray(480,640);
  for (i=0;i<480;i++) for (j=0;j<640;j++) img(i,j) = (short)(i-j);
  if (img(25,77) != -52 || img(331,224) != 107) return 2;
  Array2D<short>::destroyArray(img);
  
  return 0;
}

int test_mwt_storage_storage()
{
  int i;
  Stackable<int> *si=NULL;
  Storage< Stackable<int> >* ssi = new Storage< Stackable<int> >(1024);
  
  for (i=0;i<2277;i++)
  {
    ssi->create()->pushOnto(si);
    si->data = 11*(i%207)+(i/207);
  }
  si = Stackable<int>::mergeSort(si);
  for (i=0;i<1457;i++) ssi->destroy( Stackable<int>::pop(si) );
  for (i=0;i<4095-2277+1457;i++)
  {
    ssi->create()->pushOnto(si);
    si->data = 4095-i;
  }
  if (ssi->used_storage==NULL || ssi->used_storage->used_storage==NULL || ssi->used_storage->used_storage->used_storage==NULL) return 1;
  if (ssi->used_storage->used_storage->used_storage->used_storage != NULL) return 2;
  if (ssi->block_size-ssi->n_used != 1) return 3;
  if (ssi->defunct != NULL) return 4;
  
  ssi->destroyList(si);
  si=NULL;
  for (i=0;i<3155;i++)
  {
    ssi->create()->pushOnto(si);
    si->data = -i;
  }
  if (ssi->defunct == NULL) return 5;
  if (ssi->used_storage->used_storage->used_storage->used_storage != NULL) return 6;
 
  delete ssi;
  
  return 0;
}

int test_mwt_storage_list()
{
  ManagedList<double> *mld = new ManagedList<double>(1024);
  int i;
  double d;
  
  d = 0.1318365179;
  for (i=0;i<1536;i++)
  {
    if (d<0.4) d = 1.0-d*d;
    else d = d*d;
    switch (i&0x3)
    {
      case 0: mld->Push(d); break;
      case 1: mld->Append(d); break;
      case 2: mld->start(); mld->advance(); mld->advance(); mld->Insert(d); break;
      case 3: mld->end(); mld->retreat(); mld->retreat(); mld->Tuck(d); break;
    }
  }
  mld->mergeSort();
  d = 0;
  i = 0;
  mld->start();
  while (mld->advance())
  {
    i++;
    if (d > mld->i()) return 1;
    d = mld->i();
    if ((i&0x3)==0) mld->Backspace();
  }
  if (mld->size != (1536*3)/4) return 2;
  
  ManagedList<double> *mld2 = new ManagedList<double>(512);
  mld2->imitate(*mld);
  mld->flush();
  
  for (i=1536;i<2048;i++)
  {
    if (d<0.4) d = 1.0-d*d*d;
    else d = d*d*d;
    mld2->Push(d);
  }
  if (mld2->size != 512+(1536*3)/4) return 3;
  while (mld2->size>0){  mld->Push( mld2->h() ); mld2->Behead(); }
  delete mld2;
  mld->mergeSort();
  for (mld->end() ; mld->current!=NULL ; mld->retreat())
  {
    if (mld->current->next!=NULL)
    {
      if (mld->current->next->data < mld->i()) return 4;
    }
  }
  delete mld;
  
  return 0;
}

int test_mwt_storage_dual()
{
  DualList<short> dls(16,false);
  
  dls.Append(1);
  dls.Append(7);
  dls.Append(4);
  dls.Append(6);
  dls.Append(3);
  dls.Append(5);
  dls.Append(2);
  dls.mergeSort2();
  
  int i;
  dls.start2();
  for (i=1 ; i<=5 ; i++) { if (!dls.advance2()) return 1; }
  if (dls.i2() != 5) return 2;
  dls.goto2();
  if (!dls.retreat1()) return 3;
  if (dls.i1() != 3) return 4;
  
  dls.goto1();
  dls.Insert(8);
  i = 0;
  dls.start1();
  while (dls.advance1()) { i++; if (dls.i1()==8) break; }
  if (i!=6) return 5;
  i=0;
  dls.start2();
  while (dls.advance2()) { i++; if (dls.i2()==8) break; }
  if (i!=4) return 6;
  
  dls.use2toSort1();
  dls.start1();
  dls.start2();
  while (dls.advance1() && dls.advance2()) if (dls.i1() != dls.i2()) return 7;
  
  return 0;
}

int test_mwt_storage()
{
  return test_mwt_storage_arrays() + 100*test_mwt_storage_storage() + 10000*test_mwt_storage_list() + 1000000*test_mwt_storage_dual();
}

#ifdef UNIT_TEST_OWNER
int main(int argc,char *argv[])
{
  int i = test_mwt_storage();
  if (argc<=1 || strcmp(argv[1],"-quiet") || i) printf("MWT_Storage test result is %d\n",i);
  return i>0;
}
#endif

