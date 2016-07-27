/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 * Copyright (C) 2016 Rex A. Kerr and Calico Life Sciences
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_STORAGE
#define MWT_STORAGE


#include <stdlib.h>
#include <new>
#include <cstddef>
#include "MWT_Lists.h"


/****************************************************************
                     Array Storage Templates
****************************************************************/

/* NOTE **
        * To use Array functions, allocate them like so:
        *   Array1D<int>& a = *Array1D<int>::createArray(5);
        * and dispose of them like thus:
        *   Array1D<int>::destroyArray(a);
        * (the original pointer also works).
** NOTE */

// Array1D provides in-line storage for a one-dimensional array.
// No constructors or destructors.

template <class T> class Array1D
{
private:
  Array1D() {}
  ~Array1D() {}
public:
  int length;
  T data[1];
  
  static Array1D<T>* createArray(int n_elements)
  {
    if (n_elements<1) n_elements=1;
    Array1D<T> *array = (Array1D<T>*) new char[sizeof(int) + n_elements*sizeof(T)];
    array->length = n_elements;
    return array;
  }
  static void destroyArray(Array1D<T>* a1d) { delete[] (char*)a1d; }
  static inline void destroyArray(Array1D<T>& a1d) { destroyArray(&a1d); }
  
  inline T& operator[](int n) { return data[n]; }
  inline T* operator+(int n) { return data+n; }
};


// Array2D provides in-line storage for a two-dimensional array.
// No constructors or destructors.

template <class T> class Array2D
{
private:
  Array2D() {}
  ~Array2D() {}
public:
  int major_length,minor_length;
  T data[1];
  
  static Array2D<T>* createArray(int n_rows,int n_cols)
  {
    if (n_rows<1) n_rows=1;
    if (n_cols<1) n_cols=1;
    Array2D<T> *array = (Array2D<T>*) new char[2*sizeof(int) + n_rows*n_cols*sizeof(T)];
    array->major_length = n_rows;
    array->minor_length = n_cols;
    return array;
  }
  static void destroyArray(Array2D<T>* a2d) { delete[] (char*)a2d; }
  static inline void destroyArray(Array2D<T>& a2d) { destroyArray(&a2d); }
  
  inline T& operator()(int R,int C) { return data[C + minor_length*R]; }
};



/****************************************************************
                   Linked List Storage Template
****************************************************************/

/* NOTE **
        * Classes must have a next data member or be wrapped by Stackable or Listable
** NOTE **
        * Only meant to be used on the heap--use "new", "delete"!
** NOTE **
        * Constructors/destructors for stored items are off by default
** NOTE **
        * It is safe to allocate a whole bunch of things with Storage and then
        * just delete the whole storage--if you've asked for it to use
        * destructors, it will clean everything up properly.
** NOTE */

// Storage class provides memory for arbitrary lists of classes
template <class T> class Storage
{
public:
  Storage<T> *used_storage;
  T *defunct;
  T *unused;
  int block_size;
  int n_used;
  bool ctor_dtor;
  
  Storage(int n_items,bool use_ctor_dtor=false) : used_storage(NULL),defunct(NULL),n_used(0),ctor_dtor(use_ctor_dtor)
  {
    if (n_items==0) { block_size=0; unused=NULL; }
    else
    {
      block_size = (n_items<1) ? 1 : n_items;
      unused = (T*) new char[block_size*sizeof(T)];
    }
  }
  ~Storage()
  {
    if (unused)
    {
      if (ctor_dtor)
      {
        int i;
        ptrdiff_t delta;
        T* t;
        char destructorable[ n_used ];
        for (i=0;i<n_used;i++) destructorable[i] = 1;
        for (t=defunct;t!=NULL;t=t->next)
        {
          defunct=t->next;
          delta = (ptrdiff_t)t - (ptrdiff_t)unused;
          if (delta < 0 || delta >= (ptrdiff_t)(sizeof(T)*block_size))  // Out of range
          {
            if (used_storage)
            {
              t->next = used_storage->defunct;
              used_storage->defunct = t;
            }
          }
          else destructorable[ delta / sizeof(T) ] = 0;  // In range, we would have destructored it already
        }
        for (i=0;i<n_used;i++) if (destructorable[i]) unused[i].~T();
      }
      delete[] (char*)unused;
      unused = NULL;
    }
    if (used_storage) { delete used_storage; used_storage=NULL; }
  }
  
  // Create an item from the store, calling default constructor if requested
  T* create(bool use_ctor)
  {
    T *t;
    if (defunct)
    {
      t = defunct;
      defunct = defunct->next;
    }
    else if (n_used < block_size) t = unused + (n_used++);
    else
    {
      Storage<T> *extra_storage = new Storage<T>(block_size);
      extra_storage->used_storage = used_storage;
      used_storage = extra_storage;
      extra_storage->n_used = n_used;
      t = extra_storage->unused;
      extra_storage->unused = unused;
      unused = t;
      n_used = 1;
    }
    t->next = NULL;
    if (use_ctor) { t = new(t) T(); }
    return t;
  }
  inline T* create() { return create(ctor_dtor); }
  
  // Free memory.  Don't try to free NULL pointers!
  void destroy(T *dead)
  {
    if (ctor_dtor) { dead->~T(); }
    dead->next = defunct;
    defunct = dead;
  }
  void destroyList(T*& dead_head , T*& dead_tail)
  {
    if (ctor_dtor)
    {
      T *t = dead_head;
      while (t!=NULL)
      {
        dead_head = dead_head->next;
        t->~T();
        t->next = defunct;
        defunct = t;
        t = dead_head;
      }
      dead_tail = NULL;
    }
    else if (dead_head != NULL)
    {
      dead_tail->next = defunct;
      defunct = dead_head;
      dead_head = dead_tail = NULL;
    }
  }
  void destroyList(T*& dead_head)
  {
    if (ctor_dtor)
    {
      T *t = dead_head;
      while (t!=NULL)
      {
        dead_head = dead_head->next;
        t->~T();
        t->next = defunct;
        defunct = t;
        t = dead_head;
      }
    }
    else if (dead_head != NULL)
    {
      T *dead_tail;
      for (dead_tail=dead_head ; dead_tail->next !=NULL ; dead_tail=dead_tail->next) { }
      destroyList(dead_head,dead_tail);
    }
  }
};



/****************************************************************
                   Managed List Template
****************************************************************/

/* NOTE **
        * ManagedList can either manage its own memory or let you do it.
        * The methods you call will depend on who is doing memory management!
        * If it is doing memory management, it needs an appropriate Storage class,
        * which it may own, or which may be a common storage.
        * If it is doing memory management, add items using the "const T& item"
        * methods--use others ONLY if they've been alloced from the same Storage!
        * Also, use Delete and Backspace to destroy things in the list.
        * If you are doing memory management, avoid the "const T& item" methods
        * for insertion and use removeL, removeR, pop, and drop to shrink the list.
** NOTE **
        * There are four ways to add items to the list
        * (1) push--adds to the beginning of the list
        * (2) insert--adds after a specified item in the list (or the current iterator)
        * (3) tuck--adds before a specified item in the list (or the current iterator)
        * (4) append--adds to the end of the list
        * There are also four corresponding ways to remove items
        *   - pop, removeR, removeL, and drop, if you are managing memory
        *   - Behead, Delete, Backspace, and Truncate, if it is managing memory
** NOTE **
        * Iterators going forwards should be started with start() and advance() should
        * be called before they're used--if advance returns false, the list is over.
        * Iterators going backwards should be started with end() and used if non-NULL.
        * retreat() will take the iterator backwards, returning false when it leaves a
        * NULL value.  If no pointer is given, the current data member will be used
        * (this is often overwritten so don't trust it!).  You can directly access the
        * struct stored in the head of the list with h(), in the tail with t(), and in
        * the current iterator with i() (but make sure current is not NULL!).  With
        * your own iterators of type Listable<T>* ptr, the data is at ptr->data.
** NOTE **
        * The = operator does not actually copy the list--it transfers ownership!
        * After that, the list on the RHS will be empty.
        * If you want to create a new deep copy of the list, use imitate().
** NOTE **
        * It is hard to get ManagedList to work nicely with classes that do their own
        * memory allocation/deallocation if MangaedList is handling memory.
        * Be mindful of whether Storage calls constructors/destructors, and turn on
        * use_dtor for ManagedList, and then things should be okay.
** NOTE **
        * If you create lists using SubordinateList or SubordinateCopy
        * actions, you must not use the subordinates after you've deallocated
        * the original if the original uses its own allocator.  If not,
        * make sure that all subordinates are finished before destroying
        * the allocator.
** NOTE */

// Doubly linked list that wraps a type and gives an iterator.
template <class T> class ManagedList
{
public:
  enum CreationType { DefaultConstructor , SubordinateList , SubordinateCopy };

  Listable<T> *list_head;
  Listable<T> *list_tail;
  Listable<T> *current;
  int size;
  bool owns_allocator,my_xtor;
  Storage< Listable<T> > *allocator;
  
  ManagedList() : list_head(NULL),list_tail(NULL),current(NULL),size(0),owns_allocator(false),my_xtor(false),allocator(NULL) { }
  ManagedList( Storage< Listable<T> >* use_this_allocator )
    : list_head(NULL),list_tail(NULL),current(NULL),size(0),owns_allocator(false),my_xtor(false)
  {
    allocator = use_this_allocator;
  }
  ManagedList(int allocator_chunk_size,bool use_cdtor=false)
    : list_head(NULL),list_tail(NULL),current(NULL),size(0),owns_allocator(true),my_xtor(use_cdtor)
  {
    // No need to have allocator do constructor stuff if we're handling everything (we can keep track better)
    allocator = new Storage< Listable<T> >(allocator_chunk_size,false);
  }
  ManagedList( const ManagedList<T>& existing , CreationType action )
    : list_head(NULL),list_tail(NULL),current(NULL),size(0),owns_allocator(false)
  {
    switch (action)
    {
      case SubordinateList:
      case SubordinateCopy:
        my_xtor = existing.my_xtor;
        allocator = existing.allocator;
        if (action==SubordinateCopy) imitate(existing);
        break;
      case DefaultConstructor:
      default:
        my_xtor = false;
        allocator = NULL;
        break;
    }
  }
  ~ManagedList() // Only deletes node memory if it's from the allocator!
  {
    if (allocator)
    {
      if (my_xtor) for (current=list_head;current!=NULL;current=current->next) current->data.~T();
      allocator->destroyList(list_head,list_tail);
      if (owns_allocator) delete allocator;
      allocator = NULL;
    }
  }
  
  // Throw away current lists--either we didn't alloc them or someone else has them now.
  inline void blank()  
  {
    list_head = list_tail = current = NULL;
    size = 0;
  }
  
  // Recycle current lists--we must have made them with the allocator
  inline void flush()
  {
    if (my_xtor) for (current=list_head;current!=NULL;current=current->next) current->data.~T();
    allocator->destroyList(list_head,list_tail);
    current = NULL;
    size = 0;
  }
  
  // make with an argument assumes that we're using an allocator--don't call it if we're not!
  Listable<T>* make(const T& item)  
  {
    Listable<T> *list_element;
    list_element = allocator->create(false);            // Don't call constructor
    list_element = new(list_element) Listable<T>(item); // Because we call it here
    return list_element;
  }
  
  // Transfer ownership of list!  RHS forgets anything it previously had.
  // LHS does no cleanup--call destructor or somesuch first; flush() doesnt' reclaim memory.
  ManagedList<T>& operator=(ManagedList<T>& ml)
  {
    list_head = ml.list_head;
    list_tail = ml.list_tail;
    current = ml.current;
    size = ml.size;
    owns_allocator = ml.owns_allocator;
    my_xtor = ml.my_xtor;
    allocator = ml.allocator;
    ml.blank();
    ml.owns_allocator = false;
    return *this;
  }
  // Deep copy of list!  Must be using own allocator
  ManagedList<T>& imitate(const ManagedList<T>& ml)
  {
    Listable<T> *iterator;
    flush();
    for (iterator=ml.list_head ; iterator!=NULL ; iterator=iterator->next)
    {
      make(iterator->data)->appendTo(list_head,list_tail);
    }
    size = ml.size;
    return *this;
  }
  
  // List addition/removal methods:
  // Caps mean we are handling memory, lower case means memory is provided or needs to be thought about
  // (lower case = lower level operation, kinda)
  
  // Push adds stuff onto the beginning of the list
  inline void push(Listable<T>* item) { item->pushOnto(list_head,list_tail); size++; }
  inline void Push(const T& item) { make(item)->pushOnto(list_head,list_tail); size++; }
  inline void push(ManagedList<T>& other_list)  // Careful about mixing memory here!
  {
    if (other_list.size==0) return;
    if (size==0)
    {
      list_head = other_list.head;
      list_tail = other_list.tail;
      size = other_list.size;
    }
    else
    {
      list_head->prev = other_list.list_tail;
      other_list.list_tail->next = list_head;
      list_head = other_list.head;
      size += other_list.size;
    }
    other_list.blank();
  }
  inline T* Push() // Allocates memory, doesn't call constructor--use placement new!
  {
    Listable<T>* item = allocator->create(false);
    item->pushOnto(list_head,list_tail); size++;
    return &(item->data);
  }
  
  // Append adds stuff onto the tail
  inline void append(Listable<T>* item) { item->appendTo(list_head,list_tail); size++; }
  inline void Append(const T& item) { make(item)->appendTo(list_head,list_tail); size++; }
  inline void append(ManagedList<T>& other_list)
  {
    if (other_list.size==0) return;
    if (size==0)
    {
      list_head = other_list.head;
      list_tail = other_list.tail;
      size = other_list.size;
    }
    else
    {
      list_tail->next = other_list.list_head;
      other_list.list_head->prev = list_tail;
      list_tail = other_list.list_tail;
      size += other_list.size;
    }
    other_list.blank();
  }
  inline T* Append() // Allocates memory, doesn't call constructor--use placement new!
  {
    Listable<T>* item = allocator->create(false);
    item->appendTo(list_head,list_tail); size++;
    return &(item->data);
  } 
  
  // Insert adds stuff after iterator (or before head, if iterator is NULL)
  inline void insert(Listable<T>*& iterator,Listable<T>* item)
  {
    if (iterator==NULL) push(item);
    else
    {
      item->addAfter(iterator);
      if (iterator==list_tail) list_tail=item;
      size++;
    }
  }
  inline void Insert(Listable<T>*& iterator,const T& item) { insert(iterator,make(item)); }
  void insert(Listable<T>*& iterator,ManagedList<T>& other_list)
  {
    if (other_list.size==0) return;
    if (size==0)
    {
      list_head = other_list.list_head;
      list_tail = other_list.list_tail;
      size = other_list.size;
    }
    else if (iterator==NULL) push(other_list);
    else
    {
      other_list.list_head->prev = iterator;
      other_list.list_tail->next = iterator->next;
      if (iterator->next==NULL) list_tail = other_list.list_tail;
      else iterator->next->prev = other_list.list_tail;
      iterator->next=other_list.list_head;
      size += other_list.size;
    }
    other_list.blank();
  }
  inline T* Insert(Listable<T>*& iterator)
  {
    Listable<T>* item = allocator->create(false);
    insert(iterator,item);
    return &(item->data);
  }
  inline void insert(Listable<T>* item) { insert(current,item); }
  inline void Insert(const T& item) { Insert(current,item); }
  inline void insert(ManagedList<T>& other_list) { insert(current,other_list); }
  inline T* Insert() { return Insert(current); }
  
  // Tuck adds stuff before iterator (or on tail if iterator is NULL)
  inline void tuck(Listable<T>*& iterator,Listable<T>* item)
  {
    if (iterator==NULL) append(item);
    else if (iterator==list_head) push(item);
    else
    {
      item->addBefore(iterator);
      size++;
    }
  }
  inline void Tuck(Listable<T>*& iterator,const T& item) { tuck(iterator,make(item)); }
  void tuck(Listable<T>*& iterator,ManagedList<T>& other_list)
  {
    if (other_list.size==0) return;
    if (size==0)
    {
      list_head = other_list.list_head;
      list_tail = other_list.list_tail;
      size = other_list.size;
    }
    else if (iterator==NULL) append(other_list);
    else if (iterator==list_head) push(other_list);
    else
    {
      other_list.list_tail->next = iterator;
      other_list.list_head->prev = iterator->prev;
      iterator->prev->next = other_list.list_head;
      iterator->prev = other_list.list_tail;
      size += other_list.size;
    }
    other_list.blank();
  }
  inline T* Tuck(Listable<T>*& iterator)
  {
    Listable<T>* item = allocator->create(false);
    tuck(iterator,item);
    return &(item->data);
  }
  inline void tuck(Listable<T>* item) { tuck(current,item); }
  inline void Tuck(const T& item) { Tuck(current,item); }
  inline void tuck(ManagedList<T>& other_list) { tuck(current,other_list); }
  inline T* Tuck() { return Tuck(current); }
  
  // Remove pulls the an element out of the list and advances
  inline Listable<T>* removeR(Listable<T>*& iterator)
  {
    if (iterator==NULL) return NULL;
    size--;
    Listable<T>* temp = iterator;
    iterator = temp->extractR(list_head,list_tail);
    return temp;
  }
  inline Listable<T>* removeR() { return removeR(current); }
  // Same thing, except we retreat
  inline Listable<T>* removeL(Listable<T>*& iterator)
  {
    if (iterator==NULL) return NULL;
    size--;
    Listable<T>* temp = iterator;
    iterator = temp->extractL(list_head,list_tail);
    return temp;
  }
  inline Listable<T>* removeL() { return removeL(current); }
  // Pop and drop remove an elements from the ends of the list
  inline Listable<T>* pop() { return removeR(list_head); }
  inline Listable<T>* drop() { return removeL(list_tail); }  
  
  // Behead destroys the first item in the list
  inline void Behead()
  {
    if (list_head==NULL) return;
    if (current==list_head) current=current->next;
    size--;
    Listable<T>* temp = list_head;
    list_head->extractR(list_head,list_tail);
    if (my_xtor) temp->data.~T();
    if (allocator!=NULL) allocator->destroy(temp);
  }
  // Truncate gets the last one
  inline void Truncate()
  {
    if (list_tail==NULL) return;
    if (current==list_tail) current=current->prev;
    size--;
    Listable<T>* temp = list_tail;
    list_tail->extractL(list_head,list_tail);
    if (my_xtor) temp->data.~T();
    if (allocator!=NULL) allocator->destroy(temp);
  }
  // Delete destroys the iterator and moves it forwards (but not off the end)
  inline void Delete(Listable<T>*& iterator)
  {
    if (iterator==NULL) return;
    size--;
    Listable<T>* temp = iterator;
    iterator = temp->extractR(list_head,list_tail);
    if (my_xtor) temp->data.~T();
    if (allocator!=NULL) allocator->destroy(temp);
  }
  inline void Delete() { Delete(current); }
  // Backspace destroys the iterator and moves it backwards (possibly off the end)
  inline void Backspace(Listable<T>*& iterator)
  {
    if (iterator==NULL) return;
    size--;
    Listable<T>* temp = iterator;
    if (iterator->prev!=NULL) iterator = temp->extractL(list_head,list_tail);
    else { temp->extractL(list_head,list_tail); iterator=NULL; }
    if (my_xtor) temp->data.~T();
    if (allocator!=NULL) allocator->destroy(temp);
  }
  inline void Backspace() { Backspace(current); }
  // Destroy destroys the given iterator and NULLs it
  inline void Destroy(Listable<T>*& iterator)
  {
    if (iterator==NULL) return;
    size--;
    iterator->extract(list_head,list_tail);
    if (my_xtor) iterator->data.~T();
    if (allocator!=NULL) allocator->destroy(iterator);
    iterator = NULL;
  }
  
  // Iterator functions move a pointer around (current by default)
  inline void start(Listable<T>*& iterator) const { iterator=NULL; }
  inline void start() { current=NULL; }
  inline bool advance(Listable<T>*& iterator) const
  {
    if (iterator==list_tail) return false;
    if (iterator==NULL) iterator=list_head;
    else iterator=iterator->next;
    return true;
  }
  inline bool advance() { return advance(current); }
  inline bool retreat(Listable<T>*& iterator) const
  {
    if (list_head==NULL || iterator==NULL) return false;
    iterator=iterator->prev;
    return (iterator==NULL) ? false : true;
  }
  inline bool retreat() { return retreat(current); }
  inline void end(Listable<T>*& iterator) const { iterator=list_tail; }
  inline void end() { current=list_tail; }
  
  // Direct access to data items--use only when they make sense (e.g. current after one advance!)
  inline T& h() { return list_head->data; }
  inline T& t() { return list_tail->data; }
  inline T& i() { return current->data; }
  inline T& ii(Listable<T>* lt) { return lt->data; }  // Syntactic sugar for separate iterator
  
  // Merge sorts using < or function pointer
  inline void mergeSort() { Listable<T>::mergeSort(list_head,list_tail); }
  inline void mergeSort(int (*lessThan)(T*,T*,void*) , void *data) { Listable<T>::mergeSort(list_head,list_tail,lessThan,data); }
};



/****************************************************************
                    Dual List Templates
****************************************************************/

// Lets you have two managed lists referring to the same data
// (but you have to manage both by hand carefully).
template <class T> class Dualled;
template <class T> class Dualable  // This one holds the data and points to the other
{
public:
  T value;
  Listable< Dualled<T> >* dual;
  Dualable() : dual(NULL) { }
  ~Dualable() { }
  bool operator<(const Dualable<T>& d) { return data<d.data; }
  inline T& data() { return value; }
};
template <class T> class Dualled  // This one doesn't hold the data and points to the one that does
{
public:
  Listable< Dualable<T> >* dual;
  Dualled() : dual(NULL) { }
  ~Dualled() { }
  bool operator<(const Dualled<T>& d) { return dual->data.value < d.dual->data.value; }
  inline T& data() { return dual->data.value; }
};


// A class to carefully manage one set of data referred to by two managed lists
template <class T> class DualList
{
public:
  ManagedList< Dualable<T> > primary;
  ManagedList< Dualled<T> > secondary;
  
  // Constructors just wrap the two lists
  DualList() { }
  DualList(int chunk_size,bool use_xtor) : primary(chunk_size,use_xtor),secondary(chunk_size,use_xtor) { }
  DualList( Storage< Listable< Dualable<T> > >* store1 , Storage< Listable< Dualled<T> > >* store2) : primary(store1),secondary(store2) { }
  ~DualList() { }
  
  // Getting rid of everything is easy
  void flush() { primary.flush(); secondary.flush(); }
  
  // Data access requires accessors for each iterator
  inline T& h1() { return primary.h().data(); }
  inline T& h2() { return secondary.h().data(); }
  inline T& i1() { return primary.i().data(); }
  inline T& i2() { return secondary.i().data(); }
  inline T& t1() { return primary.t().data(); }
  inline T& t2() { return secondary.t().data(); }
  
  // Moving around is not too hard, but we need to move separately
  inline void start1() { primary.start(); }
  inline void start2() { secondary.start(); }
  inline bool advance1() { return primary.advance(); }
  inline bool advance2() { return secondary.advance(); }
  inline bool retreat1() { return primary.retreat(); }
  inline bool retreat2() { return secondary.retreat(); }
  inline void end1() { primary.end(); }
  inline void end2() { secondary.end(); }
  inline Listable< Dualable<T> >* getIterator1() { return primary.current; }
  inline Listable< Dualled<T> >* getIterator2() { return secondary.current; }
  inline void setIterator1(Listable< Dualable<T> >* ld1) { primary.current=ld1; }
  inline void setIterator2(Listable< Dualled<T> >* ld2) { secondary.current=ld2; }
  
  // Sometimes you want both iterators to be at the same place
  void goto1() { if (primary.current==NULL) secondary.current=NULL; else secondary.current = primary.i().dual; }
  void goto2() { if (secondary.current==NULL) primary.current=NULL; else primary.current = secondary.i().dual; }
  
  // Adding objects wraps the two lists, but here care is needed
  T* Push()
  {
    Dualable<T>* d1 = primary.Push();
    Dualled<T>* d2 = secondary.Push();
    d2->dual = primary.list_head;
    d1->dual = secondary.list_head;
    return &d1->value;
  }
  T* Append()
  {
    Dualable<T>* d1 = primary.Append();
    Dualled<T>* d2 = secondary.Append();
    d2->dual = primary.list_tail;
    d1->dual = secondary.list_tail;
    return &d1->value;
  }
  T* Insert()
  {
    Dualable<T>* d1;
    Dualled<T>* d2;
    if (primary.current==NULL) d1 = primary.Push();
    else d1 = primary.Insert();
    if (secondary.current==NULL) d2 = secondary.Push();
    else d2 = secondary.Insert();
    if (primary.current==NULL) d2->dual = primary.list_head;
    else d2->dual = primary.current->next;
    if (secondary.current==NULL) d1->dual = secondary.list_head;
    else d1->dual = secondary.current->next;
    return &d1->value;
  }
  inline void Push(const T& t) { T* tp = Push(); *tp = t; }
  inline void Append(const T& t) { T* tp = Append(); *tp = t; }
  inline void Insert(const T& t) { T* tp = Insert(); *tp = t; }
  
  // Removing objects is tricky, and there are separate versions depending on the iterator you want
  // Also, this can clobber the iterator of the other list, so take care (use a gotoN or startN or endN)
  void Behead1()
  {
    secondary.Destroy( primary.h().dual );
    primary.Behead();
  }
  void Behead2()
  {
    primary.Destroy( secondary.h().dual );
    secondary.Behead();
  }
  void Truncate1()
  {
    secondary.Destroy( primary.t().dual );
    primary.Truncate();
  }
  void Truncate2()
  {
    primary.Destroy( secondary.t().dual );
    secondary.Truncate();
  }
  void Backspace1()
  {
    secondary.Destroy( primary.i().dual );
    primary.Backspace();
  }
  void Backspace2()
  {
    primary.Destroy( secondary.i().dual );
    secondary.Backspace();
  }
  void Delete1()
  {
    secondary.Destroy( primary.i().dual );
    primary.Delete();
  }
  void Delete2()
  {
    primary.Destroy( secondary.i().dual );
    secondary.Delete();
  }
  
  // Standard merge sorts
  inline void mergeSort1() { primary.mergeSort(); }
  inline void mergeSort1(int (*lessThan)(Dualable<T>*,Dualable<T>*,void*) , void* data) { primary.mergeSort(lessThan,data); }
  inline void mergeSort2() { secondary.mergeSort(); }
  inline void mergeSort2(int (*lessThan)(Dualled<T>*,Dualled<T>*,void*) , void* data) { secondary.mergeSort(lessThan,data); }
  // Re-point secondary list so it follows primary list order
  void use1toSort2()
  {
    if (primary.size==0) return;
    primary.start();
    secondary.start();
    while (primary.advance() && secondary.advance())
    {
      secondary.i().dual = primary.current;
    }
  }
  // Relink primary list so it follows secondary list order
  void use2toSort1()
  {
    if (secondary.size==0) return;
    primary.list_head = secondary.h().dual;
    primary.list_tail = secondary.t().dual;
    secondary.start();
    while (secondary.advance())
    {
      if (secondary.current->prev==NULL) secondary.i().dual->prev = NULL;
      else secondary.i().dual->prev = secondary.current->prev->data.dual;
      if (secondary.current->next==NULL) secondary.i().dual->next = NULL;
      else secondary.i().dual->next = secondary.current->next->data.dual;
    }
  }
};



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_storage_arrays();
int test_mwt_storage_storage();
int test_mwt_storage_list();
int test_mwt_storage_dual();
int test_mwt_storage();


#endif

