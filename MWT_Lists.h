/* Copyright (C) 2007 - 2010 Howard Hughes Medical Institute
 * Copyright (C) 2007 - 2010 Rex A. Kerr, Nicholas A. Swierczek
 *
 * This file is a part of the Multi-Worm Tracker and is distributed under the
 * terms of the GNU Lesser General Public Licence version 2.1 (LGPL 2.1).
 * For details, see file LICENSE, or http://www.gnu.org/licences
 */
 
#ifndef MWT_LISTS
#define MWT_LISTS


#include <stdlib.h>


/****************************************************************
                   Singly Linked List Template
****************************************************************/

/* NOTE **
        * Stackable only supports NULL-terminated lists, not circular ones.
** NOTE **
        * Mergesort is a high-performance implementation: qsort is slower for
        * all but long random lists (and then is less than 2x faster)
** NOTE **
        * If you keep track of both ends of the list, always use the head/tail
        * versions.  They work properly if head and tail are NULL.  Otherwise,
        * never use the tail versions.  If you want to append to a list, keep
        * track of both head and tail.  If you want to remove from the end
        * of a list, use the Listable class instead!
** NOTE */

// Class to wrap another struct/class into a singly linked list
template <class Node> class Stackable
{
public:
  Stackable<Node> *next;
  Node data;
  
  Stackable() { }
  Stackable(const Node& node_data) : next(NULL),data(node_data) { }
  
  inline void pushOnto(Stackable<Node>*& head) { next=head; head=this; }
  inline void pushOnto(Stackable<Node>*& head , Stackable<Node>*& tail) { pushOnto(head); if (tail==NULL) tail=this; }
  inline void appendTo(Stackable<Node>*& head , Stackable<Node>*& tail)
  {
    next=NULL;
    if (head==NULL) head = this;
    else tail->next = this;
    tail=this;
  }
  
  static inline void pushList(Stackable<Node>*& head , Stackable<Node> *list)
  {
    if (list==NULL) return;
    else if (head==NULL) head=list;
    else
    {
      Stackable<Node> *tmp=list;
      for ( ; tmp->next!=NULL ; tmp=tmp->next) { }
      tmp->next=head;
      head=list;
    }
  }
  static inline void pushList(Stackable<Node>*& head , Stackable<Node>*& tail , Stackable<Node>* list)
  {
    pushList(head,list);
    if (tail==NULL && list!=NULL) for (tail=head ; tail->next!=NULL; tail=tail->next) { }
  }
  static inline void appendList(Stackable<Node>*& head , Stackable<Node>*& tail , Stackable<Node> *list)
  {
    if (list==NULL) return;
    else if (head==NULL) head=tail=list;
    else
    {
      Stackable<Node> *tmp=list;
      for ( ; tmp->next!=NULL ; tmp=tmp->next) { }
      tail->next=list;
      tail=tmp;
    }
  }

  static inline Stackable<Node>* pop(Stackable<Node>*& head)
  {
    if (head==NULL) return NULL;
    Stackable<Node> *temp = head;
    head=head->next;
    temp->next=NULL;
    return temp;
  }
  static inline Stackable<Node>* pop(Stackable<Node>*& head , Stackable<Node>*& tail)
  {
    Stackable<Node> *temp = pop(head);
    if (head==NULL) tail=NULL;
    return temp;
  }
  
  // This version of mergeSort does comparison with < (on data field)
  static Stackable<Node>* mergeSort( Stackable<Node> *head )
  {
    Stackable<Node> *subhead[32],*subtail[32];  /* Allow merging of lists up to 2^32 long */
    Stackable<Node> *result,*left,*right,*tail,*old_tail;
    int depth,max_depth;
    
    max_depth = 0;
    subhead[0] = subtail[0] = NULL;
    result = NULL;
    while (head!=NULL)
    {
      // Sort up to 4 at a time by hand
      if (head->next!=NULL)
      {
        if (head->next->data < head->data)
        {
          tail=head;
          result=head->next;
          head=result->next;
          result->next=tail;
        }
        else
        {
          result=head;
          tail=head->next;
          head=tail->next;
        }
        tail->next=NULL;
        if (head!=NULL)
        {
          if (head->data < result->data)
          {
            left=result;
            result=head;
            head=head->next;
            result->next=left;
          }
          else if (head->data < tail->data)
          {
            left=head;
            head=head->next;
            result->next=left;
            left->next=tail;
          }
          else
          {
            left=tail;
            tail=head;
            head=head->next;
            left->next=tail;
            tail->next=NULL;
          }
          if (head!=NULL)
          {
            if (head->data < left->data)
            {
              if (head->data < result->data)
              {
                left=result;
                result=head;
                head=head->next;
                result->next=left;
              }
              else
              {
                right=left;
                left=head;
                head=head->next;
                result->next=left;
                left->next=right;
              }
            }
            else
            {
              if (head->data < tail->data)
              {
                right=head;
                head=head->next;
                left->next=right;
                right->next=tail;
              }
              else
              {
                right=tail;
                tail=head;
                head=head->next;
                right->next=tail;
                tail->next=NULL;
              }
            }
          }
        }
      }
      else
      {
        result=tail=head;
        head=NULL;
      }
      
      // Now merge
      for (depth=0 ; depth<max_depth; depth++) // Merge
      {
        if (subhead[depth]==NULL)  // Don't have anything to merge
        {
          if (head!=NULL) break;                   // Grab more if we can
          while (subhead[depth]==NULL) depth++;    // Otherwise we need to finish merging all the way down
        }
        if (!(result->data < subtail[depth]->data))  // No merging needed, just put lists together
        {
          subtail[depth]->next=result;
          result=subhead[depth];
          subtail[depth] = subhead[depth] = NULL;
        }
        else // Actually need to merge
        {
          left=subhead[depth];
          right=result;
          old_tail = tail;
          if (right->data < left->data)
          {
            result = tail = right;
            right=right->next;
          }
          else
          {
            result = tail = left;
            left=left->next;
          }
          while ( left!=NULL && right!=NULL )
          {
            if (right->data < left->data)
            {
              tail->next = right;
              right = right->next;
            }
            else
            {
              tail->next = left;
              left = left->next;
            }
            tail=tail->next;
            tail->next=NULL;
          }
          if (left != NULL)
          {
            tail->next=left;
            tail=subtail[depth];
          }
          else // right != NULL
          {
            tail->next = right;
            tail=old_tail;
          }
          subhead[depth] = subtail[depth] = NULL;
        }
      }
      
      // If we can get more, we need to store what we've done so far
      if (head!=NULL)
      {
        if (depth==max_depth)
        {
          max_depth++;
          subhead[max_depth] = subtail[max_depth] = NULL;
        }
        subhead[depth] = result;
        subtail[depth] = tail;
      }
    }
    return result;
  }
  
  // This version of mergeSort does comparison with a function
  static Stackable<Node>* mergeSort( Stackable<Node> *head , int (*lessThan)(const Node&,const Node&,void*) , void *data )
  {
    Stackable<Node> *subhead[32],*subtail[32];  /* Allow merging of lists up to 2^32 long */
    Stackable<Node> *result,*left,*right,*tail,*old_tail;
    int depth,max_depth;
    
    max_depth = 0;
    subhead[0] = subtail[0] = NULL;
    result = NULL;
    while (head!=NULL)
    {
      // Sort up to 4 at a time by hand
      if (head->next!=NULL)
      {
        if ((*lessThan)(head->next->data,head->data,data))
        {
          tail=head;
          result=head->next;
          head=result->next;
          result->next=tail;
        }
        else
        {
          result=head;
          tail=head->next;
          head=tail->next;
        }
        tail->next=NULL;
        if (head!=NULL)
        {
          if ((*lessThan)(head->data,result->data,data))
          {
            left=result;
            result=head;
            head=head->next;
            result->next=left;
          }
          else if ((*lessThan)(head->data,tail->data,data))
          {
            left=head;
            head=head->next;
            result->next=left;
            left->next=tail;
          }
          else
          {
            left=tail;
            tail=head;
            head=head->next;
            left->next=tail;
            tail->next=NULL;
          }
          if (head!=NULL)
          {
            if ((*lessThan)(head->data,left->data,data))
            {
              if ((*lessThan)(head->data,result->data,data))
              {
                left=result;
                result=head;
                head=head->next;
                result->next=left;
              }
              else
              {
                right=left;
                left=head;
                head=head->next;
                result->next=left;
                left->next=right;
              }
            }
            else
            {
              if ((*lessThan)(head->data,tail->data,data))
              {
                right=head;
                head=head->next;
                left->next=right;
                right->next=tail;
              }
              else
              {
                right=tail;
                tail=head;
                head=head->next;
                right->next=tail;
                tail->next=NULL;
              }
            }
          }
        }
      }
      else
      {
        result=tail=head;
        head=NULL;
      }
      
      // Now merge
      for (depth=0 ; depth<max_depth; depth++) // Merge
      {
        if (subhead[depth]==NULL)  // Don't have anything to merge
        {
          if (head!=NULL) break;                   // Grab more if we can
          while (subhead[depth]==NULL) depth++;    // Otherwise we need to finish merging all the way down
        }
        if (!(*lessThan)(result->data,subtail[depth]->data,data))  // No merging needed, just put lists together
        {
          subtail[depth]->next=result;
          result=subhead[depth];
          subtail[depth] = subhead[depth] = NULL;
        }
        else // Actually need to merge
        {
          left=subhead[depth];
          right=result;
          old_tail = tail;
          if ((*lessThan)(right->data,left->data,data))
          {
            result = tail = right;
            right=right->next;
          }
          else
          {
            result = tail = left;
            left=left->next;
          }
          while ( left!=NULL && right!=NULL )
          {
            if ((*lessThan)(right->data,left->data,data))
            {
              tail->next = right;
              right = right->next;
            }
            else
            {
              tail->next = left;
              left = left->next;
            }
            tail=tail->next;
            tail->next=NULL;
          }
          if (left != NULL)
          {
            tail->next=left;
            tail=subtail[depth];
          }
          else // right != NULL
          {
            tail->next = right;
            tail=old_tail;
          }
          subhead[depth] = subtail[depth] = NULL;
        }
      }
      
      // If we can get more, we need to store what we've done so far
      if (head!=NULL)
      {
        if (depth==max_depth)
        {
          max_depth++;
          subhead[max_depth] = subtail[max_depth] = NULL;
        }
        subhead[depth] = result;
        subtail[depth] = tail;
      }
    }
    return result;
  }

  static void join(Stackable<Node>*& head , Stackable<Node>*& tail , Stackable<Node>* join_head , Stackable<Node>* join_tail)
  {
    if (join_head==NULL) return;
    if (head==NULL) { head=join_head; tail=join_tail; }
    else { tail->next=join_head; tail=join_tail; }
  }
};



/****************************************************************
                   Doubly Linked List Template
****************************************************************/

/* NOTE **
         * Listable only supports NULL-terminated lists, not circular ones.
** NOTE **
         * Listable uses the same high-performance mergesort as Stackable.
** NOTE **
         * Lists should always be used with both head and tail.  Watch out for
         * the helper member functions that don't pay attention to heads and
         * tails!  Only use these in special situations.
** NOTE */

// Class to wrap another class or struct in a doubly linked list
template <class Node> class Listable
{
public:
  Listable<Node> *next;
  Listable<Node> *prev;
  Node data;
  
  Listable() { }
  Listable(const Node& node_data) : next(NULL),prev(NULL),data(node_data) { }
  
  // Utility methods that insert into a non-NULL list and ignore head/tail
  inline void addBefore(Listable<Node> *item)
  {
    prev=item->prev;
    next=item;
    if (prev!=NULL) prev->next=this;
    item->prev=this;
  }
  inline void addAfter(Listable<Node> *item)
  {
    next=item->next;
    prev=item;
    if (next!=NULL) next->prev=this;
    item->next=this;
  }
  
  // Utility methods that remove themselves and return the next or previous
  // element, respectively; if that's NULL, they go the other way
  inline Listable<Node>* extractR()
  {
    if (prev!=NULL) prev->next=next;
    if (next!=NULL)
    {
      next->prev=prev;
      return next;
    }
    else return prev;
  }
  inline Listable<Node>* extractL()
  {
    if (next!=NULL) next->prev=prev;
    if (prev!=NULL)
    {
      prev->next=next;
      return prev;
    }
    else return next;
  }
  // Or just remove self without returning anything
  inline void extract()
  {
    if (next!=NULL) next->prev = prev;
    if (prev!=NULL) prev->next = next;
  }
  
  // Standard addition of nodes to list
  inline void pushOnto(Listable<Node>*& head , Listable<Node>*& tail)
  {
    if (head==NULL) { next=prev=NULL; head=tail=this; }
    else { addBefore(head); head=this; }
  }
  
  inline void appendTo(Listable<Node>*& head , Listable<Node>*& tail)
  {
    if (head==NULL) { next=prev=NULL; head=tail=this; }
    else { addAfter(tail); tail=this; }
  }
  
  // Removal from the list (returns node after extracted one)
  inline Listable<Node>* extractR(Listable<Node>*& head , Listable<Node>*& tail)
  {
    Listable<Node>* tmp = extractR();
    if (tmp==NULL) head=tail=NULL;
    else if (head==this) head=tmp;
    else if (tail==this) tail=tmp;
    return tmp;
  }
  inline Listable<Node>* extractR(Listable<Node>*& head)
  {
    Listable<Node>* tmp = extractR();
    if (tmp==NULL) head=NULL;
    else if (head==this) head=tmp;
    return tmp;
  }
  // Removal returning the one before
  inline Listable<Node>* extractL(Listable<Node>*& head , Listable<Node>*& tail)
  {
    Listable<Node>* tmp = extractL();
    if (tmp==NULL) head=tail=NULL;
    else if (head==this) head=tmp;
    else if (tail==this) tail=tmp;
    return tmp;
  }
  inline Listable<Node>* extractL(Listable<Node>*& head)
  {
    Listable<Node>* tmp = extractL();
    if (tmp==NULL) head=NULL;
    else if (head==this) head=tmp;
    return tmp;
  }
  // Removal without returning anything
  inline void extract(Listable<Node>*& head , Listable<Node>*& tail)
  {
    if (this==head) head=head->next;
    if (this==tail) tail=tail->prev;
    extract();
  }
  inline void extract(Listable<Node>*& head)
  {
    if (this==head) head=head->next;
    extract();
  }
  
  // Mergesort with < as comparator
  static void mergeSort( Listable<Node>*& head , Listable<Node>*& tail )
  {
    Listable<Node> *subhead[32],*subtail[32];  /* Allow merging of lists up to 2^32 long */
    Listable<Node> *result,*left,*right,*old_tail;
    int depth,max_depth;
    
    max_depth = 0;
    subhead[0] = subtail[0] = NULL;
    result = NULL;
    while (head!=NULL)
    {
      // Sort up to 4 at a time by hand
      if (head->next!=NULL)
      {
        if (head->next->data < head->data)
        {
          tail=head;
          result=head->next;
          head=result->next;
          result->next=tail;
        }
        else
        {
          result=head;
          tail=head->next;
          head=tail->next;
        }
        tail->next=NULL;
        if (head!=NULL)
        {
          if (head->data < result->data)
          {
            left=result;
            result=head;
            head=head->next;
            result->next=left;
          }
          else if (head->data < tail->data)
          {
            left=head;
            head=head->next;
            result->next=left;
            left->next=tail;
          }
          else
          {
            left=tail;
            tail=head;
            head=head->next;
            left->next=tail;
            tail->next=NULL;
          }
          if (head!=NULL)
          {
            if (head->data < left->data)
            {
              if (head->data < result->data)
              {
                left=result;
                result=head;
                head=head->next;
                result->next=left;
              }
              else
              {
                right=left;
                left=head;
                head=head->next;
                result->next=left;
                left->next=right;
              }
            }
            else
            {
              if (head->data < tail->data)
              {
                right=head;
                head=head->next;
                left->next=right;
                right->next=tail;
              }
              else
              {
                right=tail;
                tail=head;
                head=head->next;
                right->next=tail;
                tail->next=NULL;
              }
            }
          }
        }
      }
      else
      {
        result=tail=head;
        head=NULL;
      }
      
      // Now merge
      for (depth=0 ; depth<max_depth; depth++) // Merge
      {
        if (subhead[depth]==NULL)  // Don't have anything to merge
        {
          if (head!=NULL) break;                   // Grab more if we can
          while (subhead[depth]==NULL) depth++;    // Otherwise we need to finish merging all the way down
        }
        if (!(result->data < subtail[depth]->data))  // No merging needed, just put lists together
        {
          subtail[depth]->next=result;
          result=subhead[depth];
          subtail[depth] = subhead[depth] = NULL;
        }
        else // Actually need to merge
        {
          left=subhead[depth];
          right=result;
          old_tail = tail;
          if (right->data < left->data)
          {
            result = tail = right;
            right=right->next;
          }
          else
          {
            result = tail = left;
            left=left->next;
          }
          while ( left!=NULL && right!=NULL )
          {
            if (right->data < left->data)
            {
              tail->next = right;
              right = right->next;
            }
            else
            {
              tail->next = left;
              left = left->next;
            }
            tail=tail->next;
            tail->next=NULL;
          }
          if (left != NULL)
          {
            tail->next=left;
            tail=subtail[depth];
          }
          else // right != NULL
          {
            tail->next = right;
            tail=old_tail;
          }
          subhead[depth] = subtail[depth] = NULL;
        }
      }
      
      // If we can get more, we need to store what we've done so far
      if (head!=NULL)
      {
        if (depth==max_depth)
        {
          max_depth++;
          subhead[max_depth] = subtail[max_depth] = NULL;
        }
        subhead[depth] = result;
        subtail[depth] = tail;
      }
    }
    
    // Fix up doubly-linked list to be doubly-linked properly again
    head=result;
    head->prev=NULL;
    for ( tail=head,result=head->next ; result!=NULL ; tail=result,result=result->next ) { result->prev = tail; }
  }

  // Mergesort with function as comparator
  static void mergeSort(Listable<Node>*& head , Listable<Node>*& tail , int (*lessThan)(const Node&,const Node&,void*) , void *data)
  {
    Listable<Node> *subhead[32],*subtail[32];  /* Allow merging of lists up to 2^32 long */
    Listable<Node> *result,*left,*right,*old_tail;
    int depth,max_depth;
    
    max_depth = 0;
    subhead[0] = subtail[0] = NULL;
    result = NULL;
    while (head!=NULL)
    {
      // Sort up to 4 at a time by hand
      if (head->next!=NULL)
      {
        if ((*lessThan)(head->next->data,head->data,data))
        {
          tail=head;
          result=head->next;
          head=result->next;
          result->next=tail;
        }
        else
        {
          result=head;
          tail=head->next;
          head=tail->next;
        }
        tail->next=NULL;
        if (head!=NULL)
        {
          if ((*lessThan)(head->data,result->data,data))
          {
            left=result;
            result=head;
            head=head->next;
            result->next=left;
          }
          else if ((*lessThan)(head->data,tail->data,data))
          {
            left=head;
            head=head->next;
            result->next=left;
            left->next=tail;
          }
          else
          {
            left=tail;
            tail=head;
            head=head->next;
            left->next=tail;
            tail->next=NULL;
          }
          if (head!=NULL)
          {
            if ((*lessThan)(head->data,left->data,data))
            {
              if ((*lessThan)(head->data,result->data,data))
              {
                left=result;
                result=head;
                head=head->next;
                result->next=left;
              }
              else
              {
                right=left;
                left=head;
                head=head->next;
                result->next=left;
                left->next=right;
              }
            }
            else
            {
              if ((*lessThan)(head->data,tail->data,data))
              {
                right=head;
                head=head->next;
                left->next=right;
                right->next=tail;
              }
              else
              {
                right=tail;
                tail=head;
                head=head->next;
                right->next=tail;
                tail->next=NULL;
              }
            }
          }
        }
      }
      else
      {
        result=tail=head;
        head=NULL;
      }
      
      // Now merge
      for (depth=0 ; depth<max_depth; depth++) // Merge
      {
        if (subhead[depth]==NULL)  // Don't have anything to merge
        {
          if (head!=NULL) break;                   // Grab more if we can
          while (subhead[depth]==NULL) depth++;    // Otherwise we need to finish merging all the way down
        }
        if (!(*lessThan)(result->data,subtail[depth]->data,data))  // No merging needed, just put lists together
        {
          subtail[depth]->next=result;
          result=subhead[depth];
          subtail[depth] = subhead[depth] = NULL;
        }
        else // Actually need to merge
        {
          left=subhead[depth];
          right=result;
          old_tail = tail;
          if ((*lessThan)(right->data,left->data,data))
          {
            result = tail = right;
            right=right->next;
          }
          else
          {
            result = tail = left;
            left=left->next;
          }
          while ( left!=NULL && right!=NULL )
          {
            if ((*lessThan)(right->data,left->data,data))
            {
              tail->next = right;
              right = right->next;
            }
            else
            {
              tail->next = left;
              left = left->next;
            }
            tail=tail->next;
            tail->next=NULL;
          }
          if (left != NULL)
          {
            tail->next=left;
            tail=subtail[depth];
          }
          else // right != NULL
          {
            tail->next = right;
            tail=old_tail;
          }
          subhead[depth] = subtail[depth] = NULL;
        }
      }
      
      // If we can get more, we need to store what we've done so far
      if (head!=NULL)
      {
        if (depth==max_depth)
        {
          max_depth++;
          subhead[max_depth] = subtail[max_depth] = NULL;
        }
        subhead[depth] = result;
        subtail[depth] = tail;
      }
    }

    // Fix up doubly-linked list to be doubly-linked properly again
    head=result;
    head->prev=NULL;
    for ( tail=head,result=head->next ; result!=NULL ; tail=result,result=result->next ) { result->prev = tail; }
  }
  
  static void joinR(Listable<Node>*& head , Listable<Node>*& tail , Listable<Node>* join_head , Listable<Node>* join_tail)
  {
    if (join_head==NULL) return;
    if (head==NULL) { head=join_head; tail=join_tail; }
    else
    {
      tail->next = join_head;
      join_head->prev = tail;
      tail = join_tail;
    }
  }
  static void joinL(Listable<Node>*& head , Listable<Node>*& tail , Listable<Node>* join_head , Listable<Node>* join_tail)
  {
    if (join_head==NULL) return;
    if (head==NULL) { head=join_head; tail=join_tail; }
    else
    {
      join_tail->next = head;
      head->prev = join_tail;
      head = join_head;
    }
  }
};



/****************************************************************
                    Unit Test-Style Functions
****************************************************************/

int test_mwt_lists_stackable();
int test_mwt_lists_listable();
int test_mwt_lists();


#endif

