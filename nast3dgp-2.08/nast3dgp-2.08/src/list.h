/* NaSt3DGP - The Parallel 3D Navier-Stokes Solver
 * Copyright (C) 2003 Institute for Numerical Simulation
 *                    University of Bonn
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIST_INCLUDED
#define LIST_INCLUDED
#include <stdio.h>
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif 


//! unsorted list 
/*!
  this template class provides an unsorted list with read/write functionality
 */
template <class T>
class NS_List {
protected:
     //! pointer to next element
     NS_List*   next;
public:
    //! data type listed
     T       data;
     //! constructor
     NS_List() {next=0;}
     //! destructor: delete all elements
     ~NS_List() {
	  NS_List* p=next,*p1;
	  while(p) {
	       p1=p->next;p->next=0;
	       delete p;p=p1;
	  }
     }        
     //! return pointer to next element
     NS_List*   GetNext() {return next;}
     //! empty the whole list
     void    Delete() {
	  NS_List* p=next,*p1;
	  while(p) {
	       p1=p->next;p->next=0;
	       delete p;p=p1;
	  }
	  Init();
     }
     void    Init() {next=0;}
     //! insert element at first position
     void    InsertFirst(NS_List* p) {
	  p->next=next;
	  next=p;
     }
#ifndef LISTNOFILE 
     //! write list to file
     int     WriteFile(FILE* f) {
	  NS_List* p=this;
	  while(p) {
	       if(fwrite((void**)&(p->next),sizeof(NS_List*),1,f)!=1) return FALSE;
	       if(!p->data.WriteFile(f)) return FALSE;
	       p=p->next;
	  }
	  return TRUE;
     }
     //! read list from file
     int     ReadFile(FILE* f) {
	  NS_List<T>* p=this;
	  int stop=FALSE;
	  if(fread(&next,sizeof(NS_List*),1,f)!=1) return FALSE;
	  if(!data.ReadFile(f)) return FALSE;
	  if(!next) stop=TRUE;p->next=0;
	  while(!stop) {
	       p=new NS_List<T>;
	       if(fread(&(p->next),sizeof(NS_List*),1,f)!=1) return FALSE;
	       if(!p->next) stop=TRUE;
	       if(!p->data.ReadFile(f)) return FALSE;                 
	       InsertFirst(p);   
	  }
	  return TRUE;
     }
     int     SkipFile(FILE* f) {
	  NS_List<T>* p=new NS_List<T>;
	  int stop=FALSE;
	  fflush(f);
	  if(fread(&(p->next),sizeof(NS_List*),1,f)!=1) return FALSE;
	  if(!data.ReadFile(f)) return FALSE;
                if(!p->next) stop=TRUE;
                delete p;
                while(!stop) {
		     p=new NS_List<T>;
		     if(fread(&(p->next),sizeof(NS_List*),1,f)!=1) return FALSE;
		     if(!p->next) stop=TRUE;
		     if(!p->data.ReadFile(f)) return FALSE;
		     delete p;
                }
                fflush(f);
                return TRUE;
     }
#endif
};



//!  sorted list
/*!
  this template class provides a sorted list.
  It inherits functionality from NS_List.
  Data is sorted based on the > operator for type T
 */
template<class T>
class SortedList:public NS_List<T> {
public:
    void    InsertSorted(SortedList* p) {
                SortedList* q=this,*r;
                do {
                    r=q;q=q->GetNext();
                } while(q!=NULL && ((p->data)>(q->data)));
                r->next=p;p->next=q;
            }
    SortedList* Find(const T& d) {
                SortedList* q=this;
                do {
                    q=q->GetNext();
                } while(q!=NULL && !((q->data)==d));
                return q;
            }
    SortedList* GetNext() {return (SortedList*)this->next;}
};

//! item in a linked list
/*!
  this class represents an item in a double linked list.
 */
template <class T>
class LinkedListItem {
public:
    //! the data which is saved in the list item
    T *data;
    //! a pointer to the next item in the list
    LinkedListItem<T> *next;
    //! a pointer to the previous item in the list
    LinkedListItem<T> *previous;
public:
    //! create a new item with data data
    LinkedListItem(T *data) {
	this->data = data;
	next = NULL;
	previous = NULL;
    }
    //! delete item and data
    ~LinkedListItem() {
	delete data;
    }

    //! return next list item
    inline LinkedListItem<T> *Next() const {
	return next;
    }

    //! return data
    inline T* Data() const {
	return data;
    }
};

//! double linked list
/*!
   this class represents a double linked list
 */
template <class T>
class LinkedList {
public:
    //! type definition of items in this list
    typedef LinkedListItem<T> item;
protected:
    //! count of items in the list
    int size;
    //! first particle in the list
    item *first;
    //! last particle in the list
    item *last;
public:
    //! create a new empty linked list
    LinkedList() {
	first = last = NULL;
	size = 0;
    }
    //! destroy the list, all items and all data
    ~LinkedList() {
	item *iter = first, *tmp;
	while (iter != NULL) {
	    tmp = iter;
	    iter = iter->next;
	    delete tmp;
	}
    }

    //! initialize list
    void Init() {
	first = last = NULL;
	size = 0;
    }

    //! return first list item
    inline item *First() const {
	return first;
    }

    //! return size of list
    inline int Size() const {
	return size;
    }

    //! append new item to the list
    void Append(T *data) {
	item *new_item = new item(data);

	size ++;
	if (last == NULL) {
	    first = last = new_item;
	} else {
	    new_item->previous = last;
	    last->next = new_item;
	    last = new_item;
	}
    }

    //! remove item from the list and delete it
    void Remove(item *i) {
	size --;
	
	if (i->previous == NULL) {
	    first = i->next;
	} else {
	    i->previous->next = i->next;
	}

	if (i->next == NULL) {
	    last = i->previous;
	} else {
	    i->next->previous = i->previous;
	}

	delete i;
    }

    //! remove all items from the list
    void RemoveAll() {
	item *iter = first, *tmp;

	while (iter != NULL) {
	    tmp = iter;
	    iter = iter->next;
	    delete tmp;
	}

	first = last = NULL;
	size = 0;
    }
};


#endif












