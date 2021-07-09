//
// Created by Moritz on 04.06.2021.
//


#include "HeapAndStructVector.h"

#include<iostream>
#include<limits>
#include<cfloat>
#include <array>
#include <vector>
int local_arr_index(int, int, int);//forward declaration

// Constructor: Builds a heap from a given array a[] of given size
MinHeap::MinHeap(const int cap)
{
    heap_size = 0;
    capacity = cap;
    std::vector<WeightedPoint> m_harr;
    table = new int[cap];
    for(int i=0; i<cap; ++i)
    {
        table[i]=-1;
    }

}

MinHeap::~MinHeap()
{
    delete[] table;
}

//returns the index of element (x,y,z)
int MinHeap::get_heap_index(const int x, const int y, const int z) const
{
    return table[local_arr_index(x,y,z)];
}

// Inserts a new point 'point'
void MinHeap::insertKey(const WeightedPoint point)
{
    if (heap_size == capacity)
    {
        std::cout << "\nOverflow: Could not insertKey\n";
        return;
    }

    // First insert the new key at the end
    ++heap_size;
    int i = heap_size - 1;
    harr.push_back(point);
    //harr[i] = point;
    table[local_arr_index(point.m_x,point.m_y,point.m_z)]=i;

    // Fix the min heap property if it is violated
    while (i != 0 && harr[parent(i)].weight > harr[i].weight)
    {
        swap(&harr[i], &harr[parent(i)]);
        i = parent(i);
    }
}

// Decreases value of key at index 'i' to new_val.  It is assumed that
// new_val is smaller than harr[i].
void MinHeap::decreaseKey( int i, double const new_val)
{
    harr[i].weight = new_val;
    while (i != 0 && harr[parent(i)].weight > harr[i].weight)
    {
        swap(&harr[i], &harr[parent(i)]);
        i = parent(i);
    }
}

// Method to remove minimum element (or root) from min heap
WeightedPoint MinHeap::extractMin()
{
    if (heap_size <= 0)
        throw "TRYING TO EXTRACT THE MINIMUM OF AN EMPTY HEAP";
    if (heap_size == 1)
    {
        heap_size--;
        table[local_arr_index(harr[0].m_x,harr[0].m_y,harr[0].m_z)]=-1;
        return harr[0];
    }

    // Store the minimum value, and remove it from heap
    WeightedPoint root = harr[0];
    harr[0] = harr[heap_size-1];
    harr.pop_back();
    heap_size--;
    MinHeapify(0);
    table[local_arr_index(root.m_x,root.m_y,root.m_z)]=-1;

    return root;
}


// This function deletes key at index i. It first reduced value to minus
// infinite, then calls extractMin()

void MinHeap::deleteKey(int i)
{
    decreaseKey(i, DBL_MIN);
    extractMin();
}

// A recursive method to heapify a subtree with the root at given index
// This method assumes that the subtrees are already heapified
void MinHeap::MinHeapify(int i)
{
    int l = left(i);
    int r = right(i);
    int smallest = i;
    if (l < heap_size && harr[l].weight < harr[i].weight)
        smallest = l;
    if (r < heap_size && harr[r].weight < harr[smallest].weight)
        smallest = r;
    if (smallest != i)
    {
        swap(&harr[i], &harr[smallest]);
        MinHeapify(smallest);
    }
}

// A utility function to swap two elements
void MinHeap::swap(WeightedPoint *x, WeightedPoint *y)
{
    int temp_index=table[local_arr_index((*x).m_x, (*x).m_y, (*x).m_z)];
    table[local_arr_index((*x).m_x, (*x).m_y, (*x).m_z)]=table[local_arr_index((*y).m_x, (*y).m_y, (*y).m_z)];
    table[local_arr_index((*y).m_x, (*y).m_y, (*y).m_z)]=temp_index;
    WeightedPoint temp = *x;
    *x = *y;
    *y = temp;
}
