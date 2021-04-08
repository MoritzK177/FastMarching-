//
// Created by Moritz on 26.03.2021.
//

#ifndef FASTMARCHING_HEAPANDSTRUCT_H
#define FASTMARCHING_HEAPANDSTRUCT_H

#include<cfloat>
#include <string>
#include "settings.h"
#include <limits>

//A simple representation of the points with weight
struct WeightedPoint
{
    int m_x{0};
    int m_y{0};
    int m_z{0};
    //std::string status{"Far"};
    double weight{std::numeric_limits<double>::infinity()};
};



// A class for Min Heap
class MinHeap
{
    WeightedPoint* harr; // pointer to array of elements in heap
    int* table; // pointer to table of indices for a given element (table[x+y*x_grid_size+z*x_grid_size*y_grid_size]=harr[(x,y,z)] if in heap, -1 else
    int capacity; // maximum possible size of min heap
    int heap_size; // Current number of elements in min heap
public:
    // Constructor
    MinHeap(int capacity);

    //Destructor
    ~MinHeap();

    // to heapify a subtree with the root at given index
    void MinHeapify(int );

    int parent(int i) { return (i-1)/2; }

    // to get index of left child of node at index i
    int left(int i) { return (2*i + 1); }

    // to get index of right child of node at index i
    int right(int i) { return (2*i + 2); }

    //returns the index of element (x,y,z)
    int get_heap_index(int x, int y, int z);

    // to extract the root which is the minimum element
    WeightedPoint extractMin();

    // Decreases key value of key at index i to new_val
    void decreaseKey(int i, double new_val);

    // Returns the minimum key (key at root) from min heap
    WeightedPoint getMin() { return harr[0]; }

    // Deletes a key stored at index i
    void deleteKey(int i);

    // Inserts a new key 'k'
    void insertKey(WeightedPoint k);

    // Prototype of a utility function to swap two points
    void swap(WeightedPoint *x, WeightedPoint *y);
};
#endif //FASTMARCHING_HEAPANDSTRUCT_H
