#include<iostream>
#include<limits>
#include<cfloat>
#include <array>
#include <algorithm>
#include <cmath>
#include "settings.h"
#include "HeapAndStruct.h"
#include <cassert>
#include <vector>


//Helper function to return the respective array index
int arr_index(int x, int y, int z)
{
    return x+y*settings::x_grid_size+z*settings::x_grid_size*settings::y_grid_size;
}


std::vector<int> get_neighbors(int x, int y, int z)
{
//returns a vector of the coordinates of the neighbors in the format<x_1,y_1,z_1,x_2,....
    std::vector<int> res={};

    if(x>0)
    {
        res.push_back(x-1);
        res.push_back(y);
        res.push_back(z);
    }
    if(x<settings::x_grid_size-1)
    {
        res.push_back(x+1);
        res.push_back(y);
        res.push_back(z);
    }
    if(y>0)
    {
        res.push_back(x);
        res.push_back(y-1);
        res.push_back(z);
    }
    if(y<settings::y_grid_size-1)
    {
        res.push_back(x);
        res.push_back(y+1);
        res.push_back(z);
    }
    if(z>0)
    {
        res.push_back(x);
        res.push_back(y);
        res.push_back(z-1);
    }
    if(z<settings::z_grid_size-1)
    {
        res.push_back(x);
        res.push_back(y);
        res.push_back(z+1);
    }

    return res;
}

double solve_eikonal_quadratic_2d(double u_north, double u_south, double u_west, double u_east, double f)
    //solves the 2d version of the quadratic update formula
{
    double u_north_south = std::min(u_north, u_south);
    double u_west_east = std::min(u_east, u_west);

    double discriminant = pow(u_north_south+u_west_east,2)-2*(pow(u_north_south,2)+pow(u_west_east,2)-pow(settings::h/f,2));

    double res;

    //naive version:
    if(discriminant >=0)
    {
        res= (u_north_south+u_west_east+sqrt(discriminant))/2.0;

        if(res<= u_north_south || res <= u_west_east)//upwinding says that our value should only depend on the vertices with smaller values then itself if this is not the case we lower the dimension
        {
            res=std::min(u_north_south, u_west_east)+settings::h/f;
        }
    }
    else
    {
        res=std::min(u_north_south, u_west_east)+settings::h/f;
    }

    return res;
}

double solve_eikonal_quadratic_3d(double u_up, double u_down, double u_north, double u_south, double u_west, double u_east, double f)
{

    //solves the 3d version of the quadratic update formula
    double u_up_down = std::min(u_up, u_down);
    double u_north_south = std::min(u_north, u_south);
    double u_west_east = std::min(u_east, u_west);

    double discriminant = pow(u_up_down+u_north_south+u_west_east,2)-3*(pow(u_up_down,2)+pow(u_north_south,2)+pow(u_west_east,2)-pow(settings::h/f,2));
    double res;

    if(discriminant >=0)
    {
        res= (u_up_down+u_north_south+u_west_east+sqrt(discriminant))/3.0;
        if (res<= u_up_down || res <= u_north_south || res <= u_west_east )//upwinding says that our value should only depend on the vertices with smaller values then itself if this is not the case we lower the dimension
        {
            res = std::min({solve_eikonal_quadratic_2d(u_up, u_down, u_north, u_south, f),solve_eikonal_quadratic_2d(u_up, u_down, u_west, u_east, f), solve_eikonal_quadratic_2d(u_north, u_south, u_west, u_east, f)});

        }
    }
    else
    {
        res = std::min({solve_eikonal_quadratic_2d(u_up, u_down, u_north, u_south, f),solve_eikonal_quadratic_2d(u_up, u_down, u_west, u_east, f), solve_eikonal_quadratic_2d(u_north, u_south, u_west, u_east, f)});
    }
    return res;
}

std::array<double,6> prepare_quadratic(WeightedPoint currNode, double weight_array[settings::total_grid_size])//hier pass by reference
{
    //returns weight of neighboring indices in correct order for solve_eikonal_quadratic , if not in grid, gets set to std::numeric_limits<double>::infinity()
    std::array<double,6> res{};
    if(currNode.m_z==settings::z_grid_size-1)
    {
        res[0]=std::numeric_limits<double>::infinity();
    }
    else
    {
        res[0]=weight_array[arr_index(currNode.m_x,currNode.m_y, currNode.m_z+1)];
    }

    if(currNode.m_z==0)
    {
        res[1]=std::numeric_limits<double>::infinity();
    }
    else
    {
        res[1]=weight_array[arr_index(currNode.m_x,currNode.m_y, currNode.m_z-1)];
    }

    if(currNode.m_y==settings::y_grid_size-1)
    {
        res[2]=std::numeric_limits<double>::infinity();
    }
    else
    {
        res[2]=weight_array[arr_index(currNode.m_x,currNode.m_y+1, currNode.m_z)];
    }

    if(currNode.m_y==0)
    {
        res[3]=std::numeric_limits<double>::infinity();
    }
    else
    {
        res[3]=weight_array[arr_index(currNode.m_x,currNode.m_y-1, currNode.m_z)];
    }

    if(currNode.m_x==settings::x_grid_size-1)
    {
        res[4]=std::numeric_limits<double>::infinity();
    }
    else
    {
        res[4]=weight_array[arr_index(currNode.m_x+1,currNode.m_y, currNode.m_z)];
    }
    if(currNode.m_x==0)
    {
        res[5]=std::numeric_limits<double>::infinity();
    }
    else
    {
        res[5]=weight_array[arr_index(currNode.m_x-1,currNode.m_y, currNode.m_z)];
    }
    return res;
}

void initialize(MinHeap &h, bool  accepted_array[settings::total_grid_size], double speed_array[settings::total_grid_size], double weight_array[settings::total_grid_size], int &accepted_counter)
{

    //Set up the masked nodes and all other ones as far, by default the status is already "Far" and weight is std::numeric_limits<double>::infinity()
    for(int x=0; x<settings::x_grid_size; ++x)
    {
        for(int y=0; y<settings::y_grid_size; ++y)
        {
            for(int z=0; z<settings::z_grid_size; ++z)
            {
                if (accepted_array[arr_index(x,y,z)])
                {
                    ++accepted_counter;
                    weight_array[arr_index(x,y,z)]=0;
                }
            }
        }
    }
    //Set up the neighbouring nodes second
    for(int x=0; x<settings::x_grid_size; ++x)
    {
        for (int y = 0; y < settings::y_grid_size; ++y)
        {
            for (int z = 0; z < settings::z_grid_size; ++z)
            {
                if(!(accepted_array[arr_index(x,y,z)]))
                {
                    std::vector<int> neighbors = get_neighbors(x, y, z);
                    for (std::size_t i = 0; i < neighbors.size(); i += 3)//see definition of get_neighbors
                    {
                        //WeightedPoint neighbor = weight_array[arr_index(neighbors[i], neighbors[i + 1], neighbors[i + 2])];
                        if (accepted_array[arr_index(neighbors[i], neighbors[i+1], neighbors[i+2])])//if this neighbor is accepted our node is considered and we update its value
                        {
                            //std::cout<<"The neighbor"<<neighbor.m_x<<neighbor.m_y<<neighbor.m_z<<"is accepted/not accepted (1/0): "<<accepted_array[arr_index(neighbor.m_x, neighbor.m_y, neighbor.m_z)]<<"\n";
                            //std::cout<<"The current point is: "<< x<<y<<z<<"\n";

                            //updates the value of the current point
                            WeightedPoint currNode{x,y,z};
                            std::array<double,6> params = prepare_quadratic(currNode, weight_array);
                            double temp = solve_eikonal_quadratic_3d(params[0], params[1], params[2], params[3], params[4], params[5], speed_array[arr_index(x,y,z)]);
                            if (temp < weight_array[arr_index(x,y,z)])
                            {
                                weight_array[arr_index(x,y,z)]=temp;
                                if(h.get_heap_index(x,y,z)==-1)
                                {
                                    //std::cout<<"INSERTED A KEY\n";
                                    currNode.weight=temp;
                                    h.insertKey(currNode);
                                }
                                else
                                {
                                    h.decreaseKey(h.get_heap_index(x,y,z), temp);
                                }
                            }
                            //it is enough if we have one accepted neighbor
                            break;
                        }
                    }
                }
            }
        }
    }
}


void fast_marching(MinHeap &h, bool  accepted_array[settings::total_grid_size], double speed_array[settings::total_grid_size], double weight_array[settings::total_grid_size], int &accepted_counter)
{
    while(accepted_counter< settings::total_grid_size)
    {
        //std::cout<<"GOT ONE\n";
        WeightedPoint a=h.extractMin();
        //std::cout<<"DOESNT GET SENT";
        accepted_array[arr_index(a.m_x,a.m_y,a.m_z)]=true;
        ++accepted_counter;

        std::vector<int> neighbors=get_neighbors(a.m_x, a.m_y,a.m_z);
        for(std::size_t i = 0; i < neighbors.size(); i+=3)
        {
            WeightedPoint currNode{neighbors[i], neighbors[i+1], neighbors[i+2], weight_array[arr_index(neighbors[i], neighbors[i+1], neighbors[i+2])]};
            //WeightedPoint currNode = weight_array[arr_index(neighbors[i], neighbors[i+1], neighbors[i+2])];
           if (!accepted_array[arr_index(currNode.m_x, currNode.m_y, currNode.m_z)])//if the current point is not accepted we update its value
            {
                std::array<double,6> params = prepare_quadratic(currNode, weight_array);
                double temp = solve_eikonal_quadratic_3d(params[0], params[1], params[2], params[3], params[4], params[5], speed_array[arr_index(currNode.m_x, currNode.m_y, currNode.m_z)]);
                if (temp < currNode.weight)
                {
                    weight_array[arr_index(currNode.m_x, currNode.m_y, currNode.m_z)]=temp;


                    if(h.get_heap_index(currNode.m_x, currNode.m_y, currNode.m_z)==-1)//if the considered point is not in the heap, we insert it
                    {
                        currNode.weight=temp;
                        h.insertKey(currNode);
                    }
                    else//if it is we update its value(which is guaranteed to be smaller)
                    {
                        h.decreaseKey(h.get_heap_index(currNode.m_x, currNode.m_y, currNode.m_z), temp);
                    }
                }
            }
        }
    }
}

int main()
{
    //test with speed 1 and start in (0,0,0), should return distance from origin
    try {
        MinHeap h(settings::total_grid_size);
        bool *mask_array{new bool[settings::total_grid_size]{}};
        mask_array[0] = true;
        double *speed_array{new double[settings::total_grid_size]{}};
        for (int i = 0; i < settings::total_grid_size; ++i) {
            speed_array[i] = 1;
        }

        double *weight_array{new double[settings::total_grid_size]{}};
        for (int i = 0; i < settings::total_grid_size; ++i) {
            weight_array[i] = std::numeric_limits<double>::infinity();
        }
        int accepted_counter{0};

        initialize(h, mask_array, speed_array, weight_array, accepted_counter);
        fast_marching(h, mask_array, speed_array, weight_array, accepted_counter);


        std::cout<<"deviation: "<<weight_array[settings::total_grid_size-1]-sqrt(3)*(settings::x_grid_size-1);

    }
    catch (const char* exception)
    {
        std::cerr << "Error: " << exception << '\n';
    }
    return 0;
}

