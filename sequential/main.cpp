#include<iostream>
#include<limits>
#include<cfloat>
#include <array>
#include <algorithm>
#include <cmath>
#include "settings.h"

#include "HeapAndStructVector.h"
#include <cassert>
#include <vector>
#include <fstream>
#include <chrono>
#include <functional>

//Helper function to return the respective array index
int arr_index(const int x, const int y, const int z)
{
    return x+y*settings::x_grid_size+z*settings::x_grid_size*settings::y_grid_size;
}

std::vector<int> get_neighbors(const int x, const int y, const int z)
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

double solve(double speed_array[settings::total_grid_size],bool const accepted_array[settings::total_grid_size], const int x, const int y, const int z){
    //TODO maybe improve or  check compatability of old version
    double temp_res = std::numeric_limits<double>::infinity();
    double speed= speed_array[arr_index(x,y,z)];
    assert(!accepted_array[arr_index(x,y,z)]);
    //If the speed is zero(outside of original domain) we should instantly return infinity
    if(speed == 0){
        std::cout<<"ABSDHABD"<<std::endl;
        bool check = accepted_array[arr_index(x,y,z)];
        return temp_res;
    }
    int node_index = arr_index(x,y,z);
    double min_res_array[3];
    double h_array[3];

    //Check the x direction to set ψ_1 and h_array[0]
    int d{0};
    if(x > 0){
        int curr_index = arr_index(x-1,y,z);
        if(accepted_array[curr_index]){
            d=-1;
        }
    }
    if(x<settings::x_grid_size-1){
        int curr_index = arr_index(x+1,y,z);
        if(accepted_array[curr_index]){
            if(d==0){
                d=1;
            }
            else if(speed_array[curr_index]< speed_array[arr_index(x-1,y,z)]){
                d=1;
            }
        }
    }
    if(d!=0){
        min_res_array[0]=speed_array[arr_index(x+d,y,z)];
        h_array[0]=pow(settings::h,-1);
    }
    else{
        min_res_array[0]=0;
        h_array[0]=0;
    }
    //Check the y direction to set ψ_2 and h_array[1]
    d=0;
    if(y > 0){
        int curr_index = arr_index(x,y-1,z);
        if(accepted_array[curr_index]){
            d=-1;
        }
    }
    if(y<settings::y_grid_size-1){
        int curr_index = arr_index(x,y+1,z);
        if(accepted_array[curr_index]){
            if(d==0){
                d=1;
            }
            else if(speed_array[curr_index]< speed_array[arr_index(x,y-1,z)]){
                d=1;
            }
        }
    }
    if(d!=0){
        min_res_array[1]=speed_array[arr_index(x,y+d,z)];
        h_array[1]=pow(settings::h,-1);
    }
    else{
        min_res_array[1]=0;
        h_array[1]=0;
    }
    //Check the z direction to set ψ_3 and h_array[2]
    d=0;
    if(z > 0){
        int curr_index = arr_index(x,y,z-1);
        if(accepted_array[curr_index]){
            d=-1;
        }
    }
    if(z<settings::z_grid_size-1){
        int curr_index = arr_index(x,y,z+1);
        if(accepted_array[curr_index]){
            if(d==0){
                d=1;
            }
            else if(speed_array[curr_index]< speed_array[arr_index(x,y,z-1)]){
                d=1;
            }
        }
    }
    if(d!=0){
        min_res_array[2]=speed_array[arr_index(x,y,z+d)];
        h_array[2]=pow(settings::h,-1);
    }
    else{
        min_res_array[2]=0;
        h_array[2]=0;
    }


    int num_dir{0};
    for(int i=0; i<3;++i){
        if(h_array[i] >0){
            ++num_dir;
        }
    }

    while(num_dir!=0){
        double a = pow(h_array[0],2)+pow(h_array[1],2)+pow(h_array[2],2);
        double b = -2*(pow(h_array[0],2)*min_res_array[0]+pow(h_array[1],2)*min_res_array[1]+pow(h_array[2],2)*min_res_array[2]);
        //double speed=speed_funct(x,y,z);
        double c = pow(h_array[0] * min_res_array[0], 2) + pow(h_array[1] * min_res_array[1], 2) +pow(h_array[2] * min_res_array[2], 2) - pow(speed, -2);

        if((pow(b,2)-4*a*c)>=0){
            double psi_t = (-1*b+std::sqrt(pow(b,2)-4*a*c))/(2*a);
            if(min_res_array[0]< psi_t && min_res_array[1]<psi_t &&min_res_array[2]<psi_t){
                temp_res=std::min(psi_t, temp_res);
            }
        }
        //TODO hier die Warnung

        int index_to_del = std::distance(min_res_array, std::max_element(min_res_array, min_res_array + 3));
        min_res_array[index_to_del]=0;
        h_array[index_to_del]=0;
        num_dir -=1;
    }
    return temp_res;

}
void update_neighbors(MinHeap &h, double speed_array[settings::total_grid_size],bool const accepted_array[settings::total_grid_size], const int x, const int y, const int z)
{
    std::vector<int> neighbors = get_neighbors(x, y, z);

    for (std::size_t i = 0; i < neighbors.size(); i += 3){
        int x_neigh = neighbors[i];
        int y_neigh = neighbors[i+1];
        int z_neigh = neighbors[i+2];
        bool check2 = accepted_array[arr_index(x_neigh, y_neigh, z_neigh)];
        double checkspeed = speed_array[x_neigh, y_neigh, z_neigh];
        int c=0;
        if(!accepted_array[arr_index(x_neigh, y_neigh, z_neigh)]) {
            bool kek = accepted_array[arr_index(x_neigh, y_neigh, z_neigh)];
            int indexkek = h.get_heap_index(14, 16, 14);
            WeightedPoint currNode{x_neigh, y_neigh, z_neigh};
            double temp_weight = solve(std::ref(speed_array), accepted_array, x_neigh, y_neigh, z_neigh);
            if (h.get_heap_index(x_neigh, y_neigh, z_neigh) == -1) {
                //std::cout<<"INSERTED A KEY\n";
                currNode.weight = temp_weight;
                h.insertKey(currNode);
            }
            else if (temp_weight < h.weight_at_index(h.get_heap_index(x_neigh, y_neigh, z_neigh))) {
                int index = h.get_heap_index(x_neigh, y_neigh, z_neigh);
                h.decreaseKey(h.get_heap_index(x_neigh, y_neigh, z_neigh), temp_weight);
                index = h.get_heap_index(x_neigh, y_neigh, z_neigh);
                int c=0;
            }

        }
    }
}
void initialize(MinHeap &h, bool const accepted_array[settings::total_grid_size], double speed_array[settings::total_grid_size], int &accepted_counter)
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
                    speed_array[arr_index(x,y,z)]=0;
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
                if(accepted_array[arr_index(x,y,z)])
                {
                    update_neighbors(h, std::ref(speed_array), accepted_array, x,y,z);
                }
            }
        }
    }
}


void fast_marching(MinHeap &h, bool  accepted_array[settings::total_grid_size], double speed_array[settings::total_grid_size], int &accepted_counter)
{
    while(accepted_counter< settings::total_grid_size)
    {
        bool c= accepted_array[arr_index(14,16,14)];
        int index = h.get_heap_index(14, 16, 14);
        int index2 = h.get_heap_index(16,27,18);
        if(c){
            //std::cout<<accepted_counter<<std::endl;
        }
        //std::cout<<"GOT ONE\n";
        WeightedPoint a=h.extractMin();

        index = h.get_heap_index(14, 16, 14);
        //assert(accepted_array[arr_index(a.m_x,a.m_y,a.m_z)]!=true);
        //std::cout<<"DOESNT GET SENT";
        accepted_array[arr_index(a.m_x,a.m_y,a.m_z)]=true;
        speed_array[arr_index(a.m_x,a.m_y,a.m_z)] = a.weight;
        ++accepted_counter;
        update_neighbors(h, std::ref(speed_array), std::ref(accepted_array), a.m_x,a.m_y,a.m_z);
    }
    for ( int i =0; i< settings::total_grid_size;++i)
    {
        //std::cout<<speed_array[i]<<std::endl;
    }
}

bool in_barrier(int x, int y, int z)
{
    double temp_x= x*settings::h - 0.5;
    double temp_y= y*settings::h - 0.5;
    double temp_z= z*settings::h - 0.5;

    double w = 1.0/24;
    double big_r = std::sqrt(std::pow(temp_x,2)+ std::pow(temp_y,2)+ std::pow(temp_z,2));
    double small_r = std::sqrt(std::pow(temp_x,2)+ std::pow(temp_y,2));
    bool in_barrier_one = (0.15 < big_r && big_r< 0.15+w) && !(small_r<0.05 && temp_z<0);
    bool in_barrier_two = (0.25 < big_r && big_r< 0.25+w) && !(small_r<0.1 && temp_z>0);
    bool in_barrier_three = (0.35 < big_r && big_r< 0.35+w) && !(small_r<0.10 && temp_z<0);
    bool in_barrier_four = (0.45 < big_r && big_r< 0.45+w) && !(small_r<0.1 && temp_z>0);
    return in_barrier_one||in_barrier_two||in_barrier_four||in_barrier_three;
}
//speed and mask functions
double speed_funct(int x, int y, int z)
{
    //return 1+ 0.5*std::sin(20*M_PI*x)*std::sin(20*M_PI*y)*std::sin(20*M_PI*z);
    return 0.001*(pow(std::sin(x),2)+pow(std::cos(y),2)+0.1);
    //return 1.0;
    //if(in_barrier(x,y,z)) return 0.0001;
    //else return 1;
}
bool in_mask(int x, int y , int z)
{
    //cube=[15,24]^3
    bool x_cor = x<=24 && x>= 15;
    bool y_cor = y<=24 && y>= 15;
    bool z_cor = z<=24 && z>= 15;
    return x_cor && y_cor && z_cor;
    //return x==0&&y==0&&z==0;
    //bool in_ball= (pow(x*settings::h -0.5,2)+ pow(y*settings::h -0.5,2)+pow(z*settings::h -0.5,2))<=1.0/16;
    //bool in_cube = x*settings::h<=0.875&&x*settings::h >=0.625&& y*settings::h<=0.875&&y*settings::h >=0.625&& z*settings::h<=0.875&&z*settings::h >=0.625;
    //return in_ball;
    //return x== settings::x_grid_size/2&&y== settings::y_grid_size/2&&z== settings::z_grid_size/2;
}
int main()
{
    //test with speed 1 and start in a ball in center of the mesh, should return distance from origin
    try {
        int counter_1=0;
        //initialize all arrays
        //total size needed for heap for the lookup table

        MinHeap h(settings::total_grid_size);
        bool *mask_array{new bool[settings::total_grid_size]{}};
        double *speed_array{new double[settings::total_grid_size]{}};
        auto startTime = std::chrono::system_clock::now();
        for(int x=0; x<settings::x_grid_size; ++x){
            for(int y=0; y<settings::y_grid_size; ++y){
                for(int z=0; z<settings::z_grid_size; ++z){
                    mask_array[arr_index(x,y,z)]= in_mask(x,y,z);
                    speed_array[arr_index(x,y,z)]= speed_funct(x,y,z);
                    //if (mask_array[arr_index(x,y,z)]== in_mask(x,y,z)) {
                    //    ++counter_1;
                    //}
                }
            }
        }

        int accepted_counter{0};

        //initialize output
        std::ofstream myfile;
        myfile.open("distance_cube_point_1.txt");
        myfile << "Dimension information\n"<<settings::x_grid_size <<"\n"<<settings::y_grid_size<<"\n"<<settings::z_grid_size<<"\n";
        myfile << "Mask information\n";
        for ( int i =0; i< settings::total_grid_size;++i)
        {
            myfile << mask_array[i]<<"\n";
        }

        initialize(h, mask_array , std::ref(speed_array), accepted_counter);
        fast_marching(h,mask_array, std::ref(speed_array), accepted_counter);

        auto endTime = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = endTime - startTime;

        std::cout<<elapsed_seconds.count()<<std::endl;

        //print result to output
        myfile<<"Result information\n";
        for ( int i =0; i< settings::total_grid_size;++i)
        {
            myfile << speed_array[i]<<"\n";
            //std::cout<<speed_array[i]<<std::endl;
        }
        myfile.close();

    }


    catch (const char* exception)
    {
        std::cerr << "Error: " << exception << '\n';
    }
    return 0;
    /* //test with speed 1 and start in (0,0,0), should return distance from origin
     try {
         std::cout<<sizeof(WeightedPoint)<<"\n";
         //initialize all arrays
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

         //initialize output
         std::ofstream myfile;
         myfile.open("distance_cube_ver_8.txt");
         myfile << "Dimension information\n"<<settings::x_grid_size <<"\n"<<settings::y_grid_size<<"\n"<<settings::z_grid_size<<"\n";
         myfile << "Mask information\n";
         for ( int i =0; i< settings::total_grid_size;++i)
         {
             myfile << mask_array[i]<<"\n";
         }
         initialize(h, mask_array, speed_array, weight_array, accepted_counter);
         fast_marching(h, mask_array, speed_array, weight_array, accepted_counter);


         std::cout<<"deviation: "<<weight_array[settings::total_grid_size-1]-sqrt(3)*((settings::x_grid_size-1)*settings::h);
         //print result to output
         myfile<<"Result information\n";
         for ( int i =0; i< settings::total_grid_size;++i)
         {
             myfile << weight_array[i]<<"\n";
         }
         myfile.close();
     }


     catch (const char* exception)
     {
         std::cerr << "Error: " << exception << '\n';
     }
     return 0;
     */
}

