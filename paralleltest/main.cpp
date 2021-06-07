#include<iostream>
#include<limits>
#include<cfloat>
#include <array>
#include <algorithm>
#include <cmath>
#include "settings.h"
#include <mutex>
#include <thread>
#include "HeapAndStructVector.h"
#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <utility>



/* legend of the statuses:
 * 0 = FAR
 * 1 = BAND_NEW
 * 2 = BAND_OLD
 * 3 = KNOWN_FIX
 * 4 = KNOWN_NEW
 * 5 = KNOWN_OLD */

//speed and mask functions
double speed_funct(int x, int y, int z)
{
    //return std::sin(x)+2*std::sin(y)+4*std::sin(z)+8;
    //return 0.001* (pow(std::sin(x),2)+pow(std::cos(y),2)+0.1);
    return 1;
}
bool in_mask(int x, int y , int z)
{
    //cube=[2]^3
    bool x_cor = x<=1 && x>= 1;
    bool y_cor = y<=1 && y>= 1;
    bool z_cor = z<=1 && z>= 1;
    return x_cor && y_cor && z_cor;
}

//Helper function to return the respective array indices
int global_arr_index(const int x, const int y, const int z)
{
    return x+y*settings::x_global_grid_size+z*settings::x_global_grid_size*settings::y_global_grid_size;
}

int local_arr_index(const int x, const int y, const int z)
{
    return x+y*settings::x_local_grid_size+z*settings::x_local_grid_size*settings::y_local_grid_size;
}
int process_index(const int x, const int y, const int z)
{
    return x+y*settings::x_num_processes+z*settings::x_num_processes*settings::y_num_processes;
}

void initialize_subdomain(int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes])
{
    //some bools to check which point is outside of the global domain
    SubdomainData subdomain = subdomain_array[subdomain_index];

    bool at_x_lower_border{subdomain.x_offset == 0};
    bool at_x_upper_border{subdomain.x_offset == settings::x_global_grid_size-(settings::x_local_grid_size-2)};
    bool at_y_lower_border{subdomain.y_offset == 0};
    bool at_y_upper_border{subdomain.y_offset == settings::y_global_grid_size-(settings::y_local_grid_size-2)};
    bool at_z_lower_border{subdomain.z_offset == 0};
    bool at_z_upper_border{subdomain.z_offset == settings::z_global_grid_size-(settings::z_local_grid_size-2)};

   //loops to initialize the points outside of the original domain and the weights and speeds of all points
    for(int x=0; x<settings::x_local_grid_size; ++x){
        for(int y=0; y<settings::y_local_grid_size; ++y){
            for(int z=0; z<settings::z_local_grid_size; ++z){
                //set every weight to infinity and status to far
                subdomain.weight_array[local_arr_index(x,y,z)] = std::numeric_limits<double>::infinity();
                subdomain.status_array[local_arr_index(x,y,z)]='0';
                //if a point is outside of the original domain, it gets the respective status and
                if(x==0 && at_x_lower_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(y==0 && at_y_lower_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(z==0 && at_z_lower_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(x==settings::x_local_grid_size-1 && at_x_upper_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(y==settings::y_local_grid_size-1 && at_y_upper_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(z==settings::z_local_grid_size-1 && at_z_upper_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }

                //the rest of the points have a well defined value in the speed and mask function
                subdomain.speed_array[local_arr_index(x,y,z)]= speed_funct(x-(subdomain.x_offset+1), y-(subdomain.y_offset+1), z-(subdomain.z_offset+1));
                if(in_mask(x+subdomain.x_offset-1, y+subdomain.y_offset-1, z+subdomain.z_offset-1)){
                    subdomain.weight_array[local_arr_index(x,y,z)] = 0;
                    subdomain.status_array[local_arr_index(x,y,z)]='3';
                }
            }
        }
    }
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
    if(x<settings::x_local_grid_size-1)
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
    if(y<settings::y_local_grid_size-1)
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
    if(z<settings::z_local_grid_size-1)
    {
        res.push_back(x);
        res.push_back(y);
        res.push_back(z+1);
    }
    return res;
}
double solve_eikonal_quadratic_3d(const int x, const int y, const int z){
    //TODO
    return 1;
}

void update_neighbors(int subdomain_index, SubdomainData &subdomain, const int x, const int y, const int z){

    std::vector<int> neighbors = get_neighbors(x, y, z);
    int curr_node_index = local_arr_index(x, y, z);
    for (std::size_t i = 0; i < neighbors.size(); i += 3){
        int neighbor_index = local_arr_index(neighbors[i], neighbors[i+1], neighbors[i+2]);
        //WeightedPoint neighbor{neighbors[i], neighbors[i+1], neighbors[i+2], subdomain.weight_array[neighbor_index]};

        if(subdomain.status_array[curr_node_index]!='3' && subdomain.weight_array[curr_node_index]< subdomain.weight_array[neighbor_index]){
            double temp = solve_eikonal_quadratic_3d(neighbors[i], neighbors[i+1], neighbors[i+2]);
            //if the solver calculates a smaller value we update the point to BAND_NEW
            if(temp< subdomain.weight_array[neighbor_index]){
                subdomain.weight_array[neighbor_index] = temp;
                subdomain.status_array[neighbor_index] = '1';
                if(subdomain.h.get_heap_index(neighbors[i], neighbors[i+1], neighbors[i+2]) ==-1){
                    subdomain.h.insertKey(WeightedPoint{neighbors[i], neighbors[i+1], neighbors[i+2], temp});
                }
                else{
                    subdomain.h.decreaseKey(neighbor_index, temp);
                }
            }
        }
    }

}

void initialize_heap(int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes]){

    SubdomainData subdomain = subdomain_array[subdomain_index];
    //wollte hier multi-threaden. Das kreirt aber zu viele THreads und das Programm stürzt ab
    /*std::vector<std::thread> thread_vector;
    for(int x=0 ; x<settings::x_local_grid_size; ++x) {
        for (int y = 0; y < settings::y_local_grid_size; ++y) {
            for (int z = 0; z < settings::z_local_grid_size; ++z) {
                if (subdomain.status_array[local_arr_index(x, y, z)] == '3') {
                    std::thread th(update_neighbors, subdomain_index, x, y, z);
                    thread_vector.push_back(std::move(th));
                }
            }
        }
    }
    for(std::thread & th : thread_vector){
            th.join();
    }*/

    for(int x=0 ; x<settings::x_local_grid_size; ++x) {
        for (int y = 0; y < settings::y_local_grid_size; ++y) {
            for (int z = 0; z < settings::z_local_grid_size; ++z) {
                if (subdomain.status_array[local_arr_index(x, y, z)] == '3') {
                    update_neighbors(subdomain_index, subdomain, x, y, z);
                }
            }
        }
    }

}



int main() {
    SubdomainData *subdomain_array{new SubdomainData[settings::total_num_processes]{}};
    //total size needed for heap for the lookup table
    MinHeap main_heap(settings::total_num_processes);
    //setting domain offsets:
    for(int x=0; x<settings::x_num_processes; ++x){
        for(int y=0; y<settings::y_num_processes; ++y){
            for(int z=0; z<settings::z_num_processes; ++z){
                subdomain_array[process_index(x,y,z)].x_offset=x*(settings::x_local_grid_size -2);
                subdomain_array[process_index(x,y,z)].y_offset=y*(settings::y_local_grid_size -2);
                subdomain_array[process_index(x,y,z)].z_offset=z*(settings::z_local_grid_size -2);
            }
        }
    }



    std::cout<<std::thread::hardware_concurrency()<<"\n";
    std::thread *thread_array{new std::thread[settings::total_num_processes]{}};

    //INITIALIZE_INTERFACE


    for(int i=0; i<settings::total_num_processes; ++i){
        std::thread th(initialize_subdomain, i, subdomain_array);
        thread_array[i]=std::move(th);
        //th.join();
    }
    for(int i=0; i<settings::total_num_processes; ++i) {
        if(thread_array[i].joinable()) {
            //thread_array[i].detach();
            thread_array[i].join();
        }
    }
    //INITIALIZE HEAP
    for(int i=0; i<settings::total_num_processes; ++i){
        std::cout<<"1\n";
        std::thread th(initialize_heap, i, subdomain_array);
        std::cout<<"2\n";
        //thread_array[i].join();
        std::cout<<"3\n";
        thread_array[i]=std::move(th);
        std::cout<<"4\n";
    }
    for(int i=0; i<settings::total_num_processes; ++i){
        if(thread_array[i].joinable()) {
            //thread_array[i].join();
        }
        else{
            std::cout<<"AASADSD\n";
        }
    }

    delete[] thread_array;
    for(int i=0; i<settings::total_num_processes; ++i){
        for(int j=0; j<settings::total_local_grid_size; ++j){
            std::cout << "STATUS: " << subdomain_array[i].status_array[j] << std::endl;
            std::cout << "WEIGHT: " << subdomain_array[i].weight_array[j] << std::endl;
            std::cout << "SPEED: " << subdomain_array[i].speed_array[j] << "\n" << std::endl;
        }
    }





    return 0;

}

/*
 * std::mutex coutMutex;
 *thread_local std::string s("hello from ");
 * void addThreadLocal(std::string const& s2, std::array<std::string, settings::total_num_processes> &shared_data, int index){

    s+=s2;
    // protect std::cout
    std::lock_guard<std::mutex> guard(coutMutex);
    std::cout << s << " " << shared_data[index] << std::endl;
    std::cout << "&s: " << &s << std::endl;
    shared_data[index]=s;

}
int main() {
    std::array<std::string,settings::total_num_processes> shared_data;
    for( int i=0; i<settings::total_num_processes; ++i )
    {
       shared_data[i]=std::to_string(i);
    }
    std::array<std::thread, settings::total_num_processes> thread_array;

    for( int i=0; i<settings::total_num_processes; ++i )
    {
        thread_array[i] = std::thread(addThreadLocal,std::to_string(i), std::ref(shared_data), i);
    }
    std::cout<<"sjfkgbahsf"<<std::endl;
    for( int i=0; i<settings::total_num_processes; ++i )
    {
        thread_array[i].join();
    }

    for( int i=0; i<settings::total_num_processes; ++i )
    {
        std::cout<< shared_data[i]<<std::endl;
    }
    return 0;
}

*/