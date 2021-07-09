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
    return 0.001* (pow(std::sin(x),2)+pow(std::cos(y),2)+0.1);
    //return 0.001;
}
bool in_mask(int x, int y , int z)
{
    //cube=[15,24]^3
    bool x_cor = x<=24 && x>= 15;
    bool y_cor = y<=24 && y>= 15;
    bool z_cor = z<=24 && z>= 15;
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
bool at_subdomain_border(const int x,const int y,const int z){
    return(x==0 || y==0 || z==0 || x==(settings::x_local_grid_size-1) || y==(settings::y_local_grid_size-1) || z==(settings::z_local_grid_size-1));
}

bool is_in_neighbor(const int node_x, const int node_y, const int node_z, const int proc_x_diff, const int proc_y_diff, const int proc_z_diff){
    //process cant be its own neighbor
    assert(proc_x_diff !=0 ||proc_y_diff !=0 ||proc_z_diff !=0);
    int x_max = settings::x_local_grid_size-1;
    int y_max = settings::y_local_grid_size-1;
    int z_max = settings::z_local_grid_size-1;

    bool res_x{true};
    if(proc_x_diff == 1){
        res_x = (node_x == 0);
    }
    else if(proc_x_diff == -1){
        res_x = (node_x == x_max);
    }

    bool res_y{true};
    if(proc_y_diff == 1){
        res_y = (node_y == 0);
    }
    else if(proc_y_diff == -1){
        res_y = (node_y == y_max);
    }

    bool res_z{true};
    if(proc_z_diff == 1){
        res_z = (node_z == 0);
    }
    else if(proc_z_diff == -1){
        res_z = (node_z == z_max);
    }

    return (res_x && res_y && res_z);
}

std::vector<int> get_node_neighbors(const int x, const int y, const int z)
{
/* returns a vector of the coordinates of the neighboring
 * procceces if c=='p'
 * nodes if c== 'n'
 * in the format<x_1,y_1,z_1,x_2,....
 * */
    std::vector<int> res={};
    int x_max{settings::x_local_grid_size-1};
    int y_max{settings::y_local_grid_size-1};
    int z_max{settings::z_local_grid_size-1};

    if(x>0)
    {
        res.push_back(x-1);
        res.push_back(y);
        res.push_back(z);
    }
    if(x<x_max)
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
    if(y<y_max)
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
    if(z<z_max)
    {
        res.push_back(x);
        res.push_back(y);
        res.push_back(z+1);
    }
    return res;
}

std::vector<int> get_process_neighbors(const int x, const int y, const int z)
{
    std::vector<int> res={};
    int x_max{settings::x_num_processes-1};
    int y_max{settings::y_num_processes-1};
    int z_max{settings::z_num_processes-1};

    for(int i = -1; i<2 ; ++i){
        for(int j = -1; j<2 ; ++j){
            for(int k = -1; k<2 ; ++k){
                //process cant be its own neighbor
                if(i !=0 || j!=0 || k!=0){
                    if( x+i >=0 && x+i <= x_max && y+j >=0 && y+j <= y_max && z+k >=0 && z+k <= z_max){
                        res.push_back(x+i);
                        res.push_back(y+j);
                        res.push_back(z+k);
                    }
                }
            }
        }
    }
    return res;
}
int get_index(const int x, const int y, const int z, const int neighbor_x, const int neighbor_y, const int neighbor_z){
    //Assumes the input is actually a neighbor
    int res{-1};
    //get neighboring proccesses:
    std::vector<int> neighbors = get_process_neighbors(x, y, z);
    for (int i = 0; i < neighbors.size(); i += 3) {
        //TODO Can maybe be parallelized more
        if(neighbor_x == neighbors[i] && neighbor_y == neighbors[i + 1] && neighbor_z == neighbors[i + 2]){
            res = i/3;
        }
    }
    //should work if actually neighbor
    assert(res >= 0);
    return res;

}
void initialize_subdomain(int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes])
{
    //some bools to check which point is outside of the global domain
    SubdomainData &subdomain = subdomain_array[subdomain_index];

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

                //if a point is outside of the original domain, it gets the respective status maybe better to give known_fix 3 instead of 6
                if(x==0 && at_x_lower_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    //subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(y==0 && at_y_lower_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    //subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(z==0 && at_z_lower_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    //subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(x==settings::x_local_grid_size-1 && at_x_upper_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    //subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(y==settings::y_local_grid_size-1 && at_y_upper_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    //subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }
                if(z==settings::z_local_grid_size-1 && at_z_upper_border){
                    subdomain.speed_array[local_arr_index(x,y,z)]=0;
                    //subdomain.status_array[local_arr_index(x,y,z)]='6';
                    continue;
                }

                //the rest of the points have a well defined value in the speed and mask function
                subdomain.speed_array[local_arr_index(x,y,z)]= speed_funct(x + subdomain.x_offset-1, y + subdomain.y_offset-1, z+subdomain.z_offset-1);
                if(in_mask(x+subdomain.x_offset-1, y+subdomain.y_offset-1, z+subdomain.z_offset-1)){
                    subdomain.weight_array[local_arr_index(x,y,z)] = 0;
                    subdomain.status_array[local_arr_index(x,y,z)]='3';
                }
            }
        }
    }
}

double solve_eikonal_quadratic_3d(SubdomainData &subdomain, const int x, const int y, const int z){
    //TODO maybe improve or  check compatability of old version
    double temp_res = std::numeric_limits<double>::infinity();
    double speed= subdomain.speed_array[local_arr_index(x,y,z)];
    assert(speed >=0);
    //If the speed is zero(outside of original domain) we should instantly return infinity
    if(speed == 0){
        return temp_res;
    }
    int node_index = local_arr_index(x,y,z);
    double min_res_array[3];
    double h_array[3];

    //Check the x direction to set ψ_1 and h_array[0]
    int d{0};
    if(x > 0){
        int curr_index = local_arr_index(x-1,y,z);
        if((subdomain.status_array[curr_index]=='3' ||subdomain.status_array[curr_index]=='4' ||subdomain.status_array[curr_index]=='5')&& subdomain.weight_array[curr_index]< subdomain.weight_array[node_index]){
            d=-1;
        }
    }
    if(x<settings::x_local_grid_size-1){
        int curr_index = local_arr_index(x+1,y,z);
        if((subdomain.status_array[curr_index]=='3' ||subdomain.status_array[curr_index]=='4' ||subdomain.status_array[curr_index]=='5')&& subdomain.weight_array[curr_index]< subdomain.weight_array[node_index]){
            if(d==0){ // ||subdomain.weight_array[curr_index]< subdomain.weight_array[local_arr_index(x-1,y,z)]){
                d=1;
            }
            else if(subdomain.weight_array[curr_index]< subdomain.weight_array[local_arr_index(x-1,y,z)]){
                d=1;
            }
        }
    }
    if(d!=0){
        min_res_array[0]=subdomain.weight_array[local_arr_index(x+d,y,z)];
        h_array[0]=pow(settings::h,-1);
    }
    else{
        min_res_array[0]=0;
        h_array[0]=0;
    }
    //Check the y direction to set ψ_2 and h_array[1]
    d=0;
    if(y > 0){
        int curr_index = local_arr_index(x,y-1,z);
        if((subdomain.status_array[curr_index]=='3' ||subdomain.status_array[curr_index]=='4' ||subdomain.status_array[curr_index]=='5')&& subdomain.weight_array[curr_index]< subdomain.weight_array[node_index]){
            d=-1;
        }
    }
    if(y<settings::y_local_grid_size-1){
        int curr_index = local_arr_index(x,y+1,z);
        if((subdomain.status_array[curr_index]=='3' ||subdomain.status_array[curr_index]=='4' ||subdomain.status_array[curr_index]=='5')&& subdomain.weight_array[curr_index]< subdomain.weight_array[node_index]){
            if(d==0){// ||subdomain.weight_array[curr_index]< subdomain.weight_array[local_arr_index(x,y-1,z)]){
                d=1;
            }
            else if(subdomain.weight_array[curr_index]< subdomain.weight_array[local_arr_index(x,y-1,z)]){
                d=1;
            }
        }
    }
    if(d!=0){
        min_res_array[1]=subdomain.weight_array[local_arr_index(x,y+d,z)];
        h_array[1]=pow(settings::h,-1);
    }
    else{
        min_res_array[1]=0;
        h_array[1]=0;
    }
    //Check the z direction to set ψ_3 and h_array[2]
    d=0;
    if(z > 0){
        int curr_index = local_arr_index(x,y,z-1);
        if((subdomain.status_array[curr_index]=='3' ||subdomain.status_array[curr_index]=='4' ||subdomain.status_array[curr_index]=='5')&& subdomain.weight_array[curr_index]< subdomain.weight_array[node_index]){
            d=-1;
        }
    }
    if(z<settings::z_local_grid_size-1){
        int curr_index = local_arr_index(x,y,z+1);
        if((subdomain.status_array[curr_index]=='3' ||subdomain.status_array[curr_index]=='4' ||subdomain.status_array[curr_index]=='5')&& subdomain.weight_array[curr_index]< subdomain.weight_array[node_index]){
            if(d==0){// ||subdomain.weight_array[curr_index]< subdomain.weight_array[local_arr_index(x,y,z-1)]){
                d=1;
            }
            else if(subdomain.weight_array[curr_index]< subdomain.weight_array[local_arr_index(x,y,z-1)]){
                d=1;
            }
        }
    }
    if(d!=0){
        min_res_array[2]=subdomain.weight_array[local_arr_index(x,y,z+d)];
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

void update_neighbors(int subdomain_index, SubdomainData &subdomain, const int x, const int y, const int z){

    std::vector<int> neighbors = get_node_neighbors(x, y, z);
    int curr_node_index = local_arr_index(x, y, z);
    for (auto i = 0; i < neighbors.size(); i += 3){
        int neighbor_index = local_arr_index(neighbors[i], neighbors[i+1], neighbors[i+2]);
        //WeightedPoint neighbor{neighbors[i], neighbors[i+1], neighbors[i+2], subdomain.weight_array[neighbor_index]};

        //we cant update fixed points (known fix, ghost points arehandled in solve_quadratic
        if(subdomain.status_array[neighbor_index]!='3'){  // && subdomain.status_array[neighbor_index]!='6' ){
            if(subdomain.weight_array[curr_node_index]< subdomain.weight_array[neighbor_index]) {

                double temp = solve_eikonal_quadratic_3d(subdomain, neighbors[i], neighbors[i + 1], neighbors[i + 2]);
                //if the solver calculates a smaller value we update the point to BAND_NEW
                if (temp < subdomain.weight_array[neighbor_index]) {
                    subdomain.weight_array[neighbor_index] = temp;
                    subdomain.status_array[neighbor_index] = '1';

                    //int check = subdomain.h.get_heap_index(neighbors[i], neighbors[i + 1], neighbors[i + 2]);

                    if (subdomain.h.get_heap_index(neighbors[i], neighbors[i + 1], neighbors[i + 2]) == -1) {
                        subdomain.h.insertKey(WeightedPoint{neighbors[i], neighbors[i + 1], neighbors[i + 2], temp});
                    }
                    else {
                        subdomain.h.decreaseKey(subdomain.h.get_heap_index(neighbors[i], neighbors[i + 1], neighbors[i + 2]), temp);
                    }
                }
            }
        }

    }

}

void initialize_heap(int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes]){

    SubdomainData &subdomain = subdomain_array[subdomain_index];
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

void collect_overlapping_data(const int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes], std::vector <std::vector<std::vector<ExchangeData>>> &exchange_vector){
    SubdomainData &subdomain = subdomain_array[subdomain_index];
    //TODO count sollte in thread data sein
    subdomain.count_new=0;
    //recover the domains' indices:
    const int x_index = subdomain.x_offset/(settings::x_local_grid_size -2);
    const int y_index = subdomain.y_offset/(settings::y_local_grid_size -2);
    const int z_index = subdomain.z_offset/(settings::z_local_grid_size -2);
    //get neighboring proccesses:
    std::vector<int> neighbors = get_process_neighbors(x_index, y_index, z_index);
    for(int x=0; x<settings::x_local_grid_size; ++x) {
        for (int y = 0; y < settings::y_local_grid_size; ++y) {
            for (int z = 0; z < settings::z_local_grid_size; ++z) {

                char curr_node_status{subdomain.status_array[local_arr_index(x,y,z)]};
                if(at_subdomain_border(x,y,z) && (curr_node_status == '1' || curr_node_status == '4')) {
                    ++subdomain.count_new;
                    if (curr_node_status == '1') {
                        subdomain.status_array[local_arr_index(x, y, z)] = '2';
                    }
                    else {
                        subdomain.status_array[local_arr_index(x, y, z)] = '5';
                    }


                    for (std::size_t i = 0; i < neighbors.size(); i += 3) {
                        //TODO Can maybe be parallelized more
                        int neighbor_index = process_index( neighbors[i], neighbors[i + 1], neighbors[i + 2]);
                        //TODO HIER GROßER FEHLER!!!
                        //assert(is_in_neighbor(x, y, z, x_index - neighbors[i], y_index - neighbors[i + 1],z_index - neighbors[i + 2]));

                        if (is_in_neighbor(x, y, z, x_index - neighbors[i], y_index - neighbors[i + 1],z_index - neighbors[i + 2])) {
                            int process_index_in_neighbor = get_index(neighbors[i], neighbors[i + 1], neighbors[i + 2], x_index, y_index, z_index);
                            //give the GLOBAL coordinates as ExchangeData
                            ExchangeData data{x+subdomain.x_offset-1, y+subdomain.y_offset-1, z+subdomain.z_offset-1, subdomain.weight_array[local_arr_index(x, y, z)]};
                            exchange_vector[neighbor_index][process_index_in_neighbor].push_back(data);
                        }
                    }

                }
            }
        }
    }
}
void integrate_overlapping_data(const int subdomain_index, double bound_band, SubdomainData subdomain_array[settings::total_num_processes], std::vector<std::vector<ExchangeData>> &exchange_vector) {
    SubdomainData &subdomain = subdomain_array[subdomain_index];
    //TODO can be further parallelized, look at suggestions
    for(int i=0; i<exchange_vector.size();++i){
        for(int j=0; j < exchange_vector[i].size();++j){
            ExchangeData temp_data = exchange_vector[i][j];
            //hier Logikfehler behoben?!
            int temp_x = temp_data.x-subdomain.x_offset+1;
            int temp_y = temp_data.y-subdomain.y_offset+1;
            int temp_z = temp_data.z-subdomain.z_offset+1;
            int temp_index = local_arr_index(temp_x, temp_y, temp_z);
            double temp_val =  temp_data.value;

            if(temp_val<subdomain.weight_array[temp_index]) {
                subdomain.weight_array[temp_index] = temp_val;

                if (subdomain.weight_array[temp_index] > bound_band) {
                    subdomain.status_array[temp_index] = '2';
                } else {
                    subdomain.status_array[temp_index] = '5';
                }
                int heap_index = subdomain.h.get_heap_index(temp_x, temp_y, temp_z);
                if(heap_index==-1){
                    subdomain.h.insertKey(WeightedPoint{temp_x, temp_y, temp_z, temp_val});
                }
                else {
                    subdomain.h.decreaseKey(heap_index, temp_val);
                }
            }

        }
    }
}
void march_narrow_band(const int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes], double bound_band) {
    SubdomainData &subdomain = subdomain_array[subdomain_index];
    while(true){
        if(subdomain.h.get_size()==0){
            break;
        }
        WeightedPoint curr_point = subdomain.h.getMin();
        int curr_index= local_arr_index(curr_point.m_x, curr_point.m_y, curr_point.m_z);
        double value = curr_point.weight;

        if(value> bound_band){
            break;
        }
        if(subdomain.status_array[curr_index]!= '5'){
            subdomain.status_array[curr_index] = '4';
        }

        WeightedPoint comp =subdomain.h.extractMin();
        subdomain.weight_array[curr_index]= value;
        update_neighbors(subdomain_index,subdomain, curr_point.m_x, curr_point.m_y, curr_point.m_z);
    }
}

/*void test(const int subdomain_index, SubdomainData subdomain_array[settings::total_num_processes]){
    SubdomainData &subdomain = subdomain_array[subdomain_index];
    const int x_index = subdomain.x_offset/(settings::x_local_grid_size -2);
    const int y_index = subdomain.y_offset/(settings::y_local_grid_size -2);
    const int z_index = subdomain.z_offset/(settings::z_local_grid_size -2);
    std::vector<int> neighbors = get_process_neighbors(x_index, y_index, z_index);
    for (std::size_t i = 0; i < neighbors.size(); i += 3) {
        std::cout << neighbors[i]<<" "<< neighbors[i + 1]<<" "<< neighbors[i + 2]<<std::endl;
        std::cout << get_index(x_index, y_index, z_index, neighbors[i], neighbors[i + 1], neighbors[i + 2]) << std::endl;
    }
    for(int x=0; x<settings::x_local_grid_size; ++x) {
        for (int y = 0; y < settings::y_local_grid_size; ++y) {
            for (int z = 0; z < settings::z_local_grid_size; ++z) {

                for (std::size_t i = 0; i < neighbors.size(); i += 3) {
                    //TODO Can maybe be parallelized more
                    if(is_in_neighbor( x, y, z,x_index - neighbors[i], y_index - neighbors[i + 1], z_index - neighbors[i + 2])){

                        //std::cout<< x<<" " << y<<" " << z<<" " <<x_index - neighbors[i]<<" " << y_index - neighbors[i + 1]<<" " << z_index - neighbors[i + 2]<<std::endl;

                    }
                }
            }
        }
    }
}*/

int main() {
    //TODO was ist das
    double width_band{std::numeric_limits<double>::infinity()};
    double stride{2.0/40};
    //double eps{1./1000};
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
    //std::thread *thread_array{new std::thread[settings::total_num_processes]{}};
    std::vector<std::thread> thread_array;
    thread_array.resize(settings::total_num_processes);
    //INITIALIZE_INTERFACE


    for(int i=0; i<settings::total_num_processes; ++i){
        std::thread th(initialize_subdomain, i, subdomain_array);
        thread_array[i]=std::move(th);
    }
    for(int i=0; i<settings::total_num_processes; ++i) {
        if(thread_array[i].joinable()) {
            thread_array[i].join();
        }
    }
    //INITIALIZE HEAP
    for(int i=0; i<settings::total_num_processes; ++i){
        std::thread th(initialize_heap, i, subdomain_array);
        thread_array[i]=std::move(th);
    }
    for(int i=0; i<settings::total_num_processes; ++i){
        if(thread_array[i].joinable()) {
            thread_array[i].join();
        }
    }

    for(int i=0; i< settings::total_num_processes;++i) {
        //collect_overlapping_data(i, subdomain_array);
    }
    //test(0,subdomain_array);
    /*std::vector<int> neighbors = get_process_neighbors(0,0,0);
    for (std::size_t i = 0; i < neighbors.size(); i += 3) {
        std::cout<<neighbors[i]<<" "<<neighbors[i+1]<<" "<<neighbors[i+2]<<std::endl;
    }*/

    while(true){
        //TODO globale reduzierung
        double min_array[settings::total_num_processes];
        int count_array[settings::total_num_processes];
        //std::array <double,settings::total_num_processes> min_array;
        //std::array <int, settings::total_num_processes> count_array;
        for(int i=0; i<settings::total_num_processes; ++i){
            count_array[i]= subdomain_array[i].count_new;
            if(subdomain_array[i].h.get_size() > 0) {
                //std::cout<<"Process: "<< i <<" minimum: "<< subdomain_array[i].h.getMin().weight<<std::endl;
                min_array[i] = subdomain_array[i].h.getMin().weight;
            }
            else min_array[i] = width_band;
        }
        double *min_val_global = std::min_element(std::begin(min_array), std::end(min_array));
        //std::cout << "MIN VAL GLOBAL: "<< *min_val_global<<std::endl;
        int *count_global = std::max_element(std::begin(count_array), std::end(count_array));
        //std::cout << "Count GLOBAL: "<< *count_global<<std::endl;
        if((*min_val_global >= width_band)  && (*count_global ==0)){
            break;
        }
        double bound_band = std::min(*min_val_global +stride, width_band);
        //std::cout<<"min_val +stride: "<< *min_val_global +stride<< std::endl;
        //std::cout<<"bound_band: "<< bound_band<< std::endl;

        //march_band
        for(int i=0; i<settings::total_num_processes; ++i){
            std::thread th(march_narrow_band,i, subdomain_array, bound_band);
            thread_array[i]=std::move(th);
        }
        for(int i=0; i<settings::total_num_processes; ++i){
            if(thread_array[i].joinable()) {
                thread_array[i].join();
            }
        }
        //TODO besser vor schleife ziehen und dann clearen?!!
        std::vector <std::vector<std::vector<ExchangeData>>> exchange_vector;
        exchange_vector.resize(settings::total_num_processes);


        for(int i=0; i< settings::total_num_processes; ++i){
            exchange_vector[i].resize(26);
        }
        //exchange data
        for(int i=0; i<settings::total_num_processes; ++i){
            std::thread th(collect_overlapping_data, i, subdomain_array, std::ref(exchange_vector));
            thread_array[i]=std::move(th);
        }
        for(int i=0; i<settings::total_num_processes; ++i){
            if(thread_array[i].joinable()) {
                thread_array[i].join();
            }
        }
        //integrate data
        for(int i=0; i<settings::total_num_processes; ++i){
            std::thread th(integrate_overlapping_data, i,bound_band, subdomain_array, std::ref(exchange_vector[i]));
            thread_array[i]=std::move(th);
        }
        for(int i=0; i<settings::total_num_processes; ++i){
            if(thread_array[i].joinable()) {
                thread_array[i].join();
            }
        }

        //march_band
        for(int i=0; i<settings::total_num_processes; ++i){
            std::thread th(march_narrow_band,i, subdomain_array, bound_band);
            thread_array[i]=std::move(th);
        }
        for(int i=0; i<settings::total_num_processes; ++i){
            if(thread_array[i].joinable()) {
                thread_array[i].join();
            }
        }
    }
    std::cout<<"pls"<<std::endl;
    //delete[] thread_array;



    double weight_array[settings::total_global_grid_size];
    for(int x =0; x<settings::x_global_grid_size;++x){
        for(int y =0; y<settings::y_global_grid_size;++y){
            for(int z =0; z<settings::z_global_grid_size;++z){

                int process_xid = x/(settings::x_local_grid_size-2);
                int process_yid = y/(settings::y_local_grid_size-2);
                int process_zid = z/(settings::z_local_grid_size-2);
                int x_local_coord = (x % (settings::x_local_grid_size-2))+1;
                int y_local_coord = (y % (settings::y_local_grid_size-2))+1;
                int z_local_coord = (z % (settings::z_local_grid_size-2))+1;
                /*std::cout<<"X-id: "<< process_xid<<" Y-id: "<<process_yid<<" Z-id: "<<process_zid<<std::endl;
                std::cout<<"X: "<< x_local_coord<<" Y: "<<y_local_coord<<" Z: "<<z_local_coord<<std::endl;
                std::cout<<"RESULT: "<<subdomain_array[process_index(process_xid, process_yid,process_zid)].weight_array[x_local_coord, y_local_coord, z_local_coord]<<std::endl;*/
                weight_array[global_arr_index(x,y,z)]= subdomain_array[process_index(process_xid, process_yid,process_zid)].weight_array[local_arr_index(x_local_coord, y_local_coord, z_local_coord)];

            }

        }
    }

    int count{0};
    for(int i=0; i<settings::total_num_processes; ++i){
        for(int x=1; x<settings::x_local_grid_size-1; ++x){
            for(int y=1; y<settings::y_local_grid_size-1;++y){
                for(int z=1;z<settings::z_local_grid_size-1;++z){
                    int j = local_arr_index(x,y,z);
                    if(subdomain_array[i].weight_array[j]<10000) {
                        ++count;
                        //std::cout << "STATUS: " << subdomain_array[i].status_array[j] << std::endl;
                        //std::cout << "WEIGHT: " << subdomain_array[i].weight_array[j] << std::endl;
                        //std::cout << "SPEED: " << subdomain_array[i].speed_array[j] << "\n" << std::endl;
                    }
                }
            }
        }
    }
    std::cout<<count<< std::endl;

    //for(int i=0; i < settings::total_global_grid_size;++i){
    //    std::cout << "WEIGHT: " << weight_array[i] << std::endl;
    //}
    std::ofstream myfile;
    myfile.open("garbage.txt");
    myfile << "Dimension information\n"<<settings::x_global_grid_size <<"\n"<<settings::y_global_grid_size<<"\n"<<settings::z_global_grid_size<<"\n";
    myfile << "Mask information\n";
    for(int x=0; x<settings::x_global_grid_size; ++x){
        for(int y=0; y<settings::y_global_grid_size; ++y){
            for(int z=0; z<settings::z_global_grid_size; ++z){
                myfile << in_mask(x,y,z)<<"\n";
            }
        }
    }

    myfile<<"Result information\n";
    for ( int i =0; i< settings::total_global_grid_size;++i)
    {
        myfile << weight_array[i]<<"\n";
    }
    myfile.close();

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