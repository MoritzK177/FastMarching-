//
// Created by Moritz on 04.06.2021.
//

#ifndef FASTMARCHINGPARALLEL1_SETTINGS_H
#define FASTMARCHINGPARALLEL1_SETTINGS_H
namespace settings{
    //The global grid sizes:
    constexpr int x_global_grid_size {128};
    constexpr int y_global_grid_size {128};
    constexpr int z_global_grid_size {128};
    constexpr int total_global_grid_size {x_global_grid_size*y_global_grid_size*z_global_grid_size};

    //The number of processes in each direction:
    constexpr int x_num_processes {1};
    constexpr int y_num_processes {1};
    constexpr int z_num_processes {1};
    constexpr int total_num_processes {x_num_processes*y_num_processes*z_num_processes};

    //The local grid sizes, assumes that each global grid size is divisible by the respective number of processes
    constexpr int x_local_grid_size {x_global_grid_size/x_num_processes +2};
    constexpr int y_local_grid_size {y_global_grid_size/y_num_processes +2};
    constexpr int z_local_grid_size {z_global_grid_size/z_num_processes +2};
    constexpr int total_local_grid_size {x_local_grid_size*y_local_grid_size*z_local_grid_size};
    constexpr double h{1.0/127};
}
#endif //FASTMARCHINGPARALLEL1_SETTINGS_H
