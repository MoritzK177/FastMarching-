//
// Created by Moritz on 31.03.2021.
//

#ifndef FASTMARCHING_SETTINGS_H
#define FASTMARCHING_SETTINGS_H



namespace settings{
    constexpr int x_grid_size {100};
    constexpr int y_grid_size {100};
    constexpr int z_grid_size {100};
    constexpr int total_grid_size {x_grid_size*y_grid_size*z_grid_size};
    constexpr double h{1.0/99};
}



#endif //FASTMARCHING_SETTINGS_H
