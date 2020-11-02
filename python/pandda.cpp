// Copyright 2018 Global Phasing Ltd.

#include <iostream>

#include "gemmi/ccp4.hpp"
#include "gemmi/gz.hpp"  // for MaybeGzipped
#include "gemmi/neighbor.hpp"
#include "gemmi/tostr.hpp"
#include "gemmi/fourier.hpp"  // for get_f_phi_on_grid

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include "common.h"  // for normalize_index

#include "ThreadPool.h"

namespace py = pybind11;
using namespace gemmi;





Grid<float> interpolate_points(
    Grid<float> moving_map,
    Grid<float> interpolated_map, 
    std::vector<std::vector<int>> point_vec,
    std::vector<std::vector<double>> pos_vec,
    std::vector<Transform> transform_vec,
    std::vector<std::vector<double>> com_moving_vec,
    std::vector<std::vector<double>> com_reference_vec,
    bool debug
    )
{
    for (std::size_t i=0; i < point_vec.size(); i++)
    {
        // Position
        // std::cout << "Getting postion\n";
        std::vector<int> point = point_vec[i];
        /*
        std::vector<float> pos_python = pos_vec[i];
        Position pos = Position(
            pos_python[0],
            pos_python[1],
            pos_python[2]
            );
        */

        // std::cout << "Getting gemmi point\n";
        //auto point_gemmi = interpolated_map.get_point(point[0], point[1], point[2]);
        // std::cout << "Getting gemmi position\n";
        //Position pos = interpolated_map.point_to_position(point_gemmi);
        Position pos = pos_vec[i];
        // std::cout << "Getting transform\n"; 
        Transform transform = transform_vec[i];
        // std::cout << "Getting com moving\n"; 
        std::vector<double> com_moving = com_moving_vec[i];
        // std::cout << "Getting reference\n"; 
        std::vector<double> com_reference = com_reference_vec[i];

        //Subtract reference com
        // std::cout << "Subtracting reference\n"; 
        if (debug) {
            std::cout << "###################" << "\n";
            std::cout << "Point: " << point[0] << " " << point[1] << " " << point[2] << "\n";
            std::vector<double> apparent_pos = pos_vec[i];
            std::cout << "Apparent pos: " << apparent_pos[0] << " " << apparent_pos[1] << " " << apparent_pos[2] << "\n";
            std::cout << "Before subtracting: " << pos.x << " " << pos.y << " " << pos.z << "\n";
            std::cout << "com reference: " << com_reference[0] << " " << com_reference[1] << " " << com_reference[2] << "\n";
        };

        pos.x -= com_reference[0];
        pos.y -= com_reference[1];
        pos.z -= com_reference[2];

        if (debug) {
            std::cout << "After subtracting: " << pos.x << " " << pos.y << " " << pos.z << "\n";
        };

        //transform
        // std::cout << "Transforming\n"; 
        Position pos_moving = Position(transform.apply(pos));

        if (debug) {
            std::cout << "After transforming: " << pos_moving.x << " " << pos_moving.y << " " << pos_moving.z << "\n";
        };

        // add moving com
        // std::cout << "Adding moving\n"; 
        pos_moving.x += com_moving[0];
        pos_moving.y += com_moving[1];
        pos_moving.z += com_moving[2];

        if (debug) {
            std::cout << "com moving: " << com_moving[0] << " " << com_moving[1] << " " << com_moving[2] << "\n";
            std::cout << "After adding: " << pos_moving.x << " " << pos_moving.y << " " << pos_moving.z << "\n";
        };

        // fractionalise
        Fractional pos_moving_fractional = moving_map.unit_cell.fractionalize(pos_moving);

        // interpolate
        //std::cout << "interpolating: " << pos_moving_fractional.x << " " << pos_moving_fractional.y << " " << pos_moving_fractional.z << "\n"; 
        Fractional wrapped = pos_moving_fractional.wrap_to_unit();
        //std::cout << "wrapped..." << wrapped.x << " " << wrapped.y << " " << wrapped.z << "\n";
        //std::cout << "interpolating..."; 
        float interpolated_value = moving_map.interpolate_value(pos_moving_fractional);
        //std::cout << "interpolated..."; 


        // assign
        // std::cout << "Assigning\n"; 
        interpolated_map.set_value(
            point[0],
            point[1],
            point[2],
            interpolated_value
            );

        // if (debug) {
        //     std::cout << "Getting gemmi point" << point_gemmi.u << " " << point_gemmi.v << " " << point_gemmi.w << "\n";
        //     std::cout << "Getting pos" << pos.x << " " << pos.y << " " << pos.z << "\n";
        //     std::cout << "com reference" << com_reference[0] << " " << com_reference[1] << " " << com_reference[2] << "\n";

        //     std::cout << "Getting pos moving" << pos_moving.x << " " << pos_moving.y << " " << pos_moving.z << "\n";
        //     std::cout << "com moving" << com_moving[0] << " " << com_moving[1] << " " << com_moving[2] << "\n";

        //     std::cout << "Getting pos_moving_fractional" << pos_moving_fractional.x << " " << pos_moving_fractional.y << " " << pos_moving_fractional.z << "\n";
        //     std::cout << "Getting wrapped" << wrapped.x << " " << wrapped.y << " " << wrapped.z << "\n";
        //     Fractional reference_fractional = interpolated_map.unit_cell.fractionalize(pos);
        //     std::cout << "Getting native fractional" << reference_fractional.x << " " << reference_fractional.y << " " << reference_fractional.z << "\n";
        // };

    };

    return interpolated_map;


}

void add_pandda(py::module& m) {
    m.def(
        "interpolate_points", 
        &interpolate_points, 
        "Interpolates a list of points."
    );

}

