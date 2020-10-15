// Copyright 2018 Global Phasing Ltd.

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
    std::vector<std::vector<float>> pos_vec,
    std::vector<Transform> transform_vec,
    std::vector<std::vector<float>> com_moving_vec,
    std::vector<std::vector<float>> com_reference_vec
    )
{

    for (std::size_t i=0; i < point_vec.size(); i++)
    {
        // Position
        std::vector<int> point = point_vec[i];
        std::vector<float> pos_python = pos_vec[i];
        Position pos = Position(
            pos_python[0],
            pos_python[1],
            pos_python[2]
            );
        Transform transform = transform_vec[i];
        std::vector<float> com_moving = com_moving_vec[i];
        std::vector<float> com_reference = com_reference_vec[i];


        //Subtract moving com
        pos.x -= com_reference[0];
        pos.y -= com_reference[1];
        pos.z -= com_reference[2];

        //transform
        Position pos_moving = Position(transform.apply(pos));

        // add reference com
        pos_moving.x += com_moving[0];
        pos_moving.y += com_moving[1];
        pos_moving.z += com_moving[2];

        // fractionalise
        Fractional pos_moving_fractional = moving_map.unit_cell.fractionalize(pos_moving);

        // interpolate
        float interpolated_value = moving_map.interpolate_value(pos_moving_fractional);

        // assign
        interpolated_map.set_value(
            point[0],
            point[1],
            point[2],
            interpolated_value
            );

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

