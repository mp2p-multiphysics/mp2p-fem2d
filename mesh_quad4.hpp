#ifndef MESH_QUAD4
#define MESH_QUAD4
#include <unordered_map>
#include "container_typedef.hpp"

struct MeshQuad4Struct
{
    /*

    struct with points and elements of a mesh domain made of quad4 elements.

    Variables
    =========
    num_point_domain : int
        Number of points in the domain.
    point_gid_vec : VectorInt
        vector with global IDs of points in the domain.
    point_position_x_vec : VectorDouble
        vector with x-coordinates of points in the domain.
    point_position_y_vec : VectorDouble
        vector with y-coordinates of points in the domain.
    point_gid_to_did_map : MapIntInt
        map that outputs the global ID of a point given its domain ID.
    num_element_domain : int
        Number of elements in the domain.
    element_gid_vec : VectorInt
        vector with global IDs of elements in the domain.
    element_p0_gid_vec : VectorInt
        vector with global IDs of point 0 in the domain.
    element_p1_gid_vec : VectorInt
        vector with global IDs of point 1 in the domain.
    element_p2_gid_vec : VectorInt
        vector with global IDs of point 2 in the domain.
    element_p3_gid_vec : VectorInt
        vector with global IDs of point 3 in the domain.
    point_gid_to_did_map : MapIntInt
        map that outputs the global ID of an element given its domain ID.

    Notes
    ====
    The global ID is a unique label for each point.
    The domain ID applies only to a domain and is used to iterate through the vectors in this code.
    The figure below is a quad4 element transformed into local coordinates. Points 0, 1, 2, 3 are labeled.

               (local y)
                   ^
                   |
              1 ---|--- 2
              |    |    | 
        <----------+----------> (local x)     
              |    |    |
              0 ---|--- 3
                   |
                   v

    */

    // did - domain ID
    // gid - global ID
    // vectors use did as input

    // point data
    int num_point_domain = 0;
    VectorInt point_gid_vec;
    VectorDouble point_position_x_vec;
    VectorDouble point_position_y_vec;
    MapIntInt point_gid_to_did_map;

    // element data
    int num_element_domain = 0;
    VectorInt element_gid_vec;
    VectorInt element_p0_gid_vec;
    VectorInt element_p1_gid_vec;
    VectorInt element_p2_gid_vec;
    VectorInt element_p3_gid_vec;
    MapIntInt element_gid_to_did_map;

};

#endif
