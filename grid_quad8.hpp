#ifndef STRUCT_QUAD8
#define STRUCT_QUAD8
#include <vector>

struct GridQuad8Struct
{
    /**
     *  Stores positions of points in 2D grid made of Quad8 elements.
     *  
     *  Parameters
     *  ----------
     *  num_point : int
     *      Number of grid points.
     *  point_id_vec : vector<int>
     *      vector with the IDs of each grid point.
     *      The index [n] contains data for point n.
     *  point_pos_x_vec : vector<double>
     *      vector with the x-coordinates of each grid point.
     *      The index [n] contains data for point n.
     *  point_pos_y_vec : vector<double>
     *      vector with the y-coordinates of each grid point.
     *      The index [n] contains data for point n.
     *  num_element : int
     *      Number of grid elements.
     *  element_id_vec : vector<int>
     *      vector with the IDs of each grid element.
     *      The index [m] contains data for element m.
     *  element_p0_id_vec : vector<int>
     *      vector with the IDs of local point 0 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p1_id_vec : vector<int>
     *      vector with the IDs of local point 1 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p2_id_vec : vector<int>
     *      vector with the IDs of local point 2 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p3_id_vec : vector<int>
     *      vector with the IDs of local point 3 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p4_id_vec : vector<int>
     *      vector with the IDs of local point 4 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p5_id_vec : vector<int>
     *      vector with the IDs of local point 5 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p6_id_vec : vector<int>
     *      vector with the IDs of local point 6 of each grid element.
     *      The index [m] contains data for element m.
     *  element_p7_id_vec : vector<int>
     *      vector with the IDs of local point 7 of each grid element.
     *      The index [m] contains data for element m.
     * 
     *  Notes
     *  -----
     *  The local points 0-7 are positioned in the following local coordinates of each element:
     *      Point 0: (-1, -1)
     *      Point 1: (-1,  0)
     *      Point 2: (-1, +1)
     *      Point 3: ( 0, +1) 
     *      Point 4: (+1, +1)
     *      Point 5: (+1,  0)
     *      Point 6: (+1, -1)
     *      Point 7: ( 0, -1) 
     * 
    */

    // point data
    int num_point = 0;
    std::vector<int> point_id_vec;
    std::vector<double> point_pos_x_vec;
    std::vector<double> point_pos_y_vec;

    // element data
    int num_element = 0;
    std::vector<int> element_id_vec;
    std::vector<int> element_p0_id_vec;
    std::vector<int> element_p1_id_vec;
    std::vector<int> element_p2_id_vec;
    std::vector<int> element_p3_id_vec;
    std::vector<int> element_p4_id_vec;
    std::vector<int> element_p5_id_vec;
    std::vector<int> element_p6_id_vec;
    std::vector<int> element_p7_id_vec;

};

struct BoundaryConfigQuad8Struct
{
    /**
     *  Stores boundary condition configuration settings.
     *  
     *  Parameters
     *  ----------
     *  boundary_type_str : string
     *      Name of the boundary condition.
     *  boundary_parameter_vec : vector<double>
     *      vector with parameters for the boundary condition.
     * 
    */
    
    std::string boundary_type_str;
    std::vector<double> boundary_parameter_vec;

};

struct BoundaryQuad8Struct
{
    /**
     *  Stores boundary condition applied on a 2D grid made of Quad8 elements.
     *  
     *  Parameters
     *  ----------
     *  num_element_flux : int
     *      Number of element edges with a flux-type boundary condition.
     *  element_flux_id_vec : vector<int>
     *      vector with IDs of elements with a flux-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_flux_pa_loc_vec : vector<int>
     *      vector with local point As of edges with a flux-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_flux_pb_loc_vec : vector<int>
     *      vector with local point Bs of edges with a flux-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_flux_pc_loc_vec : vector<int>
     *      vector with local point Cs of edges with a flux-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_flux_config_id_vec : vector<int>
     *      vector with IDs of flux-type boundary conditions.
     *      The index [k] contains data for edge k.
     *  num_element_value : int
     *      Number of element edges with value-type boundary conditions.
     *  element_value_id_vec : vector<int>
     *      vector with IDs of elements with a value-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_value_pa_loc_vec : vector<int>
     *      vector with local point As of edges with a value-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_value_pb_loc_vec : vector<int>
     *      vector with local point Bs of edges with a value-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_value_pc_loc_vec : vector<int>
     *      vector with local point Cs of edges with a value-type boundary condition.
     *      The index [k] contains data for edge k.
     *  element_value_config_id_vec : vector<int>
     *      vector with IDs of value-type boundary conditions.
     *      The index [k] contains data for edge k.
     *  boundary_config_vec : vector<BoundaryConfigQuad4Struct>
     *      vector with data on each boundary condition.
     *      The index [l] contains data for the boundary condition of ID l.
     * 
     *  Notes
     *  -----
     *  A, B, and C are the local point numbers along the edge with an applied boundary condition.
     *      For example, if the edge is along the local coordinates (1, 1), (1, 0), and (1, -1), then A = 4, B = 5, and C = 6.
     *
    */
    
    // flux boundary condition data   
    int num_element_flux = 0;
    std::vector<int> element_flux_id_vec;
    std::vector<int> element_flux_pa_loc_vec;
    std::vector<int> element_flux_pb_loc_vec;
    std::vector<int> element_flux_pc_loc_vec;
    std::vector<int> element_flux_config_id_vec;
    
    // value boundary condition data
    int num_element_value = 0;
    std::vector<int> element_value_id_vec;
    std::vector<int> element_value_pa_loc_vec;
    std::vector<int> element_value_pb_loc_vec;
    std::vector<int> element_value_pc_loc_vec;
    std::vector<int> element_value_config_id_vec;

    // boundary condition data
    std::vector<BoundaryConfigQuad8Struct> boundary_config_vec;

};

#endif
