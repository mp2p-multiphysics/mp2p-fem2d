#ifndef BOUNDARY_UNIT
#define BOUNDARY_UNIT
#include <fstream>
#include <functional>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "container_typedef.hpp"
#include "domain_unit.hpp"
#include "variable_unit.hpp"

namespace FEM2D
{

class BoundaryUnit
{
    /*

    Boundary conditions (BC) for quad4 mesh elements.

    Variables
    =========
    domain_in : DomainQuad4
        Domain where boundary conditions are applied.
    file_in_str_in : string
        Path to CSV file with data for BCs.

    Functions
    =========
    set_boundary : void
        Assigns a BC type ID and parameters to a BC configuration ID.
    set_boundary_type : void
        Classifies BCs as essential or natural.
    update_parameter : void
        Recalculates non-constant BC parameters.

    Notes
    ====
    The input CSV file must have the following columns:
        global element ID where BC is applied
        local point ID where BC is applied (0 to 3)
        BC configuration ID

    */

    public:

    // bid - boundary ID - used to index BC data
    // bnid - boundary natural ID - used to index natural BC data
    // beid - boundary essential ID - used to index essential BC data
    // btid - boundary type ID - denotes type of BC
    // bcid - boundary config ID - denotes location of BC
    // blid - boundary local ID - local point along boundary
    // vectors use did as input

    // mesh where variable is applied
    DomainUnit* domain_ptr;

    // type of element
    // 0 - tri3; 1 - quad4; 2 - tri6; 3 - quad8
    int num_neighbor = 0;

    // boundary condition data
    // index with bid
    int num_boundary = 0;
    VectorInt boundary_egid_vec;
    VectorInt boundary_bcid_vec;
    VectorInt boundary_btid_vec;
    VectorInt2D boundary_blid_to_plid_vec;  // [bid][blid] -> plid
    VectorDouble3D boundary_blid_to_parameter_vec; // [bid][blid] -> parameter_vec

    // use for non-constant boundary conditions
    // index with bcid
    std::unordered_map<int, bool> boundary_is_parameter_constant_map;
    std::unordered_map<int, std::function<VectorDouble(double, double, VectorDouble)>> boundary_parameter_function_map;
    std::unordered_map<int, std::vector<VariableUnit*>> boundary_variable_ptr_map;

    // essential boundary condition data
    // index with beid
    int num_essential = 0;
    VectorInt essential_egid_vec;
    VectorInt essential_bcid_vec;
    VectorInt essential_btid_vec;
    VectorInt2D essential_blid_to_plid_vec;  // [beid][blid] -> plid
    VectorDouble3D essential_blid_to_parameter_vec; // [beid][blid] -> parameter_vec

    // natural boundary condition data
    // index with bnid
    int num_natural = 0;
    VectorInt natural_egid_vec;
    VectorInt natural_bcid_vec;
    VectorInt natural_btid_vec;
    VectorInt2D natural_blid_to_plid_vec;  // [bnid][blid] -> plid
    VectorDouble3D natural_blid_to_parameter_vec; // [bnid][blid] -> parameter_vec

    // file names
    std::string file_in_str;

    // functions
    void set_boundary(int boundaryconfig_id, int boundarytype_id, VectorDouble parameter_vec);
    void set_boundary(int boundaryconfig_id, int boundarytype_id, std::function<VectorDouble(double, double, VectorDouble)> parameter_function, std::vector<VariableUnit*> variable_ptr_vec);
    void set_boundary_type(VectorInt boundarytype_essential_vec, VectorInt boundarytype_natural_vec);
    void update_parameter();

    // default constructor
    BoundaryUnit() {}

    // constructor
    BoundaryUnit(DomainUnit &domain_in, std::string file_in_str_in)
    {

        // store variables
        domain_ptr = &domain_in;
        file_in_str = file_in_str_in;

        // get number of neighboring points per boundary
        switch (domain_ptr->type_element)
        {
            case 0: num_neighbor = 2; break;  // tri3
            case 1: num_neighbor = 2; break;  // quad4
            case 2: num_neighbor = 3; break;  // tri6
            case 3: num_neighbor = 3; break;  // quad8
        }

        // read input files and store boundary condition data
        read_boundary(file_in_str);

    }

    private:
    void read_boundary(std::string file_in_str);

};

void BoundaryUnit::set_boundary(int boundaryconfig_id, int boundarytype_id, VectorDouble parameter_vec)
{
    /*

    Assigns a BC type ID and parameters to a BC configuration ID.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    boundarytype_id : int
        BC type ID assigned to BC configuration ID.
    parameter_vec : VectorDouble
        Parameters assigned to BC configuration ID.

    Returns
    =======
    (none)

    Notes
    ====
    The BC configuration ID indicates the location of the BC.
    The BC type ID denotes the type of BC (e.g., Dirichlet, Neumann, etc.).
        The int corresponding to each BC type depends on the physics.

    */

    // iterate through each boundary condition
    for (int bid = 0; bid < num_boundary; bid++)
    {

        // skip if different boundary config ID
        if (boundary_bcid_vec[bid] != boundaryconfig_id)
        {
            continue;
        }

        // set boundary type and parameters
        boundary_btid_vec[bid] = boundarytype_id;
        boundary_blid_to_parameter_vec[bid] = VectorDouble2D(num_neighbor, parameter_vec);

    }

    // set boundary parameters as constant
    boundary_is_parameter_constant_map[boundaryconfig_id] = true;

}

void BoundaryUnit::set_boundary(int boundaryconfig_id, int boundarytype_id, std::function<VectorDouble(double, double, VectorDouble)> parameter_function, std::vector<VariableUnit*> variable_ptr_vec)
{
    /*

    Assigns a BC type ID and parameters to a BC configuration ID.

    Arguments
    =========
    boundaryconfig_id : int
        BC configuration ID.
    boundarytype_id : int
        BC type ID assigned to BC configuration ID.
    parameter_function : function<double, double, VectorDouble> -> VectorDouble
        Function used to compute parameter values based on variable values.
    variable_ptr_vec : vector<VariableQuad4*>
        vector of pointers to variable objects needed to compute parameter values.

    Returns
    =======
    (none)

    Notes
    ====
    The BC configuration ID indicates the location of the BC.
    The BC type ID denotes the type of BC (e.g., Dirichlet, Neumann, etc.).
        The int corresponding to each BC type depends on the physics.
    The inputs to the parameter function are the x-coordinate (double), y-coordinate (double), and vector of variable values (VectorDouble) at a specified point.
        The variable values are in the same order as the variables in variable_ptr_vec.
    The output of the parameter function is the vector of parameter values (VectorDouble).

    */

    // iterate through each boundary condition
    for (int bid = 0; bid < num_boundary; bid++)
    {

        // skip if different boundary config ID
        if (boundary_bcid_vec[bid] != boundaryconfig_id)
        {
            continue;
        }

        // set boundary type and parameters
        boundary_btid_vec[bid] = boundarytype_id;

    }

    // set boundary parameters as non-constant
    boundary_is_parameter_constant_map[boundaryconfig_id] = false;
    boundary_parameter_function_map[boundaryconfig_id] = parameter_function;
    boundary_variable_ptr_map[boundaryconfig_id] = variable_ptr_vec;

}

void BoundaryUnit::set_boundary_type(VectorInt boundarytype_essential_vec, VectorInt boundarytype_natural_vec)
{
    /*

    Classifies BCs as essential or natural.

    Arguments
    =========
    boundarytype_essential_vec : VectorInt
        vector with BC type IDs that denote essential BCs.
    boundarytype_essential_vec : VectorInt
        vector with BC type IDs that denote natural BCs.

    Returns
    =======
    (none)

    */

    // iterate through boundaries
    for (int bid = 0; bid < num_boundary; bid++)
    {

        // get boundary type
        int btid = boundary_btid_vec[bid];

        // fill up vectors for essential boundary conditions
        auto iter_essential = std::find(boundarytype_essential_vec.begin(), boundarytype_essential_vec.end(), btid);
        if (iter_essential != boundarytype_essential_vec.end())
        {
            num_essential++;
            essential_egid_vec.push_back(boundary_egid_vec[bid]);
            essential_bcid_vec.push_back(boundary_bcid_vec[bid]);
            essential_btid_vec.push_back(boundary_btid_vec[bid]);
            essential_blid_to_plid_vec.push_back(boundary_blid_to_plid_vec[bid]);
            essential_blid_to_parameter_vec.push_back(boundary_blid_to_parameter_vec[bid]);
        }

        // fill up vectors for natural boundary conditions
        auto iter_natural = std::find(boundarytype_natural_vec.begin(), boundarytype_natural_vec.end(), btid);
        if (iter_natural != boundarytype_natural_vec.end())
        {
            num_natural++;
            natural_egid_vec.push_back(boundary_egid_vec[bid]);
            natural_bcid_vec.push_back(boundary_bcid_vec[bid]);
            natural_btid_vec.push_back(boundary_btid_vec[bid]);
            natural_blid_to_plid_vec.push_back(boundary_blid_to_plid_vec[bid]);
            natural_blid_to_parameter_vec.push_back(boundary_blid_to_parameter_vec[bid]);
        }

    }

}

void BoundaryUnit::update_parameter()
{
    /*

    Recalculates non-constant BC parameters.

    Arguments
    =========
    (none)
    
    Returns
    =========
    (none)

    */

    // iterate through essential boundaries
    for (int beid = 0; beid < num_essential; beid++)
    {

        // get boundary type
        int bcid = essential_bcid_vec[beid];

        // skip if parameters are constant
        if (boundary_is_parameter_constant_map[bcid])
        {
            continue;
        }

        // get element and surrounding points
        int egid = essential_egid_vec[beid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];
        VectorInt pgid_vec = domain_ptr->element_edid_plid_to_pgid_vec[edid];

        // iterate through each boundary point
        for (int blid = 0; blid < num_neighbor; blid++)
        {

            // get boundary point
            int plid = essential_blid_to_plid_vec[beid][blid];
            int pgid = pgid_vec[blid];
            int pdid = domain_ptr->point_pgid_to_pdid_map[pgid];

            // subset location on mesh
            double position_x = domain_ptr->point_position_x_vec[pdid];
            double position_y = domain_ptr->point_position_y_vec[pdid];

            // subset value of variables
            VectorDouble value_vec;
            for (auto variable_ptr : boundary_variable_ptr_map[bcid])
            {
                double value_sub = variable_ptr->point_value_vec[pdid];
                value_vec.push_back(value_sub);
            }

            // calculate parameter value
            essential_blid_to_parameter_vec[beid][blid] = boundary_parameter_function_map[bcid](position_x, position_y, value_vec);

        }

    }

    // iterate through natural boundaries
    for (int bnid = 0; bnid < num_natural; bnid++)
    {

        // get boundary type
        int bcid = natural_bcid_vec[bnid];

        // skip if parameters are constant
        if (boundary_is_parameter_constant_map[bcid])
        {
            continue;
        }

        // get element and surrounding points
        int egid = natural_egid_vec[bnid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];
        VectorInt pgid_vec = domain_ptr->element_edid_plid_to_pgid_vec[edid];

        // iterate through each boundary point
        for (int blid = 0; blid < num_neighbor; blid++)
        {

            // get boundary point
            int plid = natural_blid_to_plid_vec[bnid][blid];
            int pgid = pgid_vec[blid];
            int pdid = domain_ptr->point_pgid_to_pdid_map[pgid];

            // subset location on mesh
            double position_x = domain_ptr->point_position_x_vec[pdid];
            double position_y = domain_ptr->point_position_y_vec[pdid];

            // subset value of variables
            VectorDouble value_vec;
            for (auto variable_ptr : boundary_variable_ptr_map[bcid])
            {
                double value_sub = variable_ptr->point_value_vec[pdid];
                value_vec.push_back(value_sub);
            }

            // calculate parameter value
            natural_blid_to_parameter_vec[bnid][blid] = boundary_parameter_function_map[bcid](position_x, position_y, value_vec);

        }

    }

}

void BoundaryUnit::read_boundary(std::string file_in_str)
{

    // read file with natural BC data
    std::ifstream file_in_stream(file_in_str);

    // initialize for iteration
    bool is_header = true;  // true while reading header
    std::string line_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_stream, line_str))
    {

        // skip header
        if (is_header)
        {
            is_header = false; // not reading header
            continue;
        }

        // count number of boundaries
        num_boundary++;

        // convert line string into stringstream
        std::stringstream line_stream(line_str);

        // initialize for iteration
        int value_num = 0;  // counts position of value
        std::string value_str;  // stores values in lines
        VectorInt plid_vec_sub;  // vector of global point IDs

        // iterate through each value
        while (std::getline(line_stream, value_str, ','))
        {

            // store values in appropriate vector
            if (value_num == 0)
            {
                boundary_bcid_vec.push_back(std::stoi(value_str));
            }
            else if (value_num == 1)
            {
                boundary_egid_vec.push_back(std::stoi(value_str));
            }
            else if (value_num < num_neighbor + 2)
            {
                plid_vec_sub.push_back(std::stoi(value_str));
            }

            // store vector of global point IDs
            if (value_num == num_neighbor + 1)
            {
                boundary_blid_to_plid_vec.push_back(plid_vec_sub);
            }

            // increment value count
            value_num++;

        }

    }

    // close point file
    file_in_stream.close();    

    // fill up vectors with preliminary values
    boundary_btid_vec = VectorInt(num_boundary, 0);
    boundary_blid_to_parameter_vec = VectorDouble3D(num_boundary, VectorDouble2D(num_neighbor, VectorDouble{}));

}

}

#endif
