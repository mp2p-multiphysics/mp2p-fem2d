#ifndef VARIABLE_TRI3
#define VARIABLE_TRI3
#include <fstream>
#include <sstream>
#include "container_typedef.hpp"
#include "mesh_tri3.hpp"

class VariableTri3
{
    /*

    Variable applied over tri3 mesh elements.

    Variables
    =========
    mesh_in : MeshTri3
        Mesh where variable value is applied.
    u_init_in : double
        Initial value of the variable.

    Functions
    =========
    output_csv : void
        Outputs a CSV file with the values of the variable.

    */

    public:

    // values in variable
    int num_point_domain = 0;  // number of points in domain
    VectorDouble point_value_vec;  // key: domain ID; value: value
    
    // mesh where variable is applied
    MeshTri3* mesh_ptr;  

    // functions
    void output_csv(std::string file_out_str);
    void output_csv(std::string file_out_base_str, int ts);

    // default constructor
    VariableTri3()
    {

    }

    // constructor
    VariableTri3(MeshTri3 &mesh_in, double u_init_in)
    {

        // store mesh
        mesh_ptr = &mesh_in;

        // get number of domain points
        num_point_domain = mesh_ptr->num_point_domain;

        // populate value vector with initial values
        for (int point_did = 0; point_did < num_point_domain; point_did++)
        {
            point_value_vec.push_back(u_init_in);
        }

    }

};

void VariableTri3::output_csv(std::string file_out_str)
{
    /*

    Outputs a CSV file with the values of the variable.

    Arguments
    =========
    file_out_str : string
        Path to CSV file.

    Returns
    =======
    (none)

    Notes
    =====
    This function is intended to be used with steady-state simulations.

    */

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "gid,position_x,position_y,value\n";
    for (int point_did = 0; point_did < num_point_domain; point_did++)
    {
        file_out_stream << mesh_ptr->point_gid_vec[point_did] << ",";
        file_out_stream << mesh_ptr->point_position_x_vec[point_did] << ",";
        file_out_stream << mesh_ptr->point_position_y_vec[point_did] << ",";
        file_out_stream << point_value_vec[point_did] << "\n";
    }

}

void VariableTri3::output_csv(std::string file_out_base_str, int ts)
{
    /*

    Outputs a CSV file with the values of the variable.

    Arguments
    =========
    file_out_base_str : string
        Path to CSV file with base file name.
    ts : int
        Timestep number.

    Returns
    =======
    (none)

    Notes
    =====
    file_out_base_str must have an asterisk '*', which will be replaced with ts.
    This function is intended to be used with transient simulations.

    */

    // split filename at '*'
    // will be replaced with timestep later
    std::vector<std::string> file_out_base_vec;
    std::stringstream file_out_base_stream(file_out_base_str);
    std::string string_sub;
    while(std::getline(file_out_base_stream, string_sub, '*'))
    {
        file_out_base_vec.push_back(string_sub);
    }

    // create output filename
    // replace '*' with timestep
    std::string file_out_str = file_out_base_vec[0];
    for (int i = 1; i < file_out_base_vec.size(); i++)
    {
        file_out_str += std::to_string(ts) + file_out_base_vec[i];
    }

    // initialize file stream
    std::ofstream file_out_stream(file_out_str);

    // write to file
    file_out_stream << "gid,position_x,position_y,value\n";
    for (int point_did = 0; point_did < num_point_domain; point_did++)
    {
        file_out_stream << mesh_ptr->point_gid_vec[point_did] << ",";
        file_out_stream << mesh_ptr->point_position_x_vec[point_did] << ",";
        file_out_stream << mesh_ptr->point_position_y_vec[point_did] << ",";
        file_out_stream << point_value_vec[point_did] << "\n";
    }

}

#endif
