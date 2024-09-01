#ifndef SCALAR_QUAD4
#define SCALAR_QUAD4
#include <vector>
#include "Eigen/Eigen"
#include "grid_quad4.hpp"

class ScalarQuad4Class
{

    public:
    
    // variables
    GridQuad4Struct gq4s;
    std::vector<double> scalar_vec;

    // functions
    void set_from_constant(double value);
    void set_from_vector(std::vector<double> value_vec, int start_id);
    void set_from_evector(Eigen::VectorXd value_evec, int start_id);

    // constructor
    ScalarQuad4Class(GridQuad4Struct &gq4s_in)
    {
        gq4s = gq4s_in;
    }

    ScalarQuad4Class(GridQuad4Struct &gq4s_in, double value)
    {
        gq4s = gq4s_in;
        set_from_constant(value);
    }

    ScalarQuad4Class(GridQuad4Struct &gq4s_in, std::vector<double> value_vec, int start_id)
    {
        gq4s = gq4s_in;
        set_from_vector(value_vec, start_id);
    }

    ScalarQuad4Class(GridQuad4Struct &gq4s_in, Eigen::VectorXd value_evec, int start_id)
    {
        gq4s = gq4s_in;
        set_from_evector(value_evec, start_id);
    }

};

void ScalarQuad4Class::set_from_constant(double value)
{

    // set each scalar_vec value to the constant
    for (int n = 0; n < gq4s.num_point; n++)
    {
        scalar_vec.push_back(value);
    }

}

void ScalarQuad4Class::set_from_vector(std::vector<double> value_vec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gq4s.num_point; n++)
    {
        scalar_vec.push_back(value_vec[start_id + n]);
    }

}

void ScalarQuad4Class::set_from_evector(Eigen::VectorXd value_evec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gq4s.num_point; n++)
    {
        scalar_vec.push_back(value_evec.coeffRef(start_id + n));
    }

}

#endif