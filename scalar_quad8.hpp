#ifndef SCALAR_QUAD8
#define SCALAR_QUAD8
#include <vector>
#include "Eigen/Eigen"
#include "grid_quad8.hpp"

class ScalarQuad8Class
{

    public:
    
    // variables
    GridQuad8Struct gq8s;
    std::vector<double> scalar_vec;

    // functions
    void set_from_constant(double value);
    void set_from_vector(std::vector<double> value_vec, int start_id);
    void set_from_evector(Eigen::VectorXd value_evec, int start_id);

    // constructor
    ScalarQuad8Class(GridQuad8Struct &gq8s_in)
    {
        gq8s = gq8s_in;
    }

    ScalarQuad8Class(GridQuad8Struct &gq8s_in, double value)
    {
        gq8s = gq8s_in;
        set_from_constant(value);
    }

    ScalarQuad8Class(GridQuad8Struct &gq8s_in, std::vector<double> value_vec, int start_id)
    {
        gq8s = gq8s_in;
        set_from_vector(value_vec, start_id);
    }

    ScalarQuad8Class(GridQuad8Struct &gq8s_in, Eigen::VectorXd value_evec, int start_id)
    {
        gq8s = gq8s_in;
        set_from_evector(value_evec, start_id);
    }

};

void ScalarQuad8Class::set_from_constant(double value)
{

    // set each scalar_vec value to the constant
    for (int n = 0; n < gq8s.num_point; n++)
    {
        scalar_vec.push_back(value);
    }

}

void ScalarQuad8Class::set_from_vector(std::vector<double> value_vec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gq8s.num_point; n++)
    {
        scalar_vec.push_back(value_vec[start_id + n]);
    }

}

void ScalarQuad8Class::set_from_evector(Eigen::VectorXd value_evec, int start_id)
{
    
    // set scalar_vec to portion of value_vec
    for (int n = 0; n < gq8s.num_point; n++)
    {
        scalar_vec.push_back(value_evec.coeffRef(start_id + n));
    }

}

#endif