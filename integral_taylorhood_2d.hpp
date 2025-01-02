#ifndef INTEGRAL_TAYLORHOOD_2D
#define INTEGRAL_TAYLORHOOD_2D
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_2d.hpp"

namespace FEM2D
{

class IntegralTaylorHood2D
{
    /*

    Taylor-Hood test function (mixed linear (N) and quadratic (M)) integrals over a domain.

    Variables
    =========
    domain_line_in : Domain2D
        Linear domain where element integrals are calculated.
    domain_quad_in : Domain2D
        Quadratic domain where element integrals are calculated.

    Functions
    =========
    evaluate_integral_Ni_derivative_Mj_x_line : void
        Calculates the integral of Ni * d(Mj)/dx over test functions in the linear domain.
    evaluate_integral_Ni_derivative_Mj_y_line : void
        Calculates the integral of Ni * d(Mj)/dy over test functions in the linear domain.
    evaluate_integral_Mi_quad : void
        Calculates the integral of Mi over test functions in the quadratic domain.
    evaluate_integral_Mi_Mj_quad : void
        Calculates the integral of Mi * Mj over test functions in the quadratic domain.
    evaluate_integral_Mi_derivative_Nj_x_quad : void
        Calculates the integral of Mi * d(Nj)/dx over test functions in the quadratic domain.
    evaluate_integral_Mi_derivative_Nj_y_quad : void
        Calculates the integral of Mi * d(Nj)/dy over test functions in the quadratic domain.
    evaluate_integral_div_Mi_dot_div_Mj_quad : void
        Calculates the integral of div(Mi) dot div(Mj) over test functions in the quadratic domain.
    evaluate_integral_Mi_Mk_derivative_Mj_x_quad : void
        Calculates the integral of Mi * Mk * d(Mj)/dx over test functions in the quadratic domain.
    evaluate_integral_Mi_Mk_derivative_Mj_y_quad : void
        Calculates the integral of Mi * Mk * d(Mj)/dy over test functions in the quadratic domain.

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the domain integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][i][j]...
        wherein edid is the domain element ID and i, j, ... are indices.

    */

    public:
    
    // domain where integral is applied
    Domain2D *domain_line_ptr;
    Domain2D *domain_quad_ptr;

    // vectors with linear test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_line_vec;
    VectorDouble3D Ni_line_vec;
    VectorDouble3D Mi_line_vec;
    VectorDouble3D derivative_Ni_x_line_vec;
    VectorDouble3D derivative_Ni_y_line_vec;
    VectorDouble3D derivative_Mi_x_line_vec;
    VectorDouble3D derivative_Mi_y_line_vec;
    VectorDouble2D jacobian_determinant_line_vec;

    // vectors with quadratic test functions
    // index as follows: [edid][integration_point][i]
    VectorDouble weight_quad_vec;
    VectorDouble3D Ni_quad_vec;
    VectorDouble3D Mi_quad_vec;
    VectorDouble3D derivative_Ni_x_quad_vec;
    VectorDouble3D derivative_Ni_y_quad_vec;
    VectorDouble3D derivative_Mi_x_quad_vec;
    VectorDouble3D derivative_Mi_y_quad_vec;
    VectorDouble2D jacobian_determinant_quad_vec;

    // vectors with linear domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble3D integral_Ni_derivative_Mj_x_line_vec;
    VectorDouble3D integral_Ni_derivative_Mj_y_line_vec;
    VectorDouble3D integral_Ni_Nj_line_vec;

    // vectors with quadratic domain integrals
    // index as follows: [edid][i][j][k]
    VectorDouble2D integral_Mi_quad_vec;
    VectorDouble3D integral_Mi_Mj_quad_vec;
    VectorDouble3D integral_Mi_derivative_Nj_x_quad_vec;
    VectorDouble3D integral_Mi_derivative_Nj_y_quad_vec;
    VectorDouble3D integral_div_Mi_dot_div_Mj_quad_vec;
    VectorDouble4D integral_Mi_Mk_derivative_Mj_x_quad_vec;
    VectorDouble4D integral_Mi_Mk_derivative_Mj_y_quad_vec;

    // functions for computing linear domain integrals
    void evaluate_integral_Ni_derivative_Mj_x_line();
    void evaluate_integral_Ni_derivative_Mj_y_line();
    void evaluate_integral_Ni_Nj_line();

    // functions for computing quadratic domain integrals
    void evaluate_integral_Mi_quad();
    void evaluate_integral_Mi_Mj_quad();
    void evaluate_integral_Mi_derivative_Nj_x_quad();
    void evaluate_integral_Mi_derivative_Nj_y_quad();
    void evaluate_integral_div_Mi_dot_div_Mj_quad();
    void evaluate_integral_Mi_Mk_derivative_Mj_x_quad();
    void evaluate_integral_Mi_Mk_derivative_Mj_y_quad();

    // default constructor
    IntegralTaylorHood2D() {}

    // constructor
    IntegralTaylorHood2D(Domain2D &domain_line_in, Domain2D &domain_quad_in)
    {
        
        // store domain and boundaries
        domain_line_ptr = &domain_line_in;
        domain_quad_ptr = &domain_quad_in;

        // evaluate test functions
        if (domain_line_ptr->type_element == 0 && domain_quad_ptr->type_element == 2)
        {
            evaluate_Mi_Ni_tri();
        }
        else if (domain_line_ptr->type_element == 1 && domain_quad_ptr->type_element == 3)
        {
            evaluate_Mi_Ni_quad();
        }

    }

    private:
    void evaluate_Mi_Ni_tri();
    void evaluate_Mi_Ni_quad();

};

void IntegralTaylorHood2D::evaluate_integral_Ni_derivative_Mj_x_line()
{
    /*

    Calculates the integral of Ni * d(Mj)/dx over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * derivative_Mi_x_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Mj_x_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Ni_derivative_Mj_y_line()
{
    /*

    Calculates the integral of Ni * d(Mj)/dy over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * derivative_Mi_y_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_derivative_Mj_y_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Ni_Nj_line()
{
    /*

    Calculates the integral of Ni * Nj over test functions in the linear domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_line_vec.size(); indx_l++) 
        {
            integral_value += weight_line_vec[indx_l] * jacobian_determinant_line_vec[edid][indx_l] * Ni_line_vec[edid][indx_l][indx_i] * Ni_line_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Ni_Nj_line_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Mi_quad()
{
    /*

    Calculates the integral of Mi over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i];
        }
        integral_part_i_vec.push_back(integral_value);
    
    }
    integral_Mi_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Mi_Mj_quad()
{
    /*

    Calculates the integral of Mi * Mj over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mj_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Mi_derivative_Nj_x_quad()
{
    /*

    Calculates the integral of Mi * d(Nj)/dx over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * derivative_Ni_x_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_derivative_Nj_x_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Mi_derivative_Nj_y_quad()
{
    /*

    Calculates the integral of Mi * d(Nj)/dy over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * derivative_Ni_y_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_derivative_Nj_y_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_div_Mi_dot_div_Mj_quad()
{
    /*

    Calculates the integral of div(Mi) dot div(Mj) over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */    

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble2D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
        
        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * (derivative_Mi_x_quad_vec[edid][indx_l][indx_i] * derivative_Mi_x_quad_vec[edid][indx_l][indx_j] + derivative_Mi_y_quad_vec[edid][indx_l][indx_i] * derivative_Mi_y_quad_vec[edid][indx_l][indx_j]);
        }
        integral_part_ij_vec.push_back(integral_value);
    
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_div_Mi_dot_div_Mj_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Mi_Mk_derivative_Mj_x_quad()
{
    /*

    Calculates the integral of Mi * Mk * d(Mj)/dx over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */ 

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
    VectorDouble integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_k] * derivative_Mi_x_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mk_derivative_Mj_x_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_integral_Mi_Mk_derivative_Mj_y_quad()
{
    /*

    Calculates the integral of Mi * Mk * d(Mj)/dy over test functions in the quadratic domain.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */ 

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++){  
    
    // iterate for each test function combination
    VectorDouble3D integral_part_i_vec;
    for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++){  
    VectorDouble2D integral_part_ij_vec;
    for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++){
    VectorDouble integral_part_ijk_vec;
    for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++){

        // iterate for each integration point
        double integral_value = 0;
        for (int indx_l = 0; indx_l < weight_quad_vec.size(); indx_l++) 
        {
            integral_value += weight_quad_vec[indx_l] * jacobian_determinant_quad_vec[edid][indx_l] * Mi_quad_vec[edid][indx_l][indx_i] * Mi_quad_vec[edid][indx_l][indx_k] * derivative_Mi_y_quad_vec[edid][indx_l][indx_j];
        }
        integral_part_ijk_vec.push_back(integral_value);
    
    }
    integral_part_ij_vec.push_back(integral_part_ijk_vec);
    }
    integral_part_i_vec.push_back(integral_part_ij_vec);
    }
    integral_Mi_Mk_derivative_Mj_y_quad_vec.push_back(integral_part_i_vec);

    }

}

void IntegralTaylorHood2D::evaluate_Mi_Ni_tri()
{

    // linear integration points

    // weights for integration
    weight_line_vec = {0.5};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    const double M_1_3 = 1./3.;
    double a_line_arr[3] = {M_1_3};
    double b_line_arr[3] = {M_1_3};

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Mi_x_part_ml_vec;
        VectorDouble2D derivative_Mi_y_part_ml_vec;

        // get points around element
        int p0_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p0_pdid = domain_line_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_line_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_line_ptr->point_pgid_to_pdid_map[p2_pgid];

        // get x values of points
        double x0 = domain_line_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_line_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_line_ptr->point_position_x_vec[p2_pdid];

        // get y values of points
        double y0 = domain_line_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_line_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_line_ptr->point_position_y_vec[p2_pdid];
 
        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 1; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_Mi_x_part_mli_vec;
            VectorDouble derivative_Mi_y_part_mli_vec;
            VectorDouble derivative_Ni_x_part_mli_vec;
            VectorDouble derivative_Ni_y_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_line_arr[indx_l];
            double b = b_line_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = x0 - x1;
            double derivative_x_b = x2 - x1;
            double derivative_y_a = y0 - y1;
            double derivative_y_b = y2 - y1;

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 3; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = a; break;
                    case 1: N = 1. - a - b; break;
                    case 2: N = b; break;
                }

                // get derivatives of test function N
                double derivative_Ni_a = 0.;
                double derivative_Ni_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Ni_a =  1.; derivative_Ni_b =  0.; break;
                    case 1: derivative_Ni_a = -1.; derivative_Ni_b = -1.; break;
                    case 2: derivative_Ni_a =  0.; derivative_Ni_b =  1.; break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_Ni_ab_vec;
                derivative_Ni_ab_vec << derivative_Ni_a, derivative_Ni_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_Ni_xy_vec = derivative_Ni_ab_vec * jacobian_inverse_mat;
                double derivative_Ni_x = derivative_Ni_xy_vec.coeffRef(0);
                double derivative_Ni_y = derivative_Ni_xy_vec.coeffRef(1);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_Ni_x_part_mli_vec.push_back(derivative_Ni_x);
                derivative_Ni_y_part_mli_vec.push_back(derivative_Ni_y);

            }

            // iterate for each test function
            for (int indx_i = 0; indx_i < 6; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case 0: M = a*(2*a - 1); break;
                    case 1: M = 4*a*(-a - b + 1); break;
                    case 2: M = (a + b - 1)*(2*a + 2*b - 1); break;
                    case 3: M = 4*b*(-a - b + 1); break;
                    case 4: M = b*(2*b - 1); break;
                    case 5: M = 4*a*b; break;
                }

                // get derivatives of test function M
                double derivative_Mi_a = 0.;
                double derivative_Mi_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Mi_a = 4*a - 1;        derivative_Mi_b = 0; break;
                    case 1: derivative_Mi_a = -8*a - 4*b + 4; derivative_Mi_b = -4*a; break;
                    case 2: derivative_Mi_a = 4*a + 4*b - 3;  derivative_Mi_b = 4*a + 4*b - 3; break;
                    case 3: derivative_Mi_a = -4*b;           derivative_Mi_b = -4*a - 8*b + 4; break;
                    case 4: derivative_Mi_a = 0;              derivative_Mi_b = 4*b - 1; break;
                    case 5: derivative_Mi_a = 4*b;            derivative_Mi_b = 4*a; break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_Mi_ab_vec;
                derivative_Mi_ab_vec << derivative_Mi_a, derivative_Mi_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_Mi_xy_vec = derivative_Mi_ab_vec * jacobian_inverse_mat;
                double derivative_Mi_x = derivative_Mi_xy_vec.coeffRef(0);
                double derivative_Mi_y = derivative_Mi_xy_vec.coeffRef(1);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_Mi_x_part_mli_vec.push_back(derivative_Mi_x);
                derivative_Mi_y_part_mli_vec.push_back(derivative_Mi_y);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Mi_x_part_ml_vec.push_back(derivative_Mi_x_part_mli_vec);
            derivative_Mi_y_part_ml_vec.push_back(derivative_Mi_y_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_line_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_line_vec.push_back(N_part_ml_vec);
        Mi_line_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_line_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_line_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Mi_x_line_vec.push_back(derivative_Mi_x_part_ml_vec);
        derivative_Mi_y_line_vec.push_back(derivative_Mi_y_part_ml_vec);
 
    }

    // quadratic integration points

    // weights for integration
    const double M_1_6 = 1./6.;
    weight_quad_vec = {M_1_6, M_1_6, M_1_6};

    // integration points
    // dimensionless coordinates if element is scaled to [0, 1] * [0, 1]
    double a_quad_arr[3] = {0.5, 0.5, 0.0};
    double b_quad_arr[3] = {0.5, 0.0, 0.5};

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_Ni_x_part_ml_vec;
        VectorDouble2D derivative_Ni_y_part_ml_vec;
        VectorDouble2D derivative_Mi_x_part_ml_vec;
        VectorDouble2D derivative_Mi_y_part_ml_vec;

        // get points around element
        int p0_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p4_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][4];
        int p5_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][5];
        int p0_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p3_pgid];
        int p4_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p4_pgid];
        int p5_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p5_pgid];

        // get x values of points
        double x0 = domain_quad_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_quad_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_quad_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_quad_ptr->point_position_x_vec[p3_pdid];
        double x4 = domain_quad_ptr->point_position_x_vec[p4_pdid];
        double x5 = domain_quad_ptr->point_position_x_vec[p5_pdid];

        // get y values of points
        double y0 = domain_quad_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_quad_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_quad_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_quad_ptr->point_position_y_vec[p3_pdid];
        double y4 = domain_quad_ptr->point_position_y_vec[p4_pdid];
        double y5 = domain_quad_ptr->point_position_y_vec[p5_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 3; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_Mi_x_part_mli_vec;
            VectorDouble derivative_Mi_y_part_mli_vec;
            VectorDouble derivative_Ni_x_part_mli_vec;
            VectorDouble derivative_Ni_y_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_quad_arr[indx_l];
            double b = b_quad_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = 4*a*x0 - 8*a*x1 + 4*a*x2 - 4*b*x1 + 4*b*x2 - 4*b*x3 + 4*b*x5 - x0 + 4*x1 - 3*x2;
            double derivative_x_b = -4*a*x1 + 4*a*x2 - 4*a*x3 + 4*a*x5 + 4*b*x2 - 8*b*x3 + 4*b*x4 - 3*x2 + 4*x3 - x4;
            double derivative_y_a = 4*a*y0 - 8*a*y1 + 4*a*y2 - 4*b*y1 + 4*b*y2 - 4*b*y3 + 4*b*y5 - y0 + 4*y1 - 3*y2;
            double derivative_y_b = -4*a*y1 + 4*a*y2 - 4*a*y3 + 4*a*y5 + 4*b*y2 - 8*b*y3 + 4*b*y4 - 3*y2 + 4*y3 - y4;

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 3; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = a; break;
                    case 1: N = 1. - a - b; break;
                    case 2: N = b; break;
                }

                // get derivatives of test function N
                double derivative_Ni_a = 0.;
                double derivative_Ni_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Ni_a =  1.; derivative_Ni_b =  0.; break;
                    case 1: derivative_Ni_a = -1.; derivative_Ni_b = -1.; break;
                    case 2: derivative_Ni_a =  0.; derivative_Ni_b =  1.; break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_Ni_ab_vec;
                derivative_Ni_ab_vec << derivative_Ni_a, derivative_Ni_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_Ni_xy_vec = derivative_Ni_ab_vec * jacobian_inverse_mat;
                double derivative_Ni_x = derivative_Ni_xy_vec.coeffRef(0);
                double derivative_Ni_y = derivative_Ni_xy_vec.coeffRef(1);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_Ni_x_part_mli_vec.push_back(derivative_Ni_x);
                derivative_Ni_y_part_mli_vec.push_back(derivative_Ni_y);

            }

            // iterate for each test function
            for (int indx_i = 0; indx_i < 6; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case 0: M = a*(2*a - 1); break;
                    case 1: M = 4*a*(-a - b + 1); break;
                    case 2: M = (a + b - 1)*(2*a + 2*b - 1); break;
                    case 3: M = 4*b*(-a - b + 1); break;
                    case 4: M = b*(2*b - 1); break;
                    case 5: M = 4*a*b; break;
                }

                // get derivatives of test function M
                double derivative_Mi_a = 0.;
                double derivative_Mi_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_Mi_a = 4*a - 1;        derivative_Mi_b = 0; break;
                    case 1: derivative_Mi_a = -8*a - 4*b + 4; derivative_Mi_b = -4*a; break;
                    case 2: derivative_Mi_a = 4*a + 4*b - 3;  derivative_Mi_b = 4*a + 4*b - 3; break;
                    case 3: derivative_Mi_a = -4*b;           derivative_Mi_b = -4*a - 8*b + 4; break;
                    case 4: derivative_Mi_a = 0;              derivative_Mi_b = 4*b - 1; break;
                    case 5: derivative_Mi_a = 4*b;            derivative_Mi_b = 4*a; break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_Mi_ab_vec;
                derivative_Mi_ab_vec << derivative_Mi_a, derivative_Mi_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_Mi_xy_vec = derivative_Mi_ab_vec * jacobian_inverse_mat;
                double derivative_Mi_x = derivative_Mi_xy_vec.coeffRef(0);
                double derivative_Mi_y = derivative_Mi_xy_vec.coeffRef(1);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_Mi_x_part_mli_vec.push_back(derivative_Mi_x);
                derivative_Mi_y_part_mli_vec.push_back(derivative_Mi_y);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_Ni_x_part_ml_vec.push_back(derivative_Ni_x_part_mli_vec);
            derivative_Ni_y_part_ml_vec.push_back(derivative_Ni_y_part_mli_vec);
            derivative_Mi_x_part_ml_vec.push_back(derivative_Mi_x_part_mli_vec);
            derivative_Mi_y_part_ml_vec.push_back(derivative_Mi_y_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_quad_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_quad_vec.push_back(N_part_ml_vec);
        Mi_quad_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_quad_vec.push_back(derivative_Ni_x_part_ml_vec);
        derivative_Ni_y_quad_vec.push_back(derivative_Ni_y_part_ml_vec);
        derivative_Mi_x_quad_vec.push_back(derivative_Mi_x_part_ml_vec);
        derivative_Mi_y_quad_vec.push_back(derivative_Mi_y_part_ml_vec);
 
    }

}

void IntegralTaylorHood2D::evaluate_Mi_Ni_quad()
{

    // linear integration points
    
    // weights for integration
    weight_line_vec = {1., 1., 1., 1.};

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1] * [-1, 1]
    const double M_1_SQRT_3 = 1./sqrt(3);
    double a_line_arr[4] = {+M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3, -M_1_SQRT_3};
    double b_line_arr[4] = {+M_1_SQRT_3, -M_1_SQRT_3, +M_1_SQRT_3, -M_1_SQRT_3};

    // iterate for each domain element
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_N_x_part_ml_vec;
        VectorDouble2D derivative_N_y_part_ml_vec;
        VectorDouble2D derivative_M_x_part_ml_vec;
        VectorDouble2D derivative_M_y_part_ml_vec;

        // get points around element
        int p0_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_line_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p0_pdid = domain_line_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_line_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_line_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_line_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get x values of points
        double x0 = domain_line_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_line_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_line_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_line_ptr->point_position_x_vec[p3_pdid];

        // get y values of points
        double y0 = domain_line_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_line_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_line_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_line_ptr->point_position_y_vec[p3_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 4; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_N_x_part_mli_vec;
            VectorDouble derivative_N_y_part_mli_vec;
            VectorDouble derivative_M_x_part_mli_vec;
            VectorDouble derivative_M_y_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_line_arr[indx_l];
            double b = b_line_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = 0.25*(b*x0 - b*x1 + b*x2 - b*x3 - x0 - x1 + x2 + x3);
            double derivative_x_b = 0.25*(a*x0 - a*x1 + a*x2 - a*x3 - x0 + x1 + x2 - x3);
            double derivative_y_a = 0.25*(b*y0 - b*y1 + b*y2 - b*y3 - y0 - y1 + y2 + y3);
            double derivative_y_b = 0.25*(a*y0 - a*y1 + a*y2 - a*y3 - y0 + y1 + y2 - y3);

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 4; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 0.25*(1. - a)*(1. - b); break;
                    case 1: N = 0.25*(1. - a)*(1. + b); break;
                    case 2: N = 0.25*(1. + a)*(1. + b); break;
                    case 3: N = 0.25*(1. + a)*(1. - b); break;
                }

                // get derivatives of test function N
                double derivative_N_a = 0.;
                double derivative_N_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_N_a = +0.25*(b - 1.); derivative_N_b = +0.25*(a - 1.); break;
                    case 1: derivative_N_a = -0.25*(b + 1.); derivative_N_b = -0.25*(a - 1.); break;
                    case 2: derivative_N_a = +0.25*(b + 1.); derivative_N_b = +0.25*(a + 1.); break;
                    case 3: derivative_N_a = -0.25*(b - 1.); derivative_N_b = -0.25*(a + 1.); break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_N_ab_vec;
                derivative_N_ab_vec << derivative_N_a, derivative_N_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_N_xy_vec = derivative_N_ab_vec * jacobian_inverse_mat;
                double derivative_N_x = derivative_N_xy_vec.coeffRef(0);
                double derivative_N_y = derivative_N_xy_vec.coeffRef(1);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_N_x_part_mli_vec.push_back(derivative_N_x);
                derivative_N_y_part_mli_vec.push_back(derivative_N_y);

            }

            // iterate for each test function
            for (int indx_i = 0; indx_i < 8; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case 0: M = 0.25*(1 - a)*(1 - b)*(-a - b - 1); break;
                    case 1: M = 0.50*(1 - a)*(1 - b*b); break;
                    case 2: M = 0.25*(1 - a)*(1 + b)*(-a + b - 1); break;
                    case 3: M = 0.50*(1 - a*a)*(1 + b); break;
                    case 4: M = 0.25*(1 + a)*(1 + b)*(+a + b - 1); break;
                    case 5: M = 0.50*(1 + a)*(1 - b*b); break;
                    case 6: M = 0.25*(1 + a)*(1 - b)*(+a - b - 1); break;
                    case 7: M = 0.50*(1 - a*a)*(1 - b); break;
                }

                // get derivatives of test function M
                double derivative_M_a = 0.;
                double derivative_M_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_M_a = 0.25*(-2*a - b)*(b - 1); derivative_M_b = 0.25*(-a - 2*b)*(a - 1); break;
                    case 1: derivative_M_a = 0.50*(b*b - 1);          derivative_M_b = b*(a - 1); break;
                    case 2: derivative_M_a = 0.25*(2*a - b)*(b + 1);  derivative_M_b = 0.25*(a - 1)*(a - 2*b); break;
                    case 3: derivative_M_a = -a*(b + 1);              derivative_M_b = 0.50*(1 - a*a); break;
                    case 4: derivative_M_a = 0.25*(2*a + b)*(b + 1);  derivative_M_b = 0.25*(a + 1)*(a + 2*b); break;
                    case 5: derivative_M_a = 0.5*(1 - b*b);           derivative_M_b = -b*(a + 1); break;
                    case 6: derivative_M_a = 0.25*(-2*a + b)*(b - 1); derivative_M_b = 0.25*(-a + 2*b)*(a + 1); break;
                    case 7: derivative_M_a = a*(b - 1);               derivative_M_b = 0.50*(a*a - 1); break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_M_ab_vec;
                derivative_M_ab_vec << derivative_M_a, derivative_M_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_M_xy_vec = derivative_M_ab_vec * jacobian_inverse_mat;
                double derivative_M_x = derivative_M_xy_vec.coeffRef(0);
                double derivative_M_y = derivative_M_xy_vec.coeffRef(1);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_M_x_part_mli_vec.push_back(derivative_M_x);
                derivative_M_y_part_mli_vec.push_back(derivative_M_y);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_N_x_part_ml_vec.push_back(derivative_N_x_part_mli_vec);
            derivative_N_y_part_ml_vec.push_back(derivative_N_y_part_mli_vec);
            derivative_M_x_part_ml_vec.push_back(derivative_M_x_part_mli_vec);
            derivative_M_y_part_ml_vec.push_back(derivative_M_y_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_line_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_line_vec.push_back(N_part_ml_vec);
        Mi_line_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_line_vec.push_back(derivative_N_x_part_ml_vec);
        derivative_Ni_y_line_vec.push_back(derivative_N_y_part_ml_vec);
        derivative_Mi_x_line_vec.push_back(derivative_M_x_part_ml_vec);
        derivative_Mi_y_line_vec.push_back(derivative_M_y_part_ml_vec);

    }

    // quadratic integration points

    // weights for integration
    const double M_25_81 = 25./81.;
    const double M_40_81 = 40./81.;
    const double M_64_81 = 64./81.;
    weight_quad_vec = {M_25_81, M_25_81, M_25_81, M_25_81, M_40_81, M_40_81, M_40_81, M_40_81, M_64_81};

    // integration points
    // dimensionless coordinates if element is scaled to [-1, 1] * [-1, 1]
    const double M_SQRT_3_SQRT_5 = sqrt(0.6);
    double a_quad_arr[9] = {+M_SQRT_3_SQRT_5, +M_SQRT_3_SQRT_5, -M_SQRT_3_SQRT_5, -M_SQRT_3_SQRT_5, +M_SQRT_3_SQRT_5, -M_SQRT_3_SQRT_5, 0., 0., 0.};
    double b_quad_arr[9] = {+M_SQRT_3_SQRT_5, -M_SQRT_3_SQRT_5, +M_SQRT_3_SQRT_5, -M_SQRT_3_SQRT_5, 0., 0., +M_SQRT_3_SQRT_5, -M_SQRT_3_SQRT_5, 0.};

    // iterate for each domain element
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // initialize
        VectorDouble jacobian_determinant_part_ml_vec;
        VectorDouble2D N_part_ml_vec;
        VectorDouble2D M_part_ml_vec;
        VectorDouble2D derivative_N_x_part_ml_vec;
        VectorDouble2D derivative_N_y_part_ml_vec;
        VectorDouble2D derivative_M_x_part_ml_vec;
        VectorDouble2D derivative_M_y_part_ml_vec;

        // get points around element
        int p0_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p4_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][4];
        int p5_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][5];
        int p6_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][6];
        int p7_pgid = domain_quad_ptr->element_edid_plid_to_pgid_vec[edid][7];
        int p0_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p3_pgid];
        int p4_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p4_pgid];
        int p5_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p5_pgid];
        int p6_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p6_pgid];
        int p7_pdid = domain_quad_ptr->point_pgid_to_pdid_map[p7_pgid];

        // get x values of points
        double x0 = domain_quad_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_quad_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_quad_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_quad_ptr->point_position_x_vec[p3_pdid];
        double x4 = domain_quad_ptr->point_position_x_vec[p4_pdid];
        double x5 = domain_quad_ptr->point_position_x_vec[p5_pdid];
        double x6 = domain_quad_ptr->point_position_x_vec[p6_pdid];
        double x7 = domain_quad_ptr->point_position_x_vec[p7_pdid];

        // get y values of points
        double y0 = domain_quad_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_quad_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_quad_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_quad_ptr->point_position_y_vec[p3_pdid];
        double y4 = domain_quad_ptr->point_position_y_vec[p4_pdid];
        double y5 = domain_quad_ptr->point_position_y_vec[p5_pdid];
        double y6 = domain_quad_ptr->point_position_y_vec[p6_pdid];
        double y7 = domain_quad_ptr->point_position_y_vec[p7_pdid];

        // iterate for each integration point (indx_l)
        for (int indx_l = 0; indx_l < 9; indx_l++)
        {

            // initialize
            VectorDouble N_part_mli_vec;
            VectorDouble M_part_mli_vec;
            VectorDouble derivative_N_x_part_mli_vec;
            VectorDouble derivative_N_y_part_mli_vec;
            VectorDouble derivative_M_x_part_mli_vec;
            VectorDouble derivative_M_y_part_mli_vec;

            // get a and b values where function is evaluated
            double a = a_quad_arr[indx_l];
            double b = b_quad_arr[indx_l];

            // get derivatives of x and y with respect to a and b
            double derivative_x_a = -a*x3*(b + 1) + a*x7*(b - 1) - x0*(a - 1)*(b - 1)*0.25 - x0*(b - 1)*(a + b + 1)*0.25 + x1*(b*b - 1)*0.5 + x2*(a - 1)*(b + 1)*0.25 + x2*(b + 1)*(a - b + 1)*0.25 + x4*(a + 1)*(b + 1)*0.25 + x4*(b + 1)*(a + b - 1)*0.25 - x5*(b*b - 1)*0.5 - x6*(a + 1)*(b - 1)*0.25 + x6*(b - 1)*(-a + b + 1)*0.25;
            double derivative_x_b = b*x1*(a - 1) - b*x5*(a + 1) - x0*(a - 1)*(b - 1)*0.25 - x0*(a - 1)*(a + b + 1)*0.25 - x2*(a - 1)*(b + 1)*0.25 + x2*(a - 1)*(a - b + 1)*0.25 - x3*(a*a - 1)*0.5 + x4*(a + 1)*(b + 1)*0.25 + x4*(a + 1)*(a + b - 1)*0.25 + x6*(a + 1)*(b - 1)*0.25 + x6*(a + 1)*(-a + b + 1)*0.25 + x7*(a*a - 1)*0.5;
            double derivative_y_a = -a*y3*(b + 1) + a*y7*(b - 1) - y0*(a - 1)*(b - 1)*0.25 - y0*(b - 1)*(a + b + 1)*0.25 + y1*(b*b - 1)*0.5 + y2*(a - 1)*(b + 1)*0.25 + y2*(b + 1)*(a - b + 1)*0.25 + y4*(a + 1)*(b + 1)*0.25 + y4*(b + 1)*(a + b - 1)*0.25 - y5*(b*b - 1)*0.5 - y6*(a + 1)*(b - 1)*0.25 + y6*(b - 1)*(-a + b + 1)*0.25;
            double derivative_y_b = b*y1*(a - 1) - b*y5*(a + 1) - y0*(a - 1)*(b - 1)*0.25 - y0*(a - 1)*(a + b + 1)*0.25 - y2*(a - 1)*(b + 1)*0.25 + y2*(a - 1)*(a - b + 1)*0.25 - y3*(a*a - 1)*0.5 + y4*(a + 1)*(b + 1)*0.25 + y4*(a + 1)*(a + b - 1)*0.25 + y6*(a + 1)*(b - 1)*0.25 + y6*(a + 1)*(-a + b + 1)*0.25 + y7*(a*a - 1)*0.5;

            // get jacobian and its inverse and determinant
            Eigen::Matrix2d jacobian_mat;
            jacobian_mat << derivative_x_a, derivative_x_b, derivative_y_a, derivative_y_b;
            Eigen::Matrix2d jacobian_inverse_mat = jacobian_mat.inverse();
            double jacobian_determinant = jacobian_mat.determinant();

            // iterate for each test function
            for (int indx_i = 0; indx_i < 4; indx_i++)
            {
        
                // get test function N
                double N = 0.;
                switch (indx_i)
                {
                    case 0: N = 0.25*(1. - a)*(1. - b); break;
                    case 1: N = 0.25*(1. - a)*(1. + b); break;
                    case 2: N = 0.25*(1. + a)*(1. + b); break;
                    case 3: N = 0.25*(1. + a)*(1. - b); break;
                }

                // get derivatives of test function N
                double derivative_N_a = 0.;
                double derivative_N_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_N_a = +0.25*(b - 1.); derivative_N_b = +0.25*(a - 1.); break;
                    case 1: derivative_N_a = -0.25*(b + 1.); derivative_N_b = -0.25*(a - 1.); break;
                    case 2: derivative_N_a = +0.25*(b + 1.); derivative_N_b = +0.25*(a + 1.); break;
                    case 3: derivative_N_a = -0.25*(b - 1.); derivative_N_b = -0.25*(a + 1.); break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_N_ab_vec;
                derivative_N_ab_vec << derivative_N_a, derivative_N_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_N_xy_vec = derivative_N_ab_vec * jacobian_inverse_mat;
                double derivative_N_x = derivative_N_xy_vec.coeffRef(0);
                double derivative_N_y = derivative_N_xy_vec.coeffRef(1);

                // store in vectors
                N_part_mli_vec.push_back(N);
                derivative_N_x_part_mli_vec.push_back(derivative_N_x);
                derivative_N_y_part_mli_vec.push_back(derivative_N_y);

            }

            // iterate for each test function
            for (int indx_i = 0; indx_i < 8; indx_i++)
            {
        
                // get test function M
                double M = 0.;
                switch (indx_i)
                {
                    case 0: M = 0.25*(1 - a)*(1 - b)*(-a - b - 1); break;
                    case 1: M = 0.50*(1 - a)*(1 - b*b); break;
                    case 2: M = 0.25*(1 - a)*(1 + b)*(-a + b - 1); break;
                    case 3: M = 0.50*(1 - a*a)*(1 + b); break;
                    case 4: M = 0.25*(1 + a)*(1 + b)*(+a + b - 1); break;
                    case 5: M = 0.50*(1 + a)*(1 - b*b); break;
                    case 6: M = 0.25*(1 + a)*(1 - b)*(+a - b - 1); break;
                    case 7: M = 0.50*(1 - a*a)*(1 - b); break;
                }

                // get derivatives of test function M
                double derivative_M_a = 0.;
                double derivative_M_b = 0.;
                switch (indx_i)
                {
                    case 0: derivative_M_a = 0.25*(-2*a - b)*(b - 1); derivative_M_b = 0.25*(-a - 2*b)*(a - 1); break;
                    case 1: derivative_M_a = 0.50*(b*b - 1);          derivative_M_b = b*(a - 1); break;
                    case 2: derivative_M_a = 0.25*(2*a - b)*(b + 1);  derivative_M_b = 0.25*(a - 1)*(a - 2*b); break;
                    case 3: derivative_M_a = -a*(b + 1);              derivative_M_b = 0.50*(1 - a*a); break;
                    case 4: derivative_M_a = 0.25*(2*a + b)*(b + 1);  derivative_M_b = 0.25*(a + 1)*(a + 2*b); break;
                    case 5: derivative_M_a = 0.5*(1 - b*b);           derivative_M_b = -b*(a + 1); break;
                    case 6: derivative_M_a = 0.25*(-2*a + b)*(b - 1); derivative_M_b = 0.25*(-a + 2*b)*(a + 1); break;
                    case 7: derivative_M_a = a*(b - 1);               derivative_M_b = 0.50*(a*a - 1); break;
                }

                // get vector with derivatives of test functions wrt a and b
                Eigen::RowVector2d derivative_M_ab_vec;
                derivative_M_ab_vec << derivative_M_a, derivative_M_b;

                // get vector with derivatives of test functions wrt x and y
                // multiply derivatives wrt a and b with inverse jacobian
                Eigen::RowVector2d derivative_M_xy_vec = derivative_M_ab_vec * jacobian_inverse_mat;
                double derivative_M_x = derivative_M_xy_vec.coeffRef(0);
                double derivative_M_y = derivative_M_xy_vec.coeffRef(1);

                // store in vectors
                M_part_mli_vec.push_back(M);
                derivative_M_x_part_mli_vec.push_back(derivative_M_x);
                derivative_M_y_part_mli_vec.push_back(derivative_M_y);

            }

            // store in vectors
            jacobian_determinant_part_ml_vec.push_back(jacobian_determinant);
            N_part_ml_vec.push_back(N_part_mli_vec);
            M_part_ml_vec.push_back(M_part_mli_vec);
            derivative_N_x_part_ml_vec.push_back(derivative_N_x_part_mli_vec);
            derivative_N_y_part_ml_vec.push_back(derivative_N_y_part_mli_vec);
            derivative_M_x_part_ml_vec.push_back(derivative_M_x_part_mli_vec);
            derivative_M_y_part_ml_vec.push_back(derivative_M_y_part_mli_vec);

        }

        // store in vectors
        jacobian_determinant_quad_vec.push_back(jacobian_determinant_part_ml_vec);
        Ni_quad_vec.push_back(N_part_ml_vec);
        Mi_quad_vec.push_back(M_part_ml_vec);
        derivative_Ni_x_quad_vec.push_back(derivative_N_x_part_ml_vec);
        derivative_Ni_y_quad_vec.push_back(derivative_N_y_part_ml_vec);
        derivative_Mi_x_quad_vec.push_back(derivative_M_x_part_ml_vec);
        derivative_Mi_y_quad_vec.push_back(derivative_M_y_part_ml_vec);

    }

}

}

#endif
