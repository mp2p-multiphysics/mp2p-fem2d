#ifndef BOUNDARYINTEGRAL_UNIT
#define BOUNDARYINTEGRAL_UNIT
#include <vector>
#include "Eigen/Eigen"
#include "boundary_unit.hpp"
#include "container_typedef.hpp"
#include "variable_unit.hpp"

namespace FEM2D
{

class BoundaryIntegralUnit
{
    /*

    Test function (N) integrals over a boundary.

    Variables
    =========
    boundary_in : BoundaryUnit
        Boundaries where boundary integrals are calculated.

    Functions
    =========
    get_boundary_key : int
        Calculates the boundary key used to access boundary integral values.
    evaluate_integral_Ni : void
        Calculates the integral of Ni at the boundary.
    evaluate_integral_Ni_Nj : void
        Calculates the integral of Ni * Nj at the boundary.

    Notes
    ====
    The calculated integrals are stored in nested maps and vectors.
    Values of the boundary integrals can be accessed from each vector using the following pattern:
        integral_vec[edid][boundary_key][i][j]...
        wherein boundary_key is an int denoting the location of the boundary.

    */

    public:
    
    // variable and boundary
    VariableUnit *variable_ptr;
    BoundaryUnit *boundary_ptr;

    // domain where boundary is applied
    DomainUnit *domain_ptr;

    // vectors with boundary test functions
    MapVector4D Ni_vec;  // [edid][boundary_key][integration_point][i]
    MapVector2D normal_x_vec;  // [edid][boundary_key]
    MapVector2D normal_y_vec;  // [edid][boundary_key]
    MapVector2D jacobian_determinant_vec;  // [edid][boundary_key]

    // vectors with boundary integrals
    // index as follows: [edid][boundary_key][i][j]
    MapVector3D integral_Ni_vec;
    MapVector4D integral_Ni_Nj_vec;

    // functions for computing boundary integrals
    int get_boundary_key(VectorInt plid_vec);
    void evaluate_integral_Ni();
    void evaluate_integral_Ni_Nj();

    // default constructor
    BoundaryIntegralUnit() {}

    // constructor
    BoundaryIntegralUnit(BoundaryUnit &boundary_in)
    {
        
        // store boundary
        boundary_ptr = &boundary_in;

        // get domain and variable is boundary is applied
        variable_ptr = boundary_ptr->variable_ptr;
        domain_ptr = variable_ptr->domain_ptr;

        // evaluate test functions
        switch (domain_ptr->type_element)
        {
            case 0:
                evaluate_Ni_tri3();
            break;
            case 1:
                evaluate_Ni_quad4();
            break;
        }

    }

    private:
    void evaluate_Ni_tri3();
    void evaluate_Ni_quad4();

};

int BoundaryIntegralUnit::get_boundary_key(VectorInt plid_vec)
{
    /*

    Calculates the boundary key used to access boundary integral values.

    Arguments
    =========
    plid_vec : VectorInt
        vector of local points along a boundary.

    Returns
    =======
    boundary_key : int
        int representing each boundary.

    */

    // sort vector
    std::sort(plid_vec.begin(), plid_vec.end());

    // initialize boundary key
    int boundary_key = 0;

    // form boundary key by combining digits
    for (int indx_i = 0; indx_i < plid_vec.size(); indx_i++)
    {
        boundary_key += plid_vec[indx_i] * pow(10, plid_vec.size() - indx_i - 1);
    }

    return boundary_key;

}

void BoundaryIntegralUnit::evaluate_integral_Ni()
{
    /*

    Calculates the integral of Ni at the boundary.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        VectorInt plid_vec = boundary_ptr->boundary_blid_to_plid_vec[bid];
        int boundary_key = get_boundary_key(plid_vec);

        // iterate for each test function combination
        VectorDouble integral_part_i_vec;
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){

            // iterate for each integration point
            double integral_value = 0;
            for (int indx_l = 0; indx_l < boundary_ptr->num_neighbor; indx_l++) 
            {
                integral_value += jacobian_determinant_vec[edid][boundary_key] * Ni_vec[edid][boundary_key][indx_l][indx_i];
            }
            integral_part_i_vec.push_back(integral_value);
        
        }
        integral_Ni_vec[edid][boundary_key] = integral_part_i_vec;

    }

}

void BoundaryIntegralUnit::evaluate_integral_Ni_Nj()
{
    /*

    Calculates the integral of Ni * Nj at the boundary.

    Arguments
    =========
    (none)

    Returns
    =========
    (none)

    */

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        VectorInt plid_vec = boundary_ptr->boundary_blid_to_plid_vec[bid];
        int boundary_key = get_boundary_key(plid_vec);

        // iterate for each test function combination
        VectorDouble2D integral_part_i_vec;
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){
        VectorDouble integral_part_ij_vec;
        for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){

            // iterate for each integration point
            double integral_value = 0;
            for (int indx_l = 0; indx_l < boundary_ptr->num_neighbor; indx_l++) 
            {
                integral_value += jacobian_determinant_vec[edid][boundary_key] * Ni_vec[edid][boundary_key][indx_l][indx_i] * Ni_vec[edid][boundary_key][indx_l][indx_j];
            }
            integral_part_ij_vec.push_back(integral_value);
        
        }
        integral_part_i_vec.push_back(integral_part_ij_vec);
        }
        integral_Ni_Nj_vec[edid][boundary_key] = integral_part_i_vec;

    }

}

void BoundaryIntegralUnit::evaluate_Ni_tri3()
{

    // initialize 2D integration points (within element)
    const double M_1_SQRT_3 = 1./sqrt(3);
    VectorDouble sample_vec = {0.5*(-M_1_SQRT_3 + 1.), 0.5*(M_1_SQRT_3 + 1.)};
    VectorDouble zero_vec = {0., 0.};

    // initialize 1D integration points (along edge)
    // use capital letters (e.g., X, A) to denote variables along an edge
    VectorDouble A_vec = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int pa_plid = boundary_ptr->boundary_blid_to_plid_vec[bid][0];
        int pb_plid = boundary_ptr->boundary_blid_to_plid_vec[bid][1];
        VectorInt plid_vec = boundary_ptr->boundary_blid_to_plid_vec[bid];
        int boundary_key = get_boundary_key(plid_vec);

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x_arr[3] = {x0, x1, x2};

        // get y values of points
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y_arr[3] = {y0, y1, y2};

        // get jacobian
        // evaluate using points along 1D edge (X; A)

        // transform 2D edge coordinates (x, y) into 1D (X)
        // get (x, y) along edge
        double xa = x_arr[pa_plid];
        double xb = x_arr[pb_plid];
        double ya = y_arr[pa_plid];
        double yb = y_arr[pb_plid];

        // let Xa = 0. Xb will be the distance from a to b
        double Xa = 0.;
        double Xb = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya));

        // get derivatives of X with respect to A
        double derivative_X_A = 0.5*(Xb - Xa);

        // get jacobian and its inverse and determinant
        double jacobian_inverse = 1./derivative_X_A;
        double jacobian_determinant = derivative_X_A;

        // determine integration points used in gaussian integration

        // get 2D integration points
        VectorDouble a_vec;
        VectorDouble b_vec;
        switch (boundary_key)
        {
            case 01: a_vec = sample_vec; b_vec = zero_vec; break;  // bottom edge (0-1)
            case 12: a_vec = zero_vec; b_vec = sample_vec; break;  // left edge (1-2)
            case 02: a_vec = sample_vec; b_vec = sample_vec; break;  // hypotenouse (2-0)
        }

        // initialize
        VectorDouble jacobian_determinant_part_l_vec;
        VectorDouble2D N_part_l_vec;

        // iterate for each integration point
        for (int indx_l = 0; indx_l < 2; indx_l++)
        {

            // initialize
            VectorDouble N_part_li_vec;

            // get test function values
            // evaluate using points in 2D element (x, y; a, b)

            // get a and b values where function is evaluated
            double a = a_vec[indx_l];
            double b = b_vec[indx_l];

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

                // store in vectors
                N_part_li_vec.push_back(N);

            }

            // store in vectors
            N_part_l_vec.push_back(N_part_li_vec);

        }

        // store in vectors
        jacobian_determinant_vec[edid][boundary_key] = jacobian_determinant;
        Ni_vec[edid][boundary_key] = N_part_l_vec;
        
        // determine outward unit normal vector

        // get vector along boundary pointing counterclocwise
        // start - stt; end - end
        double x_start = 0.; double y_start = 0.;
        double x_end = 0.; double y_end = 0.;    
        switch (boundary_key)
        {
            case 01:  // bottom edge (0-1)
                x_start = x0; y_start = y0;
                x_end = x1; y_end = y1;
            break;  
            case 12:  // left edge (1-2)
                x_start = x1; y_start = y1;
                x_end = x2; y_end = y2;
            break;  
            case 02:  // hypotenouse (2-0)
                x_start = x2; y_start = y2;
                x_end = x0; y_end = y0;
            break;  
        }

        // calculate normal to boundary vector
        // should be pointing out of element
        double boundarynormal_x = y_end - y_start;
        double boundarynormal_y = x_start - x_end;

        // calculate unit normal vector components
        double boundarynormal_mag = sqrt(boundarynormal_x*boundarynormal_x + boundarynormal_y*boundarynormal_y);
        normal_x_vec[edid][boundary_key] = boundarynormal_x/boundarynormal_mag;
        normal_y_vec[edid][boundary_key] = boundarynormal_y/boundarynormal_mag;

    }

}

void BoundaryIntegralUnit::evaluate_Ni_quad4()
{

    // initialize 2D integration points (within element)
    const double M_1_SQRT_3 = 1./sqrt(3);
    VectorDouble sample_vec = {-M_1_SQRT_3, +M_1_SQRT_3};
    VectorDouble positive_vec = {+1., +1.};
    VectorDouble negative_vec = {-1., -1.};

    // initialize 1D integration points (along edge)
    // use capital letters (e.g., X, A) to denote variables along an edge
    VectorDouble A_vec = {-M_1_SQRT_3, +M_1_SQRT_3};

    // iterate through each boundary
    for (int bid = 0; bid < boundary_ptr->num_boundary; bid++)
    {

        // get the element ID
        int egid = boundary_ptr->boundary_egid_vec[bid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get points around element
        int p0_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][0];
        int p1_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][1];
        int p2_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][2];
        int p3_pgid = domain_ptr->element_edid_plid_to_pgid_vec[edid][3];
        int p0_pdid = domain_ptr->point_pgid_to_pdid_map[p0_pgid];
        int p1_pdid = domain_ptr->point_pgid_to_pdid_map[p1_pgid];
        int p2_pdid = domain_ptr->point_pgid_to_pdid_map[p2_pgid];
        int p3_pdid = domain_ptr->point_pgid_to_pdid_map[p3_pgid];

        // get boundary key
        // use a symmetric pairing function
        // this indicates which edge the boundary is applied on
        int pa_plid = boundary_ptr->boundary_blid_to_plid_vec[bid][0];
        int pb_plid = boundary_ptr->boundary_blid_to_plid_vec[bid][1];
        VectorInt plid_vec = boundary_ptr->boundary_blid_to_plid_vec[bid];
        int boundary_key = get_boundary_key(plid_vec);

        // get x values of points
        double x0 = domain_ptr->point_position_x_vec[p0_pdid];
        double x1 = domain_ptr->point_position_x_vec[p1_pdid];
        double x2 = domain_ptr->point_position_x_vec[p2_pdid];
        double x3 = domain_ptr->point_position_x_vec[p3_pdid];
        double x_arr[4] = {x0, x1, x2, x3};

        // get y values of points
        double y0 = domain_ptr->point_position_y_vec[p0_pdid];
        double y1 = domain_ptr->point_position_y_vec[p1_pdid];
        double y2 = domain_ptr->point_position_y_vec[p2_pdid];
        double y3 = domain_ptr->point_position_y_vec[p3_pdid];
        double y_arr[4] = {y0, y1, y2, y3};

        // get jacobian
        // evaluate using points along 1D edge (X; A)

        // transform 2D edge coordinates (x, y) into 1D (X)
        // get (x, y) along edge
        double xa = x_arr[pa_plid];
        double xb = x_arr[pb_plid];
        double ya = y_arr[pa_plid];
        double yb = y_arr[pb_plid];

        // let Xa = 0. Xb will be the distance from a to b
        double Xa = 0.;
        double Xb = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya));

        // get derivatives of X with respect to A
        double derivative_X_A = 0.5*(Xb - Xa);

        // get jacobian and its inverse and determinant
        double jacobian_inverse = 1./derivative_X_A;
        double jacobian_determinant = derivative_X_A;

        // determine integration points used in gaussian integration

        // get 2D integration points
        VectorDouble a_vec;
        VectorDouble b_vec;
        switch (boundary_key)
        {
            case 01: a_vec = negative_vec; b_vec = sample_vec; break;  // left edge (0-1)
            case 12: a_vec = sample_vec; b_vec = positive_vec; break;  // top edge (1-2)
            case 23: a_vec = positive_vec; b_vec = sample_vec; break;  // right edge (2-3)
            case 03: a_vec = sample_vec; b_vec = negative_vec; break;  // bottom edge (3-0)
        }

        // initialize
        VectorDouble jacobian_determinant_part_l_vec;
        VectorDouble2D N_part_l_vec;

        // iterate for each integration point
        for (int indx_l = 0; indx_l < 2; indx_l++)
        {

            // initialize
            VectorDouble N_part_li_vec;

            // get test function values
            // evaluate using points in 2D element (x, y; a, b)

            // get a and b values where function is evaluated
            double a = a_vec[indx_l];
            double b = b_vec[indx_l];

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

                // store in vectors
                N_part_li_vec.push_back(N);

            }

            // store in vectors
            N_part_l_vec.push_back(N_part_li_vec);

        }

        // store in vectors
        jacobian_determinant_vec[edid][boundary_key] = jacobian_determinant;
        Ni_vec[edid][boundary_key] = N_part_l_vec;
        
        // determine outward unit normal vector

        // get vector along boundary pointing counterclocwise
        double x_start = 0.; double y_start = 0.;
        double x_end = 0.; double y_end = 0.;    
        switch (boundary_key)
        {
            case 01:  // left edge (0-1)
                x_start = x0; y_start = y0;
                x_end = x1; y_end = y1;
            break;  
            case 12:  // top edge (1-2)
                x_start = x1; y_start = y1;
                x_end = x2; y_end = y2;
            break;  
            case 23:  // right edge (2-3)
                x_start = x2; y_start = y2;
                x_end = x3; y_end = y3;
            break;
            case 03:  // bottom edge (3-0)
                x_start = x3; y_start = y3;
                x_end = x0; y_end = y0;
            break; 
        }

        // calculate normal to boundary vector
        // should be pointing out of element
        double boundarynormal_x = y_end - y_start;
        double boundarynormal_y = x_start - x_end;

        // calculate unit normal vector components
        double boundarynormal_mag = sqrt(boundarynormal_x*boundarynormal_x + boundarynormal_y*boundarynormal_y);
        normal_x_vec[edid][boundary_key] = boundarynormal_x/boundarynormal_mag;
        normal_y_vec[edid][boundary_key] = boundarynormal_y/boundarynormal_mag;

    }

}

}

#endif
