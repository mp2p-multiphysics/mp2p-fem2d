#ifndef PHYSICSSTEADY_NAVIERSTOKES
#define PHYSICSSTEADY_NAVIERSTOKES
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_0d.hpp"
#include "domain_1d.hpp"
#include "domain_2d.hpp"
#include "integral_taylorhood_2d.hpp"
#include "physicssteady_base.hpp"
#include "scalar_0d.hpp"
#include "scalar_1d.hpp"
#include "scalar_2d.hpp"
#include "variable_2d.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class PhysicsSteadyNavierStokes : public PhysicsSteadyBase
{
    /*

    Steady-state Navier-Stokes equation (NSE) and continuity equation.
    
    NSE x-component: rho * (u * du/dx + v * du/dy) = -dp/dx + mu * (d^2(u)/dx^2 + d^2(u)/dy^2) + rho * g_x
    NSE y-component: rho * (u * dv/dx + v * dv/dy) = -dp/dy + mu * (d^2(v)/dx^2 + d^2(v)/dy^2) + rho * g_y
    Continuity: du/dx + dv/dy = 0

    Variables
    =========
    velocity_x_in : VariableGroup
        u in the NSE and continuity equation.
        Must be defined over a quadratic domain.
    velocity_y_in : VariableGroup
        v in the NSE and continuity equation.
        Must be defined over a quadratic domain.
    pressure_in : VariableGroup
        p in the NSE.
        Must be defined over a linear domain.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_variablegroup : void
        Set variables used in this physics.
    set_domain : void
        Set scalars applied to 2D domains.
    set_boundary_velocity : void
        Set a velocity boundary condition along a 1D domain.
    set_boundary_pressure : void
        Set a pressure boundary condition along a 1D domain.

    */

    public:

    // variables
    VariableGroup* velocity_x_ptr;
    VariableGroup* velocity_y_ptr;
    VariableGroup* pressure_ptr;

    // vectors indicating dirichlet BCs
    // [pfid] -> true if dirichlet bc
    std::vector<bool> is_velocity_x_dirichlet_vec;  
    std::vector<bool> is_velocity_y_dirichlet_vec;
    std::vector<bool> is_pressure_dirichlet_vec;

    // domain objects
    std::vector<Domain2D*> domain_line_ptr_vec;
    std::vector<Domain2D*> domain_quad_ptr_vec;
    std::vector<IntegralTaylorHood2D*> integral_ptr_vec;
    std::vector<Scalar2D*> density_ptr_vec;
    std::vector<Scalar2D*> viscosity_ptr_vec;
    std::vector<Scalar2D*> force_x_ptr_vec;
    std::vector<Scalar2D*> force_y_ptr_vec;

    // boundary objects - velocity
    std::vector<Domain1D*> velocity_domain_ptr_vec;
    std::vector<Scalar1D*> velocity_x_ptr_vec;
    std::vector<Scalar1D*> velocity_y_ptr_vec;

    // boundary objects - pressure
    std::vector<Domain1D*> pressure_domain_ptr_vec;
    std::vector<Scalar1D*> pressure_ptr_vec;

    // boundary objects - pressure point
    std::vector<Domain0D*> pressure_point_domain_ptr_vec;
    std::vector<Scalar0D*> pressure_point_ptr_vec;

    // vectors of objects to update
    std::vector<Scalar1D*> scalar1d_ptr_vec;
    std::vector<Scalar2D*> scalar2d_ptr_vec;
    std::vector<VariableGroup*> variablegroup_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec);
    void set_domain(Domain2D &domain_line_in, Domain2D &domain_quad_in, IntegralTaylorHood2D &integral_in, Scalar2D &density_in, Scalar2D &viscosity_in, Scalar2D &force_x_in, Scalar2D &force_y_in);
    void set_boundary_velocity(Domain1D &domain_quad_in, Scalar1D &velocity_x_in, Scalar1D &velocity_y_in);
    void set_boundary_pressure(Domain1D &domain_quad_in, Scalar1D &pressure_in); 
    void set_boundary_pressure(Domain0D &domain_in, Scalar0D &pressure_in); 

    // getter and setter functions
    void set_start_row(int start_row_in) {start_row = start_row_in;}
    int get_start_row() {return start_row;}
    std::vector<Scalar1D*> get_scalar1d_ptr_vec() {return scalar1d_ptr_vec;}
    std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsSteadyNavierStokes() {}

    // default constructor
    PhysicsSteadyNavierStokes(VariableGroup &velocity_x_in, VariableGroup &velocity_y_in, VariableGroup &pressure_in)
    {

        // set variable groups
        velocity_x_ptr = &velocity_x_in;
        velocity_y_ptr = &velocity_y_in;
        pressure_ptr = &pressure_in;

        // add to vector of variable groups
        variablegroup_ptr_vec.push_back(&velocity_x_in);
        variablegroup_ptr_vec.push_back(&velocity_y_in);
        variablegroup_ptr_vec.push_back(&pressure_in);

        // vector indicating if dirichlet bc is applied
        is_velocity_x_dirichlet_vec = std::vector<bool> (velocity_x_in.num_point, false);
        is_velocity_y_dirichlet_vec = std::vector<bool> (velocity_y_in.num_point, false);
        is_pressure_dirichlet_vec = std::vector<bool> (pressure_in.num_point, false);

    }

    private:

    void matrix_fill_domain
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain2D *domain_line_ptr, Domain2D *domain_quad_ptr, IntegralTaylorHood2D *integral_ptr,
        Scalar2D *density_ptr, Scalar2D *viscosity_ptr, Scalar2D *force_x_ptr, Scalar2D *force_y_ptr
    );
    void matrix_fill_velocity
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain1D *domain_ptr,
        Scalar1D *velocity_x_ptr, Scalar1D *velocity_y_ptr
    );
    void matrix_fill_pressure
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain1D *domain_ptr,
        Scalar1D *pressure_ptr
    );
    void matrix_fill_pressure_point
    (
        EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
        Domain0D *domain_ptr,
        Scalar0D *pressure_ptr
    );

};

void PhysicsSteadyNavierStokes::set_domain(Domain2D &domain_line_in, Domain2D &domain_quad_in, IntegralTaylorHood2D &integral_in, Scalar2D &density_in, Scalar2D &viscosity_in, Scalar2D &force_x_in, Scalar2D &force_y_in)
{
    /*
    
    Set scalars applied to 2D domains.

    Arguments
    =========
    domain_line_in : Domain2D
        Linear domain that this physics applies to.
    domain_quad_in : Domain2D
        Quadratic domain that this physics applies to.
    integral_in : IntegralTaylorHood2D
        Test function integrals over the domains.
    density_in : Scalar2D
        rho in NSE.
    viscosity_in : Scalar2D
        mu in NSE.
    force_x_in : Scalar2D
        g_x in NSE.
    force_y_in : Scalar2D
        g_y in NSE.

    Returns
    =======
    (none)

    Notes
    =====
    Coinciding points in the linear and quadratic domains must have the same IDs.

    */

    // add to vector of domain objects
    domain_line_ptr_vec.push_back(&domain_line_in);
    domain_quad_ptr_vec.push_back(&domain_quad_in);
    integral_ptr_vec.push_back(&integral_in);
    density_ptr_vec.push_back(&density_in);
    viscosity_ptr_vec.push_back(&viscosity_in);
    force_x_ptr_vec.push_back(&force_x_in);
    force_y_ptr_vec.push_back(&force_y_in);

    // add to vector of scalar2d objects
    scalar2d_ptr_vec.push_back(&density_in);
    scalar2d_ptr_vec.push_back(&viscosity_in);
    scalar2d_ptr_vec.push_back(&force_x_in);
    scalar2d_ptr_vec.push_back(&force_y_in);

    // calculate linear domain integrals
    integral_in.evaluate_integral_Ni_derivative_Mj_x_line();
    integral_in.evaluate_integral_Ni_derivative_Mj_y_line();

    // calculate quadratic domain integrals
    integral_in.evaluate_integral_Mi_quad();
    integral_in.evaluate_integral_Mi_derivative_Nj_x_quad();
    integral_in.evaluate_integral_Mi_derivative_Nj_y_quad();
    integral_in.evaluate_integral_div_Mi_dot_div_Mj_quad();
    integral_in.evaluate_integral_Mi_Mk_derivative_Mj_x_quad();
    integral_in.evaluate_integral_Mi_Mk_derivative_Mj_y_quad();

}

void PhysicsSteadyNavierStokes::set_boundary_velocity(Domain1D &domain_quad_in, Scalar1D &velocity_x_in, Scalar1D &velocity_y_in)
{
    /*
    
    Set a velocity boundary condition along a 1D domain.

    Arguments
    =========
    domain_quad_in : Domain1D
        Quadratic domain that this boundary condition applies to.
    velocity_x_in : Scalar1D
        x velocity prescribed by the boundary condition.
    velocity_y_in : Scalar1D
        y velocity prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    velocity_domain_ptr_vec.push_back(&domain_quad_in);
    velocity_x_ptr_vec.push_back(&velocity_x_in);
    velocity_y_ptr_vec.push_back(&velocity_y_in);

    // add to vector of scalar1d objects
    scalar1d_ptr_vec.push_back(&velocity_x_in);
    scalar1d_ptr_vec.push_back(&velocity_y_in);

    // specify points with dirichlet BC
    for (auto pgid : domain_quad_in.point_pdid_to_pgid_vec)
    {
        int velx_pfid = velocity_x_ptr->point_pgid_to_pfid_map[pgid];
        int vely_pfid = velocity_y_ptr->point_pgid_to_pfid_map[pgid];
        is_velocity_x_dirichlet_vec[velx_pfid] = true;
        is_velocity_y_dirichlet_vec[vely_pfid] = true;
    }

}

void PhysicsSteadyNavierStokes::set_boundary_pressure(Domain1D &domain_line_in, Scalar1D &pressure_in)
{
    /*
    
    Set a pressure boundary condition along a 1D domain.

    Arguments
    =========
    domain_line_in : Domain1D
        Linear domain that this boundary condition applies to.
    pressure_in : Scalar1D
        Pressure prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    pressure_domain_ptr_vec.push_back(&domain_line_in);
    pressure_ptr_vec.push_back(&pressure_in);

    // add to vector of scalar1d objects
    scalar1d_ptr_vec.push_back(&pressure_in);

    // specify points with dirichlet BC
    for (auto pgid : domain_line_in.point_pdid_to_pgid_vec)
    {
        int pres_pfid = pressure_ptr->point_pgid_to_pfid_map[pgid];
        is_pressure_dirichlet_vec[pres_pfid] = true;
    }

}

void PhysicsSteadyNavierStokes::set_boundary_pressure(Domain0D &domain_in, Scalar0D &pressure_in)
{
    /*
    
    Set a pressure boundary condition along a 0D domain.

    Arguments
    =========
    domain_in : Domain1D
        Point domain that this boundary condition applies to.
    pressure_in : Scalar1D
        Pressure prescribed by the boundary condition.

    Returns
    =======
    (none)

    */

    // add to vector of boundary objects
    pressure_point_domain_ptr_vec.push_back(&domain_in);
    pressure_point_ptr_vec.push_back(&pressure_in);

    // specify points with dirichlet BC
    for (auto pgid : domain_in.point_pdid_to_pgid_vec)
    {
        int pres_pfid = pressure_ptr->point_pgid_to_pfid_map[pgid];
        is_pressure_dirichlet_vec[pres_pfid] = true;
    }

}

void PhysicsSteadyNavierStokes::matrix_fill
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec
)
{
    /*

    Fill up the matrix equation Ax = b with entries as dictated by the physics. 

    Arguments
    =========
    a_trivec : EigenTripletVector
        A in Ax = b.
    b_vec : EigenVector
        b in Ax = b.
    x_vec : EigenVector
        x in Ax = b.

    Returns
    =======
    (none)

    */

    // iterate through each domain
    for (int indx_d = 0; indx_d < domain_line_ptr_vec.size(); indx_d++)
    {

        // subset domain objects
        Domain2D *domain_line_ptr = domain_line_ptr_vec[indx_d];
        Domain2D *domain_quad_ptr = domain_quad_ptr_vec[indx_d];
        IntegralTaylorHood2D *integral_ptr = integral_ptr_vec[indx_d];
        Scalar2D *density_ptr = density_ptr_vec[indx_d];
        Scalar2D *viscosity_ptr = viscosity_ptr_vec[indx_d];
        Scalar2D *force_x_ptr = force_x_ptr_vec[indx_d];
        Scalar2D *force_y_ptr = force_y_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(a_trivec, b_vec, x_vec, domain_line_ptr, domain_quad_ptr, integral_ptr, density_ptr, viscosity_ptr, force_x_ptr, force_y_ptr);

    }

    // iterate through each dirichlet boundary
    for (int indx_d = 0; indx_d < velocity_domain_ptr_vec.size(); indx_d++)
    {
        Domain1D *domain_ptr = velocity_domain_ptr_vec[indx_d];
        Scalar1D *velocity_x_ptr = velocity_x_ptr_vec[indx_d];
        Scalar1D *velocity_y_ptr = velocity_y_ptr_vec[indx_d];
        matrix_fill_velocity(a_trivec, b_vec, x_vec, domain_ptr, velocity_x_ptr, velocity_y_ptr);
    }
    for (int indx_d = 0; indx_d < pressure_domain_ptr_vec.size(); indx_d++)
    {
        Domain1D *domain_ptr = pressure_domain_ptr_vec[indx_d];
        Scalar1D *pressure_ptr = pressure_ptr_vec[indx_d];
        matrix_fill_pressure(a_trivec, b_vec, x_vec, domain_ptr, pressure_ptr);
    }
    for (int indx_d = 0; indx_d < pressure_point_domain_ptr_vec.size(); indx_d++)
    {
        Domain0D *domain_ptr = pressure_point_domain_ptr_vec[indx_d];
        Scalar0D *pressure_ptr = pressure_point_ptr_vec[indx_d];
        matrix_fill_pressure_point(a_trivec, b_vec, x_vec, domain_ptr, pressure_ptr);
    }

}

void PhysicsSteadyNavierStokes::matrix_fill_domain
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
    Domain2D *domain_line_ptr, Domain2D *domain_quad_ptr, IntegralTaylorHood2D *integral_ptr,
    Scalar2D *density_ptr, Scalar2D *viscosity_ptr, Scalar2D *force_x_ptr, Scalar2D *force_y_ptr
)
{

    // offset the starting row depending on the equation
    int offset_nsex = 0;
    int offset_nsey = offset_nsex + velocity_x_ptr->num_point;
    int offset_cont = offset_nsey + velocity_y_ptr->num_point;

    // fill up matrix with x-component of NSE
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble den_vec = density_ptr->get_neighbor_value(edid);
        VectorDouble visc_vec = viscosity_ptr->get_neighbor_value(edid);
        VectorDouble fcex_vec = force_x_ptr->get_neighbor_value(edid);

        // get velocities of points around element
        VectorDouble velx_vec = velocity_x_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble vely_vec = velocity_y_ptr->get_neighbor_value(domain_quad_ptr, edid);

        // get group ID of points
        int egid = domain_quad_ptr->element_edid_to_egid_vec[edid];
        int edid_line = domain_line_ptr->element_egid_to_edid_map[egid];  // local element ID of linear domain
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_quad_ptr, edid);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid_line);

        // matrix indexing
        // matrix row = start_row of test function (physics) + offset + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // iterate through test functions
        // associated with matrix row
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {

            // skip if matrix row has dirichlet condition
            if (is_velocity_x_dirichlet_vec[velx_pfid_vec[indx_i]])
            {
                continue;
            }

            // calculate matrix row
            int mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];

            // iterate through velocity x trial functions
            // associated with matrix column
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = velocity_x_ptr->start_col + velx_pfid_vec[indx_j];

                // calculate velocity integrals
                double integ_Mi_velx_derv_Mj_x = 0;
                double integ_Mi_vely_derv_Mj_y = 0;
                for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++)
                {
                    integ_Mi_velx_derv_Mj_x += velx_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_x_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_vely_derv_Mj_y += vely_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_y_quad_vec[edid][indx_i][indx_j][indx_k];
                }
                
                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col,
                    den_vec[indx_i] * (integ_Mi_velx_derv_Mj_x + integ_Mi_vely_derv_Mj_y) +
                    visc_vec[indx_i] * integral_ptr->integral_div_Mi_dot_div_Mj_quad_vec[edid][indx_i][indx_j]
                ));

            }

            // iterate through pressure trial functions
            for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_j];

                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Mi_derivative_Nj_x_quad_vec[edid][indx_i][indx_j]));
            
            }

            // append to b_vec
            b_vec.coeffRef(mat_row) += den_vec[indx_i] * fcex_vec[indx_i] * integral_ptr->integral_Mi_quad_vec[edid][indx_i];

        }

    }

    // fill up matrix with y-component of NSE
    for (int edid = 0; edid < domain_quad_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble den_vec = density_ptr->get_neighbor_value(edid);
        VectorDouble visc_vec = viscosity_ptr->get_neighbor_value(edid);
        VectorDouble fcey_vec = force_y_ptr->get_neighbor_value(edid);

        // get velocities of points around element
        VectorDouble velx_vec = velocity_x_ptr->get_neighbor_value(domain_quad_ptr, edid);
        VectorDouble vely_vec = velocity_y_ptr->get_neighbor_value(domain_quad_ptr, edid);

        // get group ID of points
        int egid = domain_quad_ptr->element_edid_to_egid_vec[edid];
        int edid_line = domain_line_ptr->element_egid_to_edid_map[egid];  // local element ID of linear domain
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_quad_ptr, edid);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid_line);

        // iterate through test functions
        for (int indx_i = 0; indx_i < domain_quad_ptr->num_neighbor; indx_i++)
        {

            // skip if matrix row has dirichlet condition
            if (is_velocity_y_dirichlet_vec[vely_pfid_vec[indx_i]])
            {
                continue;
            }

            // calculate matrix row
            int mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];

            // iterate through velocity y trial functions
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = velocity_y_ptr->start_col + vely_pfid_vec[indx_j];

                // calculate velocity integrals
                double integ_Mi_velx_derv_Mj_x = 0;
                double integ_Mi_vely_derv_Mj_y = 0;
                for (int indx_k = 0; indx_k < domain_quad_ptr->num_neighbor; indx_k++)
                {
                    integ_Mi_velx_derv_Mj_x += velx_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_x_quad_vec[edid][indx_i][indx_j][indx_k];
                    integ_Mi_vely_derv_Mj_y += vely_vec[indx_k] * integral_ptr->integral_Mi_Mk_derivative_Mj_y_quad_vec[edid][indx_i][indx_j][indx_k];
                }
                
                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col,
                    den_vec[indx_i] * (integ_Mi_velx_derv_Mj_x + integ_Mi_vely_derv_Mj_y) +
                    visc_vec[indx_i] * integral_ptr->integral_div_Mi_dot_div_Mj_quad_vec[edid][indx_i][indx_j]
                ));

            }

            // iterate through pressure trial functions
            for (int indx_j = 0; indx_j < domain_line_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_j];

                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Mi_derivative_Nj_y_quad_vec[edid][indx_i][indx_j]));

            }

            // append to b_vec
            b_vec.coeffRef(mat_row) += den_vec[indx_i] * fcey_vec[indx_i] * integral_ptr->integral_Mi_quad_vec[edid][indx_i];

        }

    }

    // fill up matrix with continuity equation
    for (int edid = 0; edid < domain_line_ptr->num_element; edid++)
    {

        // get group ID of points
        int egid = domain_line_ptr->element_edid_to_egid_vec[edid];
        int edid_quad = domain_quad_ptr->element_egid_to_edid_map[egid];  // local element ID of quadratic domain
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_quad_ptr, edid_quad);
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_quad_ptr, edid_quad);
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_line_ptr, edid);

        // iterate through test functions
        for (int indx_i = 0; indx_i < domain_line_ptr->num_neighbor; indx_i++)
        {
            
            // skip if matrix row has dirichlet condition
            if (is_pressure_dirichlet_vec[pres_pfid_vec[indx_i]])
            {
                continue;
            }
            
            // calculate matrix row
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];

            // iterate through velocity x trial functions
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = velocity_x_ptr->start_col + velx_pfid_vec[indx_j];

                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Ni_derivative_Mj_x_line_vec[edid][indx_i][indx_j]));
            
            }

            // iterate through velocity y trial functions
            for (int indx_j = 0; indx_j < domain_quad_ptr->num_neighbor; indx_j++)
            {
                
                // calculate matrix column
                int mat_col = velocity_y_ptr->start_col + vely_pfid_vec[indx_j];

                // append to a_trivec
                a_trivec.push_back(EigenTriplet(mat_row, mat_col, integral_ptr->integral_Ni_derivative_Mj_y_line_vec[edid][indx_i][indx_j]));
            
            }

        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_velocity
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
    Domain1D *domain_ptr,
    Scalar1D *velx_ptr, Scalar1D *vely_ptr
)
{

    // offset the starting row depending on the equation
    int offset_nsex = 0;
    int offset_nsey = offset_nsex + velocity_x_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble velx_vec = velx_ptr->get_neighbor_value(edid);
        VectorDouble vely_vec = vely_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt velx_pfid_vec = velocity_x_ptr->get_neighbor_pfid(domain_ptr, edid);
        VectorInt vely_pfid_vec = velocity_y_ptr->get_neighbor_pfid(domain_ptr, edid);

        // apply velocity x dirichlet BC
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsex + velx_pfid_vec[indx_i];
            int mat_col = velocity_x_ptr->start_col + velx_pfid_vec[indx_i];
            a_trivec.push_back(EigenTriplet(mat_row, mat_col, 1.));
            b_vec.coeffRef(mat_row) += velx_vec[indx_i];
        }

        // apply velocity y dirichlet BC
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_nsey + vely_pfid_vec[indx_i];
            int mat_col = velocity_y_ptr->start_col + vely_pfid_vec[indx_i];
            a_trivec.push_back(EigenTriplet(mat_row, mat_col, 1.));
            b_vec.coeffRef(mat_row) += vely_vec[indx_i];
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_pressure
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
    Domain1D *domain_ptr,
    Scalar1D *pres_ptr
)
{

    // offset the starting row depending on the equation
    int offset_cont = velocity_x_ptr->num_point + velocity_y_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble pres_vec = pres_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_ptr, edid);

        // apply pressure dirichlet BC
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
            int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_i];
            a_trivec.push_back(EigenTriplet(mat_row, mat_col, 1.));
            b_vec.coeffRef(mat_row) += pres_vec[indx_i];
        }

    }

}

void PhysicsSteadyNavierStokes::matrix_fill_pressure_point
(
    EigenTripletVector &a_trivec, EigenVector &b_vec, EigenVector &x_vec,
    Domain0D *domain_ptr,
    Scalar0D *pres_ptr
)
{

    // offset the starting row depending on the equation
    int offset_cont = velocity_x_ptr->num_point + velocity_y_ptr->num_point;

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble pres_vec = pres_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pres_pfid_vec = pressure_ptr->get_neighbor_pfid(domain_ptr, edid);

        // apply pressure dirichlet BC
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + offset_cont + pres_pfid_vec[indx_i];
            int mat_col = pressure_ptr->start_col + pres_pfid_vec[indx_i];
            a_trivec.push_back(EigenTriplet(mat_row, mat_col, 1.));
            b_vec.coeffRef(mat_row) += pres_vec[indx_i];
        }

    }

}

}

#endif
