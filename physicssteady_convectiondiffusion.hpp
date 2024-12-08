#ifndef PHYSICSSTEADY_CONVECTIONDIFFUSION
#define PHYSICSSTEADY_CONVECTIONDIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "boundary_group.hpp"
#include "boundaryintegral_group.hpp"
#include "container_typedef.hpp"
#include "domain_group.hpp"
#include "domainintegral_group.hpp"
#include "physicssteady_base.hpp"
#include "scalar_group.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class PhysicsSteadyConvectionDiffusion : public PhysicsSteadyBase
{
    /*

    Single-component steady-state convection-diffusion equation.    
    
    0 = -div(-b * grad(u)) - v * grad(u) + c

    Variables
    =========
    domain_group_in : DomainGroup
        Domains where this physics is applied to.
    domainintegral_group_in : DomainIntegralGroup
        Test function integrals over the domains.
    value_group_in : VariableGroup
        u in 0 = -div(-b * grad(u)) + c.
        This will be solved for by the matrix equation.
    diffusioncoefficient_group_in : ScalarGroup
        b in 0 = -div(-b * grad(u)) + c.
    velocity_x_field_in : ScalarFieldGroup
        x-component of v in -div(-b * grad(u)) - v * grad(u) + c.
    velocity_x_field_in : ScalarFieldGroup
        y-component of v in -div(-b * grad(u)) - v * grad(u) + c.
    generationcoefficient_group_in : ScalarGroup
        c in 0 = -div(-b * grad(u)) + c.
    boundary_group_in : BoundaryGroup
        Boundary conditions applied on u.
    boundaryintegral_group_in : BoundaryIntegralGroup
        Test function integrals over the boundaries.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_start_row : void
        Sets the starting row in A and b where entries are filled up.
    get_start_row : int
        Returns the starting row.
    get_boundary_group_ptr_vec : vector<BoundaryGroup*>
        Returns the vector containing pointers to BoundaryGroup objects tied to this physics.
    get_scalar_group_ptr_vec : vector<ScalarGroup*>
        Returns the vector containing pointers to ScalarGroup objects tied to this physics.
    get_variable_group_ptr_vec : vector<VariableGroup*>
        Returns the vector containing pointers to VariableGroup objects tied to this physics.

    */

    public:

    // variables
    DomainGroup *domain_group_ptr;
    DomainIntegralGroup *domainintegral_group_ptr;
    VariableGroup *value_group_ptr;
    ScalarGroup *diffusioncoefficient_group_ptr;
    ScalarGroup *velocity_x_group_ptr;
    ScalarGroup *velocity_y_group_ptr;
    ScalarGroup *generationcoefficient_group_ptr;
    BoundaryGroup *boundary_group_ptr;
    BoundaryIntegralGroup *boundaryintegral_group_ptr;

    // vector of scalar and variable groups
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;
    std::vector<BoundaryGroup*> boundary_group_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);
    void set_start_row(int start_row_in);
    int get_start_row();
    std::vector<ScalarGroup*> get_scalar_group_ptr_vec();
    std::vector<VariableGroup*> get_variable_group_ptr_vec();
    std::vector<BoundaryGroup*> get_boundary_group_ptr_vec();

    // default constructor
    PhysicsSteadyConvectionDiffusion() {}

    // constructor
    PhysicsSteadyConvectionDiffusion
    (
        DomainGroup &domain_group_in, DomainIntegralGroup &domainintegral_group_in,
        VariableGroup &value_group_in,
        ScalarGroup &diffusioncoefficient_group_in, ScalarGroup &velocity_x_group_in, ScalarGroup &velocity_y_group_in, ScalarGroup &generationcoefficient_group_in,
        BoundaryGroup &boundary_group_in, BoundaryIntegralGroup &boundaryintegral_group_in
    )
    {
        
        // store variables
        domain_group_ptr = &domain_group_in;
        domainintegral_group_ptr = &domainintegral_group_in;
        value_group_ptr = &value_group_in;
        diffusioncoefficient_group_ptr = &diffusioncoefficient_group_in;
        velocity_x_group_ptr = &velocity_x_group_in;
        velocity_y_group_ptr = &velocity_y_group_in;
        generationcoefficient_group_ptr = &generationcoefficient_group_in;
        boundary_group_ptr = &boundary_group_in;
        boundaryintegral_group_ptr = &boundaryintegral_group_in;

        // set boundary conditions
        boundary_group_ptr->set_boundary_type({0}, {1, 2});

        // vector of scalar and variable groups 
        scalar_group_ptr_vec = {diffusioncoefficient_group_ptr, velocity_x_group_ptr, velocity_y_group_ptr, generationcoefficient_group_ptr};
        variable_group_ptr_vec = {value_group_ptr};
        boundary_group_ptr_vec = {boundary_group_ptr};

        // calculate integrals
        domainintegral_group_ptr->evaluate_integral_div_Ni_dot_div_Nj();
        domainintegral_group_ptr->evaluate_integral_Ni_derivative_Nj_x();
        domainintegral_group_ptr->evaluate_integral_Ni_derivative_Nj_y();
        domainintegral_group_ptr->evaluate_integral_Ni();
        boundaryintegral_group_ptr->evaluate_integral_Ni();
        boundaryintegral_group_ptr->evaluate_integral_Ni_Nj();

    }

    private:
    void matrix_fill_domain
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainUnit *domain_ptr, DomainIntegralUnit *domainintegral_ptr,
        ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *velocity_x_ptr, ScalarUnit *velocity_y_ptr, ScalarUnit *generationcoefficient_ptr
    );
    void matrix_fill_natural
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainUnit *domain_ptr,
        ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *velocity_x_ptr, ScalarUnit *velocity_y_ptr, ScalarUnit *generationcoefficient_ptr,
        BoundaryUnit *boundary_ptr, BoundaryIntegralUnit *boundaryintegral_ptr
    );
    void matrix_fill_essential_clear
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainUnit *domain_ptr,
        BoundaryUnit *boundary_ptr
    );
    void matrix_fill_essential
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        DomainUnit *domain_ptr,
        ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *velocity_x_ptr, ScalarUnit *velocity_y_ptr, ScalarUnit *generationcoefficient_ptr,
        BoundaryUnit *boundary_ptr, BoundaryIntegralUnit *boundaryintegral_ptr
    );

};

void PhysicsSteadyConvectionDiffusion::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec
)
{
    /*

    Fill up the matrix equation Ax = b with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : Eigen::SparseMatrix<double>
        A in Ax = b.
    b_vec : Eigen::VectorXd
        b in Ax = b.
    x_vec : Eigen::VectorXd
        x in Ax = b.

    Returns
    =======
    (none)

    */

   // iterate through each domain
   for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
   {

        // subset the domain and integrals
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        DomainIntegralUnit *domainintegral_ptr = domainintegral_group_ptr->integral_ptr_vec[indx_d];

        // subset the scalars
        ScalarUnit *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *velocity_x_ptr = velocity_x_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *velocity_y_ptr = velocity_y_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(a_mat, b_vec, x_vec, domain_ptr, domainintegral_ptr, diffusioncoefficient_ptr, velocity_x_ptr, velocity_y_ptr, generationcoefficient_ptr);

   }

   // iterate through each natural boundary
   for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
   {

        // subset the boundary and integrals
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        BoundaryUnit *boundary_ptr = boundary_group_ptr->boundary_ptr_vec[indx_d];
        BoundaryIntegralUnit *boundaryintegral_ptr = boundaryintegral_group_ptr->integral_ptr_vec[indx_d];

        // subset the scalars
        ScalarUnit *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *velocity_x_ptr = velocity_x_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *velocity_y_ptr = velocity_y_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_natural(a_mat, b_vec, x_vec, domain_ptr, diffusioncoefficient_ptr, velocity_x_ptr, velocity_y_ptr, generationcoefficient_ptr, boundary_ptr, boundaryintegral_ptr);

   }

   // clear equations with essential boundary conditions
   for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
   {

        // subset the boundary
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        BoundaryUnit *boundary_ptr = boundary_group_ptr->boundary_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_essential_clear(a_mat, b_vec, x_vec, domain_ptr, boundary_ptr);

   }

   // iterate through each essential boundary
   for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
   {

        // subset the boundary
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        BoundaryUnit *boundary_ptr = boundary_group_ptr->boundary_ptr_vec[indx_d];
        BoundaryIntegralUnit *boundaryintegral_ptr = boundaryintegral_group_ptr->integral_ptr_vec[indx_d];

        // subset the scalars
        ScalarUnit *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *velocity_x_ptr = velocity_x_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *velocity_y_ptr = velocity_y_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_essential(a_mat, b_vec, x_vec, domain_ptr, diffusioncoefficient_ptr, velocity_x_ptr, velocity_y_ptr, generationcoefficient_ptr, boundary_ptr, boundaryintegral_ptr);

   }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainUnit *domain_ptr, DomainIntegralUnit *domainintegral_ptr,
    ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *velocity_x_ptr, ScalarUnit *velocity_y_ptr, ScalarUnit *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble diffcoeff_vec = diffusioncoefficient_ptr->get_neighbor_value(edid);
        VectorDouble velx_vec = velocity_x_ptr->get_neighbor_value(edid);
        VectorDouble vely_vec = velocity_y_ptr->get_neighbor_value(edid);
        VectorDouble gencoeff_vec = generationcoefficient_ptr->get_neighbor_value(edid);

        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // get group ID of points
        // used for getting matrix rows and columns
        VectorInt pfid_vec = value_group_ptr->get_neighbor_pfid(domain_ptr, edid);

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){
        for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_group_ptr->start_col + pfid_vec[indx_j];
            a_mat.coeffRef(mat_row, mat_col) += (
                diffcoeff_vec[indx_i] * domainintegral_ptr->integral_div_Ni_dot_div_Nj_vec[edid][indx_i][indx_j] +
                velx_vec[indx_i] * domainintegral_ptr->integral_Ni_derivative_Nj_x_vec[edid][indx_i][indx_j] +
                vely_vec[indx_i] * domainintegral_ptr->integral_Ni_derivative_Nj_y_vec[edid][indx_i][indx_j]
            );
        }}

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            b_vec.coeffRef(mat_row) += gencoeff_vec[indx_i] * domainintegral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_natural
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainUnit *domain_ptr,
    ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *velocity_x_ptr, ScalarUnit *velocity_y_ptr, ScalarUnit *generationcoefficient_ptr,
    BoundaryUnit *boundary_ptr, BoundaryIntegralUnit *boundaryintegral_ptr
)
{

    // iterate for each natural boundary element
    for (int bnid = 0; bnid < boundary_ptr->num_natural; bnid++)
    {

        // get element
        int egid = boundary_ptr->natural_egid_vec[bnid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary data
        int btid = boundary_ptr->boundary_btid_vec[bnid];
        VectorInt plid_vec = boundary_ptr->natural_blid_to_plid_vec[bnid];
        VectorDouble2D parameter_mat = boundary_ptr->natural_blid_to_parameter_vec[bnid];
        
        // calculate boundary key
        int boundary_key = boundaryintegral_ptr->get_boundary_key(plid_vec);

        // get group ID of points
        // used for getting matrix rows and columns
        VectorInt pfid_vec = value_group_ptr->get_neighbor_pfid(domain_ptr, edid);

        // iterate for each boundary point
        for (int indx_b = 0; indx_b < boundary_ptr->num_neighbor; indx_b++)
        {

            // get boundary point data
            int indx_i = plid_vec[indx_b];
            VectorDouble parameter_vec = parameter_mat[indx_b];

            // apply boundary condition
            int mat_row = start_row + pfid_vec[indx_i];
            switch (btid)
            {
                case 1:  // neumann
                    b_vec.coeffRef(mat_row) += parameter_vec[0] * boundaryintegral_ptr->integral_Ni_vec[edid][boundary_key][indx_i];
                break;
                case 2:  // robin
                    b_vec.coeffRef(mat_row) += parameter_vec[0] * boundaryintegral_ptr->integral_Ni_vec[edid][boundary_key][indx_i];
                    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++)
                    {
                        int mat_col = value_group_ptr->start_col + pfid_vec[indx_i];
                        a_mat.coeffRef(mat_row, mat_col) += -parameter_vec[1] * boundaryintegral_ptr->integral_Ni_Nj_vec[edid][boundary_key][indx_i][indx_j];
                    }
                break;
            }
            
        }
        
    }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_essential_clear
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainUnit *domain_ptr,
    BoundaryUnit *boundary_ptr
)
{

    // iterate for each essential boundary element
    for (int beid = 0; beid < boundary_ptr->num_essential; beid++)
    {

        // get element
        int egid = boundary_ptr->essential_egid_vec[beid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary data
        VectorInt plid_vec = boundary_ptr->essential_blid_to_plid_vec[beid];

        // get group ID of points
        // used for getting matrix rows and columns
        VectorInt pfid_vec = value_group_ptr->get_neighbor_pfid(domain_ptr, edid);

        // iterate for each boundary point
        for (int indx_b = 0; indx_b < boundary_ptr->num_neighbor; indx_b++)
        {
            
            // get boundary point data
            int indx_i = plid_vec[indx_b];

            // clear row
            int mat_row = start_row + pfid_vec[indx_i];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;

        }

    }

}

void PhysicsSteadyConvectionDiffusion::matrix_fill_essential
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    DomainUnit *domain_ptr,
    ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *velocity_x_ptr, ScalarUnit *velocity_y_ptr, ScalarUnit *generationcoefficient_ptr,
    BoundaryUnit *boundary_ptr, BoundaryIntegralUnit *boundaryintegral_ptr
)
{

    // iterate for each essential boundary element
    for (int beid = 0; beid < boundary_ptr->num_essential; beid++)
    {

        // get element
        int egid = boundary_ptr->essential_egid_vec[beid];
        int edid = domain_ptr->element_egid_to_edid_map[egid];

        // get boundary data
        int btid = boundary_ptr->boundary_btid_vec[beid];
        VectorInt plid_vec = boundary_ptr->essential_blid_to_plid_vec[beid];
        VectorDouble2D parameter_mat = boundary_ptr->essential_blid_to_parameter_vec[beid];
        
        // calculate boundary key
        int boundary_key = boundaryintegral_ptr->get_boundary_key(plid_vec);

        // get group ID of points
        // used for getting matrix rows and columns
        VectorInt pfid_vec = value_group_ptr->get_neighbor_pfid(domain_ptr, edid);

        // iterate for each boundary point
        for (int indx_b = 0; indx_b < boundary_ptr->num_neighbor; indx_b++)
        {

            // get boundary point data
            int indx_i = plid_vec[indx_b];
            VectorDouble parameter_vec = parameter_mat[indx_b];

            // apply boundary condition
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_group_ptr->start_col + pfid_vec[indx_i];
            switch (btid)
            {
                case 0:  // dirichlet
                    a_mat.coeffRef(mat_row, mat_col) += 1.;
                    b_vec.coeffRef(mat_row) += parameter_vec[0];
                break;
            }
            
        }

    }

}

void PhysicsSteadyConvectionDiffusion::set_start_row(int start_row_in)
{
    /*

    Sets the starting row in A and b where entries are filled up.

    Arguments
    =========
    start_row_in : int
        Starting row in A and b.

    Returns
    =======
    (none)

    */

    start_row = start_row_in;

}

int PhysicsSteadyConvectionDiffusion::get_start_row()
{
    /*

    Returns the starting row.

    Arguments
    =========
    (none)

    Returns
    =======
    start_row : int
        Starting row in A and b.

    */

    return start_row;

}

std::vector<BoundaryGroup*> PhysicsSteadyConvectionDiffusion::get_boundary_group_ptr_vec()
{
    /*

    Returns the pointer to the BoundaryGroup object tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    boundary_group_ptr : BoundaryGroup*
        Pointer to BoundaryGroup object.

    */
    
    return boundary_group_ptr_vec;

}

std::vector<ScalarGroup*> PhysicsSteadyConvectionDiffusion::get_scalar_group_ptr_vec()
{
    /*

    Returns the vector containing pointers to ScalarGroup objects tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    scalar_group_ptr : vector<ScalarGroup*>
        Vector containing pointers to ScalarGroup objects.

    */
    
    return scalar_group_ptr_vec;

}

std::vector<VariableGroup*> PhysicsSteadyConvectionDiffusion::get_variable_group_ptr_vec()
{
    /*

    Returns the vector containing pointers to VariableGroup objects tied to this physics.

    Arguments
    =========
    (none)

    Returns
    =======
    variable_group_ptr : vector<VariableGroup*>
        Vector containing pointers to VariableGroup objects.

    */
    
    return variable_group_ptr_vec;

}

}

#endif
