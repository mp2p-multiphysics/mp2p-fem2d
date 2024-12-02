#ifndef PHYSICSTRANSIENT_DIFFUSION
#define PHYSICSTRANSIENT_DIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "boundary_group.hpp"
#include "container_typedef.hpp"
#include "integral_group.hpp"
#include "domain_group.hpp"
#include "physicstransient_base.hpp"
#include "scalar_group.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class PhysicsTransientDiffusion : public PhysicsTransientBase
{
    /*

    Single-component transient diffusion equation.    
    
    a * du/dt = -div(-b * grad(u)) + c

    Variables
    =========
    domain_group_in : DomainGroup
        Domains where this physics is applied to.
    boundary_group_in : BoundaryGroup
        Boundaries where this physics is applied to.
    integral_group_in : IntegralGroup
        Test function integrals that this physics uses.
    value_group_in : VariableGroup
        u in a * du/dt = -div(-b * grad(u)) + c.
        This will be solved for by the matrix equation.
    derivativecoefficient_group_ptr : ScalarGroup
        a in a * du/dt = -div(-b * grad(u)) + c.
    diffusioncoefficient_group_in : ScalarGroup
        b in a * du/dt = -div(-b * grad(u)) + c.
    generationcoefficient_group_in : ScalarGroup
        c in a * du/dt = -div(-b * grad(u)) + c.

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
    BoundaryGroup *boundary_group_ptr;
    IntegralGroup *integral_group_ptr;
    VariableGroup *value_group_ptr;
    ScalarGroup *derivativecoefficient_group_ptr;
    ScalarGroup *diffusioncoefficient_group_ptr;
    ScalarGroup *generationcoefficient_group_ptr;

    // vector of scalar and variable groups
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt
    );
    void set_start_row(int start_row_in);
    int get_start_row();
    BoundaryGroup* get_boundary_group_ptr();
    std::vector<ScalarGroup*> get_scalar_group_ptr_vec();
    std::vector<VariableGroup*> get_variable_group_ptr_vec();

    // default constructor
    PhysicsTransientDiffusion() {}

    // constructor
    PhysicsTransientDiffusion
    (
        DomainGroup &domain_group_in, BoundaryGroup &boundary_group_in, IntegralGroup &integral_group_in,
        VariableGroup &value_group_in,
        ScalarGroup &derivativecoefficient_group_in, ScalarGroup &diffusioncoefficient_group_in, ScalarGroup &generationcoefficient_group_in
    )
    {
        
        // store variables
        domain_group_ptr = &domain_group_in;
        boundary_group_ptr = &boundary_group_in;
        integral_group_ptr = &integral_group_in;
        value_group_ptr = &value_group_in;
        derivativecoefficient_group_ptr = &derivativecoefficient_group_in;
        diffusioncoefficient_group_ptr = &diffusioncoefficient_group_in;
        generationcoefficient_group_ptr = &generationcoefficient_group_in;

        // set boundary conditions
        boundary_group_ptr->set_boundary_type({0}, {1, 2});

        // vector of variable groups
        scalar_group_ptr_vec = {derivativecoefficient_group_ptr, diffusioncoefficient_group_ptr, generationcoefficient_group_ptr};
        variable_group_ptr_vec = {value_group_ptr};

        // calculate integrals
        integral_group_ptr->evaluate_integral_div_Ni_dot_div_Nj();
        integral_group_ptr->evaluate_integral_Ni_Nj();
        integral_group_ptr->evaluate_integral_Ni();
        integral_group_ptr->evaluate_integral_boundary_Ni();
        integral_group_ptr->evaluate_integral_boundary_Ni_Nj();

    }

    private:
    void matrix_fill_domain
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        DomainUnit *domain_ptr, IntegralUnit *integral_ptr,
        ScalarUnit *derivativecoefficient_ptr, ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *generationcoefficient_ptr
    );
    void matrix_fill_natural
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        DomainUnit *domain_ptr, BoundaryUnit *boundary_ptr,  IntegralUnit *integral_ptr,
        ScalarUnit *derivativecoefficient_ptr, ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *generationcoefficient_ptr
    );
    void matrix_fill_essential_clear
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        DomainUnit *domain_ptr, BoundaryUnit *boundary_ptr
    );
    void matrix_fill_essential
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
        Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
        DomainUnit *domain_ptr, BoundaryUnit *boundary_ptr,  IntegralUnit *integral_ptr,
        ScalarUnit *derivativecoefficient_ptr, ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *generationcoefficient_ptr
    );

};

void PhysicsTransientDiffusion::matrix_fill
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt
)
{
    /*

    Fill up the matrix equation Ax(t+1) = Cx(t) + d with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : Eigen::SparseMatrix<double>
        A in Ax(t+1) = Cx(t) + d.
    c_mat : Eigen::SparseMatrix<double>
        C in Ax(t+1) = Cx(t) + d.
    d_vec : Eigen::VectorXd
        d in Ax(t+1) = Cx(t) + d.
    x_vec : Eigen::VectorXd
        x(t+1) in Ax(t+1) = Cx(t) + d.
    x_last_timestep_vec : Eigen::VectorXd
        x(t) in Ax(t+1) = Cx(t) + d.
    dt : double
        Length of the timestep.

    Returns
    =======
    (none)

    */

    // iterate through each domain
    for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
    {

        // subset the domain and integrals
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        IntegralUnit *integral_ptr = integral_group_ptr->integral_ptr_vec[indx_d];

        // subset the scalars
        ScalarUnit *derivativecoefficient_ptr = derivativecoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(
            a_mat, c_mat, d_vec,
            x_vec, x_last_timestep_vec, dt,
            domain_ptr, integral_ptr,
            derivativecoefficient_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr
        );

    }

    // iterate through each natural boundary
    for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
    {

        // subset the domain and integrals
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        BoundaryUnit *boundary_ptr = boundary_group_ptr->boundary_ptr_vec[indx_d];
        IntegralUnit *integral_ptr = integral_group_ptr->integral_ptr_vec[indx_d];

        // subset the scalars
        ScalarUnit *derivativecoefficient_ptr = derivativecoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_natural(
            a_mat, c_mat, d_vec,
            x_vec, x_last_timestep_vec, dt,
            domain_ptr, boundary_ptr, integral_ptr,
            derivativecoefficient_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr
        );

    }

    // clear equations with essential boundary conditions
    for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
    {

        // subset the domain and integrals
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        BoundaryUnit *boundary_ptr = boundary_group_ptr->boundary_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_essential_clear(
            a_mat, c_mat, d_vec,
            x_vec, x_last_timestep_vec, dt,
            domain_ptr, boundary_ptr
        );

    }

    // iterate through each essential boundary
    for (int indx_d = 0; indx_d < domain_group_ptr->domain_ptr_vec.size(); indx_d++)
    {

        // subset the domain and integrals
        DomainUnit *domain_ptr = domain_group_ptr->domain_ptr_vec[indx_d];
        BoundaryUnit *boundary_ptr = boundary_group_ptr->boundary_ptr_vec[indx_d];
        IntegralUnit *integral_ptr = integral_group_ptr->integral_ptr_vec[indx_d];

        // subset the scalars
        ScalarUnit *derivativecoefficient_ptr = derivativecoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *diffusioncoefficient_ptr = diffusioncoefficient_group_ptr->scalar_ptr_vec[indx_d];
        ScalarUnit *generationcoefficient_ptr = generationcoefficient_group_ptr->scalar_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_essential(
            a_mat, c_mat, d_vec,
            x_vec, x_last_timestep_vec, dt,
            domain_ptr, boundary_ptr, integral_ptr,
            derivativecoefficient_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr
        );

    }

}

void PhysicsTransientDiffusion::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    DomainUnit *domain_ptr, IntegralUnit *integral_ptr,
    ScalarUnit *derivativecoefficient_ptr, ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble dervcoeff_vec = derivativecoefficient_ptr->get_neighbor_value(edid);
        VectorDouble diffcoeff_vec = diffusioncoefficient_ptr->get_neighbor_value(edid);
        VectorDouble gencoeff_vec = generationcoefficient_ptr->get_neighbor_value(edid);

        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // get group ID of points
        // used for getting matrix rows and columns
        VectorInt pfid_vec = value_group_ptr->get_neighbor_pfid(domain_ptr, edid);

        // calculate a_mat and c_mat coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){
        for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_group_ptr->start_col + pfid_vec[indx_j];
            a_mat.coeffRef(mat_row, mat_col) += (
                (dervcoeff_vec[indx_i]/dt) * integral_ptr->integral_Ni_Nj_vec[edid][indx_i][indx_j] +
                diffcoeff_vec[indx_i] * integral_ptr->integral_div_Ni_dot_div_Nj_vec[edid][indx_i][indx_j]
            );
            c_mat.coeffRef(mat_row, mat_col) += (dervcoeff_vec[indx_i]/dt) * integral_ptr->integral_Ni_Nj_vec[edid][indx_i][indx_j];
        }}

        // calculate d_vec coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            d_vec.coeffRef(mat_row) += gencoeff_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsTransientDiffusion::matrix_fill_natural
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    DomainUnit *domain_ptr, BoundaryUnit *boundary_ptr,  IntegralUnit *integral_ptr,
    ScalarUnit *derivativecoefficient_ptr, ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *generationcoefficient_ptr
)
{

    // iterate for each flux boundary element
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
        int boundary_key = integral_ptr->get_boundary_key(plid_vec);

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
                    d_vec.coeffRef(mat_row) += parameter_vec[0] * integral_ptr->integral_boundary_Ni_vec[edid][boundary_key][indx_i];
                break;
                case 2:  // robin
                    d_vec.coeffRef(mat_row) += parameter_vec[0] * integral_ptr->integral_boundary_Ni_vec[edid][boundary_key][indx_i];
                    for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++)
                    {
                        int mat_col = value_group_ptr->start_col + pfid_vec[indx_i];
                        a_mat.coeffRef(mat_row, mat_col) += -parameter_vec[1] * integral_ptr->integral_boundary_Ni_Nj_vec[edid][boundary_key][indx_i][indx_j];
                    }
                break;
            }
            
        }

    }

}

void PhysicsTransientDiffusion::matrix_fill_essential_clear
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    DomainUnit *domain_ptr, BoundaryUnit *boundary_ptr
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
            c_mat.row(mat_row) *= 0.;
            d_vec.coeffRef(mat_row) = 0.;

        }

    }

}

void PhysicsTransientDiffusion::matrix_fill_essential
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::SparseMatrix<double> &c_mat, Eigen::VectorXd &d_vec,
    Eigen::VectorXd &x_vec, Eigen::VectorXd &x_last_timestep_vec, double dt,
    DomainUnit *domain_ptr, BoundaryUnit *boundary_ptr,  IntegralUnit *integral_ptr,
    ScalarUnit *derivativecoefficient_ptr, ScalarUnit *diffusioncoefficient_ptr, ScalarUnit *generationcoefficient_ptr
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
        int boundary_key = integral_ptr->get_boundary_key(plid_vec);

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
                    d_vec.coeffRef(mat_row) += parameter_vec[0];
                break;
            }
            
        }

    }

}

void PhysicsTransientDiffusion::set_start_row(int start_row_in)
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

int PhysicsTransientDiffusion::get_start_row()
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

BoundaryGroup* PhysicsTransientDiffusion::get_boundary_group_ptr()
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
    
    return boundary_group_ptr;

}

std::vector<ScalarGroup*> PhysicsTransientDiffusion::get_scalar_group_ptr_vec()
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

std::vector<VariableGroup*> PhysicsTransientDiffusion::get_variable_group_ptr_vec()
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
