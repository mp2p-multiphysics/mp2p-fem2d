#ifndef PHYSICSSTEADY_DIFFUSION
#define PHYSICSSTEADY_DIFFUSION
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "domain_1d.hpp"
#include "domain_2d.hpp"
#include "integral_1d.hpp"
#include "integral_2d.hpp"
#include "physicssteady_base.hpp"
#include "scalar_1d.hpp"
#include "scalar_2d.hpp"
#include "variable_2d.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class PhysicsSteadyDiffusion : public PhysicsSteadyBase
{
    /*

    Single-component steady-state diffusion equation.    
    
    0 = -div(-b * grad(u)) + c

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
    VariableGroup* value_ptr;

    // domain objects
    std::vector<Domain2D*> domain_ptr_vec;
    std::vector<Integral2D*> integral_ptr_vec;
    std::vector<Scalar2D*> diffusioncoefficient_ptr_vec;
    std::vector<Scalar2D*> generationcoefficient_ptr_vec;

    // boundary objects - dirichlet
    std::vector<Domain1D*> dirichlet_domain_ptr_vec;
    std::vector<Scalar1D*> dirichlet_constant_ptr_vec;

    // boundary objects - neumann
    std::vector<Domain1D*> neumann_domain_ptr_vec;
    std::vector<Integral1D*> neumann_integral_ptr_vec;
    std::vector<Scalar1D*> neumann_flux_ptr_vec;

    // vectors of objects to update
    std::vector<Scalar1D*> scalar1d_ptr_vec;
    std::vector<Scalar2D*> scalar2d_ptr_vec;
    std::vector<VariableGroup*> variablegroup_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    void matrix_fill(Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec);
    void set_variablegroup(VariableGroup &value_in);
    void set_domain(Domain2D &domain_in, Integral2D &integral_in, Scalar2D &diffusioncoefficient_in, Scalar2D &generationcoefficient_in);
    void set_boundary_dirichlet(Domain1D &domain_in, Scalar1D &value_constant_in);
    void set_boundary_neumann(Domain1D &domain_in, Integral1D &integral_in, Scalar1D &value_flux_in);

    // getter and setter functions
    void set_start_row(int start_row_in) {start_row = start_row_in;}
    int get_start_row() {return start_row;}
    std::vector<Scalar1D*> get_scalar1d_ptr_vec() {return scalar1d_ptr_vec;}
    std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsSteadyDiffusion() {}

    private:

    void matrix_fill_domain
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        Domain2D *domain_ptr, Integral2D *integral_ptr,
        Scalar2D *diffusioncoefficient_ptr, Scalar2D *generationcoefficient_ptr
    );
    void matrix_fill_neumann
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        Domain1D *domain_ptr, Integral1D *integral_ptr,
        Scalar1D *value_flux_ptr
    );
    void matrix_fill_dirichlet_clear
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        Domain1D *domain_ptr
    );
    void matrix_fill_dirichlet
    (
        Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
        Domain1D *domain_ptr,
        Scalar1D *value_constant_ptr
    );

};

void PhysicsSteadyDiffusion::set_variablegroup(VariableGroup &value_in)
{
    
    // set variable groups
    value_ptr = &value_in;

    // add to vector of variable groups
    variablegroup_ptr_vec.push_back(&value_in);

}

void PhysicsSteadyDiffusion::set_domain(Domain2D &domain_in, Integral2D &integral_in, Scalar2D &diffusioncoefficient_in, Scalar2D &generationcoefficient_in)
{

    // add to vector of domain objects
    domain_ptr_vec.push_back(&domain_in);
    integral_ptr_vec.push_back(&integral_in);
    diffusioncoefficient_ptr_vec.push_back(&diffusioncoefficient_in);
    generationcoefficient_ptr_vec.push_back(&generationcoefficient_in);

    // add to vector of scalar2d objects
    scalar2d_ptr_vec.push_back(&diffusioncoefficient_in);
    scalar2d_ptr_vec.push_back(&generationcoefficient_in);

    // calculate integrals
    integral_in.evaluate_integral_Ni();
    integral_in.evaluate_integral_div_Ni_dot_div_Nj();

}

void PhysicsSteadyDiffusion::set_boundary_dirichlet(Domain1D &domain_in, Scalar1D &value_constant_in)
{

    // add to vector of boundary objects
    dirichlet_domain_ptr_vec.push_back(&domain_in);
    dirichlet_constant_ptr_vec.push_back(&value_constant_in);

    // add to vector of scalar1d objects
    scalar1d_ptr_vec.push_back(&value_constant_in);

}

void PhysicsSteadyDiffusion::set_boundary_neumann(Domain1D &domain_in, Integral1D &integral_in, Scalar1D &value_flux_in)
{

    // add to vector of boundary objects
    neumann_domain_ptr_vec.push_back(&domain_in);
    neumann_integral_ptr_vec.push_back(&integral_in);
    neumann_flux_ptr_vec.push_back(&value_flux_in);

    // add to vector of scalar1d objects
    scalar1d_ptr_vec.push_back(&value_flux_in);

    // calculate integrals
    integral_in.evaluate_integral_Ni();

}

void PhysicsSteadyDiffusion::matrix_fill
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
   for (int indx_d = 0; indx_d < domain_ptr_vec.size(); indx_d++)
   {

        // subset domain objects
        Domain2D *domain_ptr = domain_ptr_vec[indx_d];
        Integral2D *integral_ptr = integral_ptr_vec[indx_d];
        Scalar2D *diffusioncoefficient_ptr = diffusioncoefficient_ptr_vec[indx_d];
        Scalar2D *generationcoefficient_ptr = generationcoefficient_ptr_vec[indx_d];

        // fill up matrix with domain equations
        matrix_fill_domain(a_mat, b_vec, x_vec, domain_ptr, integral_ptr, diffusioncoefficient_ptr, generationcoefficient_ptr);

   }

   // iterate through each neumann boundary
   for (int indx_d = 0; indx_d < neumann_domain_ptr_vec.size(); indx_d++)
   {

        // subset domain objects
        Domain1D *domain_ptr = neumann_domain_ptr_vec[indx_d];
        Integral1D *integral_ptr = neumann_integral_ptr_vec[indx_d];
        Scalar1D *value_flux_ptr = neumann_flux_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_neumann(a_mat, b_vec, x_vec, domain_ptr, integral_ptr, value_flux_ptr);

   }

   // clear equations with dirichlet boundary conditions
   for (int indx_d = 0; indx_d < dirichlet_domain_ptr_vec.size(); indx_d++)
   {

        // subset domain objects
        Domain1D *domain_ptr = dirichlet_domain_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_dirichlet_clear(a_mat, b_vec, x_vec, domain_ptr);

   }

   // iterate through each dirichlet boundary
   for (int indx_d = 0; indx_d < dirichlet_domain_ptr_vec.size(); indx_d++)
   {

        // subset domain objects
        Domain1D *domain_ptr = dirichlet_domain_ptr_vec[indx_d];
        Scalar1D *value_constant_ptr = dirichlet_constant_ptr_vec[indx_d];

        // fill up matrix with boundary conditions
        matrix_fill_dirichlet(a_mat, b_vec, x_vec, domain_ptr, value_constant_ptr);

   }

}

void PhysicsSteadyDiffusion::matrix_fill_domain
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    Domain2D *domain_ptr, Integral2D *integral_ptr,
    Scalar2D *diffusioncoefficient_ptr, Scalar2D *generationcoefficient_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble diffcoeff_vec = diffusioncoefficient_ptr->get_neighbor_value(edid);
        VectorDouble gencoeff_vec = generationcoefficient_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // matrix indexing
        // matrix row = start_row of test function (physics) + group ID of variable
        // matrix column = start_column of variable + group ID of variable

        // calculate a_mat coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++){
        for (int indx_j = 0; indx_j < domain_ptr->num_neighbor; indx_j++){
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_ptr->start_col + pfid_vec[indx_j];
            a_mat.coeffRef(mat_row, mat_col) += diffcoeff_vec[indx_i] * integral_ptr->integral_div_Ni_dot_div_Nj_vec[edid][indx_i][indx_j];
        }}

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            b_vec.coeffRef(mat_row) += gencoeff_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsSteadyDiffusion::matrix_fill_neumann
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    Domain1D *domain_ptr, Integral1D *integral_ptr,
    Scalar1D *value_flux_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble value_flux_vec = value_flux_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // calculate b_vec coefficients
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            b_vec.coeffRef(mat_row) += value_flux_vec[indx_i] * integral_ptr->integral_Ni_vec[edid][indx_i];
        }

    }

}

void PhysicsSteadyDiffusion::matrix_fill_dirichlet_clear
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    Domain1D *domain_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            a_mat.row(mat_row) *= 0.;
            b_vec.coeffRef(mat_row) = 0.;
        }

    }

}

void PhysicsSteadyDiffusion::matrix_fill_dirichlet
(
    Eigen::SparseMatrix<double> &a_mat, Eigen::VectorXd &b_vec, Eigen::VectorXd &x_vec,
    Domain1D *domain_ptr,
    Scalar1D *value_constant_ptr
)
{

    // iterate for each domain element
    for (int edid = 0; edid < domain_ptr->num_element; edid++)
    {

        // get scalar values of points around element
        VectorDouble value_constant_vec = value_constant_ptr->get_neighbor_value(edid);

        // get group ID of points
        VectorInt pfid_vec = value_ptr->get_neighbor_pfid(domain_ptr, edid);

        // clear rows
        for (int indx_i = 0; indx_i < domain_ptr->num_neighbor; indx_i++)
        {
            int mat_row = start_row + pfid_vec[indx_i];
            int mat_col = value_ptr->start_col + pfid_vec[indx_i];
            a_mat.coeffRef(mat_row, mat_col) += 1.;
            b_vec.coeffRef(mat_row) += value_constant_vec[indx_i];
        }

    }

}

}

#endif
