#ifndef MATRIXEQUATION_STEADY
#define MATRIXEQUATION_STEADY
#include <iostream>
#include <limits>
#include <set>
#include <vector>
#include "Eigen/Eigen"
#include "physicssteady_base.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class MatrixEquationSteady
{
    /*

    Represents the matrix equation (system of linear equations) Ax = b for use in steady-state problems.

    Variables
    =========
    physics_ptr_vec_in : vector<PhysicsSteadyBase*>
        vector with pointers to PhysicsSteadyBase objects.

    Functions
    =========
    set_iteration : void
        Set parameters needed for iterations when solving the matrix equation.
    solve : void
        Solves the matrix equation.

    Notes
    =====
    In the code; A, x, and b are referred to as a_mat, x_vec, and b_vec respectively.

    */

    public:

    // vector of physics
    std::vector<PhysicsSteadyBase*> physics_ptr_vec;
    std::vector<BoundaryGroup*> boundary_group_ptr_vec;
    std::vector<ScalarGroup*> scalar_group_ptr_vec;
    std::vector<VariableGroup*> variable_group_ptr_vec;

    // matrix equation variables
    Eigen::SparseMatrix<double> a_mat;
    Eigen::VectorXd b_vec;
    Eigen::VectorXd x_vec;
    int num_equation = 0;

    // settings for iteration
    int num_iteration_max = 0;
    double residual_tol = 0.;

    // function
    void set_iteration(int num_iteration_max_in, double residual_tol_in);
    void solve(bool verbose);

    // default constructor
    MatrixEquationSteady() {};

    // constructor
    MatrixEquationSteady(std::vector<PhysicsSteadyBase*> physics_ptr_vec_in)
    {

        // store vector of pointers to physics
        physics_ptr_vec = physics_ptr_vec_in;

        // generate starting rows and columns
        // assign starting row in matrix to test functions (physics)
        // assign starting column in matrix to variables

        // initialize starting rows and columns
        int assign_start_row = 0;
        int assign_start_col = 0;

        // iterate through each physics
        for (auto physics_ptr : physics_ptr_vec)
        {

            // iterate through each variable group
            for (auto variable_group_ptr : physics_ptr->get_variable_group_ptr_vec())
            {
                
                // assign starting column to variable if none yet
                // increment assign_start_col by number of domain points
                if (variable_group_ptr->start_col == -1)
                {
                    variable_group_ptr->start_col = assign_start_col;
                    assign_start_col += variable_group_ptr->num_point;
                }

                // assign starting row to physics if none yet
                // increment assign_start_row by number of new domain points
                if (physics_ptr->get_start_row() == -1)
                {
                    physics_ptr->set_start_row(assign_start_row);
                    assign_start_row = assign_start_col;
                }

            }

        }

        // get number of linear equations (total number of domain points)
        num_equation = assign_start_col;

        // get vector of boundary groups
        // iterate through each physics and store boundary groups
        std::set<BoundaryGroup*> boundary_group_ptr_set;
        for (auto physics_ptr : physics_ptr_vec)
        {
            boundary_group_ptr_set.insert(physics_ptr->get_boundary_group_ptr());
        }
        boundary_group_ptr_vec = std::vector<BoundaryGroup*>(boundary_group_ptr_set.begin(), boundary_group_ptr_set.end());

        // get vector of scalar groups
        // iterate through each physics and store scalar groups
        std::set<ScalarGroup*> scalar_group_ptr_set;
        for (auto physics_ptr : physics_ptr_vec) {
        for (auto scalar_group_ptr : physics_ptr->get_scalar_group_ptr_vec()) {
            scalar_group_ptr_set.insert(scalar_group_ptr);
        }}
        scalar_group_ptr_vec = std::vector<ScalarGroup*>(scalar_group_ptr_set.begin(), scalar_group_ptr_set.end());

        // get vector of variable groups
        // iterate through each physics and store variable groups
        std::set<VariableGroup*> variable_group_ptr_set;
        for (auto physics_ptr : physics_ptr_vec) {
        for (auto variable_group_ptr : physics_ptr->get_variable_group_ptr_vec()) {
            variable_group_ptr_set.insert(variable_group_ptr);
        }}
        variable_group_ptr_vec = std::vector<VariableGroup*>(variable_group_ptr_set.begin(), variable_group_ptr_set.end());

        // initialize matrix equation variables
        a_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
        b_vec = Eigen::VectorXd::Zero(num_equation);
        x_vec = Eigen::VectorXd::Zero(num_equation);
        
        // populate x_vec with initial values

        // iterate through each variable group
        for (auto variable_group_ptr : variable_group_ptr_vec)
        {

            // get starting row
            // note: column in a_mat = row in x_vec
            int start_row = variable_group_ptr->start_col;

            // iterate through each variable
            for (auto variable_ptr : variable_group_ptr->variable_t3_ptr_vec)
            {

                // iterate through each global ID
                for (auto pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
                {

                    // get domain and group IDs
                    int pfid = variable_group_ptr->point_pgid_to_pfid_map[pgid];
                    int pdid = variable_ptr->domain_ptr->point_pgid_to_pdid_map[pgid];

                    // get value from variable
                    double value = variable_ptr->point_value_vec[pdid];

                    // store value in x_vec
                    int vec_row = start_row + pfid;
                    x_vec.coeffRef(vec_row) = value;

                }

            }
            for (auto variable_ptr : variable_group_ptr->variable_q4_ptr_vec)
            {

                // iterate through each global ID
                for (auto pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
                {

                    // get domain and group IDs
                    int pfid = variable_group_ptr->point_pgid_to_pfid_map[pgid];
                    int pdid = variable_ptr->domain_ptr->point_pgid_to_pdid_map[pgid];

                    // get value from variable
                    double value = variable_ptr->point_value_vec[pdid];

                    // store value in x_vec
                    int vec_row = start_row + pfid;
                    x_vec.coeffRef(vec_row) = value;

                }

            }

        }

    }

    private:

    // functions
    void iterate_solution();
    void store_solution();
    void output_solution();

};

void MatrixEquationSteady::set_iteration(int num_iteration_max_in = 500, double residual_tol_in = 1e-5)
{
    /*

    Set parameters needed for iterations when solving the matrix equation.

    Arguments
    =========
    num_iteration_max_in : int
        Maximum number of iterations.
    residual_tol_in : double
        Residual value at which iterations are stopped.

    Returns
    =======
    (none)

    */

    // set values
    num_iteration_max = num_iteration_max_in;
    residual_tol = residual_tol_in;

}

void MatrixEquationSteady::solve(bool verbose = true)
{
    /*

    Solves the matrix equation.

    Arguments
    =========
    verbose : bool
        If true, outputs the iteration number and residual.
        Defaults to true.

    Returns
    =======
    (none)

    */

    // initalize norm of residual
    double residual_norm = std::numeric_limits<double>::max();

    // iterate to convergence
    for (int it = 0; it < num_iteration_max; it++)
    {

        // store previous value of x_vec
        Eigen::VectorXd x_last_iteration_vec = x_vec;

        // perform one iteration and store x_vec values into variables and scalars
        // this automatically updates a_mat, x_vec, and b_vec
        iterate_solution();
        store_solution();

        // calculate residual using previous x_vec and current a_mat and b_vec
        residual_norm = (a_mat*x_last_iteration_vec - b_vec).norm();

        // display iteration count and norm
        if (verbose)
        {
            std::cout << "Iteration: " << it << "; Residual L2 Norm: " << residual_norm << "\n";
        }
        
        // stop if convergence is reached
        if (residual_norm < residual_tol)
        {
            break;
        }

    }

    // write output files
    output_solution();

}

void MatrixEquationSteady::iterate_solution()
{

    // update boundary parameters using most recent variable values
    for (auto boundary_group_ptr : boundary_group_ptr_vec)
    {
        boundary_group_ptr->update_parameter();
    }

    // update scalar values using most recent variable values
    for (auto scalar_group_ptr : scalar_group_ptr_vec)
    {
        scalar_group_ptr->update_value();
    }

    // reset matrices
    a_mat = Eigen::SparseMatrix<double> (num_equation, num_equation);
    b_vec = Eigen::VectorXd::Zero(num_equation);
    x_vec = Eigen::VectorXd::Zero(num_equation);

    // fill up a_mat and b_vec with each physics
    for (auto physics_ptr : physics_ptr_vec)
    {
        physics_ptr->matrix_fill(a_mat, b_vec, x_vec);
    }

    // solve the matrix equation
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(a_mat);
    solver.factorize(a_mat);
    x_vec = solver.solve(b_vec);

}

void MatrixEquationSteady::store_solution()
{

    // iterate through each variable group
    for (auto variable_group_ptr : variable_group_ptr_vec)
    {

        // get starting row
        // note: column in a_mat = row in x_vec
        int start_row = variable_group_ptr->start_col;

        // iterate through each variable
        for (auto variable_ptr : variable_group_ptr->variable_t3_ptr_vec)
        {

            // iterate through each global ID
            for (auto pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {

                // get domain and group IDs
                int pfid = variable_group_ptr->point_pgid_to_pfid_map[pgid];
                int pdid = variable_ptr->domain_ptr->point_pgid_to_pdid_map[pgid];

                // get value from x_vec
                int vec_row = start_row + pfid;
                double value = x_vec.coeffRef(vec_row);

                // store value in variable
                variable_ptr->point_value_vec[pdid] = value;

            }

        }
        for (auto variable_ptr : variable_group_ptr->variable_q4_ptr_vec)
        {

            // iterate through each global ID
            for (auto pgid : variable_ptr->domain_ptr->point_pdid_to_pgid_vec)
            {

                // get domain and group IDs
                int pfid = variable_group_ptr->point_pgid_to_pfid_map[pgid];
                int pdid = variable_ptr->domain_ptr->point_pgid_to_pdid_map[pgid];

                // get value from x_vec
                int vec_row = start_row + pfid;
                double value = x_vec.coeffRef(vec_row);

                // store value in variable
                variable_ptr->point_value_vec[pdid] = value;

            }

        }

    }

}

void MatrixEquationSteady::output_solution()
{

    // iterate through each variable group
    for (auto variable_group_ptr : variable_group_ptr_vec)
    {
        variable_group_ptr->output_csv();
    }

    // iterate through each scalar group
    for (auto scalar_group_ptr : scalar_group_ptr_vec)
    {
        scalar_group_ptr->output_csv();
    }

}

}

#endif
