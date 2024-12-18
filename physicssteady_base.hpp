#ifndef PHYSICSSTEADY_BASE
#define PHYSICSSTEADY_BASE
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "scalar_1d.hpp"
#include "scalar_2d.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class PhysicsSteadyBase
{
    /*

    Base class for steady-state physics.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_start_row : void
        Sets the starting row in A and b where entries are filled up.
    get_start_row : int
        Returns the starting row.
    get_scalar1d_ptr_vec : vector<BoundaryGroup*>
        Returns the vector containing pointers to Scalar1D objects tied to this physics.
    get_scalar2d_ptr_vec : vector<ScalarGroup*>
        Returns the vector containing pointers to Scalar2D objects tied to this physics.
    get_variablegroup_ptr_vec : vector<VariableGroup*>
        Returns the vector containing pointers to VariableGroup objects tied to this physics.

    */

    public:

    // vector of scalars and vectors
    std::vector<Scalar1D*> scalar1d_ptr_vec;
    std::vector<Scalar2D*> scalar2d_ptr_vec;
    std::vector<VariableGroup*> variablegroup_ptr_vec;

    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    virtual void matrix_fill(EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec);
    virtual void set_start_row(int start_row_in) {start_row = start_row_in;}
    virtual int get_start_row() {return start_row;}
    virtual std::vector<Scalar1D*> get_scalar1d_ptr_vec() {return scalar1d_ptr_vec;}
    virtual std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    virtual std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsSteadyBase() {}

};

void PhysicsSteadyBase::matrix_fill(EigenSparseMatrix &a_mat, EigenVector &b_vec, EigenVector &x_vec)
{
    /*

    Fill up the matrix equation Ax = b with entries as dictated by the physics. 

    Arguments
    =========
    a_mat : EigenSparseMatrix
        A in Ax = b.
    b_vec : EigenVector
        b in Ax = b.
    x_vec : EigenVector
        x in Ax = b.

    Returns
    =======
    (none)

    */

}

}

#endif
