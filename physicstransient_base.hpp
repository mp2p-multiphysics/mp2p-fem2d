#ifndef PHYSICSTRANSIENT_BASE
#define PHYSICSTRANSIENT_BASE
#include <vector>
#include "Eigen/Eigen"
#include "container_typedef.hpp"
#include "scalar_1d.hpp"
#include "scalar_2d.hpp"
#include "variable_group.hpp"

namespace FEM2D
{

class PhysicsTransientBase
{
    /*

    Base class for transient physics.

    Functions
    =========
    matrix_fill : void
        Fill up the matrix equation Ax = b with entries as dictated by the physics. 
    set_start_row : void
        Sets the starting row in A and b where entries are filled up.
    get_start_row : int
        Returns the starting row.
    get_scalar1d_ptr_vec : vector<Scalar1D*>
        Returns the vector containing pointers to Scalar1D objects tied to this physics.
    get_scalar2d_ptr_vec : vector<Scalar2D*>
        Returns the vector containing pointers to Scalar2D objects tied to this physics.
    get_variablegroup_ptr_vec : vector<VariableGroup*>
        Returns the vector containing pointers to VariableGroup objects tied to this physics.

    */

    public:

    // vector of scalar and variable groups
    std::vector<Scalar1D*> scalar1d_ptr_vec;
    std::vector<Scalar2D*> scalar2d_ptr_vec;
    std::vector<VariableGroup*> variablegroup_ptr_vec;
    
    // starting row of test functions in matrix equation
    int start_row = -1;

    // functions
    virtual void matrix_fill(
        EigenTripletVector &a_trivec, EigenTripletVector &c_trivec, EigenVector &d_vec,
        EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt
    );
    virtual void set_start_row(int start_row_in) {start_row = start_row_in;}
    virtual int get_start_row() {return start_row;}
    virtual std::vector<Scalar1D*> get_scalar1d_ptr_vec() {return scalar1d_ptr_vec;}
    virtual std::vector<Scalar2D*> get_scalar2d_ptr_vec() {return scalar2d_ptr_vec;}
    virtual std::vector<VariableGroup*> get_variablegroup_ptr_vec() {return variablegroup_ptr_vec;}

    // default constructor
    PhysicsTransientBase() {}

};

void PhysicsTransientBase::matrix_fill
(
    EigenTripletVector &a_trivec, EigenTripletVector &c_trivec, EigenVector &d_vec,
    EigenVector &x_vec, EigenVector &x_last_timestep_vec, double dt
)
{
    /*

    Fill up the matrix equation Ax(t+1) = Cx(t) + d with entries as dictated by the physics. 

    Arguments
    =========
    a_trivec : EigenTripletVector
        A in Ax(t+1) = Cx(t) + d.
    c_trivec : EigenTripletVector
        C in Ax(t+1) = Cx(t) + d.
    d_vec : EigenVector
        d in Ax(t+1) = Cx(t) + d.
    x_vec : EigenVector
        x(t+1) in Ax(t+1) = Cx(t) + d.
    x_last_timestep_vec : EigenVector
        x(t) in Ax(t+1) = Cx(t) + d.
    dt : double
        Length of the timestep.

    Returns
    =======
    (none)

    */

}

}

#endif
