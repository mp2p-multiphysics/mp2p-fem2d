#ifndef INTEGRAL_PHYSICSGROUP_VECTOR
#define INTEGRAL_PHYSICSGROUP_VECTOR
#include <vector>
#include "integral_physicsgroup.hpp"

class IntegralPhysicsGroupVector
{
    /*

    Groups IntegralPhysicsGroup objects into a column vector.
    Used in physics with systems of partial differential equations.

    Variables
    =========
    integral_physicsgroup_ptr_vec : vector<IntegralPhysicsGroup*>
        vector with pointers to IntegralPhysicsGroup objects.
    
    Functions
    =========
    get_entry : IntegralPhysicsGroup*
        Returns the pointer to the IntegralPhysicsGroup object at a given vector position.

    */

    public:

    // variables
    int num_entry = 0;
    std::vector<IntegralPhysicsGroup*> integral_physicsgroup_ptr_vec;

    // functions
    IntegralPhysicsGroup* get_entry(int vector_row);

    // default constructor
    IntegralPhysicsGroupVector()
    {

    }

    // constructor
    IntegralPhysicsGroupVector(std::vector<IntegralPhysicsGroup*> integral_physicsgroup_ptr_vec_in)
    {

        // store vector of variables
        integral_physicsgroup_ptr_vec = integral_physicsgroup_ptr_vec_in;

        // count number of variables
        num_entry = integral_physicsgroup_ptr_vec.size();

    }

};

IntegralPhysicsGroup* IntegralPhysicsGroupVector::get_entry(int vector_row)
{
    /*

    Returns the pointer to the IntegralPhysicsGroup object at a vector position.

    Arguments
    =========
    vector_row : int
        Row in column vector
    
    Returns
    =======
    integral_physicsgroup_ptr : IntegralPhysicsGroup*
        Pointer to the IntegralPhysicsGroup object at the specified row.

    */
    
    return integral_physicsgroup_ptr_vec[vector_row];
}

#endif
