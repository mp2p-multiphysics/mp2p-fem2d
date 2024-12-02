#ifndef DOMAIN_GROUP
#define DOMAIN_GROUP
#include <vector>
#include "domain_unit.hpp"

namespace FEM2D
{

class DomainGroup
{
    /*

    Groups domains that are used in the same physics.

    Variables
    =========
    domain_t3_ptr_vec_in : vector<DomainTri3*>
        vector with pointers to DomainTri3 objects.
    domain_q4_ptr_vec_in : vector<DomainQuad4*>
        vector with pointers to DomainQuad4 objects.

    */

    public:

    // vector with domains in group
    std::vector<DomainUnit*> domain_ptr_vec;

    // default constructor
    DomainGroup() {}

    // constructor
    DomainGroup(std::vector<DomainUnit*> domain_ptr_vec_in)
    {
        
        // store variables
        domain_ptr_vec = domain_ptr_vec_in;

    }

};

}

#endif
