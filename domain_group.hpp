#ifndef DOMAIN_GROUP
#define DOMAIN_GROUP
#include <vector>
#include "domain_quad4.hpp"
#include "domain_tri3.hpp"

namespace FEM2D
{

class DomainGroup
{
    /*

    Groups domains that are used in the same physics.

    Variables
    =========
    domain_l2_ptr_vec_in : vector<DomainLine2*>
        vector with pointers to DomainLine2 objects.
    
    */

    public:

    // vector with domains in group
    std::vector<DomainTri3*> domain_t3_ptr_vec;
    std::vector<DomainQuad4*> domain_q4_ptr_vec;

    // default constructor
    DomainGroup() {}

    // constructor
    DomainGroup(std::vector<DomainTri3*> domain_t3_ptr_vec_in, std::vector<DomainQuad4*> domain_q4_ptr_vec_in)
    {
        
        // store variables
        domain_t3_ptr_vec = domain_t3_ptr_vec_in;
        domain_q4_ptr_vec = domain_q4_ptr_vec_in;

    }

};

}

#endif
