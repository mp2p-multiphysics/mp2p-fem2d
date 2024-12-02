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
    domain_ptr_vec_in : vector<DomainUnit*>
        vector with pointers to DomainUnit objects.

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
