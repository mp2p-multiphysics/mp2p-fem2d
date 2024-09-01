#ifndef INPUT_BOUNDARY_QUAD4_CSV
#define INPUT_BOUNDARY_QUAD4_CSV
#include <sstream>
#include <fstream>
#include <vector>
#include "grid_quad4.hpp"

BoundaryQuad4Struct input_boundary_quad4_csv(std::string file_in_flux_str, std::string file_in_value_str, std::string file_in_config_str)
{

    // initialize struct with boundary condition (BC) data
    BoundaryQuad4Struct bq4s;

    // read file with flux BC data
    std::ifstream file_in_flux_stream(file_in_flux_str);

    // initialize for iteration
    bool is_flux_header = true;  // true while reading header
    std::string line_flux_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_flux_stream, line_flux_str))
    {

        // skip header
        if (is_flux_header)
        {
            is_flux_header = false; // not reading header
            continue;
        }

        // count number of particles
        bq4s.num_element_flux++;

        // convert line string into stringstream
        std::stringstream line_flux_stream(line_flux_str);

        // initialize for iteration
        int value_flux_num = 0;  // counts position of value
        std::string value_flux_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_flux_stream, value_flux_str, ','))
        {

            // store values in appropriate vector
            switch (value_flux_num)
            {
                case 0: bq4s.element_flux_id_vec.push_back(std::stoi(value_flux_str)); break;
                case 1: bq4s.element_flux_pa_loc_vec.push_back(std::stoi(value_flux_str)); break;
                case 2: bq4s.element_flux_pb_loc_vec.push_back(std::stoi(value_flux_str)); break;
                case 3: bq4s.element_flux_config_id_vec.push_back(std::stoi(value_flux_str)); break;
            }

            // increment value count
            value_flux_num++;

        }

    }

    // close point file
    file_in_flux_stream.close();    



    // read file with value BC data
    std::ifstream file_in_value_stream(file_in_value_str);

    // initialize for iteration
    bool is_value_header = true;  // true while reading header
    std::string line_value_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_value_stream, line_value_str))
    {

        // skip header
        if (is_value_header)
        {
            is_value_header = false; // not reading header
            continue;
        }

        // count number of particles
        bq4s.num_element_value++;

        // convert line string into stringstream
        std::stringstream line_value_stream(line_value_str);

        // initialize for iteration
        int value_value_num = 0;  // counts position of value
        std::string value_value_str;  // stores values in lines

        // iterate through each value
        while (std::getline(line_value_stream, value_value_str, ','))
        {

            // store values in appropriate vector
            switch (value_value_num)
            {
                case 0: bq4s.element_value_id_vec.push_back(std::stoi(value_value_str)); break;
                case 1: bq4s.element_value_pa_loc_vec.push_back(std::stoi(value_value_str)); break;
                case 2: bq4s.element_value_pb_loc_vec.push_back(std::stoi(value_value_str)); break;
                case 3: bq4s.element_value_config_id_vec.push_back(std::stoi(value_value_str)); break;
            }

            // increment value count
            value_value_num++;

        }

    }

    // close point file
    file_in_value_stream.close();



    // read file with BC type
    std::ifstream file_in_config_stream(file_in_config_str);

    // initialize struct with BC type data
    BoundaryConfigQuad4Struct bcq4s_sub;  // temporarily stores BoundaryTypeData for insertion into bq4s
    BoundaryConfigQuad4Struct bcq4s_new;  // new copy of BoundaryTypeData object

    // initialize for iteration
    bool is_config_header = true;  // true while reading header
    bool is_config_next = false;  // true if next entry is the BC type
    bool is_value_next = false;  // true if next entry is the BC value
    std::string line_config_str;  // stores lines in files

    // iterate for each line in the file
    while (std::getline(file_in_config_stream, line_config_str))
    {

        // skip first bracket
        if (is_config_header)
        {
            is_config_header = false; // not reading header
            continue;
        }

        // convert line string into stringstream
        std::stringstream line_config_stream(line_config_str);

        // initialize for iteration
        std::string subline_config_str;  // stores parts of lines (subline)

        // iterate through each substring
        while (std::getline(line_config_stream, subline_config_str, ' '))
        {
            
            // skip blank lines
            if (subline_config_str == "")
            {
                continue;
            }

            // store BC type
            if (is_config_next)
            {
                subline_config_str = subline_config_str.substr(1, subline_config_str.length()-3);  // remove quotes and comma
                bcq4s_sub.boundary_type_str = subline_config_str;
                is_config_next = false; // next line is not a BC type
            }

            // store BC value
            if (is_value_next)
            {
                
                // remove [ and ] from start and end
                subline_config_str = subline_config_str.substr(1, subline_config_str.length()-1);

                // convert subline string into stringstream
                std::stringstream subline_config_stream(subline_config_str);

                // initialize for iteration
                std::string value_config_str;  // stores parts of lines (subline)

                // iterate through each value
                while (std::getline(subline_config_stream, value_config_str, ','))
                {
                    bcq4s_sub.boundary_parameter_vec.push_back(std::stod(value_config_str));
                }

                // store BoundaryTypeData object
                bq4s.boundary_config_vec.push_back(bcq4s_sub);
                bcq4s_sub = bcq4s_new;  // clears bcq4s_sub

                // next line is not a BC value
                is_value_next = false;

            }

            // check if BC type or value will be given in next subline
            // "type": indicates that next entry is the BC type
            // "value": indicates that next entry is the BC value
            if (subline_config_str == "\"type\":")  
            {
                is_config_next = true;
            }
            if (subline_config_str == "\"value\":")  
            {
                is_value_next = true;
            }

        }

    }

    // close BC type file
    file_in_config_stream.close();

    return bq4s;

}

#endif
