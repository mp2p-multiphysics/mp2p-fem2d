#ifndef OUTPUT_SOLUTION_QUAD4_CSV
#define OUTPUT_SOLUTION_QUAD4_CSV
#include <vector>
#include <fstream>
#include "Eigen/Eigen"
#include "grid_quad4.hpp"
#include "scalar_quad4.hpp"

void output_solution_quad4_csv(ScalarQuad4Class &sfq4, std::string file_out_base_str)
{

        // output file name
        std::string file_out_str = file_out_base_str + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridQuad4Struct gq4s = sfq4.gq4s;

        // write to file
        file_out_stream << "id,pos_x,pos_y,value\n";
        for (int n = 0; n < gq4s.num_point; n++)
        {
            file_out_stream << gq4s.point_id_vec[n] << "," << gq4s.point_pos_x_vec[n] << "," << gq4s.point_pos_y_vec[n] << "," << sfq4.scalar_vec[n] << "\n";
        }

}

void output_solution_quad4_csv(ScalarQuad4Class &sfq4, std::string file_out_base_str, int ts)
{

        // output file name
        std::string file_out_str = file_out_base_str + std::to_string(ts) + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridQuad4Struct gq4s = sfq4.gq4s;

        // write to file
        file_out_stream << "id,pos_x,pos_y,value\n";
        for (int n = 0; n < gq4s.num_point; n++)
        {
            file_out_stream << gq4s.point_id_vec[n] << "," << gq4s.point_pos_x_vec[n] << "," << gq4s.point_pos_y_vec[n] << "," << sfq4.scalar_vec[n] << "\n";
        }

}

#endif
