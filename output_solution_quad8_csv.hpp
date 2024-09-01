#ifndef OUTPUT_SOLUTION_QUAD8_CSV
#define OUTPUT_SOLUTION_QUAD8_CSV
#include <vector>
#include <fstream>
#include "Eigen/Eigen"
#include "grid_quad8.hpp"
#include "scalar_quad8.hpp"

void output_solution_quad8_csv(ScalarQuad8Class &sfq8, std::string file_out_base_str, int ts)
{

        // output file name
        std::string file_out_str = file_out_base_str + std::to_string(ts) + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridQuad8Struct gq8s = sfq8.gq8s;

        // write to file
        file_out_stream << "id,pos_x,pos_y,value\n";
        for (int n = 0; n < gq8s.num_point; n++)
        {
            file_out_stream << gq8s.point_id_vec[n] << "," << gq8s.point_pos_x_vec[n] << "," << gq8s.point_pos_y_vec[n] << "," << sfq8.scalar_vec[n] << "\n";
        }

}

void output_solution_quad8_csv(ScalarQuad8Class &sfq8, std::string file_out_base_str)
{

        // output file name
        std::string file_out_str = file_out_base_str + ".csv";

        // initialize file stream
        std::ofstream file_out_stream(file_out_str);

        // initialize grid object
        GridQuad8Struct gq8s = sfq8.gq8s;

        // write to file
        file_out_stream << "id,pos_x,pos_y,value\n";
        for (int n = 0; n < gq8s.num_point; n++)
        {
            file_out_stream << gq8s.point_id_vec[n] << "," << gq8s.point_pos_x_vec[n] << "," << gq8s.point_pos_y_vec[n] << "," << sfq8.scalar_vec[n] << "\n";
        }

}

#endif
