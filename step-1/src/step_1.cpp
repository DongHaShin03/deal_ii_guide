#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii; 

void visualize_grid(Triangulation<2> triangulation, std::string filename)
{
    std::ofstream output_file(filename); 
    GridOut gridout; 
    gridout.write_svg(triangulation, output_file); 
}

Triangulation<2> create_2D_grid(int number_of_refinements)
{
    Triangulation<2> triangulation; 
    GridGenerator::hyper_cube(triangulation); 
    triangulation.refine_global(number_of_refinements); 

    return triangulation; 
}

Triangulation<2> create_coarse_ring(int number_of_cells)
{
    Triangulation<2> ring; 
    const Point<2> center(1, 0); 
    double in_radius = 0.5, out_radius = 1.0; 

    GridGenerator::hyper_shell(ring, center, in_radius, out_radius, number_of_cells); 
    return ring; 
}

Triangulation<2> create_refined_ring(int number_of_cells, unsigned int num_steps)
{
    Triangulation<2> ring; 
    const Point<2> center(1, 0); 
    double in_radius = 0.5, out_radius = 1.0; 

    GridGenerator::hyper_shell(ring, center, in_radius, out_radius, number_of_cells); 

    for(unsigned int step = 0; step < num_steps; ++step)
    {
        // Loop over cells
        for(const auto &cell : ring.active_cell_iterators())
        {
            // Loop over vertices of each cell
            for(const auto v : cell->vertex_indices())
            {
                const double distance_from_center = center.distance(cell->vertex(v)); 

                if(std::fabs(distance_from_center - in_radius) < 1e-16)
                {
                    cell->set_refine_flag(); 
                    break; 
                }
            }
        }
        
        ring.execute_coarsening_and_refinement(); 
    }

    return ring; 
}

int main()
{
    visualize_grid(create_2D_grid(/* refinementsz*/ 4), "../imgs/grid_1.svg"); 
    visualize_grid(create_coarse_ring(/*num_cells*/ 5), "../imgs/coarse_ring.svg"); 
    visualize_grid(create_refined_ring(/*num_cells*/ 10, /* steps */ 5), "../imgs/fine_ring.svg"); 

    return 0; 
}
