#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <fstream>

using namespace dealii; 

const int dim = 2; 

void make_grid(Triangulation<dim> &triangulation)
{
    Point<dim> center(1, 0); 
    const double inner_radius = .5, outer_radius = 1.0; 
    GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius, 5); // 5 cells ring

    // refine towards the inner ring: 
    for(unsigned int step = 0; step < 5; ++step)
    {
        for(const auto &cell : triangulation.active_cell_iterators())
        {
            for(const auto v : cell -> vertex_indices())
            {
                const double distance_from_center = center.distance(cell->vertex(v)); 
                if(std::fabs(distance_from_center - inner_radius) <= 1e-6 * inner_radius)
                {
                    cell->set_refine_flag(); 
                    break; 
                }
            }
        }
        triangulation.execute_coarsening_and_refinement(); 
    }

    std::string filename = "../mesh/mesh.gnuplot"; 
    std::ofstream mesh_file(filename); 
    GridOut().write_gnuplot(triangulation, mesh_file); 
}

void write_dof_locations(const DoFHandler<dim> &dof_handler, const std::string &filename)
{
    const std::map<types::global_dof_index, Point<2>> dof_location_map 
        = DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler); 
    
    std::ofstream dof_location_file(filename);
    DoFTools::write_gnuplot_dof_support_point_info(dof_location_file, dof_location_map); 
}

int main()
{
    Triangulation<dim> triangulation; 
    std::string img_file = "../imgs/mesh.svg"; 
    make_grid(triangulation);


    return 0; 
}