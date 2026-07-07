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

#include <iostream>
#include <fstream>

using namespace dealii; 

const int dim = 2; 

void make_grid(Triangulation<dim> &triangulation, int n_cells, int n_steps)
{
    Point<dim> center(1, 0); 
    const double inner_radius = .5, outer_radius = 1.0; 
    GridGenerator::hyper_shell(triangulation, center, inner_radius, outer_radius, n_cells); // 5 cells ring

    // refine towards the inner ring: 
    for(int step = 0; step < n_steps; ++step)
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

void visualize(Triangulation<dim> &triangulation)
{
    std::ofstream output_file("../imgs/mesh.svg"); 
    GridOut gridout; 
    gridout.write_svg(triangulation, output_file); 
}

void write_dof_locations(const DoFHandler<dim> &dof_handler, const std::string &filename)
{
    const std::map<types::global_dof_index, Point<2>> dof_location_map 
        = DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler); 
    
    std::ofstream dof_location_file(filename);
    DoFTools::write_gnuplot_dof_support_point_info(dof_location_file, dof_location_map); 
}

void distribute_dofs(DoFHandler<dim> &dof_handler, int order)
{
    const FE_Q<2> finite_element(order);  
    dof_handler.distribute_dofs(finite_element); 
    write_dof_locations(dof_handler, "dof-locations-1.gnuplot");    

    DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(), dof_handler.n_dofs()); // useful in order not to have waste of memory with low overhead
    DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern); // fill the object with the non-zero elements
    SparsityPattern sparsity_pattern; 
    sparsity_pattern.copy_from(dynamic_sparsity_pattern); // create the sparsity pattern 

    std::ofstream out("../imgs/sparsity-pattern-1.svg"); 
    sparsity_pattern.print_svg(out); 
}

  int main(int argc, char** argv)
  {

    if(argc < 3)
    {
        std::cerr << "Correct usage: ./run <int> <int> <int>\n" << std::endl; 
        return 0; 
    }

    int order = atoi(argv[1]); 
    int n_cells = atoi(argv[2]); 
    int n_steps = atoi(argv[3]); 

    Triangulation<2> triangulation;
    make_grid(triangulation, n_cells, n_steps);
  
    DoFHandler<2> dof_handler(triangulation);
  
    distribute_dofs(dof_handler, order);

    visualize(triangulation); 
  }