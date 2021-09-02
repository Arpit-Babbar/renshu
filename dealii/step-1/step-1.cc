/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2021 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 */


#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <cmath>

// From Step 49

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

using namespace dealii;

void first_grid()
{
  Triangulation<2> triangulation;

  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  std::ofstream out("grid-1.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);
  std::cout << "Grid written to grid-1.svg" << std::endl;
}

void second_grid()
{
  Triangulation<2> triangulation;
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  for (unsigned int step = 0; step < 5; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          for (const auto v : cell->vertex_indices())
            {
              const double distance_from_center =
                center.distance(cell->vertex(v));

              if (std::fabs(distance_from_center - inner_radius) <=
                  1e-6 * inner_radius)
                {
                  cell->set_refine_flag();
                  break;
                }
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }


  std::ofstream out("grid-2.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-2.svg" << std::endl;
}

// Exercise 1 - Different kind of refinement
void third_grid()
{
  Triangulation<2> triangulation;
  const Point<2> center(1, 0);
  const double   inner_radius = 0.5, outer_radius = 1.0;
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10);
  for (unsigned int step = 0; step < 10; ++step)
    {
      for (auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->center()[0] > 1.0 && cell->center()[1] > 0.0)
            {
              cell->set_refine_flag();
              break;
            }
        }

      triangulation.execute_coarsening_and_refinement();
    }

  std::ofstream out("grid-3.svg");
  GridOut       grid_out;
  grid_out.write_svg(triangulation, out);

  std::cout << "Grid written to grid-3.svg" << std::endl;
}

// Exercise 2 - Make an entirely new grid. Source -
// https://www.dealii.org/current/doxygen/deal.II/step_49.html#grid_5DemonstratingGridToolstransformpart1

template <int dim>
void print_mesh_info(const Triangulation<dim> &triangulation,
                     const std::string &       filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
  {
    std::map<types::boundary_id, unsigned int> boundary_count;
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary())
        boundary_count[face->boundary_id()]++;
    std::cout << " boundary indicators: ";
    for (const std::pair<const types::boundary_id, unsigned int> &pair :
         boundary_count)
      {
        std::cout << pair.first << "(" << pair.second << " times) ";
      }
    std::cout << std::endl;
  }
  std::ofstream out(filename);
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
  std::cout << " written to " << filename << std::endl << std::endl;
}


void fourth_grid()
{
  Triangulation<2>          triangulation;
  std::vector<unsigned int> repetitions(2);
  repetitions[0] = 14;
  repetitions[1] = 2;
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            Point<2>(-3.0,  0.0),
                                            Point<2>( 3.0,  1.0));
  // Writing a lambda that maps points to points. Source -
  // https://www.learncpp.com/cpp-tutorial/introduction-to-lambdas-anonymous-functions/
  auto trans
  {
    [](const Point<2> &in)
      {
        // return Point<2>(in[0], in[1] + std::sin(numbers::PI * in[0] / 5.0));
        return Point<2>(in[0], in[1] + in[0]*in[0]);
      }
  };
  GridTools::transform(
    trans,
    triangulation);
  print_mesh_info(triangulation, "grid-4.vtu");
}


int main()
{
  first_grid();
  second_grid();
  third_grid();
  // Comment from documentation - You can see that we used cell -> center here,
  // which is the second cell->something() we did here. In general, what you
  // can do with operations
  // of the form cell->something() is a bit difficult to find in the
  // documentation because cell is not a pointer but an iterator. The
  // functions you can call on a cell can be found in the documentation of
  // the classes TriaAccessor (which has functions that can also be called on
  // faces of cells or, more generally, all sorts of geometric objects that
  // appear in a triangulation), and CellAccessor (which adds a few
  // functions that are specific to cells).
  fourth_grid();
}
