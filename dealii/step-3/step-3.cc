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

 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 */

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <assert.h>

// Exercise - boundary function
#include <deal.II/base/function_lib.h>
// Exercise - Adding point and mean into .h5 file
#include <deal.II/base/hdf5.h>
using namespace dealii;

class Step3
{
public:
  Step3(std::string grid_indicator, unsigned int n_refinements);

  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  Triangulation<2> triangulation;
  FE_Q<2> fe;
  DoFHandler<2> dof_handler;

  SparsityPattern sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;

  std::string grid_indicator;
  unsigned int n_refinements;
};

Step3::Step3(std::string grid_indicator, unsigned int n_refinements)
    : fe(1), dof_handler(triangulation), grid_indicator(grid_indicator),
      n_refinements(n_refinements)
{}

void Step3::make_grid()
{
  if (grid_indicator == "hyper_cube")
  {
      GridGenerator::hyper_cube(triangulation, -1, 1);
      triangulation.begin_active()->face(0)->set_boundary_id(1); // 0 neumann
  }
  else if (grid_indicator == "hyper_shell")
  {
    const Point<2> center(1, 0);
    const double   inner_radius = 0.2, outer_radius = 1.0;
    GridGenerator::hyper_shell(
      triangulation, center, inner_radius, outer_radius, 5);
    for (auto &face : triangulation.active_face_iterators())
      if (sqrt(pow(face->center()(0),2)+pow(face->center()(1),2))
          < 0.2 + 1e-12)
      face->set_boundary_id(1); // 0 neumann bc
  }
  else if (grid_indicator == "simplex")
  {
    const int dim = 2;
    std::vector<Point<dim>> vertices{{0.0,0.0},{1.0,0.0},{0.0,1.0}};
    GridGenerator::simplex(triangulation, vertices);
    triangulation.begin_active()->face(0)->set_boundary_id(1); // 0 neumann
  }
  else
    assert(false && "Incorrect Grid indicator");
  triangulation.refine_global(n_refinements);
}

void Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}

void Step3::assemble_system()
{
  QGauss<2> quadrature_formula(fe.degree + 1);
  FEValues<2> fe_values(fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Boundary function
  const Point<2> direction{1.0,0.0};
  const double steepness = 1.0;
  Functions::JumpFunction<2> boundary_data(direction, steepness);

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);

    cell_matrix = 0;
    cell_rhs = 0;

    for (const unsigned int q_index : fe_values.quadrature_point_indices())
    {
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) +=
              (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
               fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
               fe_values.JxW(q_index));           // dx

      for (const unsigned int i : fe_values.dof_indices())
        cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                        1. *                                // f(x_q)
                        fe_values.JxW(q_index));            // dx
    }
    cell->get_dof_indices(local_dof_indices);

    for (const unsigned int i : fe_values.dof_indices())
      for (const unsigned int j : fe_values.dof_indices())
        system_matrix.add(local_dof_indices[i],
                          local_dof_indices[j],
                          cell_matrix(i, j));

    for (const unsigned int i : fe_values.dof_indices())
      system_rhs(local_dof_indices[i]) += cell_rhs(i);
  }

  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,               // boundary id
                                           boundary_data,
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}

void Step3::solve()
{
  SolverControl solver_control(10000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);

  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}

void Step3::output_results() const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
  // Writing hdf5 file (Not working, output is crazy)
  const std::string filename_h5 = "solution_" + std::to_string(n_refinements) + ".h5";
  DataOutBase::DataOutFilterFlags flags(true, true);
  DataOutBase::DataOutFilter data_filter(flags);
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter, filename_h5, MPI_COMM_WORLD);

  // Printing value at a point to check convergence
  std::cout << "Solution at (1/3,1/3): "
            << VectorTools::point_value(dof_handler, solution,
                                        Point<2>(1./3, 1./3))
            << std::endl;

  // Printing Mean value to check convergence
  std::cout << "Mean value: "
            << VectorTools::compute_mean_value(dof_handler,
                                               QGauss<2>(fe.degree + 1),
                                               solution, 0)
            << std::endl;
  // Adding a point value and mean value to HDF5 file(visualization incomplete)
  HDF5::File data_file(filename_h5, HDF5::File::FileAccessMode::open, MPI_COMM_WORLD);
  Vector<double> point_value(1);
  point_value[0] = VectorTools::point_value(dof_handler,
                                            solution, Point<2>(1. / 3, 1. / 3));
  data_file.write_dataset("point_value", point_value);
  Vector<double> mean_value(1);
  mean_value[0] = VectorTools::compute_mean_value(dof_handler,
                                                  QGauss<2>(fe.degree + 1),
                                                  solution, 0);
  data_file.write_dataset("mean_value", mean_value);
}

void Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}

int main(int argc, char **argv)
{
  // Initialize MPI Communicator with 1 processor
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  deallog.depth_console(2);

  if (argc != 3)
  {
    cout << "Run as ./step-3 grid_indicator n_refinements "
         << "Choices for grid_indicator are hyper_cube, hyper_shell.";
    assert(false);
  }

  std::string grid_indicator = argv[1];
  unsigned int n_refinements = atoi(argv[2]);
  cout << "Grid_indicator = "<<grid_indicator << "\n";
  Step3 laplace_problem(grid_indicator, n_refinements);
  laplace_problem.run();

  return 0;
}
