/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-01 15:43:34
 * @LastEditTime: 2024-08-29 11:36:54
 * @FilePath: /cpgrid/example/twophase/convergence_test_1/main.C
 * @Description:
 *
 */
#include <petscvec.h>
#include <petscdm.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscdmtypes.h>
#include "utils.h"
#include "petscdmcp.h"
#include "vec_func.h"
#include "system.h"
#include "dof_map.h"
#include "ghosted_array_wrapper.h"
#include "fluid.h"
#include "mat_func.h"
#include "case_fluid.h"
extern PetscErrorCode FormFunction(SNES, Vec, Vec, void *);
extern PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void *);
extern PetscErrorCode FormInitialGuess(SNES snes, Vec x, void *ctx);
extern PetscErrorCode set_para(Para &para);
PetscErrorCode compute_errorvec(Vec &numerical_solution,
                                System &system, CaseFluid &fluid,
                                PetscReal time, Vec *error_vec);
PetscErrorCode compute_variable_errornorm(CaseFluid &fluid, PetscReal time, 
                                          Vec &error_vec, System &system, NormType norm_type, 
                                          PetscReal &p_errornorm, PetscReal &s_errornorm);
PetscErrorCode fill_exact_solution(Vec &exact_solution, Mesh &mesh, DofMap &dof_map, CaseFluid &fluid, PetscReal time, int pressure_var, int saturation_var);
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  // Mesh settings
  Mesh mesh;
  std::string case_name, mesh_name;
  char casename_char[256], meshname_char[256];
  PetscCall(PetscOptionsGetString(nullptr, nullptr, "-case_name", casename_char, sizeof(casename_char), nullptr));
  case_name = std::string(casename_char);
  PetscCall(PetscOptionsGetString(nullptr, nullptr, "-mesh_name", meshname_char, sizeof(meshname_char), nullptr));
  mesh_name = std::string(meshname_char);

  std::string filepath = "./grid/" + mesh_name;

  mesh.read(filepath);
  mesh.init();

  mesh.printInfo();
  mesh.partition();
  mesh.reset_elem_id();

  // System settings
  System sys(mesh);
  sys.add_variable("p");
  sys.add_variable("s");
  sys.init();

  sys.set_name(case_name);
  Para &para = sys.para();
  PetscCall(set_para(para));
  CaseFluid fluid(para);

  // DM (Distributed Mesh) settings
  DM dm_cp;
  PetscCall(DMRegister(DMCP, DMCreate_CP));
  PetscCall(DMCreate(PETSC_COMM_WORLD, &dm_cp));
  PetscCall(DMSetType(dm_cp, DMCP));

  // Associate System with DM
  PetscErrorCode (*setSystem)(DM, System &) = nullptr;
  PetscCall(PetscObjectQueryFunction((PetscObject)dm_cp, "DMCPSetSystem_C", &setSystem));
  PetscCall((*setSystem)(dm_cp, sys));
  PetscCall(DMSetFromOptions(dm_cp));
  PetscCall(DMSetUp(dm_cp));

  // SNES (Nonlinear solver) settings
  SNES snes;
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  PetscCall(SNESSetDM(snes, dm_cp));

  // Create global vectors
  Vec x, res;
  PetscCall(DMCreateGlobalVector(dm_cp, &x));
  PetscCall(VecDuplicate(x, &res));

  // Set initial guess and function
  // PetscCall(SNESSetComputeInitialGuess(snes, FormInitialGuess, &sys));
  PetscCall(FormInitialGuess(snes, x, &sys));
  PetscCall(SNESSetFunction(snes, res, FormFunction, &sys));
  SaveVecToMatlab(x, "init_s.m", "p_s");

  Mat j;
  sys.create_mat(j);
  PetscCall(SNESSetJacobian(snes, j, j, FormJacobian, NULL));
  /*
  SNESComputeJacobian(snes, x, j, j);
  SaveMatToMatlab(j, "analytic.m", "analytic_jacobian");
  SNESComputeJacobianDefault(snes, x, j, j, PETSC_NULL);
  SaveMatToMatlab(j, "fd.m", "fd_jacobian");
  */
  // Set SNES options from command line
  PetscCall(SNESSetFromOptions(snes));
  PetscCall(SNESSetUp(snes));

  // Time stepping
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscInt n_timesteps = para.get_parameter<PetscInt>("-n_steps");
  PetscInt t_order = para.get_parameter<PetscInt>("-time_order");
  PetscScalar t_end = dt * n_timesteps;
  PetscInt t_step = 0;
  char filename[256];
  sys.set_old_solution(x);

  ierr = SNESComputeFunction(snes,x,res);CHKERRQ(ierr);
  SaveVecToMatlab(x, "x.m", "x");
  SaveVecToMatlab(sys.get_old_solution(), "x_old.m", "x_old");
  SaveVecToMatlab(res, "res.m", "res");
  // return 0;

  SNESConvergedReason reason;
  while (sys.time < t_end)
  {
    sys.time += dt;
    t_step++;
    PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Solving time step %d, time = %g, dt = %g ***\n", t_step, sys.time, dt);
    // Save old solution and solve
    if (t_order == 2)
    {
      Vec x_old = sys.get_old_solution();
      sys.set_older_solution(x_old);
    }
    sys.set_old_solution(x);
    PetscCall(SNESSolve(snes, NULL, x));
    SaveVecToMatlab(x, "num_solution.m", "p_s");
    Vec error_vec;
    PetscCall(compute_errorvec(x, sys, fluid, sys.time, &error_vec));

    PetscReal p_errornorm, s_errornorm;
    compute_variable_errornorm(fluid, sys.time, error_vec, sys, NORM_2, p_errornorm, s_errornorm);
    PetscPrintf(PETSC_COMM_WORLD, "Pressure Error L_2 Norm: %.2e\nSaturation Error L_2 Norm: %.2e\n", p_errornorm, s_errornorm);
    compute_variable_errornorm(fluid, sys.time, error_vec, sys, NORM_INFINITY, p_errornorm, s_errornorm);
    PetscPrintf(PETSC_COMM_WORLD, "Pressure Error L_inf Norm: %.2e\nSaturation Error L_inf Norm: %.2e\n", p_errornorm, s_errornorm);

    PetscCall(SNESGetConvergedReason(snes, &reason));
    if (reason < 0)
    {
      PetscCall(PetscPrintf(PetscObjectComm((PetscObject)snes), "\t\t ***************** SNES diverged with reason %d at timestep %d.\n", (int)reason, t_step));
      PetscFunctionReturn(PETSC_SUCCESS);
    }

    // Get and save the solution
    PetscCall(SNESGetSolution(snes, &x));
    if (t_step % 1 == 0)
    {
      snprintf(filename, sizeof(filename), "./output/solution_%d.m", t_step);
      PetscCall(SaveVecToMatlab(x, filename, "p_s"));
    }
  }

  // Clean up resources
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&res));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm_cp));
  PetscCall(PetscFinalize());
  return 0;
}

// Function to fill the exact solution vector
PetscErrorCode fill_exact_solution(Vec &exact_solution, Mesh &mesh, DofMap &dof_map, CaseFluid &fluid, PetscReal time, int pressure_var, int saturation_var)
{
  PetscFunctionBeginUser;
  PetscScalar *exact_values;
  PetscCall(VecGetArray(exact_solution, &exact_values));

  for (const Polyhedron &element : mesh.elem_local_range())
  {
    int pressure_dof = dof_map.dof_local_indices(element, pressure_var);
    int saturation_dof = dof_map.dof_local_indices(element, saturation_var);
    PetscReal exact_pressure, exact_saturation;
    Point centroid = element.compute_cell_centroid();
    exact_saturation = fluid.Rsw(time, centroid.x(), centroid.y());
    exact_pressure = fluid.Rpw(time, centroid.x(), centroid.y());

    exact_values[pressure_dof] = exact_pressure;
    exact_values[saturation_dof] = exact_saturation;
  }
  PetscCall(VecRestoreArray(exact_solution, &exact_values));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// Main function to compute the error vector
PetscErrorCode compute_errorvec(Vec &numerical_solution, System &system, CaseFluid &fluid, PetscReal time, Vec *error_vec)
{
  PetscFunctionBeginUser;

  Mesh &mesh = system.get_mesh();
  DofMap &dof_map = system.dof_map();

  int pressure_var = dof_map.variable_num("p");
  int saturation_var = dof_map.variable_num("s");

  Vec exact_solution;
  PetscCall(VecDuplicate(numerical_solution, &exact_solution));
  PetscCall(VecDuplicate(numerical_solution, error_vec));

  PetscCall(fill_exact_solution(exact_solution, mesh, dof_map, fluid, time, pressure_var, saturation_var));
  SaveVecToMatlab(exact_solution, "exact_solution.m", "p_s");
  PetscCall(VecWAXPY(*error_vec, -1.0, numerical_solution, exact_solution));

  PetscCall(VecDestroy(&exact_solution));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode compute_variable_errornorm(CaseFluid &fluid, PetscReal time, Vec &error_vec, System &system, NormType norm_type, PetscReal &p_errornorm, PetscReal &s_errornorm)
{
  PetscFunctionBeginUser;

  // Get mesh and degrees of freedom (DoF) map
  Mesh   &mesh    = system.get_mesh();
  DofMap &dof_map = system.dof_map();

  // Identify variable numbers for pressure and saturation
  int pressure_var   = dof_map.variable_num("p");
  int saturation_var = dof_map.variable_num("s");

  // Vectors to store global indices for pressure and saturation variables
  std::vector<PetscInt> pressure_indices, saturation_indices;

  // Loop over each local element in the mesh and extract DoF indices
  for (const Polyhedron &element : mesh.elem_local_range()) {
    // Get DoF indices for pressure and saturation
    int pressure_dof   = dof_map.dof_indices(element, pressure_var);
    int saturation_dof = dof_map.dof_indices(element, saturation_var);

    pressure_indices.push_back(pressure_dof);
    saturation_indices.push_back(saturation_dof);
  }

  // Create IS (index sets) for pressure and saturation variables
  IS pressure_is, saturation_is;
  PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, pressure_indices.size(), pressure_indices.data(), PETSC_COPY_VALUES, &pressure_is));
  PetscCall(ISCreateGeneral(PETSC_COMM_WORLD, saturation_indices.size(), saturation_indices.data(), PETSC_COPY_VALUES, &saturation_is));

  // Compute exact solution
  Vec exact_solution;
  PetscCall(VecDuplicate(error_vec, &exact_solution));
  PetscCall(fill_exact_solution(exact_solution, mesh, dof_map, fluid, time, pressure_var, saturation_var));
  SaveVecToMatlab(exact_solution, "exact_solution.m", "p_s");

  // Extract subvectors for pressure and saturation error
  Vec p_error_vec, s_error_vec;
  PetscCall(VecGetSubVector(error_vec, pressure_is, &p_error_vec));
  PetscCall(VecGetSubVector(error_vec, saturation_is, &s_error_vec));
  SaveVecToMatlab(p_error_vec, "p_error_vec.m", "p_s");
  SaveVecToMatlab(s_error_vec, "s_error_vec.m", "p_s");

  Vec p_exact_vec, s_exact_vec;
  PetscCall(VecGetSubVector(exact_solution, pressure_is, &p_exact_vec));
  PetscCall(VecGetSubVector(exact_solution, saturation_is, &s_exact_vec));

  // Compute the desired norm for the pressure error vector
  PetscCall(VecNorm(p_error_vec, norm_type, &p_errornorm));
  PetscCall(VecNorm(s_error_vec, norm_type, &s_errornorm));

  // Compute the desired norm for the pressure exact vector
  PetscReal p_exact_norm, s_exact_norm;
  PetscCall(VecNorm(p_exact_vec, norm_type, &p_exact_norm));
  PetscCall(VecNorm(s_exact_vec, norm_type, &s_exact_norm));

  // Compute relative error
  p_errornorm = p_errornorm / p_exact_norm;
  s_errornorm = s_errornorm / s_exact_norm;

  // Output results
  // PetscPrintf(PETSC_COMM_WORLD, "Pressure Error Norm: %.2e\n", p_errornorm);
  // PetscPrintf(PETSC_COMM_WORLD, "Saturation Error Norm: %.2e\n", s_errornorm);

  // Restore the subvectors after use
  PetscCall(VecRestoreSubVector(error_vec, pressure_is, &p_error_vec));
  PetscCall(VecRestoreSubVector(error_vec, saturation_is, &s_error_vec));

  // Destroy the index sets
  PetscCall(ISDestroy(&pressure_is));
  PetscCall(ISDestroy(&saturation_is));

  PetscFunctionReturn(PETSC_SUCCESS);
}
