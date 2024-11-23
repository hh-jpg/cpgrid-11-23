/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-06 16:21:56
 * @LastEditTime: 2024-06-14 10:16:05
 * @FilePath: /cpgrid/src/grid.C
 * @Description:
 *
 */
#include "grid.h"
#include "compare.h"

PetscErrorCode initGrid(Grid **grid, const PetscInt nx, const PetscInt ny, const PetscInt nz)
{
  PetscErrorCode ierr; // To handle error codes from PETSc functions
  Grid *gridtemp;
  PetscFunctionBegin;
  // Allocate memory for the Grid structure
  gridtemp = (Grid *)malloc(sizeof(Grid));
  gridtemp->size.nx = nx;
  gridtemp->size.ny = ny;
  gridtemp->size.nz = nz;
  // Initialize the DMs
  PetscCall(initializeDMs(gridtemp, nx, ny, nz));
  // Initialize vectors
  PetscCall(initializeVectors(gridtemp));
  // Assign the newly created Grid to the output pointer
  *grid = gridtemp;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode initializeDMs(Grid *grid, PetscInt nx, PetscInt ny, PetscInt nz)
{
  PetscFunctionBegin;
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                         nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &grid->dmdas.da1));
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                         nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 3, 1, NULL, NULL, NULL, &grid->dmdas.da3));
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                         nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 6, 1, NULL, NULL, NULL, &grid->dmdas.da6));
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                         nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 18, 1, NULL, NULL, NULL, &grid->dmdas.da18));
  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                         nx, ny, nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 24, 1, NULL, NULL, NULL, &grid->dmdas.da24));
  PetscCall(DMSetFromOptions(grid->dmdas.da1));
  PetscCall(DMSetUp(grid->dmdas.da1));
  PetscCall(DMSetFromOptions(grid->dmdas.da3));
  PetscCall(DMSetUp(grid->dmdas.da3));
  PetscCall(DMSetFromOptions(grid->dmdas.da6));
  PetscCall(DMSetUp(grid->dmdas.da6));
  PetscCall(DMSetFromOptions(grid->dmdas.da18));
  PetscCall(DMSetUp(grid->dmdas.da18));
  PetscCall(DMSetFromOptions(grid->dmdas.da24));
  PetscCall(DMSetUp(grid->dmdas.da24));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode LoadVectorFromBinary(const char *path, DM da, Vec *vec)
{
  PetscViewer viewer;
  Vec g_vector;
  PetscCall(DMCreateGlobalVector(da, &g_vector));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path, FILE_MODE_READ, &viewer));
  PetscCall(VecLoad(g_vector, viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(DMGlobalToLocalBegin(da, g_vector, INSERT_VALUES, *vec));
  PetscCall(DMGlobalToLocalEnd(da, g_vector, INSERT_VALUES, *vec));
  PetscCall(VecDestroy(&g_vector));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CreateLocalVectors(Grid *grid)
{
  PetscFunctionBegin;
  PetscCall(DMCreateLocalVector(grid->dmdas.da6, &grid->vectors.face_areas));
  int g_size, l_size;
  VecGetSize(grid->vectors.face_areas, &g_size);
  VecGetLocalSize(grid->vectors.face_areas, &l_size);
  PetscPrintf(PETSC_COMM_WORLD, "gsize is %d, local size %d \n", g_size, l_size);
 
   Vec globalVec;
  PetscCall(DMCreateGlobalVector(grid->dmdas.da6, &globalVec));
  // 获取全局和本地大小
  PetscCall(VecGetSize(globalVec, &g_size));
  PetscCall(VecGetLocalSize(globalVec, &l_size));

  // 打印信息
  PetscPrintf(PETSC_COMM_WORLD, "Global size is %d, local size is %d\n", g_size, l_size);

  PetscCall(DMCreateLocalVector(grid->dmdas.da18, &grid->vectors.face_normals));
  PetscCall(DMCreateLocalVector(grid->dmdas.da18, &grid->vectors.face_centroids));
  PetscCall(DMCreateLocalVector(grid->dmdas.da3, &grid->vectors.cell_centroid));
  PetscCall(DMCreateLocalVector(grid->dmdas.da1, &grid->vectors.cell_volume));
  PetscCall(DMCreateLocalVector(grid->dmdas.da24, &grid->vectors.nodes));
  PetscCall(DMCreateLocalVector(grid->dmdas.da6, &grid->vectors.neighbor_indices));
  PetscCall(DMCreateLocalVector(grid->dmdas.da1, &grid->vectors.perm));
  PetscCall(DMCreateLocalVector(grid->dmdas.da6, &grid->vectors.htrans));
  PetscCall(DMCreateLocalVector(grid->dmdas.da6, &grid->vectors.trans));
  PetscCall(DMCreateLocalVector(grid->dmdas.da6, &grid->vectors.efftrans));
  PetscCall(DMCreateLocalVector(grid->dmdas.da1, &grid->vectors.totmob));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode initializeVectors(Grid *grid)
{
  PetscFunctionBegin;
  PetscCall(CreateLocalVectors(grid));
  char path_nodes[100] = "./node.bin";
  PetscCall(PetscOptionsGetString(NULL, NULL, "-nodes_file", path_nodes, sizeof path_nodes, NULL));
  PetscCall(LoadVectorFromBinary(path_nodes, grid->dmdas.da24, &grid->vectors.nodes));

  // 创建全局向量
  Vec global_face_areas_vec, global_face_normals_vec, global_face_centroids_vec;
  Vec global_cell_centroid_vec, global_cell_volume_vec;

  PetscCall(DMCreateGlobalVector(grid->dmdas.da6, &global_face_areas_vec));
  PetscCall(DMCreateGlobalVector(grid->dmdas.da18, &global_face_normals_vec));
  PetscCall(DMCreateGlobalVector(grid->dmdas.da18, &global_face_centroids_vec));
  PetscCall(DMCreateGlobalVector(grid->dmdas.da3, &global_cell_centroid_vec));
  PetscCall(DMCreateGlobalVector(grid->dmdas.da1, &global_cell_volume_vec));

  HexNodes ***cell_nodes_array;
  PetscCall(DMDAVecGetArray(grid->dmdas.da24, grid->vectors.nodes, &cell_nodes_array));

  HexFaces ***face_areas_array;
  HexFaceNormals ***face_normals_array;
  HexFaceCentroids ***face_centroids_array;
  Point ***cell_centroid_array;
  PetscScalar ***cell_volume_array;

  PetscCall(DMDAVecGetArray(grid->dmdas.da6, global_face_areas_vec, &face_areas_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da18, global_face_normals_vec, &face_normals_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da18, global_face_centroids_vec, &face_centroids_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da3, global_cell_centroid_vec, &cell_centroid_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da1, global_cell_volume_vec, &cell_volume_array));

  // Process local array data
  int xl, yl, zl, nxl, nyl, nzl, mx, my, mz;
  PetscCall(DMDAGetCorners(grid->dmdas.da6, &xl, &yl, &zl, &nxl, &nyl, &nzl));
  PetscCall(DMDAGetInfo(grid->dmdas.da6, 0, &mx, &my, &mz, 0, 0, 0, 0, 0, 0, 0, 0, 0));
  for (int k = zl; k < zl + nzl; k++)
  {
    for (int j = yl; j < yl + nyl; j++)
    {
      for (int i = xl; i < xl + nxl; i++)
      {
        HexNodes current_cell_nodes = cell_nodes_array[k][j][i];
        face_areas_array[k][j][i] = faceareas_setup(current_cell_nodes);
        face_normals_array[k][j][i] = facenormals_setup(current_cell_nodes);
        face_centroids_array[k][j][i] = facecentroids_setup(current_cell_nodes);
        CellCentroidVol cv = calculateHexCellCV(current_cell_nodes);
        cell_centroid_array[k][j][i] = cv.centroid;
        cell_volume_array[k][j][i] = cv.volume;
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da6, global_face_areas_vec, &face_areas_array));
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da18, global_face_normals_vec, &face_normals_array));
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da18, global_face_centroids_vec, &face_centroids_array));
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da3, global_cell_centroid_vec, &cell_centroid_array));
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da1, global_cell_volume_vec, &cell_volume_array));
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da24, grid->vectors.nodes, &cell_nodes_array));

  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da6, global_face_areas_vec, INSERT_VALUES, grid->vectors.face_areas));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da6, global_face_areas_vec, INSERT_VALUES, grid->vectors.face_areas));
  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da18, global_face_normals_vec, INSERT_VALUES, grid->vectors.face_normals));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da18, global_face_normals_vec, INSERT_VALUES, grid->vectors.face_normals));
  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da18, global_face_centroids_vec, INSERT_VALUES, grid->vectors.face_centroids));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da18, global_face_centroids_vec, INSERT_VALUES, grid->vectors.face_centroids));
  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da3, global_cell_centroid_vec, INSERT_VALUES, grid->vectors.cell_centroid));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da3, global_cell_centroid_vec, INSERT_VALUES, grid->vectors.cell_centroid));
  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da1, global_cell_volume_vec, INSERT_VALUES, grid->vectors.cell_volume));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da1, global_cell_volume_vec, INSERT_VALUES, grid->vectors.cell_volume));

  // 销毁全局向量
  PetscCall(VecDestroy(&global_face_areas_vec));
  PetscCall(VecDestroy(&global_face_normals_vec));
  PetscCall(VecDestroy(&global_face_centroids_vec));
  PetscCall(VecDestroy(&global_cell_centroid_vec));
  PetscCall(VecDestroy(&global_cell_volume_vec));
  // PetscCall(DMDAVecRestoreArray(grid->dmdas.da24, grid->vectors.nodes, &cell_nodes_array));
  PetscCall(find_neighbor(grid));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode cleanupGrid(Grid *grid)
{
  PetscFunctionBegin;
  // 销毁 Vec 对象
  if (grid->vectors.face_areas)
  {
    VecDestroy(&grid->vectors.face_areas);
  }
  if (grid->vectors.face_normals)
  {
    VecDestroy(&grid->vectors.face_normals);
  }
  if (grid->vectors.face_centroids)
  {
    VecDestroy(&grid->vectors.face_centroids);
  }
  if (grid->vectors.cell_centroid)
  {
    VecDestroy(&grid->vectors.cell_centroid);
  }
  if (grid->vectors.cell_volume)
  {
    VecDestroy(&grid->vectors.cell_volume);
  }
  if (grid->vectors.nodes)
  {
    VecDestroy(&grid->vectors.nodes);
  }
  if (grid->vectors.neighbor_indices)
  {
    VecDestroy(&grid->vectors.neighbor_indices);
  }
  if (grid->vectors.perm)
  {
    VecDestroy(&grid->vectors.perm);
  }
  if (grid->vectors.trans)
  {
    VecDestroy(&grid->vectors.trans);
  }
  if (grid->vectors.htrans)
  {
    VecDestroy(&grid->vectors.htrans);
  }
  if (grid->vectors.efftrans)
  {
    VecDestroy(&grid->vectors.efftrans);
  }
  if (grid->vectors.totmob)
  {
    VecDestroy(&grid->vectors.totmob);
  }
  // 销毁 DM 对象
  if (grid->dmdas.da1)
  {
    DMDestroy(&grid->dmdas.da1);
  }
  if (grid->dmdas.da3)
  {
    DMDestroy(&grid->dmdas.da3);
  }
  if (grid->dmdas.da6)
  {
    DMDestroy(&grid->dmdas.da6);
  }
  if (grid->dmdas.da18)
  {
    DMDestroy(&grid->dmdas.da18);
  }
  if (grid->dmdas.da24)
  {
    DMDestroy(&grid->dmdas.da24);
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
