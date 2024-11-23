/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-20 15:17:03
 * @LastEditTime: 2024-08-06 15:14:38
 * @FilePath: /cpgrid/src/trans.C
 * @Description:
 *
 */
#include "grid.h"
#include "utils.h"

PetscErrorCode tpfa_htrans_compute(Grid *grid)
{
  PetscFunctionBeginUser;
  HexFaces ***faces_trans_array;
  HexFaceNormals ***faces_normals_array;
  HexFaces ***faces_areas_array;
  HexFaceCentroids ***faces_centroids_array;
  Point ***cell_centroid_array;
  PetscReal ***perm_array;
  Vec global_htrans;
  PetscCall(DMDAVecGetArray(grid->dmdas.da3, grid->vectors.cell_centroid, &cell_centroid_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da18, grid->vectors.face_normals, &faces_normals_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, grid->vectors.face_areas, &faces_areas_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da18, grid->vectors.face_normals, &faces_normals_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da18, grid->vectors.face_centroids, &faces_centroids_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da1, grid->vectors.perm, &perm_array));

  PetscCall(DMCreateGlobalVector(grid->dmdas.da6, &global_htrans));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, global_htrans, &faces_trans_array));
  // Process local array data
  int xl, yl, zl, nxl, nyl, nzl;
  PetscCall(DMDAGetCorners(grid->dmdas.da6, &xl, &yl, &zl, &nxl, &nyl, &nzl));
  for (int k = zl; k < zl + nzl; k++)
  {
    for (int j = yl; j < yl + nyl; j++)
    {
      for (int i = xl; i < xl + nxl; i++)
      {
        for (int f = 0; f < 6; f++)
        {
          Point cc2fc = faces_centroids_array[k][j][i].centroids[f] - cell_centroid_array[k][j][i];
          faces_trans_array[k][j][i].face[f] = faces_areas_array[k][j][i].face[f] * perm_array[k][j][i] *
                                               cc2fc.dot(faces_normals_array[k][j][i].normals[f]) / (cc2fc.two_norm() * cc2fc.two_norm());
        }
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da6, global_htrans, &faces_trans_array));
  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da6, global_htrans, INSERT_VALUES, grid->vectors.htrans));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da6, global_htrans, INSERT_VALUES, grid->vectors.htrans));
  PetscCall(VecDestroy(&global_htrans));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode tpfa_trans_compute(Grid *grid)
{
  PetscFunctionBeginUser;
  HexFaces ***htrans_array, ***trans_array, ***neighbor_indices_array;
  Vec global_trans;

  PetscCall(DMDAVecGetArray(grid->dmdas.da6, grid->vectors.htrans, &htrans_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, grid->vectors.neighbor_indices, &neighbor_indices_array));

  PetscCall(DMCreateGlobalVector(grid->dmdas.da6, &global_trans));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, global_trans, &trans_array));

  int xl, yl, zl, nxl, nyl, nzl;
  int ni, nj, nk, nf;
  PetscCall(DMDAGetCorners(grid->dmdas.da6, &xl, &yl, &zl, &nxl, &nyl, &nzl));
  PetscScalar cht, nht;
  for (int k = zl; k < zl + nzl; k++)
  {
    for (int j = yl; j < yl + nyl; j++)
    {
      for (int i = xl; i < xl + nxl; i++)
      {
        for (int f = 0; f < 6; f++)
        {
          PetscInt neighbor_id = PetscInt(neighbor_indices_array[k][j][i].face[f]);
          if (neighbor_id == -1)
          {
            trans_array[k][j][i].face[f] = htrans_array[k][j][i].face[f];
          }
          else
          {
            index_to_ijk(neighbor_id, grid->size, ni, nj, nk);
            nf = neighbor_face[f];
            cht = htrans_array[k][j][i].face[f];
            nht = htrans_array[nk][nj][ni].face[nf];
            trans_array[k][j][i].face[f] = half_harmonicMean(cht, nht);
          }
        }
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da6, global_trans, &trans_array));
  PetscCall(DMGlobalToLocal(grid->dmdas.da6, global_trans, INSERT_VALUES, grid->vectors.trans));
  PetscCall(VecDestroy(&global_trans));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode tpfa_efftrans_compute(Grid *grid)
{
  PetscFunctionBeginUser;
  HexFaces ***htrans_array, ***neighbor_indices_array, ***efftrans_array;
  PetscReal ***totmob_array;
  Vec global_efftrans;
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, grid->vectors.htrans, &htrans_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da1, grid->vectors.totmob, &totmob_array));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, grid->vectors.neighbor_indices, &neighbor_indices_array));

  PetscCall(DMCreateGlobalVector(grid->dmdas.da6, &global_efftrans));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, global_efftrans, &efftrans_array));

  // Process local array data
  int xl, yl, zl, nxl, nyl, nzl, ni, nj, nk, nf;
  PetscScalar c_efhtrans, n_efhtrans;
  PetscCall(DMDAGetCorners(grid->dmdas.da6, &xl, &yl, &zl, &nxl, &nyl, &nzl));
  for (int k = zl; k < zl + nzl; k++)
  {
    for (int j = yl; j < yl + nyl; j++)
    {
      for (int i = xl; i < xl + nxl; i++)
      {
        for (int f = 0; f < 6; f++)
        {
          PetscInt neighbor_id = PetscInt(neighbor_indices_array[k][j][i].face[f]);
          if (neighbor_id == -1)
          {
            efftrans_array[k][j][i].face[f] = totmob_array[k][j][i] * htrans_array[k][j][i].face[f];
          }
          else
          {
            index_to_ijk(neighbor_id, grid->size, ni, nj, nk);
            nf = neighbor_face[f];
            c_efhtrans = htrans_array[k][j][i].face[f] * totmob_array[k][j][i];
            n_efhtrans = htrans_array[nk][nj][ni].face[nf] * totmob_array[nk][nj][ni];
            efftrans_array[k][j][i].face[f] = half_harmonicMean(c_efhtrans, n_efhtrans);
          }
        }
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da6, global_efftrans, &efftrans_array));
  PetscCall(DMGlobalToLocal(grid->dmdas.da6, global_efftrans, INSERT_VALUES, grid->vectors.efftrans));
  PetscCall(VecDestroy(&global_efftrans));
  PetscFunctionReturn(PETSC_SUCCESS);
}