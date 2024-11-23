/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-08 17:04:59
 * @LastEditTime: 2024-06-07 15:46:24
 * @FilePath: /cpgrid/include/grid.h
 * @Description:
 *
 */
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscdmtypes.h>
#include "geo.h"
#include <algorithm>
#include <map>
#include <array>
#include <utility>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <iostream>
#include "compare.h"
#include <unordered_set>
#include "config.h"

#ifndef GRID_H
#define GRID_H

typedef struct
{
  int nx, ny, nz;
} GridSize;

typedef struct
{
  DM da1;  // 有1个自由度的DMDA
  DM da3;  // 有3个自由度的DMDA
  DM da6;  // 有6个自由度的DMDA
  DM da18; // 有18个自由度的DMDA
  DM da24; // 有24个自由度的DMDA
} GridDMDAs;

typedef struct
{
  Vec face_areas;
  Vec face_normals;
  Vec face_centroids;
  Vec cell_centroid;
  Vec cell_volume;
  Vec nodes;
  Vec neighbor_indices;
  Vec perm;
  Vec htrans; // half transmissibilities
  Vec trans;
  Vec efftrans;
  Vec totmob;
} GridVectors;

typedef struct
{
  GridSize size;       // 网格尺寸
  GridDMDAs dmdas;     // 不同类型的 DMDA 对象
  GridVectors vectors; // 向量变量
} Grid;


PetscErrorCode cleanupGrid(Grid *grid);

PetscErrorCode initGrid(Grid **grid, const PetscInt nx, const PetscInt ny, const PetscInt nz);

PetscErrorCode LoadVectorFromBinary(const char *path, DM da, Vec *vec);

PetscErrorCode CreateLocalVectors(Grid *grid);

PetscErrorCode initializeVectors(Grid *grid);

PetscErrorCode initializeDMs(Grid *gridtemp, PetscInt nx, PetscInt ny, PetscInt nz);

#endif // POINT_H
