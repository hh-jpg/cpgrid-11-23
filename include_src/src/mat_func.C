/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-07 16:16:20
 * @LastEditTime: 2024-08-07 16:16:22
 * @FilePath: /cpgrid/src/mat_func.C
 * @Description: 
 * 
 */
#include "mat_func.h"
PetscErrorCode SaveMatToMatlab(Mat j, const char *filename, const char *varname)
{
  PetscViewer viewer;
  PetscFunctionBegin;
  PetscCall(PetscObjectSetName((PetscObject)j, varname));
  PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer));
  PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
  PetscCall(MatView(j, viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}