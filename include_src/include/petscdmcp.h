/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-24 09:33:54
 * @LastEditTime: 2024-08-02 18:22:38
 * @FilePath: /cpgrid/include/petscdmcp.h
 * @Description:
 *
 */
#ifndef PETSCDMCP_H
#define PETSCDMCP_H
#include <petscsys.h> // Include the basic PETSc system header
#include <petscdm.h>
#include "system.h"
#include "config.h"

PETSC_EXTERN PetscErrorCode DMCPSetSystem_CP(DM, System &);
PETSC_EXTERN PetscErrorCode DMCPGetSystem_CP(DM, System *&);
#define DMCP "cp"
PETSC_EXTERN PetscErrorCode DMCreate_CP(DM);
#endif // #ifdef PETSCDMCP
