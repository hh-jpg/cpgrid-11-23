/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-13 16:14:00
 * @LastEditTime: 2024-08-28 10:17:55
 * @FilePath: /cpgrid/include/petsc_utils.h
 * @Description: 
 * 
 */
#include <petscversion.h>
#include "config.h"

#define PETSC_VERSION_LT(MAJOR, MINOR, SUBMINOR) \
    (PETSC_VERSION_LT_(MAJOR, MINOR, SUBMINOR, PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR))

#define PETSC_VERSION_LT_(REQ_MAJOR, REQ_MINOR, REQ_SUBMINOR, CUR_MAJOR, CUR_MINOR, CUR_SUBMINOR) \
    ((CUR_MAJOR < REQ_MAJOR) || \
    (CUR_MAJOR == REQ_MAJOR && CUR_MINOR < REQ_MINOR) || \
    (CUR_MAJOR == REQ_MAJOR && CUR_MINOR == REQ_MINOR && CUR_SUBMINOR < REQ_SUBMINOR))

// Example usage:
// #if PETSC_VERSION_LT(3, 12, 0)
//     // Code for PETSc versions lower than 3.12.0
// #endif
