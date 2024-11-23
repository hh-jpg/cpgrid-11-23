/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-02 13:49:04
 * @LastEditTime: 2024-08-28 10:18:41
 * @FilePath: /cpgrid/include/vec_func.h
 * @Description:
 *
 */
#ifndef VEC_FUNC
#define VEC_FUNC
#include <petscsys.h>
#include <petscvec.h>
#include <vector>
#include <petscmat.h>
#include "config.h"
PetscErrorCode globalvec_to_local(Vec & global_vec, const std::vector<unsigned int> &ghosted_index, Vec *local_vec);
PetscErrorCode localvec_to_global(Vec& local_vec, Vec *global_vec);
PetscErrorCode SaveVecToMatlab(Vec& vec, const char *filename, const char *vecname);
#endif
