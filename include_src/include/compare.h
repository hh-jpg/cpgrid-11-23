/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-31 09:30:32
 * @LastEditTime: 2024-08-28 10:15:09
 * @FilePath: /cpgrid/include/compare.h
 * @Description: some function used to compare results to mrst
 *
 */
#ifndef COMPARE_H
#define COMPARE_H
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscdmtypes.h>
#include "geo.h"
#include "config.h"
#include <algorithm>
#include <map>
#include <array>
#include <utility>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <iostream>

PetscErrorCode DataSaveASCII(Vec x, char *filename);
void TwoVecsToCSV(Vec vec1, Vec vec2, const std::string &filename);
PetscErrorCode SaveVecToMatlab(DM dm, Vec vec, const char *filename, const char *varname);
PetscErrorCode CompareVecDifference(Vec x, Vec y, PetscReal *diff_norm);

#endif