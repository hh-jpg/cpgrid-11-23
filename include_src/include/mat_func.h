#ifndef MAT_FUNC
#define MAT_FUNC
#include <petscsys.h>
#include <petscmat.h>
#include "config.h"
PetscErrorCode SaveMatToMatlab(Mat j, const char *filename, const char *varname);
#endif
