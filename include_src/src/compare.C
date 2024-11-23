/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-31 09:33:25
 * @LastEditTime: 2024-06-21 16:20:34
 * @FilePath: /cpgrid/src/compare.C
 * @Description:
 *
 */
#include "compare.h"

PetscErrorCode DataSaveASCII(Vec x, char *filename)
{
  PetscErrorCode ierr;
  PetscViewer viewer;

  PetscFunctionBegin;
  // ierr = PetscViewerHDF5Open(comm,filename,FILE_MODE_WRITE,&viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
  CHKERRQ(ierr);
  // format: PETSC_VIEWER_ASCII_SYMMODU, PETSC_VIEWER_ASCII_MATLAB, PETSC_VIEWER_ASCII_COMMON
  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_SYMMODU);
  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);
  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode CompareVecDifference(Vec x, Vec y, PetscReal *diff_norm)
{
  PetscFunctionBegin;
  Vec diff;
  // Ensure x and y are of the same size
  PetscInt size_x, size_y;
  PetscCall(VecGetSize(x, &size_x));
  PetscCall(VecGetSize(y, &size_y));
  if (size_x != size_y)
  {
    PetscPrintf(PETSC_COMM_WORLD, "size_x is %d, size_y is %d", size_x, size_y);
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Vectors must be of the same size");
  }
  // Create a vector to hold the difference
  PetscCall(VecDuplicate(x, &diff));

  // Compute the difference: diff = x - y
  PetscCall(VecWAXPY(diff, -1.0, y, x));

  // Compute the L2 norm of the difference
  PetscCall(VecNorm(diff, NORM_2, diff_norm));

  // 获取向量的数据指针
  const PetscScalar *array;
  PetscCall(VecGetArrayRead(diff, &array));
  // 找到最大绝对值及其位置
  double maxVal = -1.0;
  PetscInt maxIndex, n;
  PetscCall(VecGetSize(diff, &n));
  for (PetscInt i = 0; i < n; ++i)
  {
    if (PetscAbsScalar(array[i]) > maxVal)
    {
      maxVal = PetscAbsScalar(array[i]);
      maxIndex = i;
    }
  }
  PetscCall(VecRestoreArrayRead(diff, &array));
  PetscPrintf(PETSC_COMM_WORLD, "Max absolute value in vector: %g at index %d\n", maxVal, maxIndex);

  // Clean up
  PetscCall(VecDestroy(&diff));
  PetscFunctionReturn(PETSC_SUCCESS);
}

void TwoVecsToCSV(Vec vec1, Vec vec2, const std::string &filename)
{
  PetscInt size1, size2;
  const PetscScalar *array1, *array2;

  // 获取 Vec 的大小
  VecGetSize(vec1, &size1);
  VecGetSize(vec2, &size2);

  // 确保两个 Vec 的大小相同
  if (size1 != size2)
  {
    std::cerr << "Error: Vec sizes do not match." << std::endl;
    return;
  }

  // 获取 Vec 的数据指针
  VecGetArrayRead(vec1, &array1);
  VecGetArrayRead(vec2, &array2);

  // 创建输出文件流
  std::ofstream outfile(filename + ".csv");

  // 检查文件是否成功打开
  if (!outfile.is_open())
  {
    std::cerr << "Error opening file: " << filename << std::endl;
    VecRestoreArrayRead(vec1, &array1);
    VecRestoreArrayRead(vec2, &array2);
    return;
  }
  std::string s1 = "mrst_";
  std::string s2 = "cpgrid_";
  outfile << (s1 + filename) << "," << "  " << "," << "  " << "," << (s2 + filename) << std::endl;
  // 写入数据到 CSV 文件，每个值占一行
  for (PetscInt i = 0; i < size1; ++i)
  {
    outfile << std::scientific << std::setprecision(10) << array1[i] << "," << "  " << "," << "  " << "," << array2[i] << std::endl;
  }

  // 释放 Vec 的数据指针
  VecRestoreArrayRead(vec1, &array1);
  VecRestoreArrayRead(vec2, &array2);

  // 关闭文件
  outfile.close();
}

PetscErrorCode SaveVecToMatlab(DM dm, Vec vec, const char *filename, const char *varname)
{
  PetscViewer viewer;
  PetscFunctionBegin;
  PetscCall(PetscObjectSetName((PetscObject)vec, varname));
  PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer));
  PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
  PetscCall(VecView(vec, viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}
