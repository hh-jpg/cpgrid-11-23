/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-22 15:17:41
 * @LastEditTime: 2024-08-24 15:21:06
 * @FilePath: /cpgrid/src/utils.C
 * @Description: 
 * 
 */
#include <petscsys.h>
#include <petscvec.h>
#include "utils.h"
#include "system.h"
#include "fluid.h"
#include <stdexcept> // 为了使用 std::runtime_error

void parallel_only(MPI_Comm comm)
{
  int size;
  MPI_Comm_size(comm, &size);

  if (size <= 1)
  {
    throw std::runtime_error("This function must be run in a parallel environment.");
  }
}
