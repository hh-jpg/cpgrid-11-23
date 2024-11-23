/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-17 14:45:47
 * @LastEditTime: 2024-08-04 10:54:26
 * @FilePath: /cpgrid/include/utils.h
 * @Description:
 *
 */
#ifndef UTILS_H
#define UTILS_H
#include <petscsys.h>
#include <iostream>
#include <vector>  // 为了使用 std::vector
#include <mpi.h>   // 为了使用 MPI_Comm
#include "config.h"

typedef MPI_Comm communicator;

enum DofArrangementType
{
  PEV, // Processor-Elem-Variable,
  PVE  // Processor-Variable-Elem
};
enum FacePosition
{
  UNVALID_POSITION,
  LEFT,
  RIGHT,
  BACK,
  FRONT,
  TOP,
  BOTTOM
};

template <typename T>
inline T half_harmonicMean(T x, T y)
{
  if (x * y <= 0)
    return 0.0;
  return x * y / (x + y);
}

template <typename T>
void printVector(const std::vector<T> &vec)
{
  for (const T &element : vec)
  {
    std::cout << element << "\n";
  }
  std::cout << std::endl;
}

 
// 声明函数
void parallel_only(MPI_Comm comm);

#endif // UTILS_H
