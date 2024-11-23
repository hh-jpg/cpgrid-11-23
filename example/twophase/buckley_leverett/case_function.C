/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-03 15:34:45
 * @LastEditTime: 2024-08-04 16:39:03
 * @FilePath: /cpgrid/example/twophase/buckley_leverett/case_function.C
 * @Description:
 *
 */
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <cstring>
#include <petscsys.h>

double wellindex(double K, double d1, double d2, double d3)
{
  double k11 = K;
  double k22 = K;
  double k21 = k22 / k11;
  double k12 = k11 / k22;
  double re1 = 2 * 0.14 * std::sqrt(d1 * d1 * std::sqrt(k21) + d2 * d2 * std::sqrt(k12));
  double re2 = std::pow(k21, 0.25) + std::pow(k12, 0.25);
  double re = re1 / re2;
  double rw = 0.1;
  double sk = 0;
  double WI = 2 * M_PI * d3 * std::sqrt(k11 * k22) / (std::log(re / rw) + sk);
  return WI;
}
