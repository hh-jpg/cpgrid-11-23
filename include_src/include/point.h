/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-13 15:00:54
 * @LastEditTime: 2024-06-19 15:35:26
 * @FilePath: /cpgrid/include/point.h
 * @Description:
 *
 */
// Point.h
#ifndef POINT_H
#define POINT_H
#include <iostream>
#include <iostream>
#include <cmath>
#include <array>
#include "config.h"

class Point
{
public:
  Point(double x = 0, double y = 0, double z = 0);
  Point operator+(const Point &rhs) const;
  Point operator-(const Point &rhs) const;
  Point operator*(const Point &rhs) const;
  Point operator/(const Point &rhs) const;

  Point& operator+=(const Point &rhs);
  Point& operator-=(const Point &rhs);
  Point& operator*=(const Point &rhs);
  Point& operator/=(const Point &rhs);

  Point operator+(const double &rhs) const;
  Point operator-(const double &rhs) const;
  Point operator*(const double &rhs) const;
  Point operator/(const double &rhs) const;

  Point& operator+=(const double &rhs);
  Point& operator-=(const double &rhs);
  Point& operator*=(const double &rhs);
  Point& operator/=(const double &rhs);

  // 其他成员函数
  double two_norm() const;
  double dot(const Point &rhs) const;
  double x() const;
  double y() const;
  double z() const;
private:
  std::array<double, 3> _coords;
};

inline std::ostream &operator<<(std::ostream &os, const Point &point)
{
  os << "(" << point.x() << ", " << point.y() << ", " << point.z() << ")";
  return os;
}

#endif // POINT_H
