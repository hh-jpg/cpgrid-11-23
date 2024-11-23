/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-13 15:01:54
 * @LastEditTime: 2024-06-19 15:38:27
 * @FilePath: /cpgrid/src/point.C
 * @Description:
 *
 */
#include "point.h"
#include <stdexcept>
#include <cmath>

// 构造函数
Point::Point(double x, double y, double z)
{
  _coords[0] = x;
  _coords[1] = y;
  _coords[2] = z;
}

// 操作符重载
Point Point::operator+(const Point &rhs) const
{
  return Point(_coords[0] + rhs.x(), _coords[1] + rhs.y(), _coords[2] + rhs.z());
}

Point Point::operator-(const Point &rhs) const
{
  return Point(_coords[0] - rhs.x(), _coords[1] - rhs.y(), _coords[2] - rhs.z());
}

Point Point::operator*(const Point &rhs) const
{
  return Point(_coords[0] * rhs.x(), _coords[1] * rhs.y(), _coords[2] * rhs.z());
}

Point Point::operator/(const Point &rhs) const
{
  if (rhs.x() == 0.0 || rhs.y() == 0.0 || rhs.z() == 0.0)
  {
    throw std::invalid_argument("Division by zero");
  }
  return Point(_coords[0] / rhs.x(), _coords[1] / rhs.y(), _coords[2] / rhs.z());
}

Point &Point::operator+=(const Point &rhs)
{
  _coords[0] += rhs.x();
  _coords[1] += rhs.y();
  _coords[2] += rhs.z();
  return *this;
}

Point &Point::operator-=(const Point &rhs)
{
  _coords[0] -= rhs.x();
  _coords[1] -= rhs.y();
  _coords[2] -= rhs.z();
  return *this;
}

Point &Point::operator*=(const Point &rhs)
{
  _coords[0] *= rhs.x();
  _coords[1] *= rhs.y();
  _coords[2] *= rhs.z();
  return *this;
}

Point &Point::operator/=(const Point &rhs)
{
  if (rhs.x() == 0.0 || rhs.y() == 0.0 || rhs.z() == 0.0)
  {
    throw std::invalid_argument("Division by zero");
  }
  _coords[0] /= rhs.x();
  _coords[1] /= rhs.y();
  _coords[2] /= rhs.z();
  return *this;
}

Point Point::operator+(const double &rhs) const
{
  return Point(_coords[0] + rhs, _coords[1] + rhs, _coords[2] + rhs);
}

Point Point::operator-(const double &rhs) const
{
  return Point(_coords[0] - rhs, _coords[1] - rhs, _coords[2] - rhs);
}

Point Point::operator*(const double &rhs) const
{
  return Point(_coords[0] * rhs, _coords[1] * rhs, _coords[2] * rhs);
}

Point Point::operator/(const double &rhs) const
{
  if (rhs == 0.0)
  {
    throw std::invalid_argument("Division by zero");
  }
  return Point(_coords[0] / rhs, _coords[1] / rhs, _coords[2] / rhs);
}

Point &Point::operator+=(const double &rhs)
{
  for (int i = 0; i < 3; ++i)
  {
    _coords[i] += rhs;
  }
  return *this;
}

Point &Point::operator-=(const double &rhs)
{
  for (int i = 0; i < 3; ++i)
  {
    _coords[i] -= rhs;
  }
  return *this;
}

Point &Point::operator*=(const double &rhs)
{
  for (int i = 0; i < 3; ++i)
  {
    _coords[i] *= rhs;
  }
  return *this;
}

Point &Point::operator/=(const double &rhs)
{
  if (rhs == 0.0)
  {
    throw std::invalid_argument("Division by zero");
  }
  for (int i = 0; i < 3; ++i)
  {
    _coords[i] /= rhs;
  }
  return *this;
}

// 其他成员函数
double Point::two_norm() const
{
  return std::sqrt(_coords[0] * _coords[0] + _coords[1] * _coords[1] + _coords[2] * _coords[2]);
}

double Point::dot(const Point &rhs) const
{
  return _coords[0] * rhs.x() + _coords[1] * rhs.y() + _coords[2] * rhs.z();
}

double Point::x() const { return _coords[0]; }
double Point::y() const { return _coords[1]; }
double Point::z() const { return _coords[2]; }
