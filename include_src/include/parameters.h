/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-15 09:35:23
 * @LastEditTime: 2024-08-23 15:07:28
 * @FilePath: /cpgrid/include/parameters.h
 * @Description:
 *
 */

#ifndef PARA_H
#define PARA_H

#include <iostream>
#include <string>
#include <petscsys.h> // 包括PETSc系统头文件
#include <unordered_map>
#include <any>
#include "config.h"

class Para
{
public:
  // 默认构造函数
  Para() = default;

  // 成员函数
  template <typename T>
  PetscErrorCode add_parameter(const char *option_name, T value)
  {
    PetscBool set;
    PetscErrorCode ierr;
    if constexpr (std::is_same_v<T, PetscInt>)
    {
      ierr = PetscOptionsGetInt(nullptr, nullptr, option_name, &value, &set);
      CHKERRQ(ierr);
    }
    else if constexpr (std::is_same_v<T, double> || std::is_same_v<T, PetscReal> || std::is_same_v<T, PetscScalar>)
    {
      ierr = PetscOptionsGetReal(nullptr, nullptr, option_name, &value, &set);
      CHKERRQ(ierr);
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
      char buffer[256];
      ierr = PetscOptionsGetString(nullptr, nullptr, option_name, buffer, sizeof(buffer), &set);
      CHKERRQ(ierr);
      if (set)
      {
        value = std::string(buffer);
      }
    }
    else
    {
      throw std::runtime_error("Unsupported type for parameter.");
    }

    if (set)
    {
      _parameters.insert({std::string(option_name), std::any(value)});
    }
    return 0;
  }

  template <typename T>
  T get_parameter(const char *para_name) const
  {
    std::string name(para_name);
    auto it = _parameters.find(name);

    // 检查是否找到对应的参数
    if (it != _parameters.end())
    {
      try
      {
        return std::any_cast<T>(it->second);
      }
      catch (const std::bad_any_cast &e)
      {
        throw std::runtime_error("Parameter found but with incorrect type for parameter: " + name);
   ;
      }
    }
    else
    {
      throw std::runtime_error("Parameter not found " + name);
    }
  }

  PetscErrorCode print_parameters() const;

private:
  std::unordered_map<std::string, std::any> _parameters;
};

#endif