/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-03 13:56:48
 * @LastEditTime: 2024-08-03 10:44:45
 * @FilePath: /cpgrid/src/variables.C
 * @Description:
 *
 */
#include "variables.h"

void Variables::add_variable(const std::string &name)
{
  if (nameToIndex.find(name) == nameToIndex.end())
  {
    int index = varnames.size();
    varnames.push_back(name);
    nameToIndex[name] = index;
  }
}

std::string Variables::variable_name(int index) const
{
  if (index >= 0 && index < varnames.size())
  {
    return varnames[index];
  }
  return "";
}

int Variables::size() const
{
  return varnames.size();
}

int Variables::variable_num(std::string_view var_name) const
{
  std::string name(var_name);
  auto it = nameToIndex.find(name);
  if (it != nameToIndex.end())
  {
    return it->second;
  }
  return -1; // 表示未找到
}
