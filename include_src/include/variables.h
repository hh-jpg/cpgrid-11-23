/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-01 09:45:20
 * @LastEditTime: 2024-08-28 10:18:38
 * @FilePath: /cpgrid/include/variables.h
 * @Description:
 *
 */
#ifndef VARIABLE_H
#define VARIABLE_H
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "config.h"
class Variables
{
public:
  Variables() = default;

  // 添加变量
  void add_variable(const std::string &name);

  // 通过序号获取名字
  std::string variable_name(int index) const;

  // 返回变量总数
  int size() const;

  // 通过名字获取序号
  int variable_num(std::string_view var_name) const;

private:
  std::vector<std::string> varnames;                // 顺序存储变量名
  std::map<std::string, int> nameToIndex; // 名字到序号的映射
};

#endif