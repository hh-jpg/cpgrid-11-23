/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-03 11:21:05
 * @LastEditTime: 2024-08-07 14:20:07
 * @FilePath: /cpgrid/src/ghosted_array_wrapper.C
 * @Description: 
 * 
 */

#include "ghosted_array_wrapper.h"

GhostedArrayWrapper::GhostedArrayWrapper(Vec x_local, DofMap &dof_map)
    : x_local(x_local), ghosted_array(nullptr), dofmap(dof_map)
{
  PetscErrorCode ierr = VecGetArrayRead(x_local, &ghosted_array);
  if (ierr)
  {
    throw std::runtime_error("Failed to get array read");
  }
}

GhostedArrayWrapper::~GhostedArrayWrapper()
{
  PetscErrorCode ierr = VecRestoreArrayRead(x_local, &ghosted_array);
  if (ierr)
  {
    std::cerr << "Failed to restore array read" << std::endl;
  }
}

const PetscScalar &GhostedArrayWrapper::operator[](PetscInt global_index) const
{
  int local_start = dofmap.first_dof();
  int local_end = dofmap.end_dof();
  std::unordered_map<int, int> &ghosted_indices = dofmap.ghosted_to_local_map();
  // 检查索引是否在局部索引范围内
  if (global_index >= local_start && global_index < local_end)
  {
    return ghosted_array[global_index - local_start];
  }
  // 检查索引是否在幽灵索引范围内
  auto it = ghosted_indices.find(global_index);
  if (it != ghosted_indices.end())
  {
    return ghosted_array[it->second];
  }
  // 如果索引不在任何已知范围内，抛出异常或进行适当处理
  std::cout << "In processor id " << dofmap.processor_id() << " invalid global index " << global_index << std::endl; 
  throw std::out_of_range("Invalid global index");
}

void GhostedArrayWrapper::print(int num_elements) const
{
  for (int i = 0; i < num_elements; ++i)
  {
    std::cout << ghosted_array[i] << " ";
  }
  std::cout << std::endl;
}
