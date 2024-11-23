/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-09 09:46:19
 * @LastEditTime: 2024-08-21 10:00:05
 * @FilePath: /cpgrid/src/sparsity.C
 * @Description:
 *
 */
#include "sparsity.h"
#include "petscmat.h"
void Sparsity::prepare_sparsity()
{
  Mesh &mesh = _dof_map.get_mesh();
  PetscInt n_vars = static_cast<PetscInt>(_dof_map.n_vars());
  int n_local_dofs = _dof_map.n_local_dofs();
  _d_nnz.resize(n_local_dofs, 0);
  _o_nnz.resize(n_local_dofs, 0);

  for (const Polyhedron &elem : mesh.elem_local_range())
  {
    for (int i = 0; i < n_vars; i++)
    {
      int local_index = _dof_map.dof_local_indices(elem, i);
      _d_nnz[local_index] += n_vars;
      for (const Face &face : elem.get_faces())
      {
        Polyhedron *neighbor = face.neighbor();
        if (neighbor != nullptr)
        {
          int neighbordof_global = _dof_map.dof_indices(*neighbor, 0);
          if (_dof_map.is_in_local_range(neighbordof_global))
          {
            _d_nnz[local_index] += n_vars;
          }
          else
          {
            _o_nnz[local_index] += n_vars;
          }
        }
      }
    }
  }
}

PetscErrorCode Sparsity::creatmat(Mat &mat)
{
  PetscFunctionBegin;

  // 创建矩阵对象
  PetscCall(MatCreate(_dof_map.comm(), &mat));

  // 获取全局和局部尺寸，以及块大小
  PetscInt m_global = static_cast<PetscInt>(_dof_map.n_dofs());
  PetscInt n_global = static_cast<PetscInt>(_dof_map.n_dofs());
  PetscInt m_local = static_cast<PetscInt>(_dof_map.n_local_dofs());
  PetscInt n_local = static_cast<PetscInt>(_dof_map.n_local_dofs());
  PetscInt blocksize = static_cast<PetscInt>(_dof_map.n_vars());

  // 设置矩阵尺寸和块大小
  PetscCall(MatSetSizes(mat, m_local, n_local, m_global, n_global));
  PetscCall(MatSetBlockSize(mat, blocksize));

  // 设置矩阵类型和选项
  PetscCall(MatSetType(mat, MATBAIJ)); // 自动选择 seqbaij 或 mpibaij
  PetscCall(MatSetOptionsPrefix(mat, ""));
  PetscCall(MatSetFromOptions(mat));
  // 设置预分配内存
  PetscCall(MatSeqAIJSetPreallocation(mat, 0, _d_nnz.empty() ? nullptr : _d_nnz.data()));
  PetscCall(MatMPIAIJSetPreallocation(mat, 0, _d_nnz.empty() ? nullptr : _d_nnz.data(), 0, _o_nnz.empty() ? nullptr : _o_nnz.data()));

  // 设置选项，防止分配新非零条目
  PetscCall(MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));

  // 检查矩阵的局部大小
  PetscInt m_l, n_l;
  PetscCall(MatGetLocalSize(mat, &m_l, &n_l));
  // 如果矩阵局部大小不为零，清零所有条目
  if (n_l)
  {
    PetscCall(MatZeroEntries(mat));
  }

  // 正常返回
  PetscFunctionReturn(PETSC_SUCCESS);
}
