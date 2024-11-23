/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-28 15:09:39
 * @LastEditTime: 2024-07-23 11:22:26
 * @FilePath: /cpgrid/include/sparsity_pattern.h
 * @Description:
 *
 */

/* This defines the sparsity pattern, or graph, of a sparse matrix.
 * The format is quite simple -- the global indices of the nonzero entries
 * in each row are packed into a vector.  The global indices (i,j) of the
 * nth nonzero entry of row i are given by j = sparsity_pattern[i][n];
 */

#include <vector>
#include <map>
#include "dof_map.h"
#include "polyhedron.h"
#include "config.h"
class SparsityPattern
{
public:
  void sorted_connected_dofs(const Polyhedron *elem,
                             std::vector<int> &dofs_vi,
                             unsigned int vi);

private:
  std::vector<std::vector<int>> Row;
  std::map<int, std::vector<std::vector<int>>> NonlocalGraph;
  std::vector<std::vector<std::vector<int>>> Graph;
  const DofMap &dof_map_in;
};