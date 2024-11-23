/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-25 17:40:30
 * @LastEditTime: 2024-08-28 10:17:46
 * @FilePath: /cpgrid/include/metis_partitioner.h
 * @Description:
 *
 */
#ifndef METIS_PARTITIONER_H
#define METIS_PARTITIONER_H

#include "metis_csr_graph.h"
#include "unstruct_mesh.h"
#include "csr_graph_visualizer.h"
#include "config.h"
#include <iostream>
#include <stdexcept>
#include <vector>
class Mesh;
class METISPartitioner
{
public:
  void partition(Mesh &mesh, METIS_CSR_Graph& csr_graph, const unsigned int n_pieces);
  
private:
  void single_partition(Mesh & mesh);
  void assignElementsToProcessors(Mesh &mesh,
                                  const std::vector<idx_t> &part);
  
};

#endif