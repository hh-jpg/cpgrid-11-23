/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-26 14:08:36
 * @LastEditTime: 2024-08-10 15:16:32
 * @FilePath: /cpgrid/src/metis_partitioner.C
 * @Description:
 *
 */
#include "metis.h"
#include "metis_partitioner.h"
#include "metis_csr_graph.h"
#include <stdexcept>
#include "utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "elem_iterator.h"
void METISPartitioner::partition(Mesh &mesh, METIS_CSR_Graph &csr_graph, const unsigned int n_pieces)
{
  int num_elem = mesh.getNumberOfPolyhedrons();
  std::vector<idx_t> part(num_elem);

  if (mesh.processor_id() == 0)
  {
    idx_t ncon = 1;
    std::vector<idx_t> vwgt(num_elem, 1);
    idx_t nparts = static_cast<idx_t>(n_pieces);
    idx_t edgecut = 0;
    if (n_pieces == 1)
    {
      single_partition(mesh);
      return;
    }
    if (n_pieces <= 8)
    {
      METIS_PartGraphRecursive(&num_elem, &ncon, csr_graph.offsets.data(), csr_graph.vals.data(),
                               vwgt.data(), nullptr, nullptr, &nparts, nullptr, nullptr, nullptr, &edgecut, part.data());
    }
    else
    {
      METIS_PartGraphKway(&num_elem, &ncon, csr_graph.offsets.data(), csr_graph.vals.data(),
                          vwgt.data(), nullptr, nullptr, &nparts, nullptr, nullptr, nullptr, &edgecut, part.data());
    }
  }

  MPI_Bcast(part.data(), num_elem, MPI_INT, 0, PETSC_COMM_WORLD);
  assignElementsToProcessors(mesh, part);
}

void METISPartitioner::assignElementsToProcessors(Mesh &mesh, const std::vector<idx_t> &part)
{
  int num_elem = mesh.getNumberOfPolyhedrons();
  for (int i = 0; i < num_elem; i++)
  {
    Polyhedron &current_elem = mesh.getPolyhedron(i);
    if (i != current_elem.id())
    {
      throw std::logic_error("i and elem id are not equal!");
    }
    current_elem.processor_id() = part[i];
    if (part[i] == mesh.processor_id())
    {
      mesh.n_local_elem()++;
    }
  }
}

void METISPartitioner::single_partition(Mesh &mesh)
{
  int num_elem = mesh.getNumberOfPolyhedrons();
  for (int i = 0; i < num_elem; i++)
  {
    Polyhedron &current_elem = mesh.getPolyhedron(i);
    if (i != current_elem.id())
    {
      throw std::logic_error("i and elem id are not equal!");
    }
    current_elem.processor_id() = 0;
  }
    mesh.n_local_elem() = num_elem;
}