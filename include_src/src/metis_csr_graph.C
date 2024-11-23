/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-28 10:05:00
 * @LastEditTime: 2024-06-28 10:28:19
 * @FilePath: /cpgrid/src/metis_csr_graph.C
 * @Description:
 *
 */

#include "metis_csr_graph.h"
#include "unstruct_mesh.h"
#include "metis.h"

METIS_CSR_Graph CSRGraphFactory::createCSRGraph(Mesh &mesh)
{
  METIS_CSR_Graph csr_graph;
  prepareForCSR(mesh, csr_graph);
  return csr_graph;
}

void CSRGraphFactory::prepareForCSR(Mesh &mesh, METIS_CSR_Graph &csr_graph)
{
  std::vector<idx_t> &offsets = csr_graph.offsets;
  std::vector<idx_t> &vals = csr_graph.vals;
  int num_elem = mesh.getNumberOfPolyhedrons();
  offsets.resize(num_elem + 1, 0);
  for (int i = 0; i < num_elem; i++)
  {
    idx_t num_neighbor = 0;
    Polyhedron elem = mesh.getPolyhedron(i);
    std::vector<Face> &faces = elem.get_faces();
    for (Face &face : faces)
    {
      if (face.neighbor() != nullptr)
      {
        num_neighbor++;
        Polyhedron *neighbor = face.neighbor();
        idx_t nd = static_cast<idx_t>(neighbor->id());
        vals.push_back(nd);
      }
    }
    offsets[i + 1] = offsets[i] + num_neighbor;
  }
}
