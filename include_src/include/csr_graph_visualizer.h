/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-27 17:23:58
 * @LastEditTime: 2024-06-28 10:37:48
 * @FilePath: /cpgrid/include/csr_graph_visualizer.h
 * @Description:
 *
 */
#ifndef CSRGRAPHVISUALIZER_H
#define CSRGRAPHVISUALIZER_H
#include <fstream>
#include <iostream>
#include <vector>
#include "metis_csr_graph.h"
#include "metis.h"
#include "petscsys.h"
#include "config.h"
#include "unstruct_mesh.h"
class CSRGraphVisualizer
{
public:
  static void visualize(const METIS_CSR_Graph &graph, const std::string &filename)
  {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::ofstream outfile(filename + ".dot");
      outfile << "graph G {" << std::endl;
      int n = graph.offsets.size() - 1;
      for (int i = 0; i < n; ++i)
      {
        for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; ++j)
        {
          int neighbor = graph.vals[j];
          if (i < neighbor)
          {
            outfile << "    " << i << " -- " << neighbor << ";" << std::endl;
          }
        }
      }
      outfile << "}" << std::endl;
      outfile.close();
      std::system(("dot -Tpng " + filename + ".dot -o " + filename + ".png").c_str());
      std::cout << "Graph visualization saved to " << filename << ".png" << std::endl;
    }
  }
};
#endif