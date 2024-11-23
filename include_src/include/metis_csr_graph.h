/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-25 16:52:26
 * @LastEditTime: 2024-08-28 10:17:08
 * @FilePath: /cpgrid/include/metis_csr_graph.h
 * @Description: 
 * 
 */
/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-25 16:52:26
 * @LastEditTime: 2024-06-28 10:28:41
 * @FilePath: /cpgrid/include/metis_csr_graph.h
 * @Description:
 *
 */

#ifndef METIS_CSR_GRAPH_H
#define METIS_CSR_GRAPH_H

#include <vector>
#include "metis.h"
#include "unstruct_mesh.h"
#include "config.h"
class METIS_CSR_Graph
{
public:
    std::vector<idx_t> offsets; // CSR格式中的偏移数组
    std::vector<idx_t> vals;    // CSR格式中的邻接数组

    // 默认构造函数
    METIS_CSR_Graph() = default;

    // 构造函数
    METIS_CSR_Graph(const std::vector<idx_t> &offsets, const std::vector<idx_t> &vals)
        : offsets(offsets), vals(vals) {}

    // 获取顶点数量
    idx_t numVertices() const
    {
        return static_cast<idx_t>(offsets.size() - 1);
    }

    // 获取边数量
    idx_t numEdges() const
    {
        return static_cast<idx_t>(vals.size());
    }

    // 清空图
    void clear()
    {
        offsets.clear();
        vals.clear();
    }
};

class CSRGraphFactory
{
public:
    static METIS_CSR_Graph createCSRGraph(Mesh &mesh);

private:
    static void prepareForCSR(Mesh &mesh, METIS_CSR_Graph & csr_graph);
};

#endif // METIS_CSR_GRAPH_H
