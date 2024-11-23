/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-17 15:03:04
 * @LastEditTime: 2024-08-28 10:18:23
 * @FilePath: /cpgrid/include/trans.h
 * @Description:
 *
 */
#ifndef TRANS_H
#define TRANS_H
#include "grid.h"
#include "config.h"
// compute the trans
// v_{i, k} \approx A_{i, k} \mathbf{K}_i \frac{\left(p_i-\pi_{i, k}\right) \vec{c}_{i, k}}
//                 {\left|\vec{c}_{i, k}\right|^2} \cdot \vec{n}_{i, k}
//                 = T_{i, k}\left(p_i-\pi_{i, k}\right)
PetscErrorCode tpfa_htrans_compute(Grid *grid);
PetscErrorCode tpfa_trans_compute(Grid *grid);
PetscErrorCode tpfa_efftrans_compute(Grid *grid);
#endif