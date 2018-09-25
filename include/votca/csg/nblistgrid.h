/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _NBLISTGRID_H
#define	_NBLISTGRID_H

#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>
#include "nblist.h"
#include <vector>

namespace votca { namespace csg {
using namespace votca::tools;

class NBListGrid
    : public NBList
{
public:
    void Generate(BeadList &list1, BeadList &list2, bool do_exclusions = true);
    void Generate(BeadList &list, bool do_exclusions = true);

protected:
    struct cell_t {
        BeadList _beads;
        std::vector<cell_t*> _neighbours;
    };

    vec _box_a, _box_b, _box_c;
    vec _norm_a, _norm_b, _norm_c;
    int _box_Na, _box_Nb, _box_Nc;

    std::vector<cell_t> _grid;
    Topology *_top;

    void InitializeGrid(const matrix &box);
    
    cell_t &getCell(const vec &r);
    cell_t &getCell(const int &a, const int &b, const int &c);

    void TestBead(cell_t &cell, std::shared_ptr<Bead> bead);
    void TestCell(cell_t &cell, std::shared_ptr<Bead> bead);
};

inline NBListGrid::cell_t &NBListGrid::getCell(const int &a, const int &b, const int &c)
{
    return _grid[a + _box_Na*b + _box_Na*_box_Nb*c];
}

}}

#endif	/* _NBLISTGRID_H */

