/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#pragma once
#ifndef VOTCA_CSG_GROWRITER_PRIV_H
#define VOTCA_CSG_GROWRITER_PRIV_H

namespace votca {
namespace csg {

template <class Topology_T>
double GROWriter<Topology_T>::formatVelocity_(const double &velocity) {
  return converter_.convert(Topology_T::velocity_unit, velocity_unit) *
         velocity;
}

template <class Topology_T>
int GROWriter<Topology_T>::formatId_(const int &id) {
  return (id + 1) % 100000;
}

template <class Topology_T>
double GROWriter<Topology_T>::formatDistance_(const double &distance) {
  return converter_.convert(Topology_T::distance_unit, distance_unit) *
         distance;
}

template <class Topology_T>
void GROWriter<Topology_T>::Open(std::string file, bool bAppend) {
  if (_out != nullptr) {
    throw std::runtime_error(
        "Cannot open file until you have closed the previously opend gro "
        "file.");
  }
  _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

template <class Topology_T>
void GROWriter<Topology_T>::Close() {
  fclose(_out);
  _out = nullptr;
}

template <class Topology_T>
void GROWriter<Topology_T>::Write(boost::any conf_any) {
  if (typeid(Topology_T *) != conf_any.type()) {
    throw std::runtime_error(
        "Error Cannot read topology using growriter write, incorrect topology "
        "type provided.");
  }
  Topology_T &top = *boost::any_cast<Topology_T *>(conf_any);
  char format[100];
  int i, resnr, l, vpr;

  fprintf(_out, "%s\n", "what a nice title");
  fprintf(_out, "%5d\n", static_cast<int>(top.BeadCount()));

  bool v = top.HasVel();
  int pr = 3;  // precision of writeout, given by the spec

  l = pr + 5;
  vpr = pr + 1;
  if (v)
    sprintf(format, "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n", l, pr,
            l, pr, l, pr, l, vpr, l, vpr, l, vpr);
  else
    sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);

  for (i = 0; static_cast<size_t>(i) < top.BeadCount(); i++) {
    resnr = top.getBead(i)->getResidueId();
    std::string resname = top.getBead(i)->getResidueType();
    std::string atomname = top.getBead(i)->getType();

    fprintf(_out, "%5d%-5.5s%5.5s%5d", formatId_(resnr), resname.c_str(),
            atomname.c_str(), formatId_(i));
    // next fprintf uses built format std::string
    Eigen::Vector3d r = top.getBead(i)->getPos();

    if (v) {
      Eigen::Vector3d vv = top.getBead(i)->getVel();
      fprintf(_out, format, formatDistance_(r.x()), formatDistance_(r.y()),
              formatDistance_(r.z()), formatVelocity_(vv.x()),
              formatVelocity_(vv.y()), formatVelocity_(vv.z()));
    } else {
      fprintf(_out, format, formatDistance_(r.x()), formatDistance_(r.y()),
              formatDistance_(r.z()));
    }
  }

  // write the boy
  Eigen::Matrix3d box = top.getBox();

  if (pr < 5) pr = 5;
  l = pr + 5;

  if (box(0, 1) || box(0, 2) || box(1, 0) || box(1, 2) || box(2, 0) ||
      box(2, 1)) {
    sprintf(format,
            "%%%d.%df%%%d.%df%%%d.%df"
            "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
            l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr);
    fprintf(_out, format, formatDistance_(box(0, 0)),
            formatDistance_(box(1, 1)), formatDistance_(box(2, 2)),
            formatDistance_(box(1, 0)), formatDistance_(box(2, 0)),
            formatDistance_(box(0, 1)), formatDistance_(box(2, 1)),
            formatDistance_(box(0, 2)), formatDistance_(box(1, 2)));
  } else {
    sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);
    fprintf(_out, format, formatDistance_(box(0, 0)),
            formatDistance_(box(1, 1)), formatDistance_(box(2, 2)));
  }
  fflush(_out);
}

}  // namespace csg
}  // namespace votca
#endif  // VOTCA_CSG_GROWRITER_PRIV_H
