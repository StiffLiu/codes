/*
 pip install pybind11
 c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) order_statistcs_pybinding.cpp -o my_qlib$(python3-config --extension-suffix)
 */

#include "order_statistics_tree.h"
#include <pybind11/pybind11.h>

using OSTreed = OSTree<double>;
namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(my_qlib, m) {
  py::class_<OSTreed >(m, "OSTreed")
      .def(py::init<double, size_t, size_t>(),
         "default"_a=NAN, "mincnt"_a=1, "maxcnt"_a=0)
      .def("find_by_rank", &OSTreed::find_by_order)
      .def("find_by_percentile", &OSTreed::find_by_percentile)
      .def("rank_of_key", &OSTreed::order_of_key)
      .def("percentile_of_key", &OSTreed::percentile_of_key)
      .def("push", &OSTreed::push)
      .def("pop", &OSTreed::pop)
      .def("clear", &OSTreed::clear)
      .def("size", &OSTreed::size)
      .def("__len__", &OSTreed::size);
}
