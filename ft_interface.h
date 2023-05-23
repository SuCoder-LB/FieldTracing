//
// Created by 苏立彪 on 2023/5/14.
//

#ifndef FIELDTRACING__FT_INTERFACE_H_
#define FIELDTRACING__FT_INTERFACE_H_

#include "ft_basic_types.h"

class TraceMesh;
//#include "tracing/mesh_type.h"



namespace field_tracing{

int FieldTracing(TraceMesh &trace_mesh, std::vector<std::array<double, 3>> &field,
                std::vector<std::array<int,2>>&features, ChartData &chart_data);

}
#endif //FIELDTRACING__FT_INTERFACE_H_
