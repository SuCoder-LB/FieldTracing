//
// Created by 苏立彪 on 2023/5/18.
//

#include "ft_basic_types.h"
#include "ft_interface.h"

//n_nodes n_subsides n_patches  //n_nodes个点，n_subsides个分区子边，n_patches个分区
//接下来n_nodes行，每行四个数
//node_id x y z        //点的编号，坐标x y z
//接下来n_subsides组数据，每组数据第1行3个数
//start end mid       //分区子边起始点编号，终止点编号（上面的分区节点编号），该条分区子边内部的点个数
//每组接下来mid行
//x y z                //中间节点依序的坐标x y z
//接下来n_patches组数据，每组数据第1行 n_sides+1 个数据
//n_sides n_subsides_1 ... n_subsides_n
//接下来n_subside1行以此类推
//subside_id reversed  //表示分区子边id，在该分区中，该条子边是否逆置

void SavePatches(const std::string &filename, field_tracing::ChartData &chart_data) {
  FILE *fp = fopen(filename.c_str(), "r");
  if (fp == nullptr)return;
  fprintf(fp, "%zu %zu %zu\n", chart_data.nodes.size(), chart_data.subsides.size(), chart_data.charts.size());

  for (auto &node : chart_data.nodes) {
    fprintf(fp, "%zd %lf %lf %lf\n", node.id, node.x, node.y, node.z);
  }

  for (auto subside : chart_data.subsides) {
    int n_mid = static_cast<int>(subside.vertices.size()) - 2;
    fprintf(fp, "%zd %zd %zu\n", subside.start, subside.end, subside.vertices.size() - 2);

    for (int j = 0; j < n_mid; ++j) {
      auto &node = subside.vertices[j + 2];
      fprintf(fp, "%lf% lf %lf\n", node.x, node.y, node.z);
    }
  }

  for (auto &chart : chart_data.charts) {
    fprintf(fp, "%zu ", chart.sides.size());
    for (auto &side : chart.sides) {
      fprintf(fp, "%zu ", side.subsides.size());
    }
    fprintf(fp, "\n");
    for (auto &side : chart.sides) {
      for (int k = 0; k < side.subsides.size(); ++k) {
        fprintf(fp, "%lu %d\n", side.subsides[k], static_cast<int>(side.reversed_subside[k]));
      }
    }
  }
  fclose(fp);
}

int LoadField(const std::string &file_name, std::vector<std::array<double, 3>> &field) {
  FILE *f = fopen(file_name.c_str(), "rt");
  if (!f) return -1;
  int num;
  fscanf(f, "%d\n", &num);
  for (int i = 0; i < num; i++) {
    double dirX, dirY, dirZ;
    fscanf(f, "%lf %lf %lf \n", &dirX, &dirY, &dirZ);
    field.push_back({dirX, dirY, dirZ});
  }
  fclose(f);
  return 0;
}

int LoadFeatures(const std::string &file_name, std::vector<std::array<int, 2>> &features) {
  FILE *f = fopen(file_name.c_str(), "rt");
  if (!f) return -1;
  int Num;
  fscanf(f, "%d\n", &Num);
  for (size_t i = 0; i < (size_t) Num; i++) {
    int FIndex, EIndex;
    fscanf(f, "%d,%d\n", &FIndex, &EIndex);
    assert(FIndex >= 0);
    assert(EIndex >= 0);
    assert(EIndex < 4);
    features.push_back({FIndex, EIndex});
  }
  fclose(f);
  return 0;
}

int main() {
  std::string mesh_file_name;

  std::cout << "Loading Remeshed M:" << mesh_file_name.c_str() << std::endl;

  std::string field_file_name;
  std::cout << "Loading Rosy Field:" << field_file_name.c_str() << std::endl;

  std::string sharp_file_name;
  std::cout << "Loading Sharp F:" << sharp_file_name.c_str() << std::endl;

  std::string chart_name;

  TraceMesh trace_mesh;
  std::vector<std::array<double, 3>> field;
  std::vector<std::array<int, 2>> features;
  field_tracing::ChartData chart_data;

  //Mesh load
  printf("Loading the mesh \n");
  bool loadedMesh = trace_mesh.LoadMesh(mesh_file_name);
  assert(loadedMesh);

  //Field load
  bool loadedField = LoadField(field_file_name, field);
  assert(loadedField);

  //Sharp load
  bool loadedFeatures = LoadFeatures(sharp_file_name, features);
  assert(loadedFeatures);

  field_tracing::FieldTracing(trace_mesh, field, features, chart_data);

  SavePatches(chart_name, chart_data);

  return 0;
}