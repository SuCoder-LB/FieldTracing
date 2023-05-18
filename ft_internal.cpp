//
// Created by 苏立彪 on 2023/5/14.
//

#include "ft_interface.h"

namespace field_tracing {

#ifndef MAXITERATIONS
#define MAXITERATIONS 100000
#endif

//It works just on triangle meshes
ChartData computeChartData(TraceMesh &mesh, std::vector<std::vector<size_t>> &corners) {
  typedef std::map<std::pair<size_t, size_t>, std::pair<int, int>> EdgeLabelMap;
  typedef std::map<std::pair<size_t, size_t>, int> EdgeSubSideMap;

  ChartData chartData;

  if (mesh.face.size() == 0)
    return chartData;

  vcg::tri::UpdateTopology<TraceMesh>::FaceFace(mesh);

  //Region growing algorithm for getting charts
  std::vector<std::vector<size_t>> faces;
  std::vector<std::vector<size_t>> border_faces;

  {
    std::set<int> labels;
    std::vector<Chart> &charts = chartData.charts;

    for (size_t fId = 0; fId < mesh.face.size(); fId++) {
      if (!mesh.face[fId].IsD() && mesh.face[fId].Q() >= 0)
        labels.insert(mesh.face[fId].Q());
    }

    int maxChartLabel = *labels.rbegin();
    charts.resize(maxChartLabel + 1);

    std::vector<bool> visited(mesh.face.size(), false);

    for (size_t i = 0; i < mesh.face.size(); i++) {
      if (!mesh.face[i].IsD() && !visited[i]) {
        //Region growing to get chart
        std::stack<size_t> stack;
        stack.push(i);

        Chart &chart = charts[mesh.face[i].Q()];

        int label = mesh.face[i].Q();

        if (label >= 0) {
          std::set<size_t> borderFacesSet;
          do {
            size_t fId = stack.top();
            stack.pop();
            assert(mesh.face[fId].Q() == label);

            if (!visited[fId]) {
              faces[label].push_back(fId);
              typename TraceMesh::FaceType *currentFacePointer = &mesh.face[fId];
              vcg::face::Pos<typename TraceMesh::FaceType> pos(currentFacePointer, 0);
              for (int k = 0; k < 3; k++) {
                pos.FlipF();
                size_t adjFace = vcg::tri::Index(mesh, pos.F());

                //Saving border faces
                if (currentFacePointer == pos.F() || mesh.face[adjFace].Q() != label) {
                  borderFacesSet.insert(fId);
                } else {
                  if (!visited[adjFace]) {
                    stack.push(adjFace);
                  }
                }
                pos.FlipF();

                //Next edge
                pos.FlipV();
                pos.FlipE();
              }

              visited[fId] = true;
            }
          } while (!stack.empty());

          assert(borderFacesSet.size() >= 1);
          std::copy(borderFacesSet.begin(), borderFacesSet.end(), std::back_inserter(border_faces[label]));
        }
      }
    }
  }

  //TODO SPLIT IN FUNCTIONS
  EdgeSubSideMap edgeSubSideMap;
  std::vector<size_t> vertexNextMap(mesh.vert.size());

  for (int pId = 0; pId < chartData.charts.size(); ++pId) {

    Chart &chart = chartData.charts[pId];

    if (faces[pId].size() == 0)
      continue;

    std::unordered_set<size_t> cornerSet(corners[pId].begin(), corners[pId].end());

#ifndef NDEBUG
    if (cornerSet.size() < 3 || cornerSet.size() > 6) {
      std::cout << "Warning 3: Given as input for " << pId << ": " << cornerSet.size() << " sides." << std::endl;
    }
#endif

    EdgeLabelMap edgeLabelMap;

    std::set<size_t> remainingVertices;

    //Fill edge map and next vertex map
    for (const size_t &fId : border_faces[pId]) {
      typename TraceMesh::FaceType *currentFacePointer = &mesh.face[fId];
      vcg::face::Pos<typename TraceMesh::FaceType> pos(currentFacePointer, 0);

      for (int k = 0; k < currentFacePointer->VN(); k++) {
        pos.FlipF();
        size_t adjFace = vcg::tri::Index(mesh, pos.F());
        int adjLabel = mesh.face[adjFace].Q();

        bool isBorderEdge = false;
        int adjChartLabel = -2;

        if (currentFacePointer == pos.F()) {
          adjChartLabel = -1;
          isBorderEdge = true;
        } else if (adjLabel != pId) {
          adjChartLabel = adjLabel;
          isBorderEdge = true;
        }
        pos.FlipF();

        //For each border edge
        if (isBorderEdge) {
          assert(adjChartLabel > -2);

          typename TraceMesh::VertexType *vStart = pos.V();
          pos.FlipV();
          typename TraceMesh::VertexType *vEnd = pos.V();
          pos.FlipV();

          size_t vStartId = vcg::tri::Index(mesh, vStart);
          size_t vEndId = vcg::tri::Index(mesh, vEnd);

          std::pair<size_t, size_t> edge(vStartId, vEndId);
          if (edge.first > edge.second) {
            std::swap(edge.first, edge.second);
          }

          edgeLabelMap.insert(std::make_pair(edge, std::make_pair(pId, adjChartLabel)));
          vertexNextMap[vStartId] = vEndId;

          remainingVertices.insert(vStartId);
          remainingVertices.insert(vEndId);
        }

        pos.FlipV();
        pos.FlipE();
      }
    }

    do {
      //Find first label
      size_t vStartId;
      size_t vCurrentId;
      size_t vNextId;

      //Corner detection variables
      typename TraceMesh::CoordType lastEdgeVec;
      bool isCorner = false;

      vCurrentId = *remainingVertices.begin();

//      std::vector<size_t> nextConfiguration = findVertexChainPath(vCurrentId, vertexNextMap);
//      vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];
      vNextId = vertexNextMap[vCurrentId];

      //Get last edge vector
      lastEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
      lastEdgeVec.Normalize();

//            std::pair<size_t, size_t> startEdge(vCurrentId, vNextId);
//            if (startEdge.first > startEdge.second) {
//                std::swap(startEdge.first, startEdge.second);
//            }

      int currentLabel;

      //Iterate in the borders to get the first corner
      vStartId = vCurrentId;
      size_t firstCornerIterations = 0;
      do {
        //Next border edge
        vCurrentId = vertexNextMap[vCurrentId];
        vNextId = vertexNextMap[vCurrentId];

        typename TraceMesh::CoordType currentEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
        currentEdgeVec.Normalize();

        //Check if it is a corner
        isCorner = cornerSet.find(vCurrentId) != cornerSet.end();

        lastEdgeVec = currentEdgeVec;

        firstCornerIterations++;
      } while (!isCorner && vCurrentId != vStartId && firstCornerIterations < MAXITERATIONS);

#ifndef NDEBUG
      if (firstCornerIterations >= MAXITERATIONS) {
        std::cout << "Error: error iterating! Cannot find the first corner or get back to the start vertex."
                  << std::endl;
      }
#endif

#ifndef NDEBUG
      if (vCurrentId == vStartId) {
        std::cout << "Warning 1: input mesh is not well-defined: no corners!" << std::endl;
      }
#endif
      vStartId = vCurrentId;
      vNextId = vertexNextMap[vCurrentId];

      ChartSide currentSide;
      size_t chartSideId = 0;

      int adjChartLabel;
      do {
        size_t subsideId = chartData.subsides.size();
        ChartSubside currentSubSide;

        //Get edge
        std::pair<size_t, size_t> edge(vCurrentId, vNextId);
        if (edge.first > edge.second) {
          std::swap(edge.first, edge.second);
        }
        //Get current label on the other side
        const std::pair<int, int> &currentEdgeLabels = edgeLabelMap.at(edge);
        assert(currentEdgeLabels.first == pId || currentEdgeLabels.second == pId);
        adjChartLabel = currentEdgeLabels.first == pId ? currentEdgeLabels.second : currentEdgeLabels.first;

        std::unordered_set<size_t> cornerSetAdj;
        if (adjChartLabel >= 0)
          cornerSetAdj.insert(corners[adjChartLabel].begin(), corners[adjChartLabel].end());

        double length = 0;

        bool newSubSide = false;

        bool firstIteration = true;
        isCorner = false;
        bool isAdjCorner = false;
        size_t vSubSideStartId = vCurrentId;
        size_t iterations = 0;
        do {
          typename TraceMesh::CoordType currentEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
          currentEdgeVec.Normalize();

          std::pair<size_t, size_t> edge(vCurrentId, vNextId);
          if (edge.first > edge.second) {
            std::swap(edge.first, edge.second);
          }

          //Check if it is a corner
          if (!firstIteration) {
            isCorner = cornerSet.find(vCurrentId) != cornerSet.end();
            isAdjCorner = cornerSetAdj.find(vCurrentId) != cornerSetAdj.end();
          }

          //Get current label on the other subside
          const std::pair<int, int> &currentEdgeLabels = edgeLabelMap.at(edge);
          assert(currentEdgeLabels.first == pId || currentEdgeLabels.second == pId);
          currentLabel = currentEdgeLabels.first == pId ? currentEdgeLabels.second : currentEdgeLabels.first;

          if (!isCorner && !isAdjCorner && currentLabel == adjChartLabel) {
            EdgeSubSideMap::iterator findIt = edgeSubSideMap.find(edge);

            //If the subside has already been processed
            if (findIt == edgeSubSideMap.end()) {
              currentSubSide.vertices.push_back({mesh.vert[vCurrentId].P()[0],
                                                 mesh.vert[vCurrentId].P()[1],
                                                 mesh.vert[vCurrentId].P()[2],
                                                 vCurrentId});

              length += (mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P()).Norm();

              edgeSubSideMap.insert(std::make_pair(edge, subsideId));

              newSubSide = true;
            } else if (firstIteration) {
              subsideId = findIt->second;
            }
            firstIteration = false;

            remainingVertices.erase(vCurrentId);

            //Next border edge
            vCurrentId = vertexNextMap[vCurrentId];
            vNextId = vertexNextMap[vCurrentId];

            lastEdgeVec = currentEdgeVec;
          }
        } while (!isCorner && !isAdjCorner && currentLabel == adjChartLabel && vCurrentId != vSubSideStartId
            && iterations < MAXITERATIONS);
#ifndef NDEBUG
        if (iterations >= MAXITERATIONS) {
          std::cout << "Error: error iterating! Cannot find a corner or get back to the start vertex." << std::endl;
        }
#endif
#ifndef NDEBUG
        if (vCurrentId == vSubSideStartId) {
          std::cout << "Warning 2: input mesh is not well-defined: single border chart with no corners!" << std::endl;
        }
#endif

        //True if the subside is reversed (from the last to the first vertex)
        bool reversed;

        if (newSubSide) {
          //Add last vertex
          currentSubSide.vertices.push_back({mesh.vert[vCurrentId].P()[0],
                                             mesh.vert[vCurrentId].P()[1],
                                             mesh.vert[vCurrentId].P()[2]});

          assert(currentSubSide.vertices.size() >= 2);

          chartData.subsides.push_back(currentSubSide);

          reversed = false;
        } else {
          assert(currentSubSide.vertices.size() == 0);

          reversed = true;
        }

        currentSide.subsides.push_back(subsideId);
        currentSide.reversed_subside.push_back(reversed);

        if (isCorner) {
          chart.sides.push_back(currentSide);
          currentSide = ChartSide();
          chartSideId++;
        }

      } while (vCurrentId != vStartId);

#ifndef NDEBUG
      if (!isCorner) {
        std::cout << "Warning 4: Chart has no final corner!" << std::endl;
      }
#endif

    } while (!remainingVertices.empty());

#ifndef NDEBUG
    if (chart.sides.size() < 3 || chart.sides.size() > 6) {
      std::cout << "Warning 3: Chart " << pId << " has " << chart.sides.size() << " sides." << std::endl;
    }
#endif
  }

  return chartData;
}

void GetAllData(PatchTracer<TraceMesh> &PTr, TraceMesh &trace_mesh,
                std::vector<std::vector<size_t>> &corners) {
  typedef typename TraceMesh::CoordType CoordType;

  std::vector<std::pair<CoordType, CoordType> > SharpCoords;
  PTr.Mesh().GetSharpCoordPairs(SharpCoords);

  //copy the mesh
  vcg::tri::Append<TraceMesh, TraceMesh>::Mesh(trace_mesh, PTr.Mesh());
  std::vector<size_t> SharpCorners;
  PTr.getCornerSharp(SharpCorners);
  std::set<CoordType> SharpCornerPos;
  for (size_t i = 0; i < SharpCorners.size(); i++)
    SharpCornerPos.insert(PTr.Mesh().vert[SharpCorners[i]].P());

  for (size_t i = 0; i < trace_mesh.face.size(); i++)
    trace_mesh.face[i].Q() = PTr.FacePartitions[i];

  //merge across sharp features
  vcg::tri::Clean<TraceMesh>::RemoveDuplicateVertex(trace_mesh);
  vcg::tri::Clean<TraceMesh>::RemoveUnreferencedVertex(trace_mesh);
  vcg::tri::Allocator<TraceMesh>::CompactEveryVector(trace_mesh);
  trace_mesh.UpdateAttributes();

  //then remove the non disk-like
  std::set<size_t> ToErasePartitions;
  PTr.GetTopologicallyNotOKPartitionsIndex(ToErasePartitions);
  std::vector<std::vector<size_t> > PatchCorners;

  //add the patches of faces that has non-manifoldness (if occours)
  int num = vcg::tri::Clean<TraceMesh>::CountNonManifoldEdgeFF(trace_mesh, true);
  if (num > 0) {
    for (size_t i = 0; i < trace_mesh.face.size(); i++) {
      if (!trace_mesh.face[i].IsS())continue;
      size_t CurrP = trace_mesh.face[i].Q();
      ToErasePartitions.insert(CurrP);
    }
  }
  bool has_erased = false;
  if (ToErasePartitions.size() > 0) {
    std::cout << ToErasePartitions.size() << "- NON DISK OR NON MANIF PATCH ERASED- " << std::endl;
    has_erased = true;
    //first remove the faces
    size_t NumPart = 0;
    for (size_t i = 0; i < trace_mesh.face.size(); i++) {
      size_t IndexP = trace_mesh.face[i].Q();
      NumPart = std::max(NumPart, IndexP);
      if (ToErasePartitions.count(IndexP) == 0)continue;
      vcg::tri::Allocator<TraceMesh>::DeleteFace(trace_mesh, trace_mesh.face[i]);
    }

    //then remap the partitions
    std::vector<int> PartitionMap(NumPart + 1, -1);
    size_t currP = 0;
    for (size_t i = 0; i < PartitionMap.size(); i++) {
      if (ToErasePartitions.count(i) > 0)
        continue;//in this case do not consider that patch

      //otherwise add the patch corners and the index too
      PatchCorners.push_back(PTr.PartitionCorners[i]);
      PartitionMap[i] = currP;
      currP++;
    }

    //final check
    assert((currP + ToErasePartitions.size() - 1) == NumPart);

    //remap faces partitions
    for (size_t i = 0; i < trace_mesh.face.size(); i++) {
      if (trace_mesh.face[i].IsD())continue;
      int currP = trace_mesh.face[i].Q();
      assert(currP >= 0);
      assert(currP < PartitionMap.size());
      int newPIdx = PartitionMap[currP];
      assert(newPIdx >= 0);
      trace_mesh.face[i].Q() = newPIdx;
    }
  } else
    PatchCorners = PTr.PartitionCorners;

  //merge vertices
  vcg::tri::Clean<TraceMesh>::RemoveDuplicateVertex(trace_mesh);
  vcg::tri::Clean<TraceMesh>::RemoveUnreferencedVertex(trace_mesh);
  vcg::tri::Allocator<TraceMesh>::CompactEveryVector(trace_mesh);
  trace_mesh.UpdateAttributes();
  trace_mesh.UpdateFromCoordPairs(SharpCoords, false);

  //save vert pos
  std::map<CoordType, size_t> VertMap;
  for (size_t i = 0; i < trace_mesh.vert.size(); i++)
    VertMap[trace_mesh.vert[i].P()] = i;

  //update sharp vertices
  trace_mesh.SharpCorners.clear();
  for (size_t i = 0; i < trace_mesh.vert.size(); i++)
    if (SharpCornerPos.count(trace_mesh.vert[i].P()) > 0)
      trace_mesh.SharpCorners.push_back(i);

  PTr.WriteInfo();

  corners.clear();
  corners.resize(PatchCorners.size());
//  std::vector<std::vector<size_t>>F(PatchCorners.size());

  for (size_t i = 0; i < PatchCorners.size(); i++) {
    //check uniqueness
    std::set<size_t> TestCorn;
    for (size_t j = 0; j < PatchCorners[i].size(); j++) {
      if (PTr.CheckQuadrangulationLimits) {
        assert(PatchCorners[i].size() >= MIN_ADMITTIBLE);
        assert(PatchCorners[i].size() <= MAX_ADMITTIBLE);
      }
      int IndexV = PatchCorners[i][j];
      assert(IndexV < PTr.Mesh().vert.size());
      assert(IndexV >= 0);
      CoordType CornerPos = PTr.Mesh().vert[IndexV].P();

      assert(VertMap.count(CornerPos) > 0);
      corners[i].push_back(VertMap[CornerPos]);
      TestCorn.insert(VertMap[CornerPos]);
    }
    if (TestCorn.size() != PatchCorners[i].size()) {
      std::cout << "WARNING DOUBLE VERT: " << TestCorn.size() << "!=" << PatchCorners[i].size() << std::endl;
    }
  }

}

int LoadField(TraceMesh &trace_mesh, std::vector<std::array<double, 3>> &field) {
  assert(field.size() == trace_mesh.fn);

  for (int i = 0; i < field.size(); i++) {
    trace_mesh.face[i].PD1() = {field[i][0], field[i][1], field[i][2]};
    trace_mesh.face[i].PD2() = trace_mesh.face[i].PD1() ^ trace_mesh.face[i].N();
    trace_mesh.face[i].PD1().Normalize();
    trace_mesh.face[i].PD2().Normalize();
  }
  vcg::tri::CrossField<TraceMesh>::OrientDirectionFaceCoherently(trace_mesh);
  vcg::tri::CrossField<TraceMesh>::UpdateSingularByCross(trace_mesh, true);
  return 0;
}

int LoadFeature(TraceMesh &trace_mesh, std::vector<std::array<int, 2>> &features) {
  trace_mesh.SharpFeatures.clear();
  for (size_t i = 0; i < features.size(); i++) {
    int FIndex = features[i][0], EIndex = features[i][1];
    assert(FIndex >= 0);
    assert(FIndex < (int) trace_mesh.face.size());
    assert(EIndex >= 0);
    assert(EIndex < 4);

    trace_mesh.face[FIndex].SetFaceEdgeS(EIndex);
    trace_mesh.SharpFeatures.push_back(std::pair<size_t, size_t>(FIndex, EIndex));
  }
  return 0;
}

int FieldTracing(TraceMesh &trace_mesh, std::vector<std::array<double, 3>> &field,
                std::vector<std::array<int,2>>&features, ChartData &chart_data) {

  trace_mesh.UpdateAttributes();
  if (!LoadField(trace_mesh, field)) {
    return -1;
  }
  trace_mesh.UpdateAttributes();

  if(!LoadFeature(trace_mesh,features)){
    return -2;
  }
  trace_mesh.UpdateAttributes();

  trace_mesh.SolveGeometricIssues();
  trace_mesh.UpdateSharpFeaturesFromSelection();

  //preprocessing mesh
  PreProcessMesh(trace_mesh);

  //initializing graph
  VertexFieldGraph<TraceMesh> VGraph(trace_mesh);
  VGraph.InitGraph(false);

  //INIT TRACER
  typedef PatchTracer<TraceMesh> TracerType;
  TracerType PTr(VGraph);
  TraceMesh::ScalarType Drift = 100;
  bool add_only_needed = true;
  bool final_removal = true;
  bool meta_mesh_collapse = true;
  bool force_split = false;
  PTr.sample_ratio = 0.01;
  PTr.CClarkability = 1;
  PTr.split_on_removal = true;
  PTr.away_from_singular = true;
  PTr.match_valence = true;
  PTr.check_quality_functor = false;
  PTr.MinVal = 3;
  PTr.MaxVal = 5;
  PTr.Concave_Need = 1;
  PTr.PrioMode = PrioModBorder;

  //TRACING
  PTr.InitTracer(Drift, false);
  RecursiveProcess<TracerType>(PTr, Drift, add_only_needed, final_removal, true,
                               meta_mesh_collapse, force_split, true, false);
  PTr.SmoothPatches();

  TraceMesh trace_mesh_1;
  std::vector<std::vector<size_t>> corners;

  // 将ptr里被各种修改过的网格和分区的节点转换出来
  GetAllData(PTr, trace_mesh_1, corners);

  chart_data = computeChartData(trace_mesh_1, corners);
}

}