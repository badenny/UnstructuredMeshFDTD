#include <iostream>
#include <vector>
#include <array>
#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include "EdgeHash.h"

using Edge = std::array<size_t, 2>;
using Hex = std::array<size_t, 8>;
/* Returns the file name of a unit cube mesh with n subdivisions in each direction.
 */
std::string meshUnitCube(const int n) {
  const double sideLength = 1.;
  using namespace gmsh::model::geo;
  using DArray = std::array<double, 3>;
  gmsh::initialize();
  //bottom square and top square corners of cube
  //bottom square has indices 0,1,2,3 clockwise; top has 4,5,6,7 clockwise
  const std::array<DArray, 8> corners {DArray{0, 0, 0}, DArray{1, 0, 0}, DArray{1, 1, 0}, DArray{0, 1, 0},
                                       DArray{0, 0, 1}, DArray{1, 0, 1}, DArray{1, 1, 1}, DArray{0, 1, 1}};
  const double dx = sideLength / n;
  std::vector<int> geoPoints;
  //loop to add corners as gmsh points
  for (auto i = 0; i < corners.size(); ++i) {
    const auto &xyz = corners[i];
    geoPoints.push_back(addPoint(xyz[0], xyz[1], xyz[2], dx, i));
  }
  auto addSquare = [&](const std::array<int, 5> indices){
    std::vector<int> lines;
    for(auto i = 0; i < 4; ++i){
      const auto ind = indices[i];
      const auto nextInd = indices[i+1];
      //std::cout << geoPoints[ind] << "," << geoPoints[nextInd] << ": " ;
      const auto line = addLine(geoPoints[ind], geoPoints[nextInd]);
      mesh::setTransfiniteCurve(line, n);
      //std::cout << "line=" << line << std::endl;
      lines.push_back(line);
    }
    const auto loop = addCurveLoop(lines);
    //std::cout << "loop=" << loop << '\n';
    //return addPlaneSurface(std::vector<int>{loop});
    return loop;
  };
  std::vector<std::array<int, 5>> squareLoopsIndices;
  squareLoopsIndices.push_back(std::array<int, 5>{0, 1, 2, 3, 0}); //xy, z=0 square
//  squareLoopsIndices.push_back(std::array<int, 5>{4, 5, 6, 7, 4}); //xy, z=1 square
//  squareLoopsIndices.push_back(std::array<int, 5>{0, 1, 5, 4, 0}); //xz, y=0 square
//  squareLoopsIndices.push_back(std::array<int, 5>{3, 2, 6, 7, 3}); //xz, y=1 square
//  squareLoopsIndices.push_back(std::array<int, 5>{0, 3, 7, 4, 0}); //yz, x=0 square
//  squareLoopsIndices.push_back(std::array<int, 5>{1, 2, 6, 5, 1}); //yz, x=1 square

  //std::vector<int> squareSurfaces;
  synchronize();
  std::vector<int> squareLoops;
  std::vector<int> squareSurfaces;
  for(const auto &indicesArray : squareLoopsIndices){
    //squareSurfaces.push_back(addSquare(indicesArray));
    const auto loop = addSquare(indicesArray);
    squareLoops.push_back(loop);
    squareSurfaces.push_back(addPlaneSurface({loop}));
  }
  synchronize();
  mesh::setTransfiniteSurface(1);
//  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2); // or 3
//  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1); //use 1 to force quads, use 2 to force hexahedra
//  gmsh::model::mesh::generate(2);
//  gmsh::model::mesh::recombine();
//  gmsh::model::mesh::refine();
  synchronize();
  std::vector<std::pair<int, int> > extrusion;
  const std::vector<int> & numElements {n};
  const std::vector<double> & heights {1};
  extrude({std::pair<int,int>{2, squareSurfaces[0]}}, 0, 0, 1, extrusion, numElements, heights, true);
//  extrude({std::pair<int,int>{2, squareSurfaces[0]}}, 0, 0, 1, extrusion);
//  for(auto ex : extrusion){
//    std::cout << "my extrusion " << ex.first << ',' << ex.second << '\n';
//  }

//  for(auto [myDim, myTag] : extrusion){
//    if(myDim == 1) mesh::setTransfiniteCurve(myTag, dx);
//    if(myDim == 2) mesh::setTransfiniteSurface(myTag);
//    //if(myDim == 3) mesh::setTransfiniteVolume(myTag);
//  }
  //mesh::setTransfiniteVolume(1);

//  int surfaceLoop = addSurfaceLoop(squareSurfaces);
//  addVolume({surfaceLoop});
  synchronize();
  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 3); // or 3
  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 2); //use 1 to force quads, use 2 to force hexahedra
  gmsh::model::mesh::generate(3);
//  gmsh::model::mesh::recombine();
//  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 2); //use 1 to force quads, use 2 to force hexahedra
//  gmsh::model::mesh::refine();

  //std::string meshFormat{".vtk"};
  std::string meshFormat{".msh"};
  //std::string meshFormat{".m"};
  auto fileName = "UnitCubeUnstructured" + std::to_string(n) + meshFormat;
  //auto fileName = "UnitSquareStructured" + std::to_string(n) + meshFormat;
  //gmsh::write("UnitSquare" + std::to_string(n) + ".msh");
  gmsh::write(fileName);
  //gmsh::write("UnitSquareStructured" + std::to_string(n) + ".m");
  //gmsh::write("UnitSquareUnstructured" + std::to_string(n) + ".m");
  gmsh::finalize();
  return fileName;
}

/*
                           2
3----------2         3----------2          +----------+
|\         |\        |\         |\         |\         |\
| \        | \       | \11      | \10      | \        | \
|  \       |  \     3|  \      1|  \       |  \    6  |  \
|   7------+---6     |   7------+---6      |   7------+---6
|   |      |   |     |   | 0    |   |      |   |      |   |
0---+------1   |     0---+------1   |      +---+------+   |
 \  |       \  |      \  |       \  |       \ 7|       \  |5
  \ |        \ |      8\ |       9\ |        \ |        \ |
   \|         \|        \|         \|         \|         \|
    4----------5         4----------5          4----------5
                                                     4
Use low index to high index ordering.
Returns local nodes of edge for a given edge index e.g. edge index 8 consists of nodes (0,4).
 */
Edge edgeIndexToNodeIndex(const size_t i){
  Edge edge;
  switch(i){
    case 0:
      edge = {0, 1};
      break;
    case 1:
      edge = {1, 2};
      break;
    case 2:
      edge = {2, 3};
      break;
    case 3:
      edge = {0, 3};
      break;
    case 4:
      edge = {4, 5};
      break;
    case 5:
      edge = {5, 6};
      break;
    case 6:
      edge = {6, 7};
      break;
    case 7:
      edge = {4, 7};
      break;
    case 8:
      edge = {0, 4};
      break;
    case 9:
      edge = {1, 5};
      break;
    case 10:
      edge = {2, 6};
      break;
    case 11:
      edge = {3, 7};
      break;
    default:
      edge = {0, 0};
      std::cout << "Invalid edge index " << i << " in edgeIndexToNodeIndex()" << '\n';
      break;
  }
  return edge;
}

/*
                           2
3----------2         3----------2          +----------+
|\         |\        |\         |\         |\         |\
| \        | \       | \11      | \10      | \        | \
|  \       |  \     3|  \      1|  \       |  \    6  |  \
|   7------+---6     |   7------+---6      |   7------+---6
|   |      |   |     |   | 0    |   |      |   |      |   |
0---+------1   |     0---+------1   |      +---+------+   |
 \  |       \  |      \  |       \  |       \  |7      \  |5
  \ |        \ |      8\ |       9\ |        \ |        \ |
   \|         \|        \|         \|         \|         \|
    4----------5         4----------5          4----------5
                                                     4
Use low index to high index ordering.
Given a local edge index, getLocalParallelEdge returns the three local edge indices which are parallel to it
E.g. the sets of local parallel edges are {0, 2, 4, 6} and {1, 3, 5, 7} and {8, 9, 10, 11}
 */
std::vector<size_t> getLocalParallelEdge(size_t edgeIndex){
  std::vector<size_t> result;
  if(!(0 <= edgeIndex && edgeIndex < 12)){
    std::cout << "Edge index passed in get getLocalParallelEdge needs to be in [0, 12)." << std::endl;
    return result;
  }
  std::vector<std::vector<size_t>> parallelEdges = {{0, 2, 4, 6},
                                                    {1, 3, 5, 7},
                                                    {8, 9, 10, 11}};
  //loop over sets of parallel edges and find set that contains edgeIndex
  for(auto pEdges : parallelEdges){
    auto it = std::find(pEdges.begin(), pEdges.end(), edgeIndex);
    //if the set pEdges contains edgeIndex, remove current index and set result to the remaining three edges
    if(it != std::end(pEdges)){
      pEdges.erase(it);
      result = pEdges;
      break;
    }
  }
  return result;
}

int localEdgeDirection(size_t edgeIndex){
  //only two cases when low to high pointing doesn't give parallel edges
  if (edgeIndex == 2 || edgeIndex == 6){
    return -1;
  }
  else{
    return 1;
  }
}

/*
 * Function to rearrange hex nodes given the index of the new desired origin
 * For example, if we want the new origin to be 5, the reordering is
 * {5, 4,
 */
std::array<size_t, 8> rearrangeNodes(size_t originIndex){
  switch(originIndex){
    case 0:
      return {0, 1, 2, 3, 4, 5, 6, 7};
    case 1:
      return {1, 5, 6, 2, 0, 4, 7, 3};
    case 2:
      return {2, 3, 0, 1, 6, 7, 4, 5};
    case 3:
      return {3, 7, 4, 0, 2, 6, 5, 1};
    case 4:
      return {4, 0, 3, 7, 5, 1, 2, 6};
    case 5:
      return {5, 4, 7, 6, 1, 0, 3, 2};
    case 6:
      return {6, 7, 3, 2, 5, 4, 0, 1};
    case 7:
      //return {7, 3, 2, 6, 4, 0, 1, 5};
      return {7, 6, 5, 4, 3, 2, 1, 0};
    default:
      std::cout << "Invalid value in rearrangeNodes()" << std::endl;
      return {0, 1, 2, 3, 4, 5, 6, 7};
  }
}

void orientEdges(const size_t globalEdgeIndex, const size_t localEdgeIndex, const size_t hexIndex,
                 const int desiredOrientation, const std::vector<std::vector<size_t>> &edge2hex,
                 const std::vector<std::array<size_t, 12>> &hex2edge, const std::vector<Edge> &edge2node,
                 const std::vector<std::array<size_t, 8>> &hex2node,
                 std::vector<int> &globalOrientation, std::vector<std::array<int, 12>> &localOrientation) {
  const auto globalOri = globalOrientation[globalEdgeIndex];
  const auto localOri = localOrientation[hexIndex][localEdgeIndex];
  //std::cout << "edge=(" << edge2node[globalEdgeIndex][0] << "," << edge2node[globalEdgeIndex][1]
  //          <<  ") local=" << localOri << ", global=" << globalOri << '\n';
  int newDesiredOrientation = 0;
  //Edge is already oriented globally and locally
  if(globalOri != 0 && localOri !=0)
    return;
  //Globally oriented but not locally oriented. This is the case when jumping across hexahedra.
  else if(globalOri != 0 && localOri == 0){
    const auto [localNode0, localNode1] = edgeIndexToNodeIndex(localEdgeIndex);
    std::array<size_t, 2> globalEdgeNodes = {hex2node[hexIndex][localNode0], hex2node[hexIndex][localNode1]};
    //check global orientation against local
//    const auto [globalNode0, globalNode1] = edge2node[globalEdgeIndex];
    //if the local orientation were positive, this is the relative orientation.
//    int initialRelativeSigma = (globalNode0 < globalNode1) ? 1 : -1;
    //if the local orientation aligned with the global, this is the relative orientation.
    int initialLocalSigma = localEdgeDirection(localEdgeIndex);
    int initialRelativeSigma = (globalEdgeNodes[0] < globalEdgeNodes[1]) ? 1 : -1;
    //correct this by the global ori found in previous hex
    int relativeSigma = initialLocalSigma * globalOri * initialRelativeSigma;
    newDesiredOrientation = relativeSigma;
    //make local edge point in same direction as global edge.
    localOrientation[hexIndex][localEdgeIndex] = newDesiredOrientation;
  }
  //no global or local orientation. Case when orienting within same hex.
  else{
    //set local orientation to desired orientation
    localOrientation[hexIndex][localEdgeIndex] = desiredOrientation;
    //check global orientation against local
    //const auto [globalNode0, globalNode1] = edge2node[globalEdgeIndex];
    const auto [localNode0, localNode1] = edgeIndexToNodeIndex(localEdgeIndex);

    std::array<size_t, 2> globalEdgeNodes = {hex2node[hexIndex][localNode0], hex2node[hexIndex][localNode1]};
//    std::cout <<  localNode0 << "," << localNode1 << ": ";
//    std::cout <<  globalEdgeNodes[0]-1 << "," << globalEdgeNodes[1]-1 << "\n";
    int initialLocalSigma = localEdgeDirection(localEdgeIndex);
    //if the local orientation were positive and low to high local index, this is the relative orientation.
    int initialRelativeSigma = (globalEdgeNodes[0] < globalEdgeNodes[1]) ? 1 : -1;
    //then correct by the desiredOrientation and the local
    int relativeSigma = initialLocalSigma * desiredOrientation * initialRelativeSigma;
    //if relative sigma is positive, then the local orientation and global agree
    globalOrientation[globalEdgeIndex] = relativeSigma;
    newDesiredOrientation = desiredOrientation;
  }
  //orient the current global edge in all of its adjacent hexahedra.
  for(auto h : edge2hex[globalEdgeIndex]){
    auto newHexEdges = hex2edge[h];
    //get new local edge index
    auto newLocalEdgeIndex = -1;
    for(auto i = 0; i < 12; ++i){
      if(newHexEdges[i] == globalEdgeIndex){
        newLocalEdgeIndex = i;
        break;
      }
    }
    orientEdges(globalEdgeIndex, newLocalEdgeIndex, h, newDesiredOrientation, edge2hex, hex2edge,
                edge2node, hex2node, globalOrientation, localOrientation);
  }
  //orient edges within hexahedron
  std::vector<size_t> localParallelEdges = getLocalParallelEdge(localEdgeIndex);
  for(auto p : localParallelEdges){
    size_t newGlobalEdgeIndex = hex2edge[hexIndex][p];
    orientEdges(newGlobalEdgeIndex, p, hexIndex, newDesiredOrientation, edge2hex, hex2edge,
                edge2node, hex2node, globalOrientation, localOrientation);
  }
}

/* Function that returns local edges and their direction stemming from a local node in a hexahedron.
 * The format returned is {edge0, edge1, edge2, dir0, dir1, dir2} where dir=+1 means the edge starts
 * at the node. dir=-1 means it ends at the node
 *                             8
 * For example, if we have 0----->4 the returned array will contain 8 as an edge and +1 as the direction of
 * that edge.
 */
std::array<int, 6> edgesStemmingFromNode(const size_t i){
  switch (i) {
    case 0:
      return {0, 3, 8,
              1, 1, 1};
    case 1:
      return {0, 1, 9,
              -1, 1, 1};
    case 2:
      return {1, 2, 10,
              -1, 1, 1};
    case 3:
      return {2, 3, 11,
              -1, -1, 1};
    case 4:
      return {4, 7, 8,
              1, 1, -1};
    case 5:
      return {4, 5, 9,
              -1, 1, -1};
    case 6:
      return {5, 6, 10,
              -1, 1, -1};
    case 7:
      return {6, 7, 11,
              -1, -1, -1};
    default:
      std::cout << "Error in edgeStemmingFromNode(i): Function needs i in [0,12)." << std::endl;
      return {-1, -1, -1, -1, -1, -1};
  }
}

int main() {
  namespace mesh = gmsh::model::mesh;
  namespace geo = gmsh::model::geo;
  auto fileName = meshUnitCube(1);
  gmsh::initialize();
  gmsh::open(fileName);
  std::vector<int> elementTypes;
  std::vector<std::vector<size_t>> elementTags, elementNodeTags;
  mesh::getElements(elementTypes, elementTags, elementNodeTags, 3, -1);
  //std::cout << elementTypes.size() << '\n';
  //std::cout << elementTags[0].size() << '\n';
  //std::cout << elementNodeTags[0].size() / 8 << " " << elementNodeTags.size() % 8;
  auto &myTags = elementTags[0];
  auto &hexNodes = elementNodeTags[0];

  std::vector<Hex> hex2node;
  //loop to populate hex2node
  for(auto i = 0; i < myTags.size(); ++i) {
    auto globalIndex = i * 8;
    Hex hex;
    //Populate local hex from nodes in global list. see Gmsh documentation on getElements
    for (auto j = 0; j < 8; ++j) {
      hex[j] = hexNodes[globalIndex + j];
    }
    hex2node.push_back(hex);
  }
  auto numHex = hex2node.size();
  //Populate edge2node and edgeMap
  std::vector<Edge> edge2node;
  std::unordered_map<Edge, size_t, EdgeHash> edgeMap;
  std::vector<std::vector<size_t>> edge2hex;
  std::vector<std::array<size_t, 12>> hex2edge(numHex);
  size_t edgeIndex = 0;
  for(size_t hexIndex = 0; hexIndex < numHex; ++hexIndex){
    auto hex = hex2node[hexIndex];
    //loop over local edges to add edges to edgeMap
    std::array<size_t, 12> localHexEdges;
    for(auto k = 0; k < 12; ++k){
      auto localEdge = edgeIndexToNodeIndex(k);
      Edge globalEdge = {hex[localEdge[0]], hex[localEdge[1]]};
      std::sort(globalEdge.begin(), globalEdge.end());
      //try to insert global edge in edgeMap
      auto [it, isNewEdge] = edgeMap.emplace(globalEdge, edgeIndex);
      //if successful it is new, update edge2node and edge2hex to reflect this
      if(isNewEdge){
        edge2node.push_back(globalEdge);
        std::vector<size_t> hexVector = {hexIndex}; //have to construct this because edge2hex[edgeIndex] has no vector yet
        edge2hex.push_back(hexVector);
        localHexEdges[k] = edgeIndex;
        ++edgeIndex;
      }
      //otherwise it is not new. This means our current hex is a part of an existing edge
      //update edge2hex to reflect this
      else{
        size_t repeatedEdgeIndex = edgeMap[globalEdge];
        edge2hex[repeatedEdgeIndex].push_back(hexIndex); //add current hex as adjacent to repeated edge
        localHexEdges[k] = repeatedEdgeIndex;
      }
    }
    hex2edge[hexIndex] = localHexEdges;
  }

  const auto numEdges = edge2node.size();
//  //Orientation set-up
  std::array<int, 12> zeros;
  zeros.fill(0);
  std::vector<std::array<int, 12>> localOrientation(numHex, zeros);
  std::vector<int> globalOrientation(numEdges, 0);
  orientEdges(10,1,2,-1,edge2hex,hex2edge,edge2node,hex2node,globalOrientation,localOrientation);
  for(auto i =0; i< numHex; ++i){
    std::array<size_t,12> myEdges = hex2edge[i];
    for(auto j = 0; j < 12; ++j){
      size_t myGlobalEdge = myEdges[j];
      orientEdges(myGlobalEdge, j, i, 1, edge2hex, hex2edge, edge2node, hex2node, globalOrientation,
                  localOrientation);
    }
  }

  for(auto i = 0; i < numEdges; ++i){
    const auto myNode = edge2node[i];
    if(globalOrientation[i] != 0){
      //std::cout << "edge=" << "(" << myNode[0]-1 << "," << myNode[1]-1 << ") Ori=" << globalOrientation[i] << '\n';
    }
  }

  std::vector<std::array<double, 3>> dirs;
  std::ofstream f;
  f.open("orient.vtk");
  f << "# vtk DataFile Version 3.0" << '\n';
  f << "vtk output\n";
  f << "ASCII" << '\n';
  f << "DATASET UNSTRUCTURED_GRID" << '\n';
  f << "POINTS " << numEdges << " double\n";
  for(const auto e : edge2node){
    std::array<double, 3> c = {0, 0, 0};
    std::array<double, 3> dir = {0, 0, 0};
    for(auto i = 0; i < 2; ++i){
      size_t n = e[i];
      std::vector<double> coord;
      std::vector<double> paraCoord;
      int d, t;
      gmsh::model::mesh::getNode(n, coord, paraCoord, d, t);

      for(auto j = 0; j < 3; ++j){
        c[j] += 0.5 * coord[j];
        dir[j] = (i==1) ? dir[j] + coord[j] : dir[j] - coord[j];
      }
    }
    dirs.push_back(dir);
    f << c[0] << " " << c[1] << " " << c[2] << '\n';
  }
//  f << "CELLS 0 0\n";
//  f << "CELL_TYPES 0\n";
  f << "POINT_DATA " << numEdges << "\n";
  f << "VECTORS edges double\n";
  for(auto k =0; k < numEdges; ++k){
    auto dir = dirs[k];
    auto go = globalOrientation[k];
        f << dir[0] * go << " " << dir[1] * go << " " << dir[2] * go << '\n';
  }
  f.close();

  std::vector<size_t> localOriginIndex(numHex);
  for(auto hexIndex =0; hexIndex < numHex; ++hexIndex){
    //loop over each node searching for origin
    for(auto i = 0; i < 8; ++i){
      const std::array<int, 6> edgeStemmingFromNode = edgesStemmingFromNode(i);
      size_t originIndex;
      bool isOrigin = true;
      //we know the three edges and their direction stemming from each node
      for(auto j = 0; j < 3; ++j){
        auto myLocalEdge = edgeStemmingFromNode[j];
        auto myDefaultDir = edgeStemmingFromNode[j+3];
        int myLocalOrientation = localOrientation[hexIndex][myLocalEdge];
        //int myTrueDir = myDefaultDir * localOrientation[hexIndex][myLocalEdge];
        //int myTrueDir = myEdgeGlobalOri * myDefaultDir * myLocalOrientation;
        int myTrueDir = myDefaultDir * myLocalOrientation * localEdgeDirection(myLocalEdge);
        isOrigin = isOrigin && (myTrueDir == 1); //if myTrueDir != 1 this is chain is false, it is not the origin
      }
      if(isOrigin){
        originIndex = i;
        localOriginIndex[hexIndex] = originIndex;
        //localOriginIndex[hexIndex] = 0;
        //break;
      }
    }
  }

  //plot origins of each hex
  //std::ofstream f;
  f.open("Origins.vtk");
  f << "# vtk DataFile Version 3.0" << '\n';
  f << "vtk output\n";
  f << "ASCII" << '\n';
  f << "DATASET UNSTRUCTURED_GRID" << '\n';
  f << "POINTS " << numHex << " double\n";
  for(auto hexIndex = 0; hexIndex < numHex; ++hexIndex){
    auto myHexGlobalNodes = hex2node[hexIndex];
    auto myLocalIndex = localOriginIndex[hexIndex];
    auto myGlobalNode = myHexGlobalNodes[myLocalIndex];
    std::vector<double> coord;
    std::vector<double> paraCoord;
    int dim, tag;
    gmsh::model::mesh::getNode(myGlobalNode, coord, paraCoord, dim, tag);
    f << coord[0] << " " << coord[1] << " " << coord[2] << '\n';
  }
  f.close();

  //Create new Hex to node vector that has the order of global nodes arranged to match the unit square
  //i.e. change ordering to match consistent orientation.
  std::vector<Hex> hex2nodeOriented(numHex, Hex{0, 0, 0, 0, 0, 0, 0, 0});
  for(auto h = 0; h < numHex; ++h){
    const auto &myHex = hex2node[h];
    const auto originIndex = localOriginIndex[h];
    auto rearrangedIndices = rearrangeNodes(originIndex);
    for(auto nodeIndex = 0; nodeIndex < 8; ++nodeIndex){
      auto newNodeIndex = rearrangedIndices[nodeIndex];
      auto newGlobalNodeIndex = hex2node[h][newNodeIndex];
      hex2nodeOriented[h][nodeIndex] = newGlobalNodeIndex; //moves global index to correct place to orient it.
    }
  }
  //Comparing adjustments of hex to node with visit visualization.
  for(auto h = 0; h < numHex; ++h){
    auto myNodes = hex2nodeOriented[h];
    auto myNodesOriginal = hex2node[h];
    int myOrigin = hex2node[h][localOriginIndex[h]];
    std::cout << "hexOrigin=" << myOrigin -1 <<   std::endl;
    for(int i = 0; i < 8; ++i){
      int ind = myNodes[i];
      std::cout << ind - 1 << ",";
    }
    std::cout << "\n";
    for(int i = 0; i < 8; ++i){
      int ind = myNodesOriginal[i];
      std::cout << ind - 1 << ",";
    }
    std::cout << "\n";
  }
  return 0;
}

