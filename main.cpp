#include <iostream>
#include <vector>
#include <array>
#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include "EdgeHash.h"
#include "Eigen/Dense"
#include <execution>

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
//      mesh::setTransfiniteCurve(line, n);
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
//  gmsh::option::setNumber("Mesh.Smoothing", 100);
  synchronize();
//  mesh::setTransfiniteSurface(1);
//  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2); // or 3
//  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1); //use 1 to force quads, use 2 to force hexahedra
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::recombine();
//  gmsh::model::mesh::refine();
  synchronize();
  std::vector<std::pair<int, int> > extrusion;
  const std::vector<int> & numElements {n};
//  const std::vector<int> & numElements {2};
  const std::vector<double> & heights {1};
  extrude({std::pair<int,int>{2, squareSurfaces[0]}}, 0, 0, 1, extrusion, numElements, heights, true);
//  extrude({std::pair<int,int>{2, squareSurfaces[0]}}, 0, 0, 1, extrusion, numElements, heights, false);
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
//  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 3); // or 3
//  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 2); //use 1 to force quads, use 2 to force hexahedra

//Use the below to make
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::recombine();
  gmsh::model::mesh::generate(3);
//  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 2); //use 1 to force quads, use 2 to force hexahedra
  size_t numSubDivs = 0;
  for(auto s = 0; s < numSubDivs; ++s){
    synchronize();
    gmsh::model::mesh::refine();
  }

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

std::string meshUnitSquare(const int n, const double sideLength) {
  using namespace gmsh::model::geo;
  using DArray = std::array<double, 2>;
  gmsh::initialize();
  //const double sideLength = 1.0;
  const std::array<DArray, 4> corners{DArray{0., 0.}, DArray{sideLength, 0.},
                                      DArray{sideLength, sideLength}, DArray{0., sideLength}};
  const double dx = sideLength / n;
  std::vector<int> geoPoints;
  //loop to add corners as gmsh points
  for (const auto &xy: corners) {
    geoPoints.push_back(addPoint(xy[0], xy[1], 0, dx));
    //geoPoints.push_back(addPoint(xy[0], xy[1], 0));
  }
  std::vector<int> geoLines;
  //connect gmsh points as gmsh lines
  for (auto i = 0; i < geoPoints.size() - 1; ++i) {
    const auto line = addLine(geoPoints[i], geoPoints[i + 1]);
    geoLines.push_back(line);
    //const double scale = 1.+5e-2;
    const double scale = 1.1;
    mesh::setTransfiniteCurve(line, n+1);
//    if(i==0) mesh::setTransfiniteCurve(line, n+1, "Linear", scale);
//    if(i==2) mesh::setTransfiniteCurve(line, n+1, "Linear", -scale);
  }
  geoLines.push_back(addLine(geoPoints[3], geoPoints[0]));
  int diag = addLine(geoPoints[0], geoPoints[2]);
  auto loop = addCurveLoop(geoLines);
  auto surface = addPlaneSurface(std::vector<int>{loop});
  synchronize();
  std::vector<int> diagAsVector {diag};
  gmsh::model::mesh::embed(1, diagAsVector, 2, surface);
  synchronize();

  mesh::setTransfiniteSurface(surface);

  //gmsh::option::setNumber("Mesh.Algorithm", 9);
  //gmsh::option::setNumber("Mesh.Smoothing", 10);
  //const double phi = 10. * M_PI / 180.;
  //const double phi = 45. * M_PI / 180.;
//  const double phi = 0;
//  gmsh::model::geo::rotate({{2, surface},{1,loop},{0,geoPoints[0]},{0,geoPoints[1]},{0,geoPoints[2]},{0,geoPoints[3]}}, 0, 0, 0, 0, 0, 1, phi);
  synchronize();

  //gmsh::model::mesh::generate(2);
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::recombine();
  synchronize();
  std::vector<std::pair<int, int> > extrusion;
  const std::vector<int> & numElements {n+5};
//  const std::vector<int> & numElements {n+5};
  const std::vector<double> & heights {1};
  extrude({std::pair<int,int>{2, surface}}, 0, 0, 1, extrusion, numElements, heights, true);
//  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 3); // or 3
  synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::recombine();
  //std::string meshFormat{".vtk"};
  gmsh::model::mesh::generate(3);
  std::string meshFormat{".msh"};
  //std::string meshFormat{".m"};

  auto fileName = "UnitSquareUnstructured" + std::to_string(n) + meshFormat;
  //auto fileName = "UnitSquareStructured" + std::to_string(n) + meshFormat;
  //gmsh::write("UnitSquare" + std::to_string(n) + ".msh");
  gmsh::write(fileName);
  //gmsh::write("UnitSquareStructured" + std::to_string(n) + ".m");
  gmsh::write("UnitSquareUnstructured" + std::to_string(n) + ".m");
  gmsh::finalize();
  return fileName;
}

class UnitCube {
public:
  const size_t mIndex;
  const Hex mNodes;
  const std::vector<double> mNodeCoords;
  std::array<double, 6> mE {0, 0, 0, 0, 0, 0};
  std::array<double, 6> mD {0, 0, 0, 0, 0, 0};
  std::array<double, 12> mH {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::array<double, 12> mB {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::array<size_t, 12> mGlobalEdges {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::array<size_t, 6> mGlobalFaces {0, 0, 0, 0, 0, 0};
  Eigen::Matrix3d mMuPhysical = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d mEpsPhysical = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d mEpsCompInv = Eigen::Matrix3d::Identity();
  std::array<size_t, 6> mFaceBoundaryFlag {0, 0, 0, 0, 0, 0};
  UnitCube(size_t ind, Hex nodes, std::vector<double> &nodeCoords, std::array<size_t, 12> globalEdges) :
    mIndex(ind),
    mNodes(nodes),
    mNodeCoords(nodeCoords),
    mGlobalEdges(globalEdges){
  }
  UnitCube(size_t ind, Hex nodes, std::vector<double> &nodeCoords, std::array<size_t, 12> globalEdges,
           Eigen::Matrix3d muBar, Eigen::Matrix3d epsBar) :
    mIndex(ind),
    mNodes(nodes),
    mNodeCoords(nodeCoords),
    mGlobalEdges(globalEdges){
    mMuPhysical = muBar;
    mEpsPhysical = epsBar;
    auto jac = jacobian({.5, .5, .5});
    mEpsCompInv = jac.transpose() * epsBar.inverse() * jac / jac.determinant();
  }
  UnitCube(size_t ind, Hex nodes, std::vector<double> &nodeCoords, std::array<size_t, 12> globalEdges,
           std::array<size_t, 6> globalFaces, Eigen::Matrix3d muBar, Eigen::Matrix3d epsBar) :
      mIndex(ind),
      mNodes(nodes),
      mNodeCoords(nodeCoords),
      mGlobalFaces(globalFaces),
      mGlobalEdges(globalEdges){
    mMuPhysical = muBar;
    mEpsPhysical = epsBar;
    auto jac = jacobian({.5, .5, .5});
    mEpsCompInv = jac.transpose() * epsBar.inverse() * jac / jac.determinant();
  }

  /*Unit square to hex map.
   * See https://math.stackexchange.com/questions/2265255/mapping-a-3d-point-inside-a-hexahedron-to-a-unit-cube
   */
  Eigen::Vector3d mapToHex(const Eigen::Vector3d unitCubeCoord){
    const auto x= unitCubeCoord[0];
    const auto y = unitCubeCoord[1];
    const auto z = unitCubeCoord[2];
    std::vector<Eigen::Vector3d> p; //corners
    for(auto i = 0; i < mNodeCoords.size(); i+=3){
      p.push_back(Eigen::Vector3d{mNodeCoords[i], mNodeCoords[i+1], mNodeCoords[i+2]});
    }
    Eigen::Vector3d result {0, 0, 0};
    result += (1 - x) * (1 - y) * (1 - z) * p[0];
    result += (1 - x) * y * (1 - z) * p[1];
    result += (1 - x) * y * z * p[2];
    result += (1 - x) * (1 - y) * z * p[3];
    result += x * (1 - y) * (1 - z) * p[4];
    result += x * y * (1 - z) * p[5];
    result += x * y * z * p[6];
    result += x * (1 - y) * z * p[7];
    return result;
  }

  Eigen::Matrix3d jacobian(const Eigen::Vector3d unitCubeCoord){
    const auto x= unitCubeCoord[0];
    const auto y = unitCubeCoord[1];
    const auto z = unitCubeCoord[2];
    std::vector<Eigen::Vector3d> p; //corners
    for(auto i = 0; i < mNodeCoords.size(); i+=3){
      p.push_back(Eigen::Vector3d{mNodeCoords[i], mNodeCoords[i+1], mNodeCoords[i+2]});
    }
    Eigen::Vector3d partialX {0, 0, 0};
    partialX += -(1 - y) * (1 - z) * p[0];
    partialX += -y * (1 - z) * p[1];
    partialX += -y * z * p[2];
    partialX += -(1 - y) * z * p[3];
    partialX += (1 - y) * (1 - z) * p[4];
    partialX += y * (1 - z) * p[5];
    partialX += y * z * p[6];
    partialX += (1 - y) * z * p[7];

    Eigen::Vector3d partialY {0, 0, 0};
    partialY += -(1 - x) * (1 - z) * p[0];
    partialY += (1 - x) * (1 - z) * p[1];
    partialY += (1 - x) * z * p[2];
    partialY += -(1 - x) * z * p[3];
    partialY += -x * (1 - z) * p[4];
    partialY += x * (1 - z) * p[5];
    partialY += x * z * p[6];
    partialY += -x * z * p[7];

    Eigen::Vector3d partialZ {0, 0, 0};
    partialZ += -(1 - x) * (1 - y) * p[0];
    partialZ += -(1 - x) * y * p[1];
    partialZ += (1 - x) * y * p[2];
    partialZ += (1 - x) * (1 - y) * p[3];
    partialZ += -x * (1 - y) * p[4];
    partialZ += -x * y * p[5];
    partialZ += x * y * p[6];
    partialZ += x * (1 - y) * p[7];

    Eigen::Matrix3d result;
    result << partialX, partialY, partialZ;
    return result;
  }

/*
  edges                  back faces         front faces
       2
 +----------+         z+----------+       z+----------+
 |\         |\         |\         |         \         |\
 | \11      | \10      | \        |          \    5   | \
3|  \    6  |1 \       |  \  0    |           \       |  \
 |   +------+---+      |   +                   +------+---+
 |   | 0    |   |      | 1 |      |            |      | 4 |
 +---+------+   |      +---+------+y           |      +y  |
  \ 7|       \  |5      \  |       \           |    3  \  |
  8\ |       9\ |        \ |   2    \          |        \ |
    \|         \|         \|         \         |         \|
     +----------+         x+----------+       x+----------+
          4
*/
  std::array<size_t, 2> facesAdjacentToEdge(size_t edgeIndex){
    switch(edgeIndex){
      case 0:
        return {0, 2};
      case 1:
        return {0, 4};
      case 2:
        return {0, 5};
      case 3:
        return {0, 1};
      case 4:
        return {2, 3};
      case 5:
        return {3, 4};
      case 6:
        return {3, 5};
      case 7:
        return {1, 3};
      case 8:
        return {1, 2};
      case 9:
        return {2, 4};
      case 10:
        return {4, 5};
      case 11:
        return {1, 5};
      default:
        std::cout << "Bad index passed to facesAdjacentToEdge()" << std::endl;
        return {0, 0};
    }
  }



  std::array<size_t, 4> edgesAdjacentToFace(size_t faceIndex){
    switch(faceIndex){
      case 0:
        return {0, 1, 2, 3};
      case 1:
        return {3, 7, 8, 11};
      case 2:
        return {0, 4, 8, 9};
      case 3:
        return {4, 5, 6, 7};
      case 4:
        return {1, 5, 9, 10};
      case 5:
        return {2, 6, 10, 11};
      default:
        std::cout << "Bad Index=" << faceIndex << " in edgesAdjacentToFace()" << std::endl;
        return {0, 0, 0, 0};
    }
  }

  std::array<int, 4> edgesOrientationToFace(size_t faceIndex){
    switch(faceIndex){
      case 0:
        return {1, 1, -1, -1};
      case 1:
        return {1, -1, -1, 1};
      case 2:
        return {-1, 1, 1, -1};
      case 3:
        return {1, 1, -1, -1};
      case 4:
        return {1, -1, -1, 1};
      case 5:
        return {-1, 1, -1, 1};
      default:
        std::cout << "Bad Index=" << faceIndex << " in edgesOrientationToFace()" << std::endl;
        return {0, 0, 0, 0};
    }
  }


    std::array<int, 2> facesOrientationToEdge(size_t edgeIndex){
    switch(edgeIndex){
      case 0:
        return {1, -1};
      case 1:
        return {1, 1};
      case 2:
        return {-1, -1};
//        return {1, 1};
      case 3:
        return {-1, 1};
      case 4:
        return {1, 1};
      case 5:
        return {1, -1};
      case 6:
        return {-1, 1};
      case 7:
        return {-1, -1};
      case 8:
        return {-1, 1};
      case 9:
        return {-1, -1};
      case 10:
        return {1, -1};
      case 11:
        return {1, 1};
      default:
        std::cout << "Bad index passed to facesOrientationToEdge()" << std::endl;
        return {0, 0};
    }
  }

  Eigen::Vector3d faceCoord(const size_t faceIndex){
    Eigen::Vector3d coord {.5, .5, .5};
    //for faceIndex = 0, 1, 2 the cases take for {0, .5, .5}, {.5, 0, .5}, and {.5, .5, 0}
    if(faceIndex < 3){
      coord[faceIndex] = 0;
    }
    //for faceIndex = 0, 1, 2 the cases take for {1, .5, .5}, {.5, 1, .5}, and {.5, .5, 1}
    else if(3 <= faceIndex && faceIndex < 6){
      coord[faceIndex-3] = 1;
    }
    else{
      std::cout << "Bad index in faceCoord()" << std::endl;
      coord = {0, 0, 0};
    }
    return coord;
  }
  Eigen::Vector3d faceNormal(const size_t faceIndex){
    Eigen::Vector3d normal {0, 0, 0};
    //for faceIndex = 0, 1, 2 the cases take for {0, .5, .5}, {.5, 0, .5}, and {.5, .5, 0}
    if(faceIndex < 3){
      normal[faceIndex] = 1;
    }
      //for faceIndex = 0, 1, 2 the cases take for {1, .5, .5}, {.5, 1, .5}, and {.5, .5, 1}
    else if(3 <= faceIndex && faceIndex < 6){
      normal[faceIndex - 3] = 1;
    }
    else{
      std::cout << "Bad index in faceNormal()" << std::endl;
      normal = {0, 0, 0};
    }
    return normal;
  }
  Eigen::Vector3d physicalFaceAreaWeighted(size_t faceIndex){
    Eigen::Matrix3d jac = faceJacobian(faceIndex);
    return jac.transpose().inverse() * faceNormal(faceIndex) * jac.determinant();
  }
  double physicalFaceArea(size_t faceIndex){
    return physicalFaceAreaWeighted(faceIndex).norm();
  }
  /*
   * Inputs: edgeIndex - local edge index in hex
   * Returns the local coordinates of an edge in the unit square given a local edge index
   * See the edge index convention in facesAdjacentToEdge() function
  */
  Eigen::Vector3d edgeCoord(const size_t edgeIndex){
    switch(edgeIndex){
      case(0):
        return {0, .5, 0};
      case(1):
        return {0, 1, .5};
      case(2):
        return {0, .5, 1};
      case(3):
        return {0, 0, .5};
      case(4):
        return {1, .5, 0};
      case(5):
        return {1, 1, .5};
      case(6):
        return {1, .5, 1};
      case(7):
        return {1, 0, .5};
      case(8):
        return {.5, 0, 0};
      case(9):
        return {.5, 1, 0};
      case(10):
        return {.5, 1, 1};
      case(11):
        return {.5, 0, 1};
      default:
        std::cout << "Bad index in edgeCoord()" << std::endl;
        return {0, 0, 0};
    }
  }

  Eigen::Vector3d edgePhysicalCoord(size_t edgeIndex){
    return mapToHex(edgeCoord(edgeIndex));
  }

  Eigen::Vector3d dualFaceCoord(const size_t edgeIndex) {
    switch (edgeIndex) {
      case (0):
        return {.25, .5, .25};
      case (1):
        return {.25, .75, .5};
      case (2):
        return {.25, .5, .75};
      case (3):
        return {.25, .25, .5};
      case (4):
        return {.75, .5, .25};
      case (5):
        return {.75, .75, .5};
      case (6):
        return {.75, .5, .75};
      case (7):
        return {.75, .25, .5};
      case (8):
        return {.5, .25, .25};
      case (9):
        return {.5, .75, .25};
      case (10):
        return {.5, .75, .75};
      case (11):
        return {.5, .25, .75};
      default:
        std::cout << "Bad index in edgeCoord()" << std::endl;
        return {0, 0, 0};
    }
  }

  Eigen::Vector3d edgeTangential(const size_t edgeIndex){
    bool isX = (edgeIndex == 8 || edgeIndex == 9 || edgeIndex == 10 ||edgeIndex == 11);
    bool isY = (edgeIndex == 0 || edgeIndex == 2 || edgeIndex == 4 ||edgeIndex == 6);
    bool isZ = (edgeIndex == 1 || edgeIndex == 3 || edgeIndex == 5 ||edgeIndex == 7);
    if(isX)
      return {1, 0, 0};
    else if(isY)
      return {0, 1, 0};
    else if(isZ)
      return {0, 0, 1};
    else{
      std::cout << "Bad index in edgeTangential()" << std::endl;
      return {0, 0, 0};
    }
  }

  std::array<Eigen::Vector3d,2> edgeNormals(size_t edgeIndex){
    bool isX = (edgeIndex == 8 || edgeIndex == 9 || edgeIndex == 10 ||edgeIndex == 11);
    bool isY = (edgeIndex == 0 || edgeIndex == 2 || edgeIndex == 4 ||edgeIndex == 6);
    bool isZ = (edgeIndex == 1 || edgeIndex == 3 || edgeIndex == 5 ||edgeIndex == 7);
    Eigen::Vector3d x {1,0,0};
    Eigen::Vector3d y {0,1,0};
    Eigen::Vector3d z {0,0,1};
    if(isX)
      return {y,z};
    else if(isY)
      return {z,x};
    else if(isZ)
      return {x,y};
    else{
      std::cout << "Bad index in edgeNormals()" << std::endl;
      return {x, x};
    }

  };

  int edgeTangentialAsIndex(const size_t edgeIndex){
    bool isX = (edgeIndex == 8 || edgeIndex == 9 || edgeIndex == 10 ||edgeIndex == 11);
    bool isY = (edgeIndex == 0 || edgeIndex == 2 || edgeIndex == 4 ||edgeIndex == 6);
    bool isZ = (edgeIndex == 1 || edgeIndex == 3 || edgeIndex == 5 ||edgeIndex == 7);
    if(isX)
      return 0;
    else if(isY)
      return 1;
    else if(isZ)
      return 2;
    else{
      std::cout << "Bad index in edgeTangentialAsIndex()" << std::endl;
      return -1;
    }
  }

  /*
   * input: local face index
   * returns index of opposite face in unit cube
   */
  int getOppositeFace(const int faceIndex){
    return (faceIndex + 3) % 6;
//    int oppositeIndex = (static_cast<int>(faceIndex + 3) % 6);
//    return static_cast<size_t>(oppositeIndex);
  }

  double getDAtCirculationLocation(const int faceIndex){
    int oppFace = getOppositeFace(faceIndex);
    return 0.75 * mD[faceIndex] + 0.25 * mD[oppFace];
  }

  void setMaterials(Eigen::Matrix3d muBar, Eigen::Matrix3d epsBar){
    mMuPhysical = muBar;
    mEpsPhysical = epsBar;
    auto jac = jacobian({.5, .5, .5});
    mEpsCompInv = jac.transpose() * epsBar.inverse() * jac / jac.determinant();
  }

  Eigen::Matrix3d getMuCompInv(size_t edgeIndex){
    auto myEdgeCoord = edgeCoord(edgeIndex);
    auto jac = jacobian(myEdgeCoord);
    auto muInv = jac.transpose() * mMuPhysical.inverse() * jac / jac.determinant();
    return muInv;
  }
  Eigen::Matrix3d getMuCompInvAtCenter(){
    auto jac = jacobian({.5,.5,.5});
    auto muInv = jac.transpose() * mMuPhysical.inverse() * jac / jac.determinant();
    return muInv;
  }

  Eigen::Matrix3d getEpsCompInv(size_t faceIndex){
    auto jac = circulationLocationJacobian(faceIndex);
    auto epsInv = jac.transpose() * mEpsPhysical.inverse() * jac / jac.determinant();
    return epsInv;
  }

  Eigen::Matrix3d getMuCompInvAtDualFace(size_t edgeIndex){
    auto jac = dualFaceJacobian(edgeIndex) ;
    auto muInv = jac.transpose() * mMuPhysical.inverse() * jac / jac.determinant();
    return muInv;
  }

  Eigen::Matrix3d getMuComp(size_t edgeIndex){
    auto myEdgeCoord = edgeCoord(edgeIndex);
    auto jac = jacobian(myEdgeCoord);
    auto muComp = jac.inverse() * mMuPhysical * jac.inverse().transpose() * jac.determinant();
    return muComp;
  }

  Eigen::Vector3d facePhysicalCoordinate(size_t faceIndex){
    auto localCoord = faceCoord(faceIndex);
    return mapToHex(localCoord);
  }

  Eigen::Vector3d hexCenterPhysicalCoordinate(){
    Eigen::Vector3d localCoord {.5, .5, .5};
    return mapToHex(localCoord);
  }

  Eigen::Matrix3d edgeJacobian(size_t edgeIndex){
    Eigen::Vector3d coord = edgeCoord(edgeIndex);
    return jacobian(coord);
  }

  Eigen::Matrix3d dualFaceJacobian(size_t edgeIndex){
    Eigen::Vector3d coord = dualFaceCoord(edgeIndex);
    return jacobian(coord);
  }

  double dualFacePhysicalArea(size_t edgeIndex){
    auto phyCenter = hexCenterPhysicalCoordinate();
    auto phyEdgeCoord = mapToHex(edgeCoord(edgeIndex));
    auto adjFaceInd = facesAdjacentToEdge(edgeIndex);
    std::vector<Eigen::Vector3d> faceCoords;
    for(auto faceInd : adjFaceInd){
      auto globalFaceCoord = mapToHex(faceCoord(faceInd));
      faceCoords.push_back(globalFaceCoord);
    }
    Eigen::Vector3d sideFromEdge0 = faceCoords[0]-phyEdgeCoord;
    Eigen::Vector3d sideFromEdge1 = faceCoords[1]-phyEdgeCoord;
    Eigen::Vector3d sideFromCenter0 = faceCoords[0]-phyCenter;
    Eigen::Vector3d sideFromCenter1 = faceCoords[1]-phyCenter;
    return 0.5 * (sideFromEdge0.cross(sideFromEdge1).norm() + sideFromCenter0.cross(sideFromCenter1).norm());
  }

  Eigen::Matrix3d faceJacobian(size_t faceIndex){
    Eigen::Vector3d coord = faceCoord(faceIndex);
    return jacobian(coord);
  }

  Eigen::Matrix3d circulationLocationJacobian(size_t faceIndex){
    auto oppositeFace = getOppositeFace(faceIndex);
    Eigen::Vector3d coord = faceCoord(faceIndex);
    Eigen::Vector3d oppositeCoord = faceCoord(oppositeFace);
    Eigen::Vector3d circulationCoord = .75 * coord + 0.25 * oppositeCoord;
    return jacobian(circulationCoord);
  }

  Eigen::Vector3d getDAtCenter(){
    Eigen::Vector3d dCenter {0, 0, 0};
    for (int f = 0; f < 6; ++f) {
      int dir = f % 3;
      dCenter[dir] += 0.5 * mD[f];
    }
    return dCenter;
  }

  Eigen::Vector3d getDVectorAtCirculationLocation(int faceIndex){
    int dir = faceIndex % 3;
    double dScale = getDAtCirculationLocation(faceIndex);
    auto dVecCirc = getDAtCenter();
    dVecCirc[dir] = dScale;
    return dVecCirc;
  }
};

class UnitEdge{
public:
  const size_t mIndex;
  const Edge mNodes;
  const std::vector<double> mNodeCoords;
  const std::vector<size_t> mAdjacentHex {};
  const std::vector<size_t> mIndexInHex {};
  const int mOrientation;
  double mH = 0;
  double mB = 0;
  double mBFlux = 0;
  double mMu = 1;
  int mBoundaryFlag = 0;

  UnitEdge(size_t ind, Edge nodes, std::vector<double> &nodeCoords, std::vector<size_t> &adjacentHex,
           std::vector<size_t> &indexInHex, int orientation) :
    mIndex(ind),
    mNodes(nodes),
    mNodeCoords(nodeCoords),
    mAdjacentHex(adjacentHex),
    mIndexInHex(indexInHex),
    mOrientation(orientation){
  }

  Eigen::Vector3d asUnitVector(){
    Eigen::Vector3d vec {mNodeCoords[3]-mNodeCoords[0],
                         mNodeCoords[4]-mNodeCoords[1],
                         mNodeCoords[5]-mNodeCoords[2]};
    return mOrientation * vec / vec.norm();
  }

  double length(){
    Eigen::Vector3d vec {mNodeCoords[3]-mNodeCoords[0],
                         mNodeCoords[4]-mNodeCoords[1],
                         mNodeCoords[5]-mNodeCoords[2]};
    return vec.norm();
  }

  Eigen::Vector3d physicalCoordinate(){
    Eigen::Vector3d vec {mNodeCoords[3]+mNodeCoords[0],
                         mNodeCoords[4]+mNodeCoords[1],
                         mNodeCoords[5]+mNodeCoords[2]};
    return 0.5 * vec;
  }

};

class UnitFace {
public:
  const size_t mIndex;
//  const std::array<size_t, 4> mNodes;
//  const std::vector<double> mNodeCoords;
//  const std::vector<size_t> mAdjacentHex{};
//  const std::vector<size_t> mIndexInHex{};
  std::vector<std::pair<size_t, size_t>> mAdjHexAndIndexInHex;
//  const int mOrientation;
  double mDFlux = 0;
  UnitFace(size_t index, std::vector<std::pair<size_t, size_t>> &adjHexAndIndexInHex) :
             mIndex(index),
             mAdjHexAndIndexInHex(adjHexAndIndexInHex){}
};

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

std::array<size_t, 4> faceIndexToNodeIndex(size_t faceIndex){
  switch(faceIndex){
    case 0:
      return {0,1,2,3};
    case 1:
      return {0,3,4,7};
    case 2:
      return {0,1,4,5};
    case 3:
      return{4,5,6,7};
    case 4:
      return {1,2,5,6};
    case 5:
      return {2,3,6,7};
    default:
      std::cout << "Bad index in faceIndexToNodeIndex()\n";
      return {0,0,0,0};
  }
}

void parseFaces(std::vector<std::array<size_t, 4>> &quad2node, const std::vector<Hex> &hex2node, std::vector<std::array<size_t, 6>> &hex2face,
                std::vector<std::vector<std::pair<size_t, size_t>>> &adjHexAndIndexInHex){
  std::map<std::array<size_t,4>, size_t> faceMap;
  size_t globalFaceIndex = 0;
  for(size_t h = 0; h < hex2node.size(); ++h){
    const auto &cube = hex2node[h];
    for(size_t f = 0; f < 6; ++f){
      //create face as an ordered tuple of node indices
      std::array<size_t, 4> localFaceNodes = faceIndexToNodeIndex(f);
      std::array<size_t, 4> globalFaceNodes;
      //fill global face nodes array using hex2node
      for(auto i = 0; i < 4; ++i){
        size_t localNode = localFaceNodes[i];
        globalFaceNodes[i] = cube[localNode];
      }
      std::sort(globalFaceNodes.begin(), globalFaceNodes.end());
      quad2node.push_back(globalFaceNodes);
      auto [it, wasInserted] = faceMap.emplace(globalFaceNodes, globalFaceIndex);
      if(wasInserted){//if face doesn't exist yet
        std::vector<std::pair<size_t, size_t>> adjHexFaceInd {std::make_pair(h, f)};
        adjHexAndIndexInHex.push_back(adjHexFaceInd);
//        face2hex.push_back(std::vector<size_t>{h});
//        indexInHex.push_back(std::vector<size_t>{f});
        hex2face[h][f] = globalFaceIndex;
        ++globalFaceIndex;
      }
      else{ //if face already exists
        const size_t existingFaceIndex = faceMap[globalFaceNodes];
        adjHexAndIndexInHex[existingFaceIndex].push_back(std::make_pair(h, f));
//        face2hex[existingFaceIndex].push_back(h);
//        indexInHex[existingFaceIndex].push_back(f);
        hex2face[h][f] = existingFaceIndex;
      }
    }
  }
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

//std::vector<size_t> localParallelEdgeCounterClockwise(size_t edgeIndex){
//  std::vector<size_t> result;
//  if(!(0 <= edgeIndex && edgeIndex < 12)){
//    std::cout << "Edge index passed in get getLocalParallelEdge needs to be in [0, 12)." << std::endl;
//    return result;
//  }
//  std::vector<std::vector<size_t>> parallelEdges = {{0, 2, 6, 4},
//                                                    {1, 3, 7, 5},
//                                                    {8, 9, 10, 11}};
//  //loop over sets of parallel edges and find set that contains edgeIndex
//  for(auto pEdges : parallelEdges){
//    auto it = std::find(pEdges.begin(), pEdges.end(), edgeIndex);
//    //if the set pEdges contains edgeIndex, remove current index and set result to the remaining three edges
//    if(it != std::end(pEdges)){
//      std::rotate(pEdges.begin(), it, pEdges.end());
//      result = pEdges;
//      break;
//    }
//  }
//  return result;
//}

std::array<size_t, 4> localParallelEdgeCounterClockwise(size_t edgeIndex){
  switch(edgeIndex) {
    case 0:
      return {0, 2, 6, 4};
    case 1:
      return {1, 3, 7, 5};
    case 2:
      return {2, 6, 4, 0};
    case 3:
      return {3, 7, 5, 1};
    case 4:
      return {4, 0, 2, 6};
    case 5:
      return {5, 1, 3, 7};
    case 6:
      return {6, 4, 0, 2};
    case 7:
      return {7, 5, 1, 3};
    case 8:
      return {8, 9, 10, 11};
    case 9:
      return {9, 10, 11, 8};
    case 10:
      return {10, 11, 8, 9};
    case 11:
      return {11, 8, 9, 10};
    default:
      std::cout << "Edge index passed in get getLocalParallelEdgeCounterClockwise needs to be in [0, 12)."
                << std::endl;
      return {0,0,0,0};
  }
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
      return {1, 2, 3, 0,  5, 6, 7, 4}; //working
//      return {1, 0, 4, 5,  2, 3, 7, 6};
//      return {1, 5, 6, 2, 0, 4, 7, 3};
    case 2:
//      return {2, 6, 7, 3,  1, 5, 4, 0};
//      return {2, 1, 5, 6,  3, 0, 4, 7};
      return {2, 3, 0, 1, 6, 7, 4, 5};
    case 3:
      return {3, 0, 1, 2,  7, 4, 5, 6}; //working
//      return {3, 7, 4, 0, 2, 6, 5, 1};
//      return {3, 2, 6, 7,  0, 1, 5, 4};
    case 4:
//      return {4, 5, 1, 0,  7, 6, 2, 3};
//      return {4, 0, 3, 7,  5, 1, 2, 6};
      return {4, 7, 6, 5,  0, 3, 2, 1}; //working
    case 5:
//      return {5, 6, 2, 1,  4, 7, 3, 0};
      return {5, 4, 7, 6, 1, 0, 3, 2}; //working
//      return {5, 1, 0, 4, 6, 2, 3, 7};
    case 6:
      return {6, 5, 4, 7,  2, 1, 0, 3}; //working
//      return {6, 7, 3, 2, 5, 4, 0, 1};
//      return {6, 2, 1, 5, 7, 3, 0, 4};
    case 7:
//      return {7, 4, 0, 3,  6, 5, 1, 2};
//      return {7, 3, 2, 6, 4, 0, 1, 5};
      return {7, 6, 5, 4, 3, 2, 1, 0}; //last
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

std::array<size_t, 4> edgesAdjacentToEdge(size_t edgeIndex) {
  switch (edgeIndex) {
    case 0:
      return {1,3,8,9};
    case 1:
      return{0,2,9,10};
    case 2:
      return {1,3,10,11};
    case 3:
      return {0,2,8,11};
    case 4:
      return {5,8,7,9};
    case 5:
      return {4,6,9,10};
    case 6:
      return {5,7,10,11};
    case 7:
      return {4,6,8,11};
    case 8:
      return {0,3,4,7};
    case 9:
      return {0,1,4,5};
    case 10:
      return {1,2,5,6};
    case 11:
      return {2,3,6,7};
    default:
      std::cout << "Bad index in edgesAdjacentToEdge()" << '\n';
      return {0, 0, 0, 0};
  }
}

void plotCirculation(const UnitEdge &edge, std::vector<UnitCube> &cubes, std::string fName){
  const auto myAdjHex = edge.mAdjacentHex;
  const auto indInAdjHex = edge.mIndexInHex;
  std::ofstream f;
  //f.open("EdgeCirculation.vtk");
  f.open(fName);
  f << "# vtk DataFile Version 3.0" << '\n';
  f << "vtk output\n";
  f << "ASCII" << '\n';
  f << "DATASET UNSTRUCTURED_GRID" << '\n';
  f << "POINTS " << 2*myAdjHex.size()+1 << " double\n";

  std::string vectorData = "";
  //loop of adj. hexahedra and write circulation locations to vtk format
  size_t nodeIndex = 0;
  for(auto i = 0; i < myAdjHex.size(); ++i){
    auto hexInd = myAdjHex[i];
    auto localEdgeInd = indInAdjHex[i];
    auto &myHex = cubes[hexInd];
    auto adjFaces = myHex.facesAdjacentToEdge(localEdgeInd);
    auto adjFacesOri = myHex.facesOrientationToEdge(localEdgeInd);
    for(auto k = 0; k < adjFaces.size(); ++k){
      auto faceIndex = adjFaces[k];
      int ori = adjFacesOri[k];
      auto faceLocalCoord = myHex.faceCoord(faceIndex);
      int oppositeInd = (static_cast<int>(faceIndex) + 3) % 6;
      auto faceLocalCoordOpposite = myHex.faceCoord(oppositeInd);
      std::cout << faceLocalCoord << std::endl;
      faceLocalCoord = 0.75 * faceLocalCoord + 0.25 * faceLocalCoordOpposite;
      std::cout << faceLocalCoord << "\n\n";
      auto globalCoord = myHex.mapToHex(faceLocalCoord);
      for(auto c : globalCoord){
        f << c << ' ';
      }
      f << '\n';
      nodeIndex++;
      auto faceNormal = myHex.faceNormal(faceIndex);
      auto faceNormalOri = faceNormal * static_cast<double>(ori);
      auto jac = myHex.jacobian(faceLocalCoord);
      auto globalNormal0 = jac.inverse().transpose() * faceNormalOri;
      auto globalNormal = globalNormal0 / globalNormal0.norm();
      for(auto c : globalNormal){
        vectorData += std::to_string(c) + " ";
      }
      vectorData += "\n";
    }
  }
  //Chunk of code to write edge coord and vector to file
  Eigen::Vector3d edgeCoord {0, 0, 0};
  Eigen::Vector3d edgeVec {0, 0, 0};
  for(auto i = 0; i < 6; i+=3){
    Eigen::Vector3d c {edge.mNodeCoords[i], edge.mNodeCoords[i+1], edge.mNodeCoords[i+2]};
    edgeCoord += 0.5 * c;
    if(i < 3){
      edgeVec -= c;
    }
    else{
      edgeVec += c;
    }
  }
  edgeVec = 10. * edgeVec * edge.mOrientation;
  for(auto i = 0; i < 3; ++i){
    f << edgeCoord[i] << " ";
    vectorData += std::to_string(edgeVec[i]) + " ";
  }
  nodeIndex++;

  f << '\n';
  f << "POINT_DATA " << 2* myAdjHex.size()+1 << "\n";
  f << "VECTORS faces double\n";
  f << vectorData;

  f.close();
}

void plotCirculation2(const UnitEdge &edge, std::vector<UnitCube> &cubes, std::string fName) {
  const auto myAdjHex = edge.mAdjacentHex;
  const auto indInAdjHex = edge.mIndexInHex;

  Eigen::Vector3d edgeCoord {0, 0, 0};
  Eigen::Vector3d edgeVecUnnormalized {edge.mNodeCoords[3]-edge.mNodeCoords[0],
                                       edge.mNodeCoords[4]-edge.mNodeCoords[1],
                                       edge.mNodeCoords[5]-edge.mNodeCoords[2]};
  Eigen::Vector3d edgeVec = edgeVecUnnormalized / edgeVecUnnormalized.norm();

  std::ofstream f;
  f.open(fName);
  f << "# vtk DataFile Version 3.0" << '\n';
  f << "vtk output\n";
  f << "ASCII" << '\n';
  f << "DATASET UNSTRUCTURED_GRID" << '\n';
  f << "POINTS " << 3 * myAdjHex.size() + 2 << " double\n";

  std::string vectorData = "";
  std::string elementData = "";
  //loop of adj. hexahedra and write circulation locations to vtk format
  size_t nodeIndex = 0;
  for (auto i = 0; i < myAdjHex.size(); ++i) {
    auto hexInd = myAdjHex[i];
    auto localEdgeInd = indInAdjHex[i];
    auto &myHex = cubes[hexInd];
    auto adjFaces = myHex.facesAdjacentToEdge(localEdgeInd);
    auto adjFacesOri = myHex.facesOrientationToEdge(localEdgeInd);
    auto myEdgeTangent = myHex.edgeTangential(localEdgeInd);
    auto myEdgeTangentCoord = myHex.edgeCoord(localEdgeInd);
    auto edgeJac = myHex.jacobian(myEdgeTangentCoord);
//    auto edgeGlobal = edgeVec.transpose() * edgeJac.inverse().transpose() * myEdgeTangent;
    auto edgeGlobal = (2*myEdgeTangent.transpose()) * (edgeJac.inverse() * edgeVec * edgeJac.determinant());
//    auto edgeGlobal = (edgeVec.transpose()) * edgeJac * myEdgeTangent / edgeJac.determinant();
//    auto edgeGlobal = 0.53033 * edgeJac * myEdgeTangent / edgeJac.determinant();
//    auto edgeGlobal = (myEdgeTangent.transpose()) * (edgeJac.transpose() * edgeVec);
//    std::cout << '\n' << "hex=" << myHex.mIndex << '\n' << edgeGlobal << '\n' << edgeVec << '\n';
//    std::cout << '\n' << "hex=" << myHex.mIndex << '\n' << edgeGlobal << '\n';
    auto centerCoord = myHex.mapToHex(Eigen::Vector3d{.5, .5, .5});
    auto centerCoordVtkInd = nodeIndex;
    for (auto c: centerCoord) {
      f << c << ' ';
    }
    f << "\n";
    ++nodeIndex;
    for (auto k = 0; k < adjFaces.size(); ++k) {
      auto faceIndex = adjFaces[k];
      int ori = adjFacesOri[k];
      auto faceLocalCoord = myHex.faceCoord(faceIndex);
      int oppositeInd = (static_cast<int>(faceIndex) + 3) % 6;
      auto faceLocalCoordOpposite = myHex.faceCoord(oppositeInd);
      auto circLocation = 0.75 * faceLocalCoord + .25 *faceLocalCoordOpposite;
      auto globalCoord = myHex.mapToHex(faceLocalCoord);
      for (auto c: globalCoord) {
        f << c << ' ';
      }
      f << '\n';
      elementData += "2 " + std::to_string(centerCoordVtkInd) + " " + std::to_string(nodeIndex) + "\n";
      nodeIndex++;
      auto faceNormal = myHex.faceNormal(faceIndex);
      auto faceNormalOri = faceNormal * static_cast<double>(ori);
      auto jac = myHex.jacobian(circLocation);
      auto globalNormal0 = jac.inverse().transpose() * faceNormalOri;
//      auto globalNormal0 = jac * faceNormalOri;
      auto globalNormal = globalNormal0 / globalNormal0.norm();
      for (auto c: globalNormal) {
        vectorData += std::to_string(c) + " ";
      }
      vectorData += "\n";
    }
  }
  //Chunk of code to write edge coord and vector to file
//  Eigen::Vector3d edgeCoord {0, 0, 0};
//  Eigen::Vector3d edgeVecUnnormalized {edge.mNodeCoords[3]-edge.mNodeCoords[0],
//                           edge.mNodeCoords[4]-edge.mNodeCoords[1],
//                           edge.mNodeCoords[5]-edge.mNodeCoords[2]};
//  Eigen::Vector3d edgeVec = edgeVecUnnormalized / edgeVecUnnormalized.norm();
  elementData += "2";
  for(auto i = 0; i < 6; i+=3){
    f << edge.mNodeCoords[i] << ' ' << edge.mNodeCoords[i+1] << ' ' << edge.mNodeCoords[i+2] << '\n';
    elementData += " " + std::to_string(nodeIndex);
    ++nodeIndex;
  }
  elementData += '\n';

  edgeVec = edgeVec * edge.mOrientation * 10.;
  for(auto i = 0; i < 3; ++i){
    //f << edgeCoord[i] << " ";
    vectorData += std::to_string(edgeVec[i]) + " ";
  }

  auto numCells = 2 * myAdjHex.size() + 1;
  f << '\n';
  f << "CELLS " << numCells << ' ' << 3 * numCells << '\n';
  f << elementData;

  f << '\n' << "CELL_TYPES " << numCells << '\n';
  for(auto i = 0; i < numCells; ++i){
    f << "3\n";
  }

  f << '\n';
  f << "CELL_DATA " << numCells << "\n";
  f << "VECTORS faces double\n";
  f << vectorData;


//  f << '\n';
//  f << "POINT_DATA " << 2* myAdjHex.size()+1 << "\n";
//  f << "VECTORS faces double\n";
//  f << vectorData;

  f.close();
}

/*
 * Input: cube whose H and B vectors you want to plot, vector of edges so that you can get H and B values
 * Returns: A pair of physical space H and B vectors at center of the Hex [H,B]
 */
std::pair<Eigen::Vector3d, Eigen::Vector3d> physicalHandBAtCenter(UnitCube &cube,
                                                                  std::vector<UnitEdge> &edges){
  Eigen::Vector3d localCoord {.5, .5, .5};
  auto jac = cube.jacobian(localCoord);
  Eigen::Vector3d hVecCenter {0, 0, 0};
  Eigen::Vector3d bVecCenter {0, 0, 0};
  //loop over local edges and get their H and B values
//  size_t e = 11;
  for(auto e = 0; e < 12; ++e)
  {
    auto localEdgeCoord = cube.edgeCoord(e);
    auto globalEdgeIndex = cube.mGlobalEdges[e];
    const auto &myEdge = edges[globalEdgeIndex];
    Eigen::Vector3d hVec = cube.edgeTangential(e) * myEdge.mH;
    Eigen::Vector3d bVec = cube.edgeTangential(e) * myEdge.mB;
//    auto jac = cube.jacobian(localEdgeCoord);
    hVecCenter += hVec * .25;
    bVecCenter += bVec * .25;
  }
  //Convert from computational space to physical space
  Eigen::Vector3d hBarVecCenter = jac.inverse().transpose() * hVecCenter;
  Eigen::Vector3d bBarVecCenter = jac * bVecCenter / jac.determinant();
//  std::cout << "hex=" << cube.mIndex << '\n';
//  std::cout << jac << '\n';
  return std::make_pair(hBarVecCenter, bBarVecCenter);
}

std::pair<Eigen::Vector3d, Eigen::Vector3d> computationalHandBAtCenter(UnitCube &cube,
                                                                  std::vector<UnitEdge> &edges){
  Eigen::Vector3d localCoord {.5, .5, .5};
  auto jac = cube.jacobian(localCoord);
  Eigen::Vector3d hVecCenter {0, 0, 0};
  Eigen::Vector3d bVecCenter {0, 0, 0};
  //loop over local edges and get their H and B values
//  size_t e = 11;
  for(auto e = 0; e < 12; ++e)
  {
    auto localEdgeCoord = cube.edgeCoord(e);
    auto globalEdgeIndex = cube.mGlobalEdges[e];
    const auto &myEdge = edges[globalEdgeIndex];
    Eigen::Vector3d hVec = cube.edgeTangential(e) * myEdge.mH;
    Eigen::Vector3d bVec = cube.edgeTangential(e) * myEdge.mB;
//    auto jac = cube.jacobian(localEdgeCoord);
    hVecCenter += hVec * .25;
    bVecCenter += bVec * .25;
  }
  return std::make_pair(hVecCenter, bVecCenter);
}


/*
 * Input: cube whose E and D vectors you want to plot, vector of edges so that you can get E and D values
 * Returns: A pair of physical space H and B vectors at center of the Hex [E,D]
 */
std::pair<Eigen::Vector3d, Eigen::Vector3d> physicalEandDAtCenter(UnitCube &cube){
  Eigen::Vector3d eVec {0, 0, 0};
  Eigen::Vector3d dVec {0, 0, 0};
  //loop over local edges and get their H and B values
  for(auto f = 0; f < 6; ++f) {
    auto faceNormal = cube.faceNormal(f);
    eVec += faceNormal * 0.5 * cube.mE[f];
    dVec += faceNormal * 0.5 * cube.mD[f];
  }
  //Convert from computational space to physical space
  Eigen::Vector3d localCoord {.5, .5, .5};
  auto jac = cube.jacobian(localCoord);
  Eigen::Vector3d eBarVec = jac.inverse().transpose() * eVec;
  Eigen::Vector3d dBarVec = jac * dVec / jac.determinant();
  return std::make_pair(eBarVec, dBarVec);
}

void plotHexFields(std::vector<double> &nodeCoords, std::vector<UnitEdge> &edges,
                   std::vector<UnitCube> &cubes, std::string fName){
  //Chunk of code to plot the nodes and hexahedra
  std::ofstream fVtk;
  fVtk.open(fName);
  fVtk << "# vtk DataFile Version 3.0" << '\n';
  fVtk << "vtk output\n";
  fVtk << "ASCII" << '\n';
  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
  fVtk << "POINTS " << nodeCoords.size() / 3 << " double\n";
  //write node coords to file
  for(auto i = 0; i < nodeCoords.size(); i += 3){
    fVtk << nodeCoords[i] << " " << nodeCoords[i + 1] << " " << nodeCoords[i + 2] << "\n";
  }
  //write hexahedra to file
  const auto numCubes = cubes.size();
  fVtk << '\n';
  fVtk << "CELLS " << numCubes << ' ' << 9 * numCubes << '\n';
  for(auto h = 0; h < numCubes; ++h){
    const UnitCube &myCube = cubes[h];
    const auto myNodes = myCube.mNodes;
    fVtk << 8 ;//<< ' ' << myNodes[0] << ' ' << myNodes[1] << '\n';
    for(auto i = 0; i < 8; ++i){
      fVtk << ' ' << myNodes[i];
    }
    fVtk << '\n';
  }
  //define cell types as hexahedra
  fVtk << '\n' << "CELL_TYPES " << numCubes << '\n';
  for(auto i = 0; i < numCubes; ++i){
    fVtk << 12 << '\n';
  }
  //Chunk of code to write B fields to file
  fVtk << '\n';
  fVtk << "CELL_DATA " << numCubes << "\n";
  fVtk << "VECTORS D_Fields double\n";
  for(auto h = 0; h < numCubes; ++h) {
    UnitCube &myCube = cubes[h];
    auto [eVec, dVec] = physicalEandDAtCenter(myCube);
//    auto [hVec, bVec] = physicalHandBAtCenter(myCube, edges);
    fVtk << dVec.transpose() << '\n';
//    fVtk << eVec.transpose() << '\n';
//    fVtk << bVec.transpose() << '\n';
//    std::cout << bVec.transpose() << "\n\n";
  }
  fVtk.close();
}

void plotHexBHFields(std::vector<double> &nodeCoords, std::vector<UnitEdge> &edges,
                   std::vector<UnitCube> &cubes, std::string fName){
  //Chunk of code to plot the nodes and hexahedra
  std::ofstream fVtk;
  fVtk.open(fName);
  fVtk << "# vtk DataFile Version 3.0" << '\n';
  fVtk << "vtk output\n";
  fVtk << "ASCII" << '\n';
  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
  fVtk << "POINTS " << nodeCoords.size() / 3 << " double\n";
  //write node coords to file
  for(auto i = 0; i < nodeCoords.size(); i += 3){
    fVtk << nodeCoords[i] << " " << nodeCoords[i + 1] << " " << nodeCoords[i + 2] << "\n";
  }
  //write hexahedra to file
  const auto numCubes = cubes.size();
  fVtk << '\n';
  fVtk << "CELLS " << numCubes << ' ' << 9 * numCubes << '\n';
  for(auto h = 0; h < numCubes; ++h){
    const UnitCube &myCube = cubes[h];
    const auto myNodes = myCube.mNodes;
    fVtk << 8 ;//<< ' ' << myNodes[0] << ' ' << myNodes[1] << '\n';
    for(auto i = 0; i < 8; ++i){
      fVtk << ' ' << myNodes[i];
    }
    fVtk << '\n';
  }
  //define cell types as hexahedra
  fVtk << '\n' << "CELL_TYPES " << numCubes << '\n';
  for(auto i = 0; i < numCubes; ++i){
    fVtk << 12 << '\n';
  }
  //Chunk of code to write B fields to file
  fVtk << '\n';
  fVtk << "CELL_DATA " << numCubes << "\n";
  fVtk << "VECTORS B_Fields double\n";
  for(auto h = 0; h < numCubes; ++h) {
    UnitCube &myCube = cubes[h];
//    auto [eVec, dVec] = physicalEandDAtCenter(myCube);
    auto [hVec, bVec] = physicalHandBAtCenter(myCube, edges);
//    fVtk << hVec.transpose() << '\n';
    fVtk << hVec.transpose() << '\n';
//    std::cout << bVec.transpose() << "\n\n";
  }
  fVtk.close();
}

void plotEdgeFields(std::vector<double> &nodeCoords, std::vector<UnitEdge> &edges,
    std::vector<UnitCube> &cubes, std::string fName){
  std::ofstream fVtk;
  fVtk.open(fName);
  fVtk << "# vtk DataFile Version 3.0" << '\n';
  fVtk << "vtk output\n";
  fVtk << "ASCII" << '\n';
  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
  fVtk << "POINTS " << nodeCoords.size() / 3 << " double\n";
  for(auto i = 0; i < nodeCoords.size(); i += 3){
    fVtk << nodeCoords[i] << " " << nodeCoords[i + 1] << " " << nodeCoords[i + 2] << "\n";
  }

  fVtk << '\n';
  fVtk << "CELLS " << edges.size() << ' ' << 3 * edges.size() << '\n';
  for(auto e = 0; e < edges.size(); ++e){
    const UnitEdge myEdge = edges[e];
    const Edge myNodes = myEdge.mNodes;
    fVtk << 2 << ' ' << myNodes[0] << ' ' << myNodes[1] << '\n';
  }

  fVtk << '\n' << "CELL_TYPES " << edges.size() << '\n';
  for(auto i = 0; i < edges.size(); ++i){
    fVtk << 3 << '\n';
  }
  fVtk << "\n CELL_DATA " << edges.size() << '\n';
  fVtk << "VECTORS H_Field double\n";
  for(auto &edge : edges){
    auto hexIndex = edge.mAdjacentHex[0];
    auto localEdgeIndex = edge.mIndexInHex[0];
    auto &hex = cubes[hexIndex];
    auto localEdgeCoord = hex.edgeCoord(localEdgeIndex);
    auto localHVec = hex.edgeTangential(localEdgeIndex) * edge.mH;
    auto jac = hex.jacobian(localEdgeCoord);
    auto physicalHVec = jac.transpose().inverse() * localHVec; //gets to physical space, but still need to project
    auto hBar = edge.asUnitVector().transpose() * physicalHVec; //we have to project onto edge to get accurate H bar.
    fVtk << edge.asUnitVector().transpose() * hBar << '\n'; //plot physical vector
  }
  fVtk.close();
}

int main() {
  namespace mesh = gmsh::model::mesh;
  namespace geo = gmsh::model::geo;
  auto fileName = meshUnitCube(5);
//  auto fileName = meshUnitSquare(10, 1);
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
//  orientEdges(10,1,2,-1,edge2hex,hex2edge,edge2node,hex2node,globalOrientation,localOrientation);
//  orientEdges(20,1,2,-1,edge2hex,hex2edge,edge2node,hex2node,globalOrientation,localOrientation);
  for(auto i =0; i< numHex; ++i){
    std::array<size_t,12> myEdges = hex2edge[i];
    for(auto j = 0; j < 12; ++j){
      size_t myGlobalEdge = myEdges[j];
      orientEdges(myGlobalEdge, j, i, 1, edge2hex, hex2edge, edge2node, hex2node, globalOrientation,
                  localOrientation);
    }
  }
//
//  for(auto i = 0; i < numEdges; ++i){
//    const auto myNode = edge2node[i];
//    if(globalOrientation[i] != 0){
//      //std::cout << "edge=" << "(" << myNode[0]-1 << "," << myNode[1]-1 << ") Ori=" << globalOrientation[i] << '\n';
//    }
//  }

////Plot Orientation of Edges
//  std::vector<std::array<double, 3>> dirs;
  std::ofstream fVtk;
//  fVtk.open("orient.vtk");
//  fVtk << "# vtk DataFile Version 3.0" << '\n';
//  fVtk << "vtk output\n";
//  fVtk << "ASCII" << '\n';
//  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
//  fVtk << "POINTS " << numEdges << " double\n";
//  for(const auto e : edge2node){
//    std::array<double, 3> c = {0, 0, 0};
//    std::array<double, 3> dir = {0, 0, 0};
//    for(auto i = 0; i < 2; ++i){
//      size_t n = e[i];
//      std::vector<double> coord;
//      std::vector<double> paraCoord;
//      int d, t;
//      gmsh::model::mesh::getNode(n, coord, paraCoord, d, t);
//
//      for(auto j = 0; j < 3; ++j){
//        c[j] += 0.5 * coord[j];
//        dir[j] = (i==1) ? dir[j] + coord[j] : dir[j] - coord[j];
//      }
//    }
//    dirs.push_back(dir);
//    fVtk << c[0] << " " << c[1] << " " << c[2] << '\n';
//  }
////  fVtk << "CELLS 0 0\n";
////  fVtk << "CELL_TYPES 0\n";
//  fVtk << "POINT_DATA " << numEdges << "\n";
//  fVtk << "VECTORS edges double\n";
//  for(auto k =0; k < numEdges; ++k){
//    auto dir = dirs[k];
//    auto go = globalOrientation[k];
//        fVtk << dir[0] * go << " " << dir[1] * go << " " << dir[2] * go << '\n';
//  }
//  fVtk.close();

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
        isOrigin = isOrigin && (myTrueDir == 1); //if myTrueDir != 1 this chain is false, it is not the origin
      }
      if(isOrigin){
        originIndex = i;
        localOriginIndex[hexIndex] = originIndex;
        //localOriginIndex[hexIndex] = 0;
        //break;
      }
    }
  }

////plot origins of each hex
////  std::ofstream fVtk;
//  fVtk.open("Origins.vtk");
//  fVtk << "# vtk DataFile Version 3.0" << '\n';
//  fVtk << "vtk output\n";
//  fVtk << "ASCII" << '\n';
//  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
//  fVtk << "POINTS " << numHex << " double\n";
//  for(auto hexIndex = 0; hexIndex < nu5mHex; ++hexIndex){
//    auto myHexGlobalNodes = hex2node[hexIndex];
//    auto myLocalIndex = localOriginIndex[hexIndex];
//    auto myGlobalNode = myHexGlobalNodes[myLocalIndex];
//    std::vector<double> coord;
//    std::vector<double> paraCoord;
//    int dim, tag;
//    gmsh::model::mesh::getNode(myGlobalNode, coord, paraCoord, dim, tag);
//    fVtk << coord[0] << " " << coord[1] << " " << coord[2] << '\n';
//  }
//  fVtk.close();

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

  //chunk of code to repopulate hex2edge or create hex2edge oriented
  for(auto h = 0; h < numHex; ++h){
    auto myEdges = &hex2edge[h];
    std::array<size_t, 12> hex2edgeOriented;
    //loop over edges, create sorted global nodes array for each edge,
    // use the edge map to get the global index
    for(auto e = 0; e < 12; ++e){
      auto localNodes = edgeIndexToNodeIndex(e);
      std::array<size_t, 2> globalNodes {hex2nodeOriented[h][localNodes[0]],
                                         hex2nodeOriented[h][localNodes[1]]};
      std::sort(globalNodes.begin(), globalNodes.end());
      const auto it = edgeMap.find(globalNodes);
      if(it == edgeMap.end()){
        std::cout << "Edge=(" << globalNodes[0] << "," << globalNodes[1]
        << ") not found in hex2edge orientation fix \n";
      }
      else{
        const auto &[globalEdgeNodes, edgeIndex] = *it;
        hex2edgeOriented[e] = edgeIndex;
      }
    }
    //update the hex2node array
    hex2edge[h] = hex2edgeOriented;
  }

  //Comparing adjustments of hex to node with visit visualization.
//  for(auto h = 0; h < numHex; ++h){
//    auto myNodes = hex2nodeOriented[h];
//    auto myNodesOriginal = hex2node[h];
//    int myOrigin = hex2node[h][localOriginIndex[h]];
//    std::cout << "hexOrigin=" << myOrigin -1 <<   std::endl;
//    for(int i = 0; i < 8; ++i){
//      int ind = myNodes[i];
//      std::cout << ind - 1 << ",";
//    }
//    std::cout << "\n";
//    for(int i = 0; i < 8; ++i){
//      int ind = myNodesOriginal[i];
//      std::cout << ind - 1 << ",";
//    }
//    std::cout << "\n";
//  }
  //return new hex2node, nodes, node2coord
  std::vector<size_t> nodeTags;
  std::vector<double> nodeCoords;
  std::vector<double> parametricCoords;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords, -1, -1, false, false);
//  for(auto nd : nodeTags){
//    std::cout << nd << ',';
//  }

  //create hex2face and face2hex
  std::vector<std::array<size_t, 6>> hex2face(numHex, {0, 0, 0, 0, 0, 0});
//  std::vector<std::vector<size_t>> face2hex;
//  std::vector<std::vector<size_t>> faceIndexInHex;
  std::vector<std::vector<std::pair<size_t, size_t>>> adjHexAndIndexInHex;
  std::vector<std::array<size_t,4>> quad2node;

  parseFaces(quad2node, hex2nodeOriented, hex2face, adjHexAndIndexInHex);
  std::vector<UnitFace> faces;
  auto numFaces = adjHexAndIndexInHex.size();
  for(auto f = 0; f < numFaces; ++f){
//    for(auto hx : face2hex[f]) std::cout <<"f=" << f << " h=" << hx <<"\n";
    faces.push_back(UnitFace(f, adjHexAndIndexInHex[f]));
  }

  //Populate UnitCube Objects
  std::vector<UnitCube> cubes;
  for(auto h = 0; h < numHex; ++h){
    const Hex myNodes = hex2nodeOriented[h];
    Hex myNodes0Index = hex2nodeOriented[h];
    std::vector<double> myNodeCoords {};
    for(auto i = 0; i < 8; ++ i){
      auto myNode = myNodes[i];
      auto myNodeIndex = myNode-1; //nodes are stored with 1-indexing; so we have to -1 to get vector index
      myNodes0Index[i] = myNodeIndex;
      for(auto c = 0; c < 3; ++c){
        auto nodeCoordIndex = myNodeIndex * 3 + c;//get exact spot in node coord vector
        myNodeCoords.push_back(nodeCoords[nodeCoordIndex]);
      }
    }
    //cubes.push_back(UnitCube(myNodes, myNodeCoords, hex2edge[h]));
    Eigen::Matrix3d id = Eigen::Matrix3d::Identity();
//    cubes.push_back(UnitCube(h, myNodes0Index, myNodeCoords, hex2edge[h], id, id));
    cubes.push_back(UnitCube(h, myNodes0Index, myNodeCoords, hex2edge[h], hex2face[h], id, id));
  }

  //get edge index in adjacent hex for each edge
  std::vector<std::vector<size_t>> edge2hexLocalIndex;
  for(auto e = 0; e < numEdges; ++e){
    auto myNodes = edge2node[e];
    auto myAdjHex = edge2hex[e];
    std::vector<size_t> myAdjHexIndex;
//    std::cout << "edge=" << e << "\n ";
    for(auto hex : myAdjHex){
      //loop over hex edges
//      std::cout << " (hex=" << hex << " ";
      for(auto i = 0; i < 12; ++i){
        auto localNodes = edgeIndexToNodeIndex(i);
        std::array<size_t, 2> globalNodes {hex2nodeOriented[hex][localNodes[0]],
                                           hex2nodeOriented[hex][localNodes[1]]};
        std::sort(globalNodes.begin(), globalNodes.end()); //sorted now
        if(myNodes == globalNodes){
          myAdjHexIndex.push_back(i);
//          std::cout << "index in=" << i << ")";
          break;
        }
      }
    }
    edge2hexLocalIndex.push_back(myAdjHexIndex);
//    std::cout << "\n";
  }

  //Populate UnitEdge objects
  std::vector<UnitEdge> edges;
  for(auto e = 0; e < numEdges; ++e){
    auto myNodes = edge2node[e];
    auto myOri = globalOrientation[e];
    auto myAdjHex = edge2hex[e];
    auto myLocalIndexInHex = edge2hexLocalIndex[e];
    std::vector<double> myNodeCoords;
    Edge myNodes0Index = myNodes; //nodes of edge with 0-index convention
    for(auto i = 0; i < 2; ++ i){
      auto myNode = myNodes[i];
      myNodes0Index[i] = myNode - 1; //store as 0-indexed
      auto myNodeIndex = myNode-1; //nodes are stored with 1-indexing; so we have to -1 to get vector index
      for(auto c = 0; c < 3; ++c){
        auto nodeCoordIndex = myNodeIndex * 3 + c;//get exact spot in node coord vector
        myNodeCoords.push_back(nodeCoords[nodeCoordIndex]);
      }
    }
//    edges.push_back(UnitEdge(myNodes, myNodeCoords, myAdjHex, myLocalIndexInHex, myOri));
    std::cout << "Edge vector creation edge index = " << e << '\n';
    edges.push_back(UnitEdge(e, myNodes0Index, myNodeCoords, myAdjHex, myLocalIndexInHex, myOri));
  }

  //Chunk to set effective mu for each edge
  for(auto &unitEdge : edges){
    auto adjHexs = unitEdge.mAdjacentHex;
    auto numAdjHexs = adjHexs.size();
    double edgeMu = 0;
    for(auto h = 0; h < adjHexs.size(); ++h){
      auto adjHexIndex = adjHexs[h];
      auto &adjHex = cubes[adjHexIndex];
      auto localEdgeIndex = unitEdge.mIndexInHex[h];
      auto localEdgeTangent = adjHex.edgeTangential(localEdgeIndex);
      double localMuTangent = localEdgeTangent.transpose() * adjHex.getMuComp(localEdgeIndex) * localEdgeTangent;
      edgeMu += localMuTangent / numAdjHexs;
    }
    unitEdge.mMu = edgeMu;
  }

//Plot Edges: Consider refactoring into own function.
  fVtk.open("Edges.vtk");
  fVtk << "# vtk DataFile Version 3.0" << '\n';
  fVtk << "vtk output\n";
  fVtk << "ASCII" << '\n';
  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
  fVtk << "POINTS " << nodeCoords.size() / 3 << " double\n";
  for(auto i = 0; i < nodeCoords.size(); i += 3){
    fVtk << nodeCoords[i] << " " << nodeCoords[i + 1] << " " << nodeCoords[i + 2] << "\n";
  }
  fVtk << '\n';
  fVtk << "CELLS " << edges.size() << ' ' << 3 * edges.size() << '\n';
  for(auto e = 0; e < edges.size(); ++e){
    const UnitEdge myEdge = edges[e];
    const Edge myNodes = myEdge.mNodes;
    fVtk << 2 << ' ' << myNodes[0] << ' ' << myNodes[1] << '\n';
  }

  fVtk << '\n' << "CELL_TYPES " << edges.size() << '\n';
  for(auto i = 0; i < edges.size(); ++i){
    fVtk << 3 << '\n';
  }
  fVtk.close();

////Plot hexahedra
//  fVtk.open("Hexahedra.vtk");
//  fVtk << "# vtk DataFile Version 3.0" << '\n';
//  fVtk << "vtk output\n";
//  fVtk << "ASCII" << '\n';
//  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
//  fVtk << "POINTS " << nodeCoords.size() / 3 << " double\n";
//  for(auto i = 0; i < nodeCoords.size(); i += 3){
//    fVtk << nodeCoords[i] << " " << nodeCoords[i + 1] << " " << nodeCoords[i + 2] << "\n";
//  }
//
//  const auto numCubes = cubes.size();
//  fVtk << '\n';
//  fVtk << "CELLS " << numCubes << ' ' << 9 * numCubes << '\n';
//  for(auto h = 0; h < numCubes; ++h){
//    const UnitCube myCube = cubes[h];
//    const auto myNodes = myCube.mNodes;
//    fVtk << 8 ;//<< ' ' << myNodes[0] << ' ' << myNodes[1] << '\n';
//    for(auto i = 0; i < 8; ++i){
//      fVtk << ' ' << myNodes[i];
//    }
//    fVtk << '\n';
//  }
//
//  fVtk << '\n' << "CELL_TYPES " << numCubes << '\n';
//  for(auto i = 0; i < numCubes; ++i){
//    fVtk << 12 << '\n';
//  }
//  fVtk.close();

//  for(auto ed : edges){
//    std::cout << "edge=" << ed.mIndex << " Adj hex=(";
//    for(auto h : ed.mAdjacentHex){
//      std::cout << h << ',';
//    }
//    std::cout << ") Index in Hex =(";
//    for(auto h : ed.mIndexInHex){
//      std::cout << h << ',';
//    }
//    std::cout << ")\n";
//  }

//  for(auto ed : edges){
//    std::string fName = "edge" + std::to_string(ed.mIndex) + ".vtk";
//    plotCirculation2(ed, cubes, fName);
//  }
//  size_t plotEdgeIndex = 44;
//  std::string fNamePlotEdgeIndex = "edge" + std::to_string(plotEdgeIndex) + ".vtk";
//  plotCirculation2(edges[plotEdgeIndex], cubes, fNamePlotEdgeIndex);

  //Loop of UnitCubes to flag faces that are boundaries
  //A face in a hex h will be a boundary if its edges do NOT share any other adjacent hexahedra besides h.
  for(auto &myCube : cubes) {
//    std::cout << "hex=" << myCube.mIndex << "\n";
    for (auto f = 0; f < 6; ++f) {
      auto localEdges = myCube.edgesAdjacentToFace(f);
      const size_t globalEdgeIndex0 = myCube.mGlobalEdges[localEdges[0]];
      std::vector<size_t> adjHexIntersection = edges[globalEdgeIndex0].mAdjacentHex;
      //find intersection of adjacent hexes
      for(auto e : localEdges){
        std::vector<size_t> tempIntersection; //holds the output of std::set_intersection
        const size_t globalEdgeIndex = myCube.mGlobalEdges[e];
        std::vector<size_t> myAdjHex = edges[globalEdgeIndex].mAdjacentHex;
        std::sort(adjHexIntersection.begin(), adjHexIntersection.end());
        std::sort(myAdjHex.begin(), myAdjHex.end());
        std::set_intersection(myAdjHex.begin(), myAdjHex.end(),
                              adjHexIntersection.begin(), adjHexIntersection.end(),
                              std::back_inserter(tempIntersection));
        adjHexIntersection = tempIntersection; //overwrite so that next edge has something to compare to
      }
      //if the face only has 1 adjacent hexahedron, then it is a boundary
      if(adjHexIntersection.size() == 1){
        //now mark the face and its edges as a boundary
        myCube.mFaceBoundaryFlag[f] = 1;
        for(auto e : localEdges){
           size_t globalEdgeIndex = myCube.mGlobalEdges[e];
           edges[globalEdgeIndex].mBoundaryFlag = 1;
        }
      }
    }
  }

  //Simple computation of dt, using min edge length
  auto itEdgeWithMinLength = std::min_element(edges.begin(), edges.end(), [](UnitEdge &a, UnitEdge &b){
    return a.length() <b.length();
  });
  auto minEdgeLength = itEdgeWithMinLength->length();
  const double dt = 0.75 * minEdgeLength / sqrt(3);
//  const double dt = .1;
  std::cout << "dt=" << dt <<std::endl;

  //Exact solution TE_
  auto exactSolH0Zero = [](Eigen::Vector3d x, double t){
    double eps = 1;
    int m = 1;
    int n = 1;
    int p = 1;
    double kx = M_PI * m;
    double ky = M_PI * n;
    double kz = M_PI * p;
    double k = sqrt(kx*kx + ky*ky + kz*kz);
    double k2 = k * k;
    double w = k / sqrt(eps);
    Eigen::Vector3d HBar {-ky*kx * cos(kx*x[0]) * sin(ky*x[1]) * sin(kz*x[2]),
                          -ky*kz * sin(kx*x[0]) * cos(ky*x[1]) * sin(kz*x[2]),
                          (k2-kz*kz) * sin(kx*x[0]) * sin(ky*x[1]) * cos(kz*x[2])};
//    HBar = HBar * sin(w*t) / (w * eps);
    HBar = -HBar * sin(w*t) / (w * eps);
    Eigen::Vector3d BBar = HBar;
    Eigen::Vector3d EBar {ky / eps * sin(kx*x[0]) * cos(ky*x[1]) * cos(kz*x[2]),
                          -kx / eps * cos(kx*x[0]) * sin(ky*x[1]) * cos(kz*x[2]),
                          0};
//    EBar = EBar * (-cos(w*t));
    EBar = EBar * cos(w*t);
    Eigen::Vector3d DBar = EBar * eps;
    return std::array<Eigen::Vector3d, 4> {HBar, BBar, EBar, DBar};
  };
//  auto exactSolH0Zero = [](Eigen::Vector3d x, double t){
//    double eps = 1;
//    int m = 1;
//    int n = 1;
//    int p = 1;
//    double kx = M_PI * m;
//    double ky = M_PI * n;
//    double kz = M_PI * p;
//    double k = sqrt(kx*kx + ky*ky + kz*kz);
//    double k2 = k * k;
//    double w = k / sqrt(eps);
//    Eigen::Vector3d HBar {ky / eps * cos(kx*x[0]) * sin(ky*x[1]) * sin(kz*x[2]),
//                          -kx / eps * sin(kx*x[0]) * cos(ky*x[1]) * sin(kz*x[2]),
//                          0};
//    HBar = HBar * sin(w*t);
//    Eigen::Vector3d BBar = HBar;
//    Eigen::Vector3d EBar {-ky*kx * sin(kx*x[0]) * cos(ky*x[1]) * cos(kz*x[2]),
//                          -ky*kz * cos(kx*x[0]) * sin(ky*x[1]) * cos(kz*x[2]),
//                          (k2-kz*kz) * cos(kx*x[0]) * cos(ky*x[1]) * sin(kz*x[2])};
//    EBar = EBar * cos(w*t) / (w * eps);
//    Eigen::Vector3d DBar = EBar;
//    return std::array<Eigen::Vector3d, 4> {HBar, BBar, EBar, DBar};
//  };

  //Set initial fields
  //Exact solution with D0=0
  auto exactSolD0Zero = [](Eigen::Vector3d x, double t){
    double eps = 1;
    int m = 1;
    int n = 1;
    int p = 1;
    double kx = M_PI * m;
    double ky = M_PI * n;
    double kz = M_PI * p;
    double k = sqrt(kx*kx + ky*ky + kz*kz);
    double k2 = k * k;
    double w = k / sqrt(eps);
    Eigen::Vector3d HBar {-ky*kx * cos(kx*x[0]) * sin(ky*x[1]) * sin(kz*x[2]),
                          -ky*kz * sin(kx*x[0]) * cos(ky*x[1]) * sin(kz*x[2]),
                          (k2-kz*kz) * sin(kx*x[0]) * sin(ky*x[1]) * cos(kz*x[2])};
    HBar = HBar * cos(w*t) / (w * eps);
    Eigen::Vector3d BBar = HBar;
    Eigen::Vector3d EBar {ky / eps * sin(kx*x[0]) * cos(ky*x[1]) * cos(kz*x[2]) * sin(w*t),
                          -kx / eps * cos(kx*x[0]) * sin(ky*x[1]) * cos(kz*x[2]) * sin(w*t),
                          0};
    Eigen::Vector3d DBar = EBar * eps;
    return std::array<Eigen::Vector3d, 4> {HBar, BBar, EBar, DBar};
  };

  //set initial fields
  for(auto &myCube : cubes){
    for(auto f = 0; f < 6; ++f){
      auto faceCoord = myCube.facePhysicalCoordinate(f);
      auto [unusedH, unusedB, unusedE, dBarVec] = exactSolH0Zero(faceCoord, -0.5 * dt);
      Eigen::Vector3d areaWeightedNormal = myCube.physicalFaceAreaWeighted(f);
      myCube.mD[f] = dBarVec.transpose() * areaWeightedNormal;
    }
    for(int f = 0; f < 6; ++f) {
      if (myCube.mFaceBoundaryFlag[f] == 0) {
        auto dVec = myCube.getDVectorAtCirculationLocation(f);
        auto epsInv = myCube.getEpsCompInv(f);
        myCube.mE[f] = myCube.faceNormal(f).transpose() * epsInv * dVec;
      }
    }
  }

  plotHexFields(nodeCoords, edges, cubes, "DFields0.vtk");
  plotHexBHFields(nodeCoords, edges, cubes, "HFields0.vtk");
  double omega = M_PI * sqrt(3);
  double tMax = 2 * M_PI / omega;
//  const int numSteps = static_cast<int>(tMax / dt);
  const int numSteps = 2;
  double t = 0;

  //FDTD Update Loop
  for(int timeStep = 0; timeStep < numSteps; ++timeStep) {
    //D-Update
    for(auto &myCube : cubes){
      //update D on each face
      for(auto f = 0; f < 6; ++f){
        if(myCube.mFaceBoundaryFlag[f]==0) {
          auto localAdjEdges = myCube.edgesAdjacentToFace(f);
          auto localAdjEdgeOri = myCube.edgesOrientationToFace(f);
          double hCirculation = 0;
          for (auto e = 0; e < localAdjEdges.size(); ++e) {
            size_t localEdge = localAdjEdges[e];
//            hCirculation += myCube.mH[localEdge] * localAdjEdgeOri[e];
            size_t globalEdgeIndex = myCube.mGlobalEdges[localEdge];
            hCirculation += edges[globalEdgeIndex].mH * localAdjEdgeOri[e];
          }
          myCube.mD[f] += dt * hCirculation;
        }
      }
      //E-update
      for(int f = 0; f < 6; ++f) {
        if (myCube.mFaceBoundaryFlag[f] == 0) {
          auto dVec = myCube.getDVectorAtCirculationLocation(f);
          auto epsInv = myCube.getEpsCompInv(f);
//          myCube.mE[f] = myCube.faceNormal(f).transpose() * myCube.mEpsCompInv * dVec; //Working
          myCube.mE[f] = myCube.faceNormal(f).transpose() * epsInv * dVec;
        }
      }
    }
    t += 0.5 * dt;

    //B-Update
    for (auto &myEdge: edges) {
      if(myEdge.mBoundaryFlag == 0) {
        const auto &adjHex = myEdge.mAdjacentHex;
        const auto &edgeIndicesInHex = myEdge.mIndexInHex;
        const size_t numAdjHex = adjHex.size();
        double circulation = 0;
        //loop over adjacent hexahedra and calculate partial circulations from each to update B on edge
        for (auto h = 0; h < numAdjHex; ++h) {
          const size_t myHexIndex = adjHex[h];
          const size_t localEdgeIndex = edgeIndicesInHex[h];
          UnitCube &myHex = cubes[myHexIndex];
          auto localEdgeDir = myHex.edgeTangential(localEdgeIndex);
          auto facesAdjToEdge = myHex.facesAdjacentToEdge(localEdgeIndex);
          auto faceOrientations = myHex.facesOrientationToEdge(localEdgeIndex);
          for (auto f = 0; f < facesAdjToEdge.size(); ++f) {
            size_t faceIndex = facesAdjToEdge[f];
            int faceOri = faceOrientations[f];
            //E is stored at circulation points, so there is no need to interpolate.
            double e = myHex.mE[faceIndex];
            //unit cube as length 0.5, Dl and e are in same dir on comp domain
            double eDotDl = e * 0.5 * faceOri;
            circulation += eDotDl;
          }
        }
        myEdge.mBFlux += -dt * 4.*circulation / numAdjHex;
//        myEdge.mBFlux += -dt *circulation;
      }
    }
//        H-Update-Working
    for (auto &myEdge: edges) {
      if(myEdge.mBoundaryFlag == 0) {
        const auto &adjHex = myEdge.mAdjacentHex;
        const auto &edgeIndicesInHex = myEdge.mIndexInHex;
        const size_t numAdjHex = adjHex.size();
        double globalEdgeScalarH = 0;
        double totalDualArea = 0;
        double totalJac = 0;
        //In each adj hex, construct a vector for B so that we can compute H = mu^-1 B
        for (auto h = 0; h < numAdjHex; ++h) {
          const size_t myHexIndex = adjHex[h];
          const size_t localEdgeIndex = edgeIndicesInHex[h];
          UnitCube &myHex = cubes[myHexIndex];
          auto localEdgeDir = myHex.edgeTangential(localEdgeIndex);
          //use local adj. edges to get dB/dt Flux as a vector
          Eigen::Vector3d localBFluxVec = localEdgeDir * myEdge.mBFlux * 9. / 16.;
          auto parEdges = localParallelEdgeCounterClockwise(localEdgeIndex);
//          Eigen::Vector3d localBFluxVec = localEdgeDir * myEdge.mBFlux;
          for(auto adjEdgeIndex = 0; adjEdgeIndex < 12; ++adjEdgeIndex){
            size_t globalAdjEdgeIndex = myHex.mGlobalEdges[adjEdgeIndex];
            Eigen::Vector3d adjEdgeDir = myHex.edgeTangential(adjEdgeIndex);
            //if the other edge is not in the same direction, then add it;
            if(adjEdgeDir != localEdgeDir){
              Eigen::Vector3d adjEdgeBFlux = myHex.edgeTangential(adjEdgeIndex) * edges[globalAdjEdgeIndex].mBFlux;
              localBFluxVec += adjEdgeBFlux * 0.25; //interpolate all other components to unit cube center
            }
            else if(adjEdgeIndex == parEdges[1] or adjEdgeIndex == parEdges[3]){
              Eigen::Vector3d adjEdgeBFlux = myHex.edgeTangential(adjEdgeIndex) * edges[globalAdjEdgeIndex].mBFlux;
              localBFluxVec += adjEdgeBFlux * 3. / 16.;
            }
            else if(adjEdgeIndex == parEdges[2]){
              Eigen::Vector3d adjEdgeBFlux = myHex.edgeTangential(adjEdgeIndex) * edges[globalAdjEdgeIndex].mBFlux;
              localBFluxVec += adjEdgeBFlux * 1. / 16.;
            }
          }
//          auto adjEdges = edgesAdjacentToEdge(localEdgeIndex);
//          for (auto adjEdgeIndex: adjEdges) {
//            size_t globalAdjEdgeIndex = myHex.mGlobalEdges[adjEdgeIndex];
//            Eigen::Vector3d adjEdgeDir = myHex.edgeTangential(adjEdgeIndex);
//            localBFluxVec += 0.5 * adjEdgeDir * edges[globalAdjEdgeIndex].mBFlux;
//          }
//          Eigen::Vector3d localHVec = myHex.getMuCompInv(localEdgeIndex) * localBFluxVec;
          Eigen::Vector3d localHVec = myHex.getMuCompInvAtCenter() * localBFluxVec;
//          Eigen::Vector3d localHVec = myHex.getMuCompInvAtDualFace(localEdgeIndex) * localBFluxVec;
          //now project H onto edge to get scalar, also we average h scalar about all adj hexes
          double localHScalar = (localEdgeDir.transpose() * localHVec);
//          globalEdgeScalarH += localHScalar * myHex.dualFacePhysicalArea(localEdgeIndex);
          globalEdgeScalarH += localHScalar;
//          globalEdgeScalarH += localHScalar * myHex.dualFaceJacobian(localEdgeIndex).determinant();
          totalJac += myHex.dualFaceJacobian(localEdgeIndex).determinant();
          totalDualArea += myHex.dualFacePhysicalArea(localEdgeIndex);
//          myHex.mH[localEdgeIndex] += dt * localHScalar;
        }
        myEdge.mH = globalEdgeScalarH / numAdjHex;
//        myEdge.mH += dt * globalEdgeScalarH / totalDualArea;
//        myEdge.mH = globalEdgeScalarH / totalDualArea;
//        myEdge.mH += dt * globalEdgeScalarH / totalJac;
      }
      else{
        myEdge.mH = 0;
      }
    }
    t += 0.5 * dt;
  }
  plotHexFields(nodeCoords, edges, cubes, "DFields1.vtk");
  plotHexBHFields(nodeCoords, edges, cubes, "HFields1.vtk");
//  plotEdgeFields(nodeCoords, edges, cubes, "HFields1.vtk");
  double runningMaxError = 0.;
  double runningMaxExact = 0.;
  double runningL2Error = 0.;
  double runningL2Exact = 0.;
  size_t maxErrorEdge = 0;
  for(auto &myEdge : edges){
    size_t adjHex = myEdge.mAdjacentHex[0];
    size_t e = myEdge.mIndexInHex[0];
    UnitCube &myCube = cubes[adjHex];
    if(myEdge.mBoundaryFlag == 0) {
      auto physicalEdgeCoord = myCube.edgePhysicalCoord(e);
      auto [hBarVec, unusedB, unusedE, unusedD] = exactSolH0Zero(physicalEdgeCoord, t);
      Eigen::Vector3d physicalEdge = myCube.edgeJacobian(e) * myCube.edgeTangential(e);
      double hExact = hBarVec.transpose() * physicalEdge; //this is in computational space
      double diff = abs(hExact - myEdge.mH);
      if(diff > runningMaxError){
        maxErrorEdge = myEdge.mIndex;
      }
      runningMaxError = std::max(runningMaxError, diff);
      runningMaxExact = std::max(runningMaxExact, hExact);
      runningL2Error += diff * diff;
      runningL2Exact += hExact * hExact;
    }
  }

  double runningMaxErrorD = 0.;
  double runningMaxExactD = 0.;
  double runningL2ErrorD = 0.;
  double runningL2ExactD = 0.;
  for(auto &myCube : cubes){
    for(auto f = 0; f<6; ++f){
      if(myCube.mFaceBoundaryFlag[f]==0) {
        auto faceCoord = myCube.facePhysicalCoordinate(f);
        auto areaWeightedNormal = myCube.physicalFaceAreaWeighted(f);
        auto [unusedH, unusedB, unusedE, dBarVec] = exactSolH0Zero(faceCoord, t - 0.5 * dt);
        double exactD = dBarVec.transpose() * areaWeightedNormal;
        double diff = abs(myCube.mD[f] - exactD);
        runningMaxErrorD = std::max(runningMaxErrorD, diff);
        runningMaxExactD = std::max(runningMaxExactD, exactD);
        runningL2ErrorD += diff * diff;
        runningL2ExactD += exactD * exactD;
      }
    }
  }

  std::cout << "H, D: rel lInf error = " << runningMaxError / runningMaxExact << ","
            << runningMaxErrorD / runningMaxExactD << '\n';
  std::cout << "H, D: rel l2 error = " << sqrt(runningL2Error / runningL2Exact)
            << "," << sqrt(runningL2ErrorD / runningL2ExactD)  << std::endl;
  std::cout << "t=" << t << ", t-max=" << tMax;
  std::cout << ", num steps = " << numSteps << std::endl;
  std::cout << "max error edge " << maxErrorEdge << '\n';

////set H-fields to error and plot
//  for(auto &edge : edges){
//    auto physicalEdgeCoord = edge.physicalCoordinate();
////    auto [hBarVec, unusedB, unusedE, unusedD] = exactSolD0Zero(physicalEdgeCoord, t-0.5*dt);
//    auto [hBarVec, unusedB, unusedE, unusedD] = exactSolH0Zero(physicalEdgeCoord, t);
//    //Piola transformation rule or H^T dl = HBar^T dl
//    double hExact = hBarVec.transpose() * edge.asUnitVector();
//    auto diff = abs(hExact*edge.length()-edge.mH);
//    edge.mH = diff;
//  }
//  plotHexBHFields(nodeCoords, edges, cubes, "HError.vtk");


//  //plot boundary faces
//  std::stringstream myBoundFaces;
//  size_t numBoundaryFaces = 0;
//  for(auto &myCube : cubes){
//    for(auto f = 0; f < 6; ++f){
//      if(myCube.mFaceBoundaryFlag[f] == 1){
//        auto myCoord = myCube.faceCoord(f);
//        auto globalCoord = myCube.mapToHex(myCoord);
//        myBoundFaces << globalCoord.transpose() << '\n';
//        numBoundaryFaces++;
//      }
//    }
//  }
//  std::stringstream myBoundEdges;
//  size_t numBoundEdges = 0;
//  for(UnitEdge &myE : edges){
//    if(myE.mBoundaryFlag ==  1){
//      myBoundEdges << myE.physicalCoordinate().transpose() << '\n';
//      numBoundEdges++;
//    }
//  }
//  fVtk.open("BoundaryFaces.vtk");
//  fVtk << "# vtk DataFile Version 3.0" << '\n';
//  fVtk << "vtk output\n";
//  fVtk << "ASCII" << '\n';
//  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
//  fVtk << "POINTS " << numBoundaryFaces << " double\n\n";
//  fVtk << myBoundFaces.str();
//  fVtk.close();
//
//  fVtk.open("BoundaryEdges.vtk");
//  fVtk << "# vtk DataFile Version 3.0" << '\n';
//  fVtk << "vtk output\n";
//  fVtk << "ASCII" << '\n';
//  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
//  fVtk << "POINTS " << numBoundEdges << " double\n\n";
//  fVtk << myBoundEdges.str();
//  fVtk.close();

////Plot max error edge
  fVtk.open("MaxErrorEdge.vtk");
  fVtk << "# vtk DataFile Version 3.0" << '\n';
  fVtk << "vtk output\n";
  fVtk << "ASCII" << '\n';
  fVtk << "DATASET UNSTRUCTURED_GRID" << '\n';
  fVtk << "POINTS " << 1 << " double\n\n";
  fVtk << edges[maxErrorEdge].physicalCoordinate().transpose() << std::endl;
  fVtk.close();
}
