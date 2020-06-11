#ifndef MESH_H
#define MESH_H

#include <QObject>


#include <common/interfaces.h>
#include <iostream>
#include <stdio.h>
#include <QtScript>
#include <vector>
#include "utilities.h"

#include <QDebug>
#include "types.h"
#include "Vertice.h"
#include "Mesher.h"


using namespace std ;
using namespace vcg ;

class Mesh{
public:
    Mesh();
    void printSignal();
    float isLeft(vcg::Point3f P0, vcg::Point3f P1,vcg::Point3f P2);
    void swap(vcg::Point3f *a, vcg::Point3f *b);
    std::vector<Point3f> sortIncreasing(std::vector<Point3f> pointVector);
    std::vector<Point3f> sortDecreasing(std::vector<Point3f> pointVector);
    std::vector<int> PointMIN_MIN_MIN_MAX(std::vector<Point3f> pointVector);
    std::vector<int> PointMAX_MIN_MAX_MAX(std::vector<Point3f> pointVector);
    std::vector<Point3f> LowerHull(std::vector<Point3f> pointVector);
    std::vector<Point3f> UpperHull(std::vector<Point3f> pointVector);
    std::vector<Point3f> ConvexHull(std::vector<Point3f> pointVector);
    void DT(MeshDocument &md);

};

#endif // MESH_H
