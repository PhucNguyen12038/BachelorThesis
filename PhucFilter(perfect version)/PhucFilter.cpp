/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "PhucFilter.h"

#include <QtScript>
#include <vector>







#define MAX 20

// When I input set of points in 2D, PHUC_FUNCTION will generate the
// convex hull of these points

// Constructor usually performs only two simple tasks of filling the two lists 
//  - typeList: with all the possible id of the filtering actions
//  - actionList with the corresponding actions. If you want to add icons to your filtering actions you can do here by construction the QActions accordingly

PhucPlugin::PhucPlugin()
{ 
        typeList << PHUC_FUNCTION;
        FilterIDType tt;
  
        foreach(tt , types()){
            actionList << new QAction(filterName(tt), this);
            if(tt == PHUC_FUNCTION){
                actionList.last()->setIcon(QIcon(":/imageFile/download.png"));
            }
        }
}

// ST() must return the very short string describing each filtering action 
// (this string is used also to define the menu entry)
QString PhucPlugin::filterName(FilterIDType filterId) const
{
  switch(filterId) {
                case PHUC_FUNCTION :  return QString("Generate Delaunay Triangulation");
		default : assert(0); 
	}
  return QString();
}

// Info() must return the longer string describing each filtering action 
// (this string is used in the About plugin dialog)
 QString PhucPlugin::filterInfo(FilterIDType filterId) const
{
  switch(filterId) {
                case PHUC_FUNCTION :  return QString("Phuc try to triangulate points.");
		default : assert(0); 
	}
	return QString("Unknown Filter");
}

// The FilterClass describes in which generic class of filters it fits. 
// This choice affect the submenu in which each filter will be placed 
// More than a single class can be choosen.
PhucPlugin::FilterClass PhucPlugin::getClass(QAction *a)
{
  switch(ID(a))
	{
                case PHUC_FUNCTION :  return MeshFilterInterface::Smoothing;
		default : assert(0); 
	}
	return MeshFilterInterface::Generic;
}

// This function define the needed parameters for each filter. Return true if the filter has some parameters
// it is called every time, so you can set the default value of parameters according to the mesh
// For each parameter you need to define, 
// - the name of the parameter, 
// - the string shown in the dialog 
// - the default value
// - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
void PhucPlugin::initParameterSet(QAction *action,MeshModel &m, RichParameterSet & parlst)
{
	 switch(ID(action))	 {
                case PHUC_FUNCTION :
             parlst.addParam(new RichFloat ("Radius",
                  1.5,
                  "Radius of ball ",
                  "Using Ball Pivot Algorithm.\n\n"
                  "If disabled the face normals will remains unchanged resulting in a visually pleasant effect."));

                        break;
											
		default : assert(0); 
	}
}
//----------Finish building the interface of function filter
//----------The below code is to declare programming function

// struct of point in 3D
struct point3D{
    float x;
    float y;
    float z;
};

// isLeft to check P2 is left of line P0P1
float isLeft(vcg::Point3f P0, vcg::Point3f P1,vcg::Point3f P2){
    return (P1.X() - P0.X())*(P2.Y() - P0.Y()) - (P2.X() - P0.X())*(P1.Y() - P0.Y());
}
// swap to swap position
void swap(vcg::Point3f *a, vcg::Point3f *b){
    vcg::Point3f t;
    t = *a;
    *a = *b;
    *b = t;
}

// sort array of point by increasing of x coordinate
// if x same - then sort by y coordinate
// After sorting vector, return a sorted vector.
std::vector<Point3f> sortIncreasing(std::vector<Point3f> pointVector){
    for(unsigned int i=0;i<pointVector.size()-1;++i){
        for(unsigned int j=0;j<pointVector.size()-1;++j){
            if(pointVector.at(j)[0] > pointVector.at(j+1)[0]){
                swap(&pointVector.at(j),&pointVector.at(j+1));
            }
            else{
                if(pointVector.at(j)[0] == pointVector.at(j+1)[0]){
                    if(pointVector.at(j)[1] > pointVector.at(j+1)[1]){
                        swap(&pointVector.at(j),&pointVector.at(j+1));
                    }
                }
            }
        }
    }
    return pointVector;

}

// the function find the point with min x min y and min x max y
// two point will be on convex hull
std::vector<int> PointMIN_MIN_MIN_MAX(std::vector<Point3f> pointVector){
    std::vector<int> TwoPoint;
    int minmin = 0,minmax;
    TwoPoint.push_back(minmin);// first point is head of sorted vector
    unsigned int i;
    for(i=1;i<pointVector.size()-1;++i){
        if(pointVector.at(i)[0] != pointVector.at(0)[0]){
            break;
        }
    }
    minmax = i-1; // second point is last point with same x coordinate
    TwoPoint.push_back(minmax);
    return TwoPoint;
}

// the function find the point with max xmin y and max x max y
// two point is on convex hull
std::vector<int> PointMAX_MIN_MAX_MAX(std::vector<Point3f> pointVector){
    std::vector<int> TwoPoint;
    int maxmax = pointVector.size()-1,maxmin;

    unsigned int i;
    for(i=maxmax-1;i>0;--i){
        if(pointVector.at(i)[0] != pointVector.at(maxmax)[0]){
            break;
        }
    }
    maxmin = i+1; // first point is last point of sorted vector
    TwoPoint.push_back(maxmin);
    TwoPoint.push_back(maxmax);// second point is head of max x in sorted vector
    return TwoPoint;
}
// Find comvex hull, first is finding lower hull
// lower hull is Pminmin to Pmaxmin
// upper hull is Pminmax to Pmaxmin
// assume that Pmaxmax = Pmaxmin because the right most is 1 point
std::vector<Point3f> LowerHull(std::vector<Point3f> pointVector){
    // find lower hull
    std::stack<Point3f> stk;
    std::vector<Point3f> LH;
    //std::vector<int> HullPoint = PointMIN_MIN_MIN_MAX(pointVector);
    std::vector<int> HullPoint1 = PointMAX_MIN_MAX_MAX(pointVector);
    //HullPoint.push_back(HullPoint1.at(0));
    //HullPoint.push_back(HullPoint1.at(1));
    stk.push(pointVector.at(0));
    int i = 1;
    //int minmax = HullPoint.at(1);
    int maxmin;
    if(HullPoint1.size() == 1){
        maxmin = pointVector.size()-1;
    }
    else{
        maxmin = pointVector.size()-2;
    }

    //int maxmax = HullPoint.at(3);
    while(++i <= maxmin){
        if(isLeft(pointVector.at(0),pointVector.at(maxmin),pointVector.at(i)) >= 0 && i < maxmin){
            continue;
        }
        while(stk.size() >= 2){
            vcg::Point3f first = stk.top();
            stk.pop();
            vcg::Point3f next = stk.top();
            stk.pop();
            if(isLeft(next,first,pointVector.at(i)) > 0){
                stk.push(next);
                stk.push(first);
                break;
            }
            else{

                stk.push(next);
            }
        }
        stk.push(pointVector.at(i));

    }
    while(!stk.empty()){
        LH.push_back(stk.top());
        stk.pop();
    }
    LH = sortIncreasing(LH);
    return LH;

}

std::vector<Point3f> UpperHull(std::vector<Point3f> pointVector){
    // find lower hull
    std::stack<Point3f> stk1;
    std::vector<Point3f> UH;
    std::vector<int> HullPoint = PointMIN_MIN_MIN_MAX(pointVector);
    std::vector<int> HullPoint1 = PointMAX_MIN_MAX_MAX(pointVector);
    HullPoint.push_back(HullPoint1.at(0));
    HullPoint.push_back(HullPoint1.at(1));

    int i = HullPoint.at(1);
    int minmax = HullPoint.at(1);
    int maxmin = HullPoint.at(2);
    int maxmax = HullPoint.at(3);
    if(maxmin != maxmax){
        stk1.push(pointVector.at(HullPoint.at(3)));// push Pmaxmax to stack
    }

    i = HullPoint.at(2);

    while(--i >= minmax){


        if(isLeft(pointVector.at(maxmax),pointVector.at(minmax),pointVector.at(i)) >=0 && i> minmax){
            //Log("%f\n",isLeft(pointVector.at(maxmax),pointVector.at(minmax),pointVector.at(i)));
            continue;

        }

        while(stk1.size()>=2){
            vcg::Point3f first = stk1.top();
            stk1.pop();
            vcg::Point3f next = stk1.top();
            stk1.pop();
            if(isLeft(next,first,pointVector.at(i)) > 0){
                stk1.push(next);
                stk1.push(first);
                break;
            }
            else{

                stk1.push(next);
            }
        }
        stk1.push(pointVector.at(i));
    }
    while(!stk1.empty()){
        UH.push_back(stk1.top());
        stk1.pop();
    }
    UH = sortIncreasing(UH);
    return UH;
}

// Internal function to check point z is internal of circumcircle of triangle a,b,c
// return true if z is internal of O(a,b,c) else false

int Internal(vcg::Point3f z, vcg::Point3f a, vcg::Point3f b,vcg::Point3f c){
    vcg::Point3f I; // I is center of triangle v,x,y;
    float D = 4*(a.X()-b.X())*(a.Y()-c.Y())-4*(a.X()-c.X())*(a.Y()-b.Y());
    float Dx=(-2)*(a.Y()-c.Y())*(b.X()*b.X()-a.X()*a.X()+b.Y()*b.Y()-a.Y()*a.Y())+2*(a.Y()-b.Y())*(c.X()*c.X()-a.X()*a.X()+c.Y()*c.Y()-a.Y()*a.Y());
    float Dy=(-2)*(a.X()-b.X())*(c.X()*c.X()-a.X()*a.X()+c.Y()*c.Y()-a.Y()*a.Y()) + 2*(a.X()-c.X())*(b.X()*b.X()-a.X()*a.X()+b.Y()*b.Y()-a.Y()*a.Y());
    if(D != 0){
        I.X() = Dx/D; // x coordinate of I
        I.Y() = Dy/D; // y coordinate of I
        //float R = (I.X()-a.X())*(I.X()-a.X()) + (I.Y()-a.Y())*(I.Y()-a.Y());
        float R = vcg::Distance(I,a); // R is distance between I and a
        float d = vcg::Distance(I,z);
        if(R>d){
            return 0; // return 0 if z is internal of O(a,b,c)
        }
        else{
            if(R<d){
                return 1;// return 1 if z is external of O(a,b,c)
            }
            else{
                return 2;// return 2 if z is on O(a,b,c)
            }
        }


    }
    else{
        return -1; // return -1 if equation has no root
    }

}

// SUCC return the point on counterclockwise of a and on convex hull
// pointVector is sorted increasing
//vcg::Point3f SUCC (std::vector<Point3f> pointVector,vcg::Point3f a ){

//}

// PRED return the point on clockwise of a and on convex hull
// pointVector is sorted increasing
//vcg::Point3f PRED (std::vector<Point3f> pointVector,vcg::Point3f a ){

//}

// Scan have some steps
// step1 : connect edge between vertex on convex hull
// step2 : scan the 4 points and check a point is internal circumcircle
// if point outside circumcircle, an edge is deleted
/*
void Scan(std::vector<Point3f> pointVector){
    vcg::Point3f first;
    // step1
    for(unsigned int i=0;i<pointVector.size()-2;++i){
        first = pointVector.at(i);
        for(unsigned int j=i+1;i<pointVector.size()-1;++i){
            // connect first to other point
        }
    }
    // finish step1.


}
*/

// when create new face, need to know how many face to add


// Function applyFilter is the main function of a program
bool PhucPlugin::applyFilter(QAction */*filter*/, MeshDocument &md, RichParameterSet & par, vcg::CallBackPos *cb)
{

        MeshModel &m=*md.mm();

        const double radius = (double)par.getFloat("Radius");
        Mesh mesh;
        mesh.DT(md);


        /*
        tri::Allocator<CMeshO>::AddVertices(m.cm,4);
        tri::Allocator<CMeshO>::AddFaces(m.cm,2);
        CMeshO::VertexPointer ivp[4];
        CMeshO::VertexIterator vi;
        unsigned int VertexCounter = 0;
        for(vi = m.cm.vert.begin();vi != m.cm.vert.end();++vi){
            ivp[VertexCounter]=&*vi;
            (*vi).P()=CMeshO::CoordType ( pointVector.at(VertexCounter)[0], pointVector.at(VertexCounter)[1], pointVector.at(VertexCounter)[2]);


        }

        CMeshO::FaceIterator fi;
        for(fi = m.cm.face.begin();fi !=m.cm.face.end();++fi){
            std::vector<int> face;
            for(int k=0;k<3;++k){
                for(int i=0;i<pointVector.size();++i){
                    if((pointVector.at(i)[0]==0.0 &&pointVector.at(i)[1]==0.0&&pointVector.at(i)[2]==0.0)){

                    }
                }
            }

            (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2];
        }
        */




        return true;

}

QString PhucPlugin::filterScriptFunctionName( FilterIDType filterID )
{
	switch(filterID) {
                case PHUC_FUNCTION :  return QString("randomVerticesDisplacement");
		default : assert(0); 
	}
	return QString();
}

Q_EXPORT_PLUGIN(PhucPlugin)
