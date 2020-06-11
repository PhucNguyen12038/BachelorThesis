#include "Mesh.h"
#include "Edge.h"
#include "Vertice.h"



Mesh::Mesh(){

}

void Mesh::printSignal(){
    qDebug()<<"hello world";
}

float Mesh::isLeft(vcg::Point3f P0, vcg::Point3f P1,vcg::Point3f P2){
    return (P1.X() - P0.X())*(P2.Y() - P0.Y()) - (P2.X() - P0.X())*(P1.Y() - P0.Y());
}

void Mesh::swap(vcg::Point3f *a, vcg::Point3f *b){
    vcg::Point3f t;
    t = *a;
    *a = *b;
    *b = t;
}
std::vector<Point3f> Mesh::sortIncreasing(std::vector<Point3f> pointVector){
    std::vector<Point3f> pointVector2 = pointVector;
    for(unsigned int i=0;i<pointVector2.size()-1;++i){
        for(unsigned int j=0;j<pointVector2.size()-1;++j){
            if(pointVector2.at(j)[0] > pointVector2.at(j+1)[0]){
                swap(&pointVector2.at(j),&pointVector2.at(j+1));
            }
            else{
                if(pointVector2.at(j)[0] == pointVector2.at(j+1)[0]){
                    if(pointVector2.at(j)[1] > pointVector2.at(j+1)[1]){
                        swap(&pointVector2.at(j),&pointVector2.at(j+1));
                    }
                }
            }
        }
    }
    return pointVector2;
}
std::vector<Point3f> Mesh::sortDecreasing(std::vector<Point3f> pointVector){
    std::vector<Point3f> pointVector2 = pointVector;
    for(unsigned int i=0;i<pointVector2.size()-1;++i){
        for(unsigned int j=0;j<pointVector2.size()-1;++j){
            if(pointVector2.at(j)[0] < pointVector2.at(j+1)[0]){
                swap(&pointVector2.at(j),&pointVector2.at(j+1));
            }
            else{
                if(pointVector2.at(j)[0] == pointVector2.at(j+1)[0]){
                    if(pointVector2.at(j)[1] < pointVector2.at(j+1)[1]){
                        swap(&pointVector2.at(j),&pointVector2.at(j+1));
                    }
                }
            }
        }
    }
    return pointVector2;
}
std::vector<int> Mesh::PointMIN_MIN_MIN_MAX(std::vector<Point3f> pointVector){
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
    if(minmax != 0){
        TwoPoint.push_back(minmax);
    }

    return TwoPoint;
}
std::vector<int> Mesh::PointMAX_MIN_MAX_MAX(std::vector<Point3f> pointVector){
    std::vector<int> TwoPoint;
    int maxmax = pointVector.size()-1,maxmin;

    unsigned int i;
    for(i=maxmax-1;i>0;--i){
        if(pointVector.at(i)[0] != pointVector.at(maxmax)[0]){
            break;
        }
    }
    maxmin = i+1; // first point is last point of sorted vector

    TwoPoint.push_back(maxmax);// second point is head of max x in sorted vector
    if(maxmin != maxmax){
        TwoPoint.push_back(maxmin);
    }
    return TwoPoint;
}
std::vector<Point3f> Mesh::LowerHull(std::vector<Point3f> pointVector){
    // find lower hull
    std::stack<Point3f> stk;
    std::vector<Point3f> LH;

    std::vector<int> HullPoint1 = PointMAX_MIN_MAX_MAX(pointVector);

    stk.push(pointVector.at(0));
    int i = 0;

    int maxmin;
    if(HullPoint1.size() == 1){
        maxmin = pointVector.size()-1;
    }
    else{
        maxmin = HullPoint1.at(0);
    }


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
std::vector<Point3f> Mesh::UpperHull(std::vector<Point3f> pointVector){
    // find lower hull
    std::stack<Point3f> stk1;
    std::vector<Point3f> UH;
    std::vector<int> HullPoint = PointMIN_MIN_MIN_MAX(pointVector);
    std::vector<int> HullPoint1 = PointMAX_MIN_MAX_MAX(pointVector);

    int minmax;
    if(HullPoint.size() == 1){
        minmax = 0;
    }
    else{
        minmax = HullPoint.at(1);
    }
    int maxmin;
    int maxmax = pointVector.size()-1;
    stk1.push(pointVector.at(maxmax));// push Pmaxmax to stack
    if(HullPoint1.size() ==1){
        maxmin = maxmax;
    }
    else{
        maxmin = HullPoint1.at(0);
    }

    int i = maxmin;

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
    UH = sortDecreasing(UH);
    return UH;
}


std::vector<Point3f> Mesh::ConvexHull(std::vector<Point3f> pointVector){
    std::vector<Point3f> LH = LowerHull(pointVector);
    std::vector<Point3f> UH = UpperHull(pointVector);

    //if(LH.at(0) == UH.at(UH.size()-1)){
    //    LH.erase(LH.begin());
    //}

    if(LH.at(LH.size()-1) == UH.at(0)){
        LH.erase(LH.end());
    }


    std::vector<Point3f> CH;
    for(unsigned int i=0;i<LH.size();++i){

        CH.push_back(LH.at(i));
    }
    for(unsigned int i=0;i<UH.size();++i){
        CH.push_back(UH.at(i));
    }

    //CH = sortIncreasing(CH);

    return CH;

}
void Mesh::DT(MeshDocument &md){
    MeshModel &m = *md.mm();
    vcg::Point3f p;
    std::vector<Point3f> pointVector; // use for storing the order of meshlab
    std::vector<Point3f> pointVector2; // use for sorting the element of points
    std::vector<Vertice*> pointList;
    Vertice* VertexArray[400][400] = {NULL};
    for(unsigned int i = 0; i< m.cm.vert.size(); i++){
            p = m.cm.vert[i].P();// assign p = input Vertice from file xyz
            pointVector.push_back(p);
    }

    pointVector2 = sortIncreasing(pointVector);
    for(unsigned int i=0;i<pointVector2.size();++i){
        Vertice *temp = new Vertice(pointVector2.at(i)[0],pointVector2.at(i)[1],pointVector2.at(i)[2]);
        pointList.push_back(temp);
        VertexArray[(int)pointVector2.at(i)[0]][(int)pointVector2.at(i)[1]] = temp;
    }

    Mesher mesher;
    mesher.GenerateDelaunayTriangle(pointList,mesher,VertexArray);

    m.cm.Clear();
    Vertice_star_list::const_iterator v;
    tri::Allocator<CMeshO>::AddVertices(m.cm,mesher.nVertices());
    tri::Allocator<CMeshO>::AddFaces(m.cm,mesher.nFacets());
    CMeshO::VertexPointer ivp[mesher.nVertices()];
    CMeshO::VertexIterator vi;
    unsigned int VertexCounter = 0;


    for(vi = m.cm.vert.begin();vi != m.cm.vert.end();++vi){
        ivp[VertexCounter]=&*vi;
        (*vi).P()=CMeshO::CoordType ( pointVector.at(VertexCounter)[0], pointVector.at(VertexCounter)[1], pointVector.at(VertexCounter)[2]);
        ++VertexCounter;
        //qDebug()<<VertexCounter;
    }


    CMeshO::FaceIterator fi2 = m.cm.face.begin();
    unsigned int faceCounter = 0;
    Facet_star_list::const_iterator fi ;

    for( fi = mesher.facets_begin(); fi != mesher.facets_end(); ++fi)
    {
        Facet *f = *fi;
        std::vector<int> VertexIndex;

        for(unsigned int i=0;i<3;++i){

            Vertice *v = f->vertice(i);

            for(unsigned int j=0;j<pointVector.size();++j){
                if((v->x() == pointVector.at(j)[0]) &&(v->y() == pointVector.at(j)[1])&&(v->z() == pointVector.at(j)[2])){
                    VertexIndex.push_back(j);

                }
            }

        }

        (*fi2).V(0)=ivp[VertexIndex.at(0)];
        (*fi2).V(1)=ivp[VertexIndex.at(1)];
        (*fi2).V(2)=ivp[VertexIndex.at(2)];
        if(faceCounter <mesher.nFacets()-1){
            ++fi2;
            ++faceCounter;
        }
        VertexIndex.clear();

    }



}


