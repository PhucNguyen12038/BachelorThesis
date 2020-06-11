/**
 * @file Mesher.cpp
 * @brief defines methods for building a surface mesh from points stored in an
 * octree these methods are declared in Mesher.h
 * @author Julie Digne julie.digne@liris.cnrs.fr
 * @date 2012/10/17
 * @copyright This file implements an algorithm possibly linked to the patent
 * US6968299B1.
 * This file is made available for the exclusive aim of serving as
 * scientific tool to verify the soundness and completeness of the
 * algorithm description. Compilation, execution and redistribution
 * of this file may violate patents rights in certain countries.
 * The situation being different for every country and changing
 * over time, it is your responsibility to determine which patent
 * rights restrictions apply to you before you compile, use,
 * modify, or redistribute this file. A patent lawyer is qualified
 * to make this determination.
 * If and only if they don't conflict with any patent terms, you
 * can benefit from the following license terms attached to this
 * file.
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANy WARRANTy; without even the implied warranty of
 * MERCHANTABILITy or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * you should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "Mesher.h"

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <sstream>
#include <queue>
#include <math.h>



using namespace std;

Mesher::Mesher()
{


    m_nfacets = 0;
    m_nvertices = 0;
    m_nEdges = 0;
}

Mesher::~Mesher()
{
    m_edge_front.clear();



    Facet_star_list::iterator fi;
    for(fi = m_facets.begin(); fi != m_facets.end(); ++fi)
    {
        delete *fi;
        *fi = NULL;
    }

    m_facets.clear();
    m_vertices.clear();
    m_nfacets = 0;
    m_nvertices = 0;
}

unsigned int Mesher::nVertices() const
{
    return m_nvertices;
}

unsigned int Mesher::nFacets() const
{
    return m_nfacets;
}

Edge* Mesher::getCommonEdge(Facet *f1, Facet *f2){
    std::vector<Vertice*> common;
    for(unsigned int i=0;i<3;++i){

        Vertice *v = f1->vertice(i);
        if(f2->hasVertice(v)){
            common.push_back(v);
            if(common.size() == 2){
                break;
            }
        }

    }
    if(common.size() < 2){
        return NULL;
    }
    else{
        Vertice *v1 = common.at(0);
        Vertice *v2 = common.at(1);
        return v1->getLinkingEdge(v2);

    }
}
// Internal function to check point z is internal of circumcircle of triangle a,b,c
// return true if z is internal of O(a,b,c) else false

int Mesher::Internal(Vertice *z, Vertice *a, Vertice *b,Vertice *c){
    // I is center of triangle v,x,y;

    float D = 4*(a->x()-b->x())*(a->y()-c->y())-4*(a->x()-c->x())*(a->y()-b->y());
    float Dx=(-2)*(a->y()-c->y())*(b->x()*b->x()-a->x()*a->x()+b->y()*b->y()-a->y()*a->y())+2*(a->y()-b->y())*(c->x()*c->x()-a->x()*a->x()+c->y()*c->y()-a->y()*a->y());
    float Dy=(-2)*(a->x()-b->x())*(c->x()*c->x()-a->x()*a->x()+c->y()*c->y()-a->y()*a->y()) + 2*(a->x()-c->x())*(b->x()*b->x()-a->x()*a->x()+b->y()*b->y()-a->y()*a->y());
    if(D != 0){
        Vertice *I = new Vertice(Dx/D,Dy/D,0);
        //qDebug()<<I->x();
        //qDebug()<<I->y();
        //qDebug()<<"center I";
        //float R = (I.x()-a.x())*(I.x()-a.x()) + (I.y()-a.y())*(I.y()-a.y());
        float R = dist2(I,a); // R is distance between I and a
        float d = dist2(I,z);
        if(R>=d){
            return 0; // return 0 if z is internal of O(a,b,c)
        }
        else{
            return 1;// return 1 if z is external of O(a,b,c)
        }

    }
    else{
        return -1; // return -1 if equation has no root
    }

}


void Mesher::addFacet(Facet* f)
{
    addVertice( f->vertice(0) );
    addVertice( f->vertice(1) );
    addVertice( f->vertice(2) );

    m_facets.push_back(f);
    m_nfacets++;
}


void Mesher::addVertice(Vertice* v)
{
    if(v->index() != -1)
        return;

    v->setIndex(m_nvertices);
    m_vertices.push_back(v);
    m_nvertices++;
}

void Mesher::setNumberOfFacet(unsigned int n){
    this->m_nfacets = n;
}


std::list< Vertice* >::const_iterator Mesher::vertices_begin() const
{
    return m_vertices.begin();
}

std::list< Vertice* >::const_iterator Mesher::vertices_end() const
{
    return m_vertices.end();
}
std::list< Facet* >::const_iterator Mesher::facets_begin() const
{
    return m_facets.begin();
}

std::list< Facet* >::const_iterator Mesher::facets_end() const
{
    return m_facets.end();
}

void fixFace(Mesher mesher){
    /*
    Facet_star_list::const_iterator fi ;

    for( fi = mesher.facets_begin(); fi != mesher.facets_end(); ++fi)
    {
        Facet *f = *fi;
        std::vector<int> VertexIndex;

        for(unsigned int i=0;i<3;++i){

            Vertice *v = f->vertice(i);
            //qDebug()<<v->x();
            //qDebug()<<v->y();
            //qDebug()<<v->z();
            for(unsigned int j=0;j<pointVector.size();++j){
                if((v->x() == pointVector.at(j)[0]) &&(v->y() == pointVector.at(j)[1])&&(v->z() == pointVector.at(j)[2])){
                    VertexIndex.push_back(j);

                }
            }

        }
    }
    */
}

void Mesher::getBorderEdgeList(std::vector<Edge*> &BorderList, Mesher &mesher){

    Facet_star_list::const_iterator fi ;
    for( fi = mesher.facets_begin(); fi != mesher.facets_end(); ++fi)
    {
        Facet *f = *fi;
        for(unsigned int i=0;i<3;++i){
            Edge* e = f->edge(i);
            if(e->getFacet2() == NULL){
                /*
               //qDebug()<<e->getSource()->x();
               //qDebug()<<e->getSource()->y();
               //qDebug()<<"1st border vertice";
               //qDebug()<<e->getTarget()->x();
               //qDebug()<<e->getTarget()->y();
               //qDebug()<<"2nd border vertice";
               */
                BorderList.push_back(e);
            }
        }
    }

}

void Mesher::FindNormalToEdge(Edge *e, double &Xdirect, double &Ydirect){
    Vertice* v4 = NULL; // v4 is the vertex combine with e to make triangle
    // v1 is source vertex of edge, v2 is target vertex of edge
    if(e->getFacet2() != NULL){
        v4 = e->getOppositeVertice();
        double XvectorLine = e->getSource()->x() - e->getTarget()->x();
        double YvectorLine = e->getSource()->y() - e->getTarget()->y();
        double XnormalOfLine = 0 - YvectorLine;
        double YnormalOfLine = XvectorLine;
        double Xvector42 = v4->x() - e->getTarget()->x();
        double Yvectir42 = v4->y() - e->getTarget()->y();
        double crossProduct2D = (XnormalOfLine*Yvectir42) - (Xvector42*YnormalOfLine);
        if(crossProduct2D > 0){
            Xdirect = XnormalOfLine;
            Ydirect = YnormalOfLine;
        }
        if(Xdirect > Ydirect){
            Xdirect = Xdirect / Xdirect;
            Ydirect = Ydirect / Xdirect;
        }
        else{
            Xdirect = Xdirect / Ydirect;
            Ydirect = Ydirect / Ydirect;
        }
        //qDebug()<<Xdirect;
        //qDebug()<<Ydirect;
        //qDebug()<<"the normal to edge\n";
    }
}
void Mesher::normalize(double &Xdirect, double &Ydirect){
    double Xtemp= Xdirect;
    double Ytemp = Ydirect;
    if(Xdirect ==0 || Ydirect ==0){
        if(Xdirect == 0){
            Ydirect = Ydirect / abs(Ydirect);
        }
        else{
            Xdirect = Xdirect / abs(Xdirect);
        }
    }
    else{
        if(abs(Xdirect) > abs(Ydirect)){
            Xdirect = (double)Xtemp / abs(Xtemp);
            Ydirect = (double)Ytemp / abs(Xtemp);
        }
        else{
            Xdirect = (double)Xtemp / abs(Ytemp);
            Ydirect = (double)Ytemp / abs(Ytemp);
        }
    }

}

void Mesher::FindVectorOfEdge(Edge* e, double &XdirectE, double &YdirectE){
    XdirectE = e->getSource()->x() - e->getTarget()->x();
    YdirectE = e->getSource()->y() - e->getTarget()->y();
    if(XdirectE > YdirectE){
        XdirectE = XdirectE / XdirectE;
        YdirectE = YdirectE / XdirectE;
    }
    else{
        XdirectE = XdirectE / YdirectE;
        YdirectE = YdirectE / YdirectE;
    }
}

std::list<Vertice*> Mesher::WooLine(Edge* e, Vertice* pointArray){

    double Xdirect;
    double Ydirect;
    // find normal vector
    FindNormalToEdge(e,Xdirect,Ydirect);
    // find vector of edge e
    double XvectorEdge;
    double YvectorEdge;
    FindVectorOfEdge(e,XvectorEdge,YvectorEdge);
    // find mid point
    Vertice midPoint = midpoint(e->getSource(),e->getTarget());
    std::list<Vertice*> candidateVertice;
    unsigned int Xcoordinate = 0, Ycoordinate=0;
    Xcoordinate = Xcoordinate + midPoint.x() + Xdirect;
    Ycoordinate = Ycoordinate + midPoint.y() + Ydirect;

    if( Xcoordinate>=0&&Ycoordinate>=0){
        Vertice* v = new Vertice((double)Xcoordinate,(double)Ycoordinate,0.0);
        //qDebug()<<Xcoordinate;
        //qDebug()<<Ycoordinate;
        //qDebug()<<"WooLine Point";
        candidateVertice.push_back(v);
    }

    unsigned int leftXcoordinate = Xcoordinate + XvectorEdge;
    unsigned int leftYcoordinate = Ycoordinate + YvectorEdge;
    if(Xcoordinate>=0&&Ycoordinate>=0){
        Vertice* v = new Vertice((double)leftXcoordinate,(double)leftXcoordinate,0.0);
        //qDebug()<<leftXcoordinate;
        //qDebug()<<leftYcoordinate;
        //qDebug()<<"WooLine Point";
        candidateVertice.push_back(v);
    }

    unsigned int rightXcoordinate = Xcoordinate - XvectorEdge;
    unsigned int rightYcoordinate = Ycoordinate - YvectorEdge;
    if(Xcoordinate>=0&&Ycoordinate>=0){
        Vertice* v = new Vertice((double)rightXcoordinate,(double)rightYcoordinate,0.0);
        //qDebug()<<rightXcoordinate;
        //qDebug()<<rightXcoordinate;
        //qDebug()<<"WooLine Point";
        candidateVertice.push_back(v);
    }

    return candidateVertice;

}

int Mesher::Angle(double xv1, double yv1, double xv2, double yv2){
    double length1 = sqrt(xv1*xv1+yv1*yv1); 
    double length2 = sqrt(xv2*xv2+yv2*yv2);
    double angle = (xv1*xv2+yv1*yv2)/(length1*length2);
    if(angle>0 && angle <=1){
        return 0; // 2 vectors are same side
    }
    else{
        if(angle<0&&angle>= -1){
            return 1; // 2 vectors are opposite side
        }
        else{
            return -1; // 2 vector are perpendicular
        }
    }

}

void Mesher::GenerateDelaunayTriangle(std::vector<Vertice*> VerticeVector, Mesher &mesher, Vertice* pointArray[400][400]){

    //----- Step 1: find seed triangle

    // the first step is finding the vertex which is closest to first vertex
    Vertice *firstVertice = VerticeVector.at(0);
    Vertice *lastVertice = VerticeVector.at(VerticeVector.size() - 1);
    Vertice* MinYMaxY = new Vertice(firstVertice->y(),lastVertice->y(),0);


    Vertice *secondVertice = VerticeVector.at(1);
    Vertice *closetVertice;
    closetVertice = secondVertice;
    double minDistance = dist2(firstVertice,secondVertice);
    // find the maximum and minimum of x,y coordinate
    for(unsigned int i=1;i<VerticeVector.size();++i){
        double distance = dist2(firstVertice,VerticeVector.at(i));
        if(minDistance > distance){
            minDistance = distance;
            closetVertice = VerticeVector.at(i);
        }
        Vertice* vv = VerticeVector.at(i);
        if(MinYMaxY->x() > vv->y()){
            MinYMaxY->setX(vv->y());
        }
        if(MinYMaxY->y() < vv->y()){
            MinYMaxY->setY(vv->y());
        }
    }



    // finish find the closet vertex with the first vertex
    // the second step is finding the third vertex such that triangle satisfy DT

    Vertice *thirdVertice = NULL;
    // counter to count the number of vertex external triangle
    // if counter = number of vertex - 3
    for(unsigned int i=0;i<VerticeVector.size();++i){
        Vertice *v = VerticeVector.at(i);
        if(v->compareWithVertex(firstVertice)==false && v->compareWithVertex(closetVertice)==false){
            if(isOnTheLine(firstVertice,closetVertice,v)==false){

                unsigned int counter = 0;
                for(unsigned int j=0;j<VerticeVector.size();++j){
                    Vertice *vNext = VerticeVector.at(j);
                    if(vNext->compareWithVertex(firstVertice)==false && vNext->compareWithVertex(closetVertice)==false && vNext->compareWithVertex(v)==false){

                        int checkInternal = isInternal(vNext,firstVertice,closetVertice,v);
                        if(checkInternal == 0){

                            break;
                        }
                        else{
                            if(checkInternal == 1){
                                counter+=1;
                            }
                        }
                    }
                }
                if(counter == VerticeVector.size()-3){
                    thirdVertice = v;
                    break;
                }
            }
        }
    }
    //finish finding the third vertex
    // if NULL so we can not find seed triangle then no facet is created
    if(thirdVertice == NULL){
        //qDebug()<<"third vertex is null";
        for(unsigned int i=0;i<VerticeVector.size();++i){
            mesher.addVertice(VerticeVector.at(i));
        }
        mesher.setNumberOfFacet(0);


    }
    else{
        // if we found seed triangle
        //qDebug()<<thirdVertice->x();
        //qDebug()<<thirdVertice->y();
        //qDebug()<<thirdVertice->z();
        //qDebug()<<"the third vertex";
        // after 2 loop for, the vertex with the first edge will make seed triangle
        // so that no points is internal circum circle of triangle
        //Vertice* choose = pointArray[1][10];
        Facet *seedTriangle = new Facet(firstVertice,closetVertice,thirdVertice);
        mesher.addFacet(seedTriangle);

        //----- Finish step 1: seed triangle is found

        //----- Step2: push edge into queue and expand triangle

        Edge* e1 = firstVertice->getLinkingEdge(thirdVertice);
        Edge* e2 = thirdVertice->getLinkingEdge(closetVertice);
        Edge* e3 = closetVertice->getLinkingEdge(firstVertice);
        std::list<Edge*> EdgeQueue;
        Edge_star_list EdgeOfMesherList;

        //WooLine(e1,pointArray);


        EdgeQueue.push_back(e1);
        EdgeQueue.push_back(e2);
        EdgeQueue.push_back(e3);
        // put edge of seed triangle to a queue, start to find new vertex
        // take out edge and find new vertex
        // new vertex is opposite with the vertex which is in facet contain edge
        // if new vertex satisfy Delaunay Trianglulation, create facet
        // continue putting new edge to queue
        // do until the queue is empty

        while(!EdgeQueue.empty()){
            Edge *e = EdgeQueue.front();
            EdgeQueue.pop_front();

            double Xdirect;
            double Ydirect;
            double XvectorLine;
            double YvectorLine;
            Vertice* v4 = NULL;
            // find normal vector
            // process edge has only one facet
            if(e->getFacet2() == NULL){
                Vertice midPoint = midpoint(e->getSource(),e->getTarget());

                //v4 is vertex combine with edge e to create facet1 of e
                v4 = e->getOppositeVertice();

                // find vector of line
                XvectorLine = e->getSource()->x() - e->getTarget()->x();
                YvectorLine = e->getSource()->y() - e->getTarget()->y();

                // find normal vector of line
                double XnormalOfLine = 0 - YvectorLine;
                double YnormalOfLine = XvectorLine;

                // find the direct of WooLine base on the angle
                // find vector4_midPoint of mid point of edge and v4
                double Xvector4_midPoint = v4->x()-midPoint.x();
                double Yvector4_midPoint = v4->y()-midPoint.y();

                // angle between vector4_midPoint and normal vector
                int angle = Angle(XnormalOfLine,YnormalOfLine,Xvector4_midPoint,Yvector4_midPoint);

                if(angle == 1){ // if 2 vector are on differrent side
                    Xdirect = XnormalOfLine;
                    Ydirect = YnormalOfLine;
                }
                else if(angle == 0){ // same side
                    Xdirect = -XnormalOfLine;
                    Ydirect = -YnormalOfLine;

                }

                //qDebug()<<"Mid point";
                //qDebug()<<midPoint.x();
                //qDebug()<<midPoint.y();

                //qDebug()<<"Vector line";
                //qDebug()<<XvectorLine;
                //qDebug()<<YvectorLine;


                // normalize the WooLine direction
                normalize(Xdirect,Ydirect);

                //qDebug()<<"Normal of line";
                //qDebug()<<Xdirect;
                //qDebug()<<Ydirect;

                // find the opposite edge of edge e
                double XoppositeSource = e->getSource()->x()+Xdirect;
                double YoppositeSource = e->getSource()->y()+Ydirect;
                double XoppositeTarget = e->getTarget()->x()+Xdirect;
                double YoppositeTarget = e->getTarget()->y()+Ydirect;

                int d=0;
                // ignore the cross line, care the horizontial and vertical line
                if((int)(e->getSource()->x()+e->getTarget()->x())%2==1
                   && (int)(e->getSource()->y()+e->getTarget()->y())%2==1){
                    // cross line
                    d=0;

                }
                else{ // not cross line
                    Vertice* srcOpp = new Vertice(XoppositeSource,YoppositeSource,0.0);
                    Vertice* tgtOpp = new Vertice(XoppositeTarget,YoppositeTarget,0.0);
                    // edgeOpp is the opposite edge of e
                    Edge* edgeOpp = new Edge(srcOpp,tgtOpp);
                    Edge_star_list::iterator EI;
                    // check if edgeOpp in queue
                    for(EI = EdgeOfMesherList.begin();EI != EdgeOfMesherList.end();++EI){
                        Edge* edgeMesh = *EI;
                        if(edgeOpp->CompareWithEdge(edgeMesh)){
                            // if edgeOpp has 2 facet, delete edge e out of queue
                            if(edgeMesh->getFacet2() != NULL){
                                d+=1;
                                break;
                            }

                            else if(edgeMesh->getFacet1() != NULL){
                                // if edgeOpp has 1 facet and v4 of edgeOpp is in
                                // edge e, then delete e out of queue
                                Vertice* oppVertex = edgeMesh->getOppositeVertice();
                                if(e->includeVertice(oppVertex)==true){
                                    d+=1;
                                    break;
                                }

                            }

                        }
                    }
                }


                if(d!=0){
                    continue;
                }

                // Apply WooLine to find candidate vertex
                std::vector<Vertice*> candidateVertice;
                double Xcoordinate, Ycoordinate;
                candidateVertice.clear();
                if(candidateVertice.size()==0){
                    Xcoordinate = midPoint.x() + Xdirect;
                    Ycoordinate = midPoint.y() + Ydirect;



                    if( Xcoordinate>=firstVertice->x()-0.5&&Ycoordinate>=MinYMaxY->x()-0.5&&
                        Xcoordinate<=lastVertice->x()+0.5&&Ycoordinate<=MinYMaxY->y()+0.5){
                        double XcoordinateUp = ceil(Xcoordinate);
                        double YcoordinateUp = ceil(Ycoordinate);

                        if( XcoordinateUp>=firstVertice->x()&&YcoordinateUp>=(int)MinYMaxY->x()&&
                            XcoordinateUp<=lastVertice->x()&&YcoordinateUp<=(int)MinYMaxY->y()){
                            Vertice* v = new Vertice((double)XcoordinateUp,(double)YcoordinateUp,0.0);

                            if(pointArray[(int)XcoordinateUp][(int)YcoordinateUp] != NULL){
                                //qDebug()<<"point scale up";
                                //qDebug()<<XcoordinateUp;
                                //qDebug()<<YcoordinateUp;
                                candidateVertice.push_back(v);
                            }
                        }
                        // point up

                        double XcoordinateDown = floor(Xcoordinate);
                        double YcoordinateDown = floor(Ycoordinate);
                        if( XcoordinateDown>=firstVertice->x()&&YcoordinateDown>=(int)MinYMaxY->x()&&
                            XcoordinateDown<=lastVertice->x()&&YcoordinateDown<=(int)MinYMaxY->y()){
                            Vertice* v1 = new Vertice((double)XcoordinateDown,(double)YcoordinateDown,0.0);

                            if(pointArray[(int)XcoordinateDown][(int)YcoordinateDown] != NULL){
                                //qDebug()<<"point scale down";
                                //qDebug()<<XcoordinateDown;
                                //qDebug()<<YcoordinateDown;
                                candidateVertice.push_back(v1);
                            }
                        }
                        // point down
                        if( XcoordinateUp>=firstVertice->x()&&YcoordinateDown>=(int)MinYMaxY->x()&&
                            XcoordinateUp<=lastVertice->x()&&YcoordinateDown<=(int)MinYMaxY->y()){
                            Vertice* v = new Vertice((double)XcoordinateUp,(double)YcoordinateDown,0.0);

                            if(pointArray[(int)XcoordinateUp][(int)YcoordinateDown] != NULL){
                                //qDebug()<<"point scale up down";
                                //qDebug()<<XcoordinateUp;
                                //qDebug()<<YcoordinateDown;
                                if(YcoordinateDown != YcoordinateUp){
                                    candidateVertice.push_back(v);
                                }

                            }
                        }
                        // point up down
                        if( XcoordinateDown>=firstVertice->x()&&YcoordinateUp>=(int)MinYMaxY->x()&&
                            XcoordinateDown<=lastVertice->x()&&YcoordinateUp<=(int)MinYMaxY->y()){
                            Vertice* v = new Vertice((double)XcoordinateDown,(double)YcoordinateUp,0.0);

                            if(pointArray[(int)XcoordinateDown][(int)YcoordinateUp] != NULL){
                                //qDebug()<<"point scale down up";
                                //qDebug()<<XcoordinateDown;
                                //qDebug()<<YcoordinateUp;
                                if(XcoordinateDown != XcoordinateUp){
                                    candidateVertice.push_back(v);
                                }

                            }
                        }
                        // point down up

                        /*
                            int leftXcoordinate = XcoordinateDown + XvectorLine;
                            int leftYcoordinate = YcoordinateDown + YvectorLine;
                            if(leftXcoordinate>=firstVertice->x()&&leftYcoordinate>=(int)MinYMaxY->x()&&
                                leftXcoordinate<=lastVertice->x()&&leftYcoordinate<=(int)MinYMaxY->y()){
                                Vertice* v = new Vertice((double)leftXcoordinate,(double)leftYcoordinate,0.0);

                                if(pointArray[leftXcoordinate][leftYcoordinate] != NULL){
                                    //qDebug()<<leftXcoordinate;
                                    //qDebug()<<leftYcoordinate;
                                    //qDebug()<<"WooLine left Point";
                                    candidateVertice.push_back(v);
                                }
                            }

                            int rightXcoordinate = XcoordinateDown - XvectorLine;
                            int rightYcoordinate = YcoordinateDown - YvectorLine;
                            if(rightXcoordinate>=firstVertice->x()&&rightYcoordinate>=(int)MinYMaxY->x()&&
                               rightXcoordinate<=lastVertice->x()&&rightYcoordinate<=(int)MinYMaxY->y()){
                                Vertice* v = new Vertice((double)rightXcoordinate,(double)rightYcoordinate,0.0);

                                if(pointArray[rightXcoordinate][rightYcoordinate] != NULL){
                                    //qDebug()<<rightXcoordinate;
                                    //qDebug()<<rightYcoordinate;
                                    //qDebug()<<"WooLine right Point";
                                    candidateVertice.push_back(v);
                                }
                            }
                            */


                    }



                }

                if(candidateVertice.size() == 1){
                    Vertice* v = candidateVertice.at(0);
                    int xCoor = (int) v->x();
                    int yCoor = (int) v->y();
                    //qDebug()<<"Candidate 1 vertex";
                    //qDebug()<<xCoor;
                    //qDebug()<<yCoor;

                    Vertice* triangleVertice = pointArray[xCoor][yCoor];
                    Facet *nf = new Facet(e,triangleVertice);
                    mesher.addFacet(nf);
                    Edge *newEdge1 = triangleVertice->getLinkingEdge(e->getSource());
                    Edge *newEdge2 = e->getTarget()->getLinkingEdge(triangleVertice);
                    EdgeQueue.push_back(newEdge1);
                    EdgeQueue.push_back(newEdge2);
                    EdgeOfMesherList.push_front(newEdge1);
                    EdgeOfMesherList.push_front(newEdge2);
                }

                else if(candidateVertice.size() > 1){
                    for(unsigned int i=0;i<candidateVertice.size();++i){
                        Vertice* v = candidateVertice.at(i);
                        unsigned int count = 0;
                        for(unsigned int j=0;j<candidateVertice.size();++j){
                            if(j != i){
                                Vertice* vtest = candidateVertice.at(j);
                                int checkInternal = isInternal(vtest,e->getSource(),e->getTarget(),v);
                                if(checkInternal == 1){
                                    count+=1;
                                }
                            }
                        }

                        if(count == candidateVertice.size() - 1){
                            int xCoor = (int) v->x();
                            int yCoor = (int) v->y();
                            //qDebug()<<"Candidate vertex";
                            //qDebug()<<xCoor;
                            //qDebug()<<yCoor;
                            Vertice* triangleVertice = pointArray[xCoor][yCoor];
                            Facet *nf = new Facet(e,triangleVertice);
                            mesher.addFacet(nf);
                            Edge *newEdge1 = triangleVertice->getLinkingEdge(e->getSource());
                            Edge *newEdge2 = e->getTarget()->getLinkingEdge(triangleVertice);
                            EdgeQueue.push_back(newEdge1);
                            EdgeQueue.push_back(newEdge2);
                            EdgeOfMesherList.push_front(newEdge1);
                            EdgeOfMesherList.push_front(newEdge2);
                            break;
                        }

                    }
                }



                // find mid point



            }

        }
        // finish taking out the edge in queue
        // all vertex is inserted into mesher
        // but not enough facets at boundary edges

        // we get the border edge because in some case, the facet is not enough
        // we make new border edge
        // first step is finding all border edge of mesh
        // second step is finding a vertex which is opposite of edge and closet to edge
        // third step is triangulating mesh

    }


    if(mesher.nVertices() != VerticeVector.size()){
        for(unsigned int i=0;i<VerticeVector.size();++i){
            mesher.addVertice(VerticeVector.at(i));
        }
    }





}











