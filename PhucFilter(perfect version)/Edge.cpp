/**
 * @file Edge.cpp
 * @author Julie Digne julie.digne@liris.cnrs.fr
 * @date 2012-10-08
 * @brief implementation of the edge methods declared in Edge.h
 *
 * @copyright This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "Edge.h"
#include "Vertice.h"
#include "Facet.h"
#include "utilities.h"

using namespace std;

Edge::Edge()
{
    m_src = NULL;
    m_tgt = NULL;
    m_facet1 = NULL;
    m_facet2 = NULL;
}

Edge::Edge(Vertice* src, Vertice* tgt)
{
    m_src = src;
    m_tgt = tgt;
    src->addAdjacentEdge(this);
    tgt->addAdjacentEdge(this);
    m_facet1 = NULL;
    m_facet2 = NULL;
    setType(1);
}

Edge::~Edge()
{
    m_src = NULL;
    m_tgt = NULL;
    m_facet1 = NULL;
    m_facet2 = NULL;
}

Vertice* Edge::getSource() const
{
    return m_src;
}

Vertice* Edge::getTarget() const
{
    return m_tgt;
}


Facet* Edge::getFacet1() const
{
    return m_facet1;
}

Facet* Edge::getFacet2() const
{
    return m_facet2;
}

void Edge::updateOrientation()
{
    // already exist edge
    // get opposite vertex of edge
    /*
    Vertice *opp = getOppositeVertice();

    double vx, vy, vz;
    // find cross product of opp and edge
    cross_product(m_tgt->x() - m_src->x(), m_tgt->y() - m_src->y(),
                  m_tgt->z() - m_src->z(), opp->x() - m_src->x(),
                  opp->y() - m_src->y(), opp->z() - m_src->z(),
                  vx, vy, vz);
    normalize(vx, vy, vz);

    double nx, ny, nz;
    nx = m_src->nx() + m_tgt->nx() + opp->nx();
    ny = m_src->ny() + m_tgt->ny() + opp->ny();
    nz = m_src->nz() + m_tgt->nz() + opp->nz();
    normalize(nx, ny, nz);

    if(   vx * nx + vy * ny + vz * nz  < 0)
    {
        Vertice *temp = m_src;
        m_src = m_tgt;
        m_tgt = temp;
    }
    */
}

// already exist an aedge. Now need to find 2 facet of edge
bool Edge::addAdjacentFacet(Facet* facet)
{
    if ((m_facet1 == facet)||(m_facet2 == facet))
        return false;

    if(m_facet1 == NULL)
    {
        m_facet1 = facet;
        updateOrientation();
        setType(1);
        return true;
    }

    if(m_facet2 == NULL)
    {
        m_facet2 = facet;
        setType(2);
        return true;
    }

    std::cout<<"Already two triangles"<<endl;
    return false;

}


bool Edge::removeAdjacentFacet(Facet* facet)
{

    if(m_facet1 == facet)
    {
        m_facet1 = NULL;
        setType(1);
        return true;
    }

    if(m_facet2 == facet)
    {
        m_facet2 = NULL;
        setType(1);
        return true;
    }

    return false;

}



bool Edge::hasVertice(Vertice *Vertice) const
{
    if((Vertice == m_src) || (Vertice == m_tgt))
        return true;
    return false;
}
bool Edge::includeVertice(Vertice *vertice) const
{
    if(m_src->x() == vertice->x() && m_src->y() == vertice->y()){
        return true;
    }
    else{
        if(m_tgt->x() == vertice->x() && m_tgt->y() == vertice->y()){
            return true;
        }
    }
    return false;
}

bool Edge::isInnerEdge() const
{
    if(m_facet2 == NULL)
        return false;
    return true;
}


int Edge::getType() const
{
    return m_type;
}

void Edge::setType(int type)
{
    m_type = type;
}


Vertice* Edge::getOppositeVertice() const
{
    if(m_facet1 == NULL)
        return NULL;
    Vertice *opp;
    for(int i = 0; i < 3; i++)
    {
        opp = m_facet1->vertice(i);
        if((opp != m_src) && (opp != m_tgt))
            return opp;
    }
    return NULL;
}

bool Edge::CompareWithEdge(Edge *testEdge){
    if(this->getSource()->x() == testEdge->getSource()->x()&&this->getSource()->y() == testEdge->getSource()->y()){
        if(this->getTarget()->x() == testEdge->getTarget()->x()&&this->getTarget()->y() == testEdge->getTarget()->y()){
            return true;
        }
    }
    else{
        if(this->getSource()->x() == testEdge->getTarget()->x()&&this->getSource()->y() == testEdge->getTarget()->y()){
            if(this->getTarget()->x() == testEdge->getSource()->x()&&this->getTarget()->y() == testEdge->getSource()->y()){
                return true;
            }
        }
    }
    return false;
}
