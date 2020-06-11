/** @file  Facet.h
 * @brief Declaration of a facet object
 * @author Julie Digne
 * @date 2012-10-08
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


#ifndef FACET_H
#define FACET_H

#include "Point.h"
#include "Vertice.h"
#include "Edge.h"

/**
 * @class Facet
 * @brief stores a triangular facet
 *
 * A facet contains pointers to its three vertices and its
 * r-circumsphere center
 */
class Facet
{
    private :
        /** @brief vertices of the facet*/
        Vertice *m_Vertice[3];


    public : //constructor+destructor

        /** @brief constructor*/
        Facet();



        /** @brief constructor from a set of vertices
         * @param v1 first Vertice
         * @param v2 second Vertice
         * @param v3 third Vertice
         */
        Facet(Vertice* v1, Vertice* v2, Vertice* v3);

        /** @brief constructor from a set of vertices
         * @param v1 first Vertice
         * @param v2 second Vertice
         * @param v3 third Vertice
         * @param ball_center center of the empty interior
         * ball incident to the three vertices
         */
        Facet(Vertice* v1, Vertice* v2, Vertice* v3, Point &ball_center);

        /** @brief constructor from an edge and a Vertice
         * prerequisite edge has at most one adjacent facet
         * @param edge edge (2 vertices to create the facet)
         * @param Vertice third Vertice to create the facet
         */
        Facet(Edge *edge, Vertice* Vertice);

        /** @brief constructor from an edge and a Vertice
         * prerequisite edge has at most one adjacent facet
         * @param edge edge
         * @param Vertice Vertice to link to the edge
         * @param ball_center center of the empty interior ball
         * incident to the three vertices
         */
        Facet(Edge *edge, Vertice* Vertice, Point &ball_center);

        /** @brief destructor*/
        ~Facet();

    public : //accessors + modifiers

        /** @brief get facet Vertice
         * @param index of the Vertice in the facet
         * @return corresponding facet
         */
        Vertice* vertice(unsigned int index) const;

        /** @brief get edge opposite to Vertice i
         * @param index index of the edge
         * @return edge
         */
        Edge* edge(unsigned int index) const;


        /** @brief test if contains Vertice
         * @param Vertice test Vertice
         * @return true if the Vertice is a Vertice of the facet
         */
        bool hasVertice(Vertice *vertice);


};

#endif
