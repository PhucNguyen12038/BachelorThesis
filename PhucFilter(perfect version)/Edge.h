/**
 * @file Edge.h
 * @author Julie Digne julie.digne@liris.cnrs.fr
 * @date 2012-10-08
 * @brief declares an edge structure (linking two vertices)
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

#ifndef EDGE_H
#define EDGE_H

#include "Vertice.h"

class Facet;

/**
* @class Edge
* @brief Stores an edge with two vertices
* 
* An edge stores pointers to its two end-vertices and its two adjacent facets
* along with its type (inner, front or border)
*/
class Edge
{
  private : //properties

    /** @brief source Vertice*/
    Vertice *m_src;

    /** @brief target Vertice*/
    Vertice *m_tgt;

    /** @brief first adjacent facet*/
    Facet *m_facet1;

    /** @brief second adjacent facet*/
    Facet *m_facet2;


    /** @brief edge type (0: border, 1: front edge, 2: inner edge)*/
    int m_type;


    public :
    //constructor+destructor
    /** @brief Default constructor*/
    Edge();

   /** @brief Edge constructor from its endpoints
    @param src source Vertice
    @param tgt target Vertice
    */
   Edge(Vertice* src, Vertice* tgt);

   /** @brief Destructor*/
   ~Edge();

  public : //accessors+modifiers

   /** @brief access source Vertice
    @return source Vertice
    */
   Vertice* getSource() const;

   /** @brief access target Vertice
   @return target Vertice
   */
   Vertice* getTarget() const;

   /** @brief recompute orientation
    * to be called after adding a first facet to the edge
    */
   void updateOrientation();

   /** @brief access facet 1
    * @return facet 1
    */
   Facet* getFacet1() const;

   /** @brief access facet 2
    @return facet 2
    */
   Facet* getFacet2() const;

   /** @brief get opposite Vertice: the Vertice of the facet adjacent to
    * the edge and opposite to it.
    * This method only makes sense for border or front edges
    * @return opposite Vertice
    */
   Vertice* getOppositeVertice() const;

   /** @brief add a facet to an edge
    * Prerequisite the facet actually contains the edge
    * @param facet the facet to add
    * @return true if the facet was successfully added, false otherwise (if the edge already had two adjacent facets)
    */
   bool addAdjacentFacet(Facet* facet);
   
   /** @brief remove a facet adjacent to an edge
    * Prerequisite the facet actually contains the edge
    * @param facet the facet to remove
    * @return true if the facet was successfully removed, false otherwise
    */
   bool removeAdjacentFacet(Facet* facet);

   /** @brief test if a Vertice is a Vertice of the edge
    * @param Vertice test Vertice
    * @return true if the Vertice is either the source or target Vertice of the edge
    */
   bool hasVertice(Vertice *Vertice) const;

   /** @brief test if an edge is an inner edge or a front/border edge
    * @return 0 the edge is a front/border edge 1 if the edge is an inner edge
    */
   bool isInnerEdge() const;

   /** @brief return edge type
    * @return type of the edge (0,1,2)
    */
    int getType() const;

    /** @brief set edge type
     * @param type (0,1,2)
     */
    void setType(int type);

    /** @brief CompareWithEdge to compare if this edge is same as the test edge
     * @brief the order of source and target vertex is not considered
     * @param an test edge
     * @return true if 2 edges is same
     * @return false if 2 edges is different
     */
    bool CompareWithEdge(Edge* testEdge);

    /** @brief includeVertice to test if a vertice is a vertice of the edge
     * @param Vertice test Vertice
     * @return true if the Vertice is either the source or target Vertice of the edge
     */
    bool includeVertice(Vertice *vertice) const;




};

#endif
