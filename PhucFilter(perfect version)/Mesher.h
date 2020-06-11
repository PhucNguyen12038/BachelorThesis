/**
 * @file Mesher.h
 * @brief declares methods for building a surface mesh from points stored in an 
 * octree
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
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */
#ifndef MESHER_H
#define MESHER_H

#include <cstdlib>
#include <cstdio>



#include "types.h"
#include "Vertice.h"
#include "Vertice.h"
#include "Facet.h"
#include "Edge.h"
#include "utilities.h"
#include "types.h"
#include <common/interfaces.h>
#include <iostream>
#include <stdio.h>

using namespace std;

/**
 * @class Mesher
 * @brief Performs the triangulation of the input points
 * 
 * This class contains all methods to perform the Ball Pivoting triangulation 
 * of the input set of vertices.
 */
class Mesher
{
    protected : //class members

        /** @brief list of active edges (edge front)*/
        Edge_star_list m_edge_front;

        /** @brief list of created triangles*/
        Facet_star_list m_facets;

        /** @brief list of vertices*/
        Vertice_star_list m_vertices;

        /** @brief number of vertices*/
        unsigned int m_nvertices;

        /** @brief number of facets*/
        unsigned int m_nfacets;

        /** @brief number of edges*/
        unsigned int m_nEdges;

    public : //constructor-destructor

        /** @brief default constructor*/
        Mesher();

        /** @brief constructor
         * @param octree octree containing the points to mesh
         * @param iterator iterator over the octree
         */


        /** @brief destructor*/
        ~Mesher();

    public : //reconstruction methods

        /** @brief do the triangulation
         * @return true if at least a triangle was created
         */
        void reconstruct();

        /** @brief do the triangulation
         * @param radii a vector of radius
         */
        void reconstruct(std::list<double> &radii);

        /** @brief parallel triangulation
         * @param radii a vector of radius
         */
        void parallelReconstruct(std::list<double> &radii);

        /** @brief fill the triangular holes that remain due to wrong
         * normal orientation
         * (post-processing methods)
         */
        void fillHoles();
        void printFacet();

    public : //accessors+modifyers

        /** @brief get the number of vertices
         * @return number of mesh vertices
         */
        unsigned int nVertices() const;

        /** @brief get the number of facets
         * @return number of mesh facets
         */
        unsigned int nFacets() const;

        /** @brief get access to the mesh vertices
         * @return begin iterator of the vertices
         */
        Vertice_star_list::const_iterator vertices_begin() const;

        /** @brief get access to the mesh vertices
         * @return end iterator of the vertices
         */
        Vertice_star_list::const_iterator vertices_end() const;

        /** @brief get access to the mesh facets
         * @return begin iterator of the facets
         */
        Facet_star_list::const_iterator facets_begin() const;

        /** @brief get access to the mesh facets
         * @return end iterator of the facets
         */
        Facet_star_list::const_iterator facets_end() const;


    public : //auxilliary methods for performing the triangulation

        /** @brief find a seed triangle
         * @return true if a seed triangle was found
         */


        /** @brief add a facet to the list of mesher facets (calls addVertice)
         * @param f input facet
         */
        void addFacet(Facet *f);

        /** @brief add a Vertice to the list of mesher vertices if not
         * already added
         * @param v input Vertice
         */
        void addVertice(Vertice *v);

        /** @brief merge with another mesher (useful for safe-thread meshing)
         * @param mesher another mesher
         */
        void merge(Mesher &mesher);

        /** @brief getCommonEdge to get the share edge between 2 facets
          * @param two facets
          */
        Edge* getCommonEdge(Facet *f1, Facet *f2);


        /** @brief Internal to check vertice z inside circum circle triangle abc
         *
         * @param test vertice z and 3 vertices of trianle abc
         * return 0 is internal 1 is external
         */
        int Internal(Vertice *z, Vertice *a, Vertice *b,Vertice *c);

        /** @brief GenerateDelaunayTriangle to triangulate points cloud
         *
         * @param a set of vertex of mesh
         * @param mesher
         */

        void GenerateDelaunayTriangle(std::vector<Vertice*> VerticeVector, Mesher &mesher, Vertice* pointArray[400][400]);

        /** @brief getBorderEdgeList to get border edges of mesh
          *
          * @param a vector of edge
          * @param mesher
          */
        void getBorderEdgeList(std::vector<Edge*> &BorderList, Mesher &mesher);

        /** @brief setNumberOfFacet to set number of face for mesher
          *
          * @param number n
          *
          */
        void setNumberOfFacet(unsigned int n);

        /** @brief WooLine to find the vertex opposite side of triangle
          *
          * @param the edge e
          * @return list of vertex
          */
        std::list<Vertice*> WooLine(Edge* e, Vertice* pointArray);

        /** @brief FindNormalToEdge to find the normal vector direct opposite side of triangle
          *
          * @param the edge e
          * @param the double Xdirect
          * @param the double Ydirect
          *
          */
        void FindNormalToEdge(Edge* e, double &Xdirect, double &Ydirect);

        /** @brief FindVectorOfEdge to find the vector of 2 vertex
          *
          * @param the edge e
          * @param the double Xdirect
          * @param the double Ydirect
          *
          */
        void FindVectorOfEdge(Edge* e, double &XdirectE, double &YdirectE);

        /** @brief Angle to find the angle of 2 vector
          *
          * @param xv1
          * @param yv1
          * @param xv2
          * @param yv2
          *
          */
        int Angle(double xv1,double yv1, double xv2, double yv2);

        /** @brief mormalize to scale a vector to unit vector
          *
          * @param x
          * @param y
          *
          */
        void normalize(double &Xdirect,double &Ydirect);

};

#endif
