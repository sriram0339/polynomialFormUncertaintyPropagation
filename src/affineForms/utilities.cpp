/* 
        Copyright (C) 2016-2019  Sriram Sankaranarayanan srirams@colorado.edu

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <https://www.gnu.org/licenses/>.

 */
#include "utilities.hpp"

bool debug = false;

typedef graph_traits<Graph>::vertex_iterator vertex_iter;
typedef graph_traits<Graph>::edge_iterator edge_iter;
/*
ostream & operator << (ostream & os, i_double const & I) {
  os.precision(10);
  os << "["<<I.lower() <<","<<I.upper()<<"]";
  return os;
}
 */

bool setsIntersect(set<int> const &a, set<int> const & b){
  set<int>::const_iterator va = a.begin(), vb = b.begin();
  
  while (va != a.end() && vb != b.end()){
    while (va != a.end() && (*va < *vb))
      ++va;
    if (va != a.end() && vb!= b.end() && *va == *vb) return true;
    while (vb != b.end() && (*vb < *va))
      ++vb;
  }

  return false;
}

void outputGraphToDot(Graph const & g, ostream & os){
  os << "graph { " << endl;
  vertex_iter vi, viEnd;
  edge_iter ei, eiEnd;
  boost::tie(vi,viEnd) = vertices(g);
  boost::tie(ei,eiEnd) = edges(g);
  // Iterate through the vertices
  for(; vi != viEnd; ++vi){
    int k = *vi;
    os << k <<" [label= \"y"<<k<<"\", shape=\"circle\"];" << endl;
  }
  // Iterate through the edges
  for (; ei != eiEnd; ++ei){
    int k = source(*ei, g);
    int j = target(*ei, g);
    os << k << "--" << j <<";" << endl;
  }
  
  os << "}" << endl;
}

int findMaximumDegree(Graph const & g){
  int maxDegree =0;
  vertex_iter vi, viEnd;
  tie(vi,viEnd) = vertices(g);
  
  for (; vi != viEnd; ++vi){
    maxDegree = MAX(maxDegree, degree(*vi,g));
  }
  return maxDegree;
}

int  getComponentsOfGraph(Graph const & g, std::vector<int> & componentIDs){
  int n = num_vertices(g);
  assert( componentIDs.size() == n);
  return connected_components(g, &componentIDs[0]);      
}


void applySubGaussianCMIInequalitiesOverRange(ostream & os, i_double range, i_double expect, double cSquared, int nSubDivs){
  double span = (range.upper() - range.lower())/nSubDivs;
  os << "{\\tiny \\[\\begin{array}{c}"<<endl;
  for (int i = 0; i <= nSubDivs; ++i){
    double v = range.lower() + i* span;
    double t = 0.0;
   
    if (v < expect.lower()) {
      // Apply a lower tail bound
      t = expect.lower() -v;
      os << "\\mathbb{P}( X \\leq " << v << ") \\leq " ;
    } else if (v > expect.upper()){
      t = v - expect.upper();
      os << "\\mathbb{P}( X \\geq  " << v << ") \\leq " ;
    } else {
      t = 0.0;
      os <<"%%"<< v << " ~~ " ;
    }

    double b = exp ( - ( t*t/ cSquared) );
    os << b << "\\\\"<<endl;
  }
  os<< "\\end{array}\\]}"<<endl;
  
}

void applyBernsteinCMIInequalityOverRange(ostream & os, i_double range, i_double expect, double var, double M, int nSubDivs){
  assert (var > 0.0);
  assert (M > 0.0);
  os << "{\\tiny \\[\\begin{array}{c}"<<endl;
  double span = (range.upper() - range.lower())/nSubDivs;
  for (int i = 0; i <= nSubDivs; ++i){
    double v = range.lower() + i* span;
    double t = 0.0;
    
    if (v < expect.lower()) {
      // Apply a lower tail bound
      t = expect.lower() -v;
      os << "\\mathbb{P}( X \\leq " << v << ") \\leq " ;
    } else if (v > expect.upper()){
      t = v - expect.upper();
      os << "\\mathbb{P}( X \\geq " << v << ") \\leq " ;
    } else {
      t = 0.0;
      os <<"%%" <<v << " ~~ " ;
    }
    double denom = 2.0 * var + 2.0/3.0 * M * t;
    double b = exp( -(t*t/denom));
   
     os << b << "\\\\"<<endl;

  }
  os<< "\\end{array}\\]}"<<endl;
}
