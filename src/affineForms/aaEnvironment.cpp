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
#include "aaEnvironment.hpp"

using namespace PolynomialForms;

bool AAEnvironment::isIndependent(int i, int j) const {
  BOOST_ASSERT( i >= 0 && i < nNoiseSymbols);
  BOOST_ASSERT( j >= 0 && j < nNoiseSymbols);
  map<int, int_list>::const_iterator vi, vj;
  vi = deps.find(i);
  vj = deps.find(j);
  BOOST_ASSERT(vi != deps.end());
  BOOST_ASSERT(vj != deps.end());
  
  int_list const & dI = vi -> second;
  int_list const & dJ = vj -> second;
  bool f1 = (dI.find(j) == dI.end());
  bool f2 = (dJ.find(i) == dJ.end());
  return (f1 && f2);
}


int AAEnvironment::makeFreshVariable(i_double const & range, i_double const & e, i_double const & m){
  int retID = this -> nNoiseSymbols;
  this -> nNoiseSymbols ++;
  int_list iList;
  iList.insert(retID);

  deps[retID] = iList;
  noise_range[retID] = range;
  noise_expect[retID] = e;
  noise_moment[retID] = m;
  return retID;
  
}

i_double AAEnvironment::getProductExpectation(int k1, int k2){
  BOOST_ASSERT( k1 >= 0 && k1 < this -> nNoiseSymbols && k2 >= 0 && k2 < this -> nNoiseSymbols);
  if (isIndependent(k1, k2)){
    return getExpect(k1) * getExpect(k2);
  } else {
    // Cut and paste from makeFreshProductVariable
    i_double cauchy_schwarz_interval0 = sqrt( noise_moment[k1] * noise_moment[k2] );
    i_double cauchy_schwarz_interval( -cauchy_schwarz_interval0.upper(), cauchy_schwarz_interval0.upper()); //hull(cauchy_schwarz_interval0, cauchy_schwarz_interval1);
    i_double rng = noise_range[k1] * noise_range[k2];
    i_double JJ =  intersect( rng, cauchy_schwarz_interval);
  
    return JJ;
  }
  
}

int AAEnvironment::makeFreshDependentVariable(i_double const & range, set<int> const & depSymbols){
  int retID = this-> nNoiseSymbols;
  this -> nNoiseSymbols++;
  // Make a list of dependent variables
  int_list iList;
  iList.insert(retID);
  BOOST_FOREACH( int var, depSymbols){
    set<int> const & dVars = deps[var];
    iList.insert(var);
    iList.insert(dVars.begin(), dVars.end());
  }
  deps[retID] = iList;
  noise_range[retID] = range;
  noise_expect[retID] = range;
  noise_moment[retID] = square(range);
  return retID;
  
  
}


int AAEnvironment::makeFreshProductVariable(int k1, int k2){
  BOOST_ASSERT( k1 >= 0 && k1 < this -> nNoiseSymbols && k2 >= 0 && k2 < this -> nNoiseSymbols);
  if (k1 > k2){
    // swap
    int tmp;
    tmp = k1;
    k1 = k2;
    k2 = tmp;
  }

  BOOST_ASSERT(k1 <= k2);

  // 1. Check if the symbol already exists.
  pair<int,int> pp(k1, k2);
  map< pair<int,int>, int>::const_iterator vi = noise_product_map.find(pp);
  if ( vi != noise_product_map.end()){
    return vi -> second;
  }
   
  // 2. Otherwise, create a fresh symbol
  
  int retID = this -> nNoiseSymbols; // This is the id of the new symbol
  this -> nNoiseSymbols ++; // Increment the number of symbols

  // create dependencies for this new symbol
  int_list iList; // Create the dependencies for the new symbols 
  iList.insert(retID); //reflexive dependence
  iList.insert(k1); // dependence on k1, k2
  iList.insert(k2);
  iList.insert(deps[k1].begin(), deps[k1].end()); // transitively close the dependences
  iList.insert(deps[k2].begin(), deps[k2].end());
  deps[retID] = iList; //set the dependencies for retID
  
  if (k1 != k2){ // Are the symbols different?
    // Noise range is given by the product of noise ranges.
    noise_range[retID] = noise_range[k1] * noise_range[k2];
    
    if (isIndependent(k1, k2)){ // Are they independent?
      noise_expect[retID] = noise_expect[k1] * noise_expect[k2];
      noise_moment[retID] = noise_moment[k1] * noise_moment[k2];
    } else { // Use Cauchy Schwarz and interval ranges both in the calculation 
      i_double cauchy_schwarz_interval0 = sqrt( noise_moment[k1] * noise_moment[k2] );
      i_double cauchy_schwarz_interval( -cauchy_schwarz_interval0.upper(), cauchy_schwarz_interval0.upper()); //hull(cauchy_schwarz_interval0, cauchy_schwarz_interval1);
     
      i_double JJ =  intersect( noise_range[retID], cauchy_schwarz_interval);
      assert(JJ.lower() <= JJ.upper());
      noise_expect[retID] = JJ;
      noise_moment[retID] = square(noise_range[retID]);
    }
    
  } else { // Are the symbols the same?
    noise_range[retID] = square(noise_range[k1]);
    noise_expect[retID] = noise_moment[k1];
    
    // Q: can we apply Cauchy Schwarz Inequality here?
    noise_moment[retID] = square(noise_range[retID]); 
  }
  // Add the pair back to the dictionary
  noise_product_map[pp] = retID;
  
  return retID;
}

void AAEnvironment::computeDependencyGraph( Graph & gRet, set<int> const & relevantNodes) const {
  BOOST_FOREACH(int k, relevantNodes){
    add_vertex(gRet);
  }
  
  BOOST_FOREACH(int k, relevantNodes){
    map<int, int_list>::const_iterator vk = deps.find(k);
    
    BOOST_ASSERT( vk != deps.end());
    int_list const & dk = vk -> second;
    
    BOOST_FOREACH(int j, relevantNodes){
      if (j >= k) continue;
      map<int, int_list>::const_iterator vj = deps.find(j);
      BOOST_ASSERT( vj != deps.end());
      int_list const & dj = vj -> second;
      
      if (setsIntersect(dk,dj)) {
	add_edge(k,j,gRet);
	
      }
    }
  }


  return;
}


ostream & operator<< (ostream & os, AAEnvironment const & aa){
  os << "Start Environment ID: " << aa.id << endl;
  // First print all the variables.
  os << "Number of noise symbols  =  " << aa.nNoiseSymbols << endl;
  for (int i = 0; i < aa.nNoiseSymbols ; ++i){
    os << "\t Variable # " << i << endl;
    os << "\t\t Range: " << aa.noise_range.at(i) << endl;
    os << "\t\t Expectation: "<< aa.noise_expect.at(i) << endl;
    os << "\t\t Moments : " << aa.noise_moment.at(i) << endl;
    os << "\t\t Dependencies: { " ;
    map<int, int_list>::const_iterator vi = aa.deps.find(i);
    BOOST_ASSERT( vi != aa.deps.end());
    BOOST_FOREACH(int j, vi -> second){
      os << j << " ";
    }
    os << "}" << endl;
  }
  os << "Ending Environment " << endl;
  return os;
}


