#ifndef D__AA_ENVIRONMENT__HPP__
#define D__AA_ENVIRONMENT__HPP__

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <boost/foreach.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/assert.hpp>
#include <string>
#include <sstream>
#include "utilities.hpp"

using namespace std;
using namespace boost;
using namespace boost::numeric;
using namespace interval_lib;

typedef std::set<int> int_list;

class AAEnvironment{
  
protected:

  int id; // An identifier for the environment -- kept  constant.
  int nNoiseSymbols; // The total number of noise symbol  
  map<int, int_list > deps;  // The dependencies between noise symbols. Keep transitively closed.
  map<int, i_double > noise_range; // The range for every noise symbol
  map<int, i_double > noise_expect; // The range for the expected value of every noise symbol
  map<int, i_double > noise_moment; // The range for the moment of every noise symbol
  map< pair<int,int>, int> noise_product_map; // If a noise symbol yi * yj is assigned to symbol yk, record it.


  i_double getValueFromMap(int k, map<int, i_double> const & m) const {
    BOOST_ASSERT( k >= 0 && k < nNoiseSymbols);
    map<int, i_double>::const_iterator vi = m.find(k);
    BOOST_ASSERT( vi != m.end());
    return vi -> second;
  }
  
public:
  AAEnvironment(int ident):id(ident), nNoiseSymbols(0) {}
  
  AAEnvironment(AAEnvironment const & a){
    // Calling the copy constructor of AAEnvironmentClass is verbotten!
    cerr << "Calling the copy constructor of the AAEnvironmentClass ist Verbotten!"<<endl;
    cerr << "Kripaya, is class ke \"copy constructor\" ko na bulayen! "<<endl;
    BOOST_ASSERT(false);
    
  }
  int getNumSymbols() const {
    return nNoiseSymbols;
  }

  int identity() const {
    return id;
  }
  string varName(int i) const {
    ostringstream ss;
    ss << "y_"<<i;
    return ss.str();
  }

  i_double getProductExpectation(int k1, int k2);
  
  bool isIndependent(int i, int j) const;
  int makeFreshVariable(i_double const & range, i_double const & e, i_double const & m);
  int makeFreshProductVariable(int k1, int k2);
  int makeFreshDependentVariable(i_double const & range, set<int> const & deps);
  void setMoments(int k, i_double e, i_double m2){
    BOOST_ASSERT( k >= 0 && k < this -> nNoiseSymbols);
    noise_expect[k] = e;
    noise_moment[k] = m2;
  }
  i_double getRange(int k) const {
    return getValueFromMap(k, noise_range);
  }

  i_double getExpect(int k) const {
    return getValueFromMap(k, noise_expect);
  }

  i_double getMoment(int k) const {
    return getValueFromMap(k, noise_moment);
  }

  void computeDependencyGraph( Graph & ret, set<int> const & relevantNodes) const;


  void computeIndependentClusters(map<int, set<int> > & res) const;


  friend ostream & operator<< (ostream & os, AAEnvironment const & aa);
  
  
};

#endif
