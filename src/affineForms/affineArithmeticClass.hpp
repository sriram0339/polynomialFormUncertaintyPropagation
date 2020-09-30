#ifndef D__AFFINE_ARITHMETIC_CLASS__
#define D__AFFINE_ARITHMETIC_CLASS__


#include <ostream>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>
#include "utilities.hpp"
#include "aaEnvironment.hpp"
using namespace std;
using namespace boost;


// Type definitions i_double has already been defined as interval<double>
typedef std::map<int, i_double > interval_map;
// ID for the constant term in the map
const int constant_idx = (-1);

class AffineArithmeticClass {

protected:
  AAEnvironment & env;
  // Contains the list of noise symbols,
  // their moments and dependencies.

  
  interval_map  m;
  // A map from integers to intervals over doubles.
  // Use the special symbol (-1) for constant


  void setMap(interval_map const & newMap);
  
public:
  // Default constructor
  AffineArithmeticClass(AAEnvironment & myEnv):env(myEnv){
  }

  // Getters and setters
  AAEnvironment & getEnv() const {
    return env;
  }
  
  i_double get (int j) const;
  
  void setCoefficient (int j, i_double const & I);

  void clear(){
    m.clear();
  }
  
  i_double  getConstant() const{
    return this -> get(constant_idx);
  }

  void setConstant( i_double const & I){
    if (I.upper() == 0.0 && I.lower() == 0.0){
      this -> m.erase(constant_idx);
    } else {
      this -> setCoefficient(constant_idx, I);
    }
  }
  

  i_double range() const;
  i_double expectation() const;
  
  // Compute the variance of this class.  This can be very expensive
  // O(n^2) in the size of the affine form. This is only called on
  // portions of an affine form given by clusters of dependent
  // variables.
  i_double secondMoment_raw() const;
  i_double variance_raw() const;

  void copy(AffineArithmeticClass const & from);
  void add_assign(AffineArithmeticClass const & op2);
  void sub_assign(AffineArithmeticClass const & op2);
  void scale_assign(i_double const & op);
  void multiply_assign(AffineArithmeticClass const & op2);
  void linearized_multiply_assign(AffineArithmeticClass const & op2);
  void square_assign(); 
  void sine_assign(AffineArithmeticClass const & op);
  void cosine_assign(AffineArithmeticClass const & op);
  void inverse_assign(AffineArithmeticClass const & op);
  void taylorFirstOrderFormAssign(AffineArithmeticClass const & h, i_double c0, i_double c1, i_double R, set<int> const & dependencies);
  void tan_assign(AffineArithmeticClass const & op);
  void specialFunctionForTumorModel_assign(AffineArithmeticClass const & op);
  void collectRelevantNoiseSymbols(std::set<int> & returnValue ) const;
  double getDenominatorBoundForConcMeasure() const;
  void clusterNoiseSymbolsForChernoffHoeffding(Graph const & depGraph, double & c1, i_double & variance, double & M) const ;
  void power_assign(AffineArithmeticClass const & op, int pow);
  void performAnalysisForCMI(ostream & os) const;
  void normalizeForm_assign();
  void compactFormClustered_assign(AffineArithmeticClass const& op);
  
  AffineArithmeticClass & operator= (AffineArithmeticClass const & op);
  AffineArithmeticClass operator- () const ;
  friend std::ostream& operator<<(std::ostream& os,  AffineArithmeticClass const & a);
  
  friend AffineArithmeticClass operator+ (AffineArithmeticClass const & a, AffineArithmeticClass const & b);
  friend AffineArithmeticClass operator+ (i_double a , AffineArithmeticClass const & b);
  friend AffineArithmeticClass operator+ (AffineArithmeticClass const & a , i_double b);
  friend AffineArithmeticClass operator- (AffineArithmeticClass const & a, AffineArithmeticClass const & b);
  friend AffineArithmeticClass operator- (i_double a , AffineArithmeticClass const & b);
  friend AffineArithmeticClass operator- (AffineArithmeticClass const & a , i_double b);
  friend AffineArithmeticClass operator* (AffineArithmeticClass const & a, AffineArithmeticClass const & b);
  friend AffineArithmeticClass operator* (i_double a , AffineArithmeticClass const & b);
  friend AffineArithmeticClass operator* (AffineArithmeticClass const & a , i_double b);

  
};

AffineArithmeticClass depRandomVariable(AAEnvironment & env,set<int> const & vars, i_double rng ,i_double expect ,i_double secondMoment);

AffineArithmeticClass randomVariable(AAEnvironment & env, i_double range, i_double mean, i_double moment2);
AffineArithmeticClass sine(AffineArithmeticClass const & a);
AffineArithmeticClass cosine(AffineArithmeticClass const & a);
AffineArithmeticClass square(AffineArithmeticClass const & a);
AffineArithmeticClass tan(AffineArithmeticClass const & a);
AffineArithmeticClass inverse(AffineArithmeticClass const & a);
AffineArithmeticClass power(AffineArithmeticClass const & a, int n);
  



#endif
