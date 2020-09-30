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
#include "affineArithmeticClass.hpp"
#include <set>
#include <fstream>
#include <cmath>

int trigonometricFunctionOrder = 2;
void AffineArithmeticClass::setMap(interval_map const & newMap ){

  m.clear();
  // Just to be safe, use the set function
  BOOST_FOREACH(interval_map::value_type const & p, newMap){
    i_double I = p.second;
    if (I.upper() == 0.0 && I.lower() == 0.0){
      // Do not add
    } else {
      this -> setCoefficient(p.first, p.second);
    }
  }
  
}

i_double AffineArithmeticClass::get(int j) const {
  BOOST_ASSERT( j >= 0 || j == constant_idx);
  BOOST_ASSERT( j < env.getNumSymbols());
  interval_map::const_iterator vj = m.find(j);
  if (vj == m.end()){
    return i_double(0.0, 0.0);
  } else {
    return vj -> second;
  }
  BOOST_ASSERT(false); // This should never be reached. 
}


void AffineArithmeticClass::setCoefficient(int j, i_double const & I){
  BOOST_ASSERT(j >= 0 || j == constant_idx);
  BOOST_ASSERT( j < env.getNumSymbols());
  if ( I.upper() == 0.0 && I.lower() == 0.0){
    m.erase(j); // erase the binding rather than explicitly store a zero.
  } else {
    m[j] = I;
  }
}

void AffineArithmeticClass::copy(AffineArithmeticClass const & from)
{
  // We will have to assume here that the environments are the same.
  BOOST_ASSERT( this -> env.identity() == from.env.identity());
  // Clear out this map
  this -> clear();
  // Next, iterate over the map for from and add entries to my own map
  BOOST_FOREACH(const interval_map::value_type & myPair, m){
    int key = myPair.first;
    i_double I = myPair.second;
    this -> setCoefficient(key, I);
  }
  
}


 
void AffineArithmeticClass::add_assign(AffineArithmeticClass const & op2){
  BOOST_ASSERT(this -> env.identity() == op2.env.identity());
  // Iterate over the map
  set<int> toErase;
  
  BOOST_FOREACH(const interval_map::value_type & myPair, m){
    int key = myPair.first;
    i_double I = myPair.second;
    i_double J = op2. get(key);
    i_double K = I + J;
    if (K.upper() == 0.0 && K.lower() == 0.0){
      toErase.insert(key);
    } else {
      this -> setCoefficient(key,K);
    }
  }

  BOOST_FOREACH(int k, toErase){
    this -> m.erase(k);
  }
  
  
  BOOST_FOREACH(const interval_map::value_type & p, op2.m){
    int key2 = p.first;
    i_double I = p.second;
    if (m.find(key2) == m.end()){
      this -> setCoefficient(key2, I);
    }
  }
  
}


void AffineArithmeticClass::sub_assign(AffineArithmeticClass const & op2){
  BOOST_ASSERT(this -> env.identity() == op2.env.identity());
  // Iterate over the map
  set<int> toErase;
  BOOST_FOREACH(const interval_map::value_type & myPair, m){
    int key = myPair.first;
    i_double I = myPair.second;
    i_double J = op2. get(key);
    i_double K = I - J;
    if (K.upper() == 0.0 && K.lower() == 0.0){
      toErase.insert(key);
    } else {
      this -> setCoefficient(key,K);
    }
  }

  BOOST_FOREACH(int k, toErase){
    this -> m.erase(k);
  }
  
  BOOST_FOREACH(const interval_map::value_type & p, op2.m){
    int key2 = p.first;
    i_double I = p.second;
    if (m.find(key2) == m.end()){
      this -> setCoefficient(key2, -I);
    }
  }
    
    
}


void AffineArithmeticClass::scale_assign(i_double const & op){
  // Iterate over the map
  if (op.upper() <= 0.0 && op.lower() >= 0.0){
    clear();
  } else {
    BOOST_FOREACH(const interval_map::value_type & myPair, m){
      int key = myPair.first;
      i_double I = myPair.second;
      i_double J = op * I;
      this -> setCoefficient(key, J);
    }
  }
}
void AffineArithmeticClass::multiply_assign(AffineArithmeticClass const & op2){
  // 1. Make forms for the linear part
  // 2. Collect nonlinear terms and their coefficients

  interval_map newMap; // Make a new map for our result.
  i_double c1 = this -> getConstant();
  i_double c2 = op2.getConstant();
  // Set the constant term for the new map
  newMap[constant_idx] = c1 * c2;
  // Now compute the linear terms
  BOOST_FOREACH(interval_map::value_type & p, m){
    int key = p.first;
    i_double val = p.second;
    if (key != constant_idx){
      i_double nVal = val * c2;
      newMap[key] = nVal;
    }
  }

  BOOST_FOREACH(interval_map::value_type const & q, op2.m){
    int key = q.first;
    i_double val = q.second;
    if (key != constant_idx){
      i_double nVal = val * c1;
      interval_map::const_iterator vi = newMap.find(key);
      if (vi == newMap.end()){
	newMap[key] = nVal;
      } else {
	i_double I = vi -> second;
	newMap[key] = I + nVal;
      }
    }
  }

  // Now newMap contains the linear part of the multiplication affine form.
  // Next, we will have to form the special variables.
  std::set< std::pair<int,int> > keyPairs;
  BOOST_FOREACH(interval_map::value_type & p, m){
    interval_map::const_iterator vi, vj;
    int key1 = p.first;
    if (key1 != constant_idx){
      BOOST_FOREACH(interval_map::value_type const & q, op2.m){
	int key2 = q.first;
	if (key2 != constant_idx){
	  pair<int,int> p1(key1, key2);
	  pair<int,int> p2(key2, key1);
	  if (keyPairs.find(p1) == keyPairs.end() && keyPairs.find(p2) == keyPairs.end()){
	    keyPairs.insert(p2);
	    keyPairs.insert(p1); 
	    int k = env.makeFreshProductVariable(key2, key1);
	    i_double nVal = p.second * q.second;
	    if (key1 != key2){
	      vi = m.find(key2);
	      if (vi != m.end()){
		i_double J = vi -> second;
		vj = op2.m.find(key1);
		if (vj != op2.m.end()){
		  i_double K = vj -> second;
		  nVal = nVal + J * K;
		}
	      }
	    }
	    interval_map::const_iterator vk = newMap.find(k);
	    if (vk == newMap.end()){
	      newMap[k] = nVal;
	    } else {
	      newMap[k] = vk -> second + nVal;
	    }
	  }
	
	}
      }
    }
  }
  setMap(newMap);
}

void AffineArithmeticClass::square_assign() {
  interval_map newMap;
  i_double c = this -> getConstant();
  // set the constant term
  newMap[constant_idx] = square(c);
  // now for the linear terms
  if (c.upper() != 0.0 || c.lower() != 0.0){
    BOOST_FOREACH(interval_map::value_type & p, m){
      // Iterate through each term
      int key = p.first;
      i_double val = p.second;
      if (key != constant_idx){
	i_double nVal = val * c;
	newMap[key] = i_double(2.0) * nVal;
      }
    }
  }
  
  // Now compute the nonlinear terms
  BOOST_FOREACH(interval_map::value_type & p, m){
    int k1 = p.first;
    i_double v1= p.second;
    interval_map::const_iterator vk;
    if (k1 != constant_idx){
      /* -- Take care of the square term --*/
      int sqTerm = env.makeFreshProductVariable(k1,k1);
      i_double sqVal = square(v1);
      vk = newMap.find(sqTerm);
      if (vk == newMap.end()){
	newMap[sqTerm] = sqVal;
      } else {
	newMap[sqTerm] = vk -> second + sqVal;
      }
      /* -- Done -- */
      
      BOOST_FOREACH(interval_map::value_type & q, m){
	int k2 = q.first;
	i_double v2 = q.second;
	if (k2 != constant_idx && k2 < k1) {
	  int nkey = env.makeFreshProductVariable(k1, k2);
	  i_double nVal =  v2 * v1;
	  nVal = i_double(2.0) * nVal;
	  vk = newMap.find(nkey);
	  if (vk == newMap.end()){
	    newMap[nkey] = nVal;
	  } else {
	    newMap[nkey] = vk -> second + nVal;
	  }
	}
      }
    }
  }

  setMap(newMap);
  
}

i_double AffineArithmeticClass::range() const{
  i_double retI(0.0);
  BOOST_FOREACH(interval_map::value_type const & p, m){
    int key = p.first;
    if (key == constant_idx){
      retI = retI + p.second;
    } else {
      i_double J = env.getRange(key);
      retI = retI + J * p.second;
    }
  }
  return retI;
}


i_double AffineArithmeticClass::expectation() const{
  i_double retI(0.0);
  BOOST_FOREACH(interval_map::value_type const & p, m){
    int key = p.first;
    if (key == constant_idx){
      retI = retI + p.second;
    } else {
      i_double J = env.getExpect(key);
      retI = retI + (J * p.second);
    }
  }
  return retI;
}


void AffineArithmeticClass::taylorFirstOrderFormAssign(AffineArithmeticClass const & h, i_double c0, i_double c1, i_double R, set<int> const & dependencies){
  // f(c+h) = c0 + h * c1 + freshVar( range(h^2/2! * R))
  this -> clear();
  add_assign(h);
  scale_assign(c1);
  i_double r = square(h.range());
  i_double rem = r * R;
  rem =  rem/i_double(2.0);
  int freshVar = env.makeFreshDependentVariable(rem,dependencies);
  this -> setCoefficient(freshVar, i_double(1.0));
  i_double c = this -> getConstant();
  c = c + c0;
  if (c.upper () != 0.0 || c.lower() != 0.0)
    this -> setConstant(c);
  return;
}



  
void AffineArithmeticClass::power_assign(AffineArithmeticClass const & op, int n){
  // (c+h)^n = c^n + n h (c^{n-1}) + freshVar( range(1/2 * h^2 * n *(n-1)* R^{n-2}))
  BOOST_ASSERT(n >= 2);
  if (n > 8) {
    set<int> dependencies;
    op. collectRelevantNoiseSymbols(dependencies);
    i_double I = op.range();
    i_double c ( median(I) );

    i_double deriv0 = pow(c,n);
    i_double deriv1 = i_double(n) * pow(c,n-1);
    i_double deriv2 = i_double(n*(n-1))* pow(I, n-2);

  
    // Make a copy of the operand
    AffineArithmeticClass h(op);
    // Subtract c from that copy
    i_double c0 = op.getConstant();
    h.setConstant(c0 - c);

    taylorFirstOrderFormAssign(h, deriv0, deriv1 ,deriv2,dependencies);
  } else {
    clear();
    this -> setConstant(1.0);
    AffineArithmeticClass tmp(op);
    while (n > 0){
      if (n%2 == 1){
	this -> multiply_assign(tmp);
      }
      n = n/2;
      if ( n > 0){
	tmp.square_assign();
      }
      std::cout << "n=" << n << "this = " << *this << " tmp =  " << tmp << std::endl;
      
    }
  }  
  
}


void AffineArithmeticClass::tan_assign(AffineArithmeticClass const & op){
  // Use the expansion
  //  tan(c+h) = tan(c) + h (1 + tan^2(c))  + freshVar( range(1/2 * h^2 * 2 * (tan(R) + tan^3(R))))

  // collect the dependencies for this noise symbol.
  set<int> dependencies;
  op. collectRelevantNoiseSymbols(dependencies);
  i_double I = op.range();
  i_double c ( median(I) );
  

  i_double tanc = tan(c);
  i_double sec2c = i_double(1.0) + square(tanc);
  i_double deriv2 = i_double(2.0)* tan(I) * ( i_double(1.0) + square(tan(I)));

  
  // Make a copy of the operand
  AffineArithmeticClass h(op);
  // Subtract c from that copy
  i_double c0 = op.getConstant();
  h.setConstant(c0 - c);

  taylorFirstOrderFormAssign(h, tanc, sec2c,deriv2,dependencies);
}

void AffineArithmeticClass::inverse_assign(AffineArithmeticClass const & op){
  // Use the expansion
  // 1/(c+h) = 1/c - h (1/c^2) + freshVar(1/2! * range(h^2 * R))
  // R = 2 * range( (c+h)^{-3})
  set<int> dependencies;
  op. collectRelevantNoiseSymbols(dependencies);
  i_double I = op.range();
  assert(! zero_in(I) ); 
  i_double c(median(I));
  AffineArithmeticClass h(op);
  i_double c0 = h.getConstant();
  h.setConstant(c0 - c);
  i_double c1 = pow(c,-1);
  i_double c2 = - square(c1);
  i_double R = i_double(2.0) * pow(I, -3);
  taylorFirstOrderFormAssign(h, c1, c2, R, dependencies);
}

void AffineArithmeticClass::specialFunctionForTumorModel_assign(AffineArithmeticClass const & op){
  // F(x) = (1 + 1/(1+x^2)) * x^2 = x^2 + (x^2/(1+x^2))
  // F'(x) = 2x + 2x/(1+x^2)^2
  // F''(x) = 2 + 2/(1+x^2)^2 - 8x^2/(1+x^2)^3 = 2 + 2(1-3x^2)/(1+x^2)^3
  
  set<int> dependencies;
  op. collectRelevantNoiseSymbols(dependencies);
  i_double I = op.range();
  i_double c(median(I));
  AffineArithmeticClass h(op);
  i_double c0 = h.getConstant();
  h.setConstant(c0 - c);
  i_double cSquare = square(c);
  i_double cSquarePlusOne = i_double(1.0) + cSquare;
  i_double deriv0 =  cSquare * ( i_double(1.0) + pow(cSquarePlusOne,-1));
  i_double deriv1 =  i_double(2) * c * ( i_double(1.0) + pow(cSquarePlusOne,-2));
  i_double Isq = square(I);
  i_double num = i_double(1.0) - i_double(3.0) * Isq;
  i_double den = i_double(1.0) + Isq;
  i_double R = i_double(2) + i_double(2.0) * num * pow(den,-3);
  taylorFirstOrderFormAssign(h, deriv0, deriv1, R, dependencies);
}


void AffineArithmeticClass::sine_assign(AffineArithmeticClass const & op){
  BOOST_ASSERT(this -> env.identity() == op.env.identity());

  // collect the dependencies for this noise symbol.
  set<int> dependencies;
  op .collectRelevantNoiseSymbols(dependencies);
 
  //1. Obtain center (c) and range R for op.
  //2. Use order 2 expansion.
  //   sin(c+h) = sin(c) + h cos(c) - freshVar(h^2/2! * sin(R))
  //   or order 3 expansion
  //   sin(c+h) = sin(c) + h cos(c) - h^2/2! * sin(c) - freshVar[h^3/3! * cos(R)]
 

  i_double I = op.range();
  i_double c ( median(I) );
 

  i_double sinC = sin(c);
  i_double cosC = cos(c);
  i_double sinR = sin(I);
  i_double cosR = cos(I);

  
  // Make a copy of the operand
  AffineArithmeticClass h(op);
  // Subtract c from that copy
  i_double c0 = op.getConstant();
  h.setConstant(c0 - c);

  
  this -> clear();
  // + h * cos(c)
  add_assign(h);
  scale_assign(cosC);


  
  if (trigonometricFunctionOrder == 2){
    i_double r = square(h.range());
    i_double rem = r * sinR;
    rem =  rem/i_double(2.0);
    int freshVar = env.makeFreshDependentVariable(rem,dependencies);
    this -> setCoefficient(freshVar, i_double(-1.0));
    // Done: hopefully, we got that correct!
  } else {
    assert(trigonometricFunctionOrder == 3);
    
    // tmp = h^2/2! * sin(c)
    AffineArithmeticClass tmp(h);
    tmp.square_assign();
    tmp.scale_assign(sinC);
    tmp.scale_assign( i_double( -0.5, -0.5));
    // - h^2/2! sin(c)
    this -> add_assign(tmp);
    // Now add the remainder term
    i_double r = pow(h.range(),3);
    i_double rem = r*cosR;
    rem = rem/i_double(6.0);
    int freshVar = env.makeFreshDependentVariable(rem,dependencies);
    this -> setCoefficient(freshVar, i_double(-1.0));
  }

  
  
  // Now add sin(c) to the constant term
  c0 = getConstant();
  c0 = c0 + sinC;
  setConstant(c0);

}


void AffineArithmeticClass::cosine_assign(AffineArithmeticClass const & op){
  BOOST_ASSERT(this -> env.identity() == op.env.identity());
  // collect the dependencies for this noise symbol.
  set<int> dependencies;
  op. collectRelevantNoiseSymbols(dependencies);
 
  
  //1. Obtain center (c) and range R for op.
  //2. Use order 2 expansion.
  //   cos(c+h) = cos(c) - h sin(c) - freshVariable(h^2/2! * cos(R))
  //3. Or Use order 3 expansion.
  //   cos(c+h) = cos(c) - h sin(c) - h^2/2! * cos(c) + freshVar(h^3/3! * sin(R))

  
  i_double I = op.range();
  i_double c ( median(I) );


  i_double sinC = sin(c);
  i_double cosC = cos(c);
  i_double sinR = sin(I);
  i_double cosR = cos(I);
  

  
  // Make a copy of the operand
  AffineArithmeticClass h(op);
  // Subtract c from that copy
  i_double c0 = op.getConstant();
  h.setConstant(c0 - c);

  
  this -> clear();
  // this = - h * sin(c)
  sub_assign(h);
  scale_assign(-sinC);

  if (trigonometricFunctionOrder == 2){
    i_double r = square(h.range());
    i_double rem = r * cosR;
    rem = rem/i_double(2.0); // rem := h^2/2 * cosR
    int freshVar = env.makeFreshDependentVariable(rem,dependencies);
    this -> setCoefficient(freshVar, i_double(-1));
  } else {
     assert(trigonometricFunctionOrder == 3);
    
    // tmp = h^2/2! * cos(c)
    AffineArithmeticClass tmp(h);
    tmp.square_assign();
    tmp.scale_assign(cosC);
    tmp.scale_assign( i_double( -0.5, -0.5));
    
    // - h^2/2! cos(c)
    this -> add_assign(tmp);
    
    // Now add the remainder term
    i_double r = pow(h.range(),3);
    i_double rem = r*sinR;
    rem = rem/i_double(6.0); // = h^3/3! * sinR
    int freshVar = env.makeFreshDependentVariable(rem,dependencies);
    this -> setCoefficient(freshVar, i_double(1.0));
  }

  
  // Now add cos(c) to the constant term
  c0 = getConstant();
  c0 = c0 + cosC;
  setConstant(c0);

  // Done: hopefully, we got that correct!
}

void AffineArithmeticClass::collectRelevantNoiseSymbols(set<int> & returnValue) const {
  BOOST_FOREACH(const interval_map::value_type & p, m){
    int key = p.first;
    if ( key == constant_idx) continue;
    else {
      i_double val = p.second;
      BOOST_ASSERT(val.upper() != 0.0 || val.lower() != 0.0);
      returnValue.insert(key);
    }
  }
}


double AffineArithmeticClass::getDenominatorBoundForConcMeasure() const {
  double retVal = 0.0;
  BOOST_FOREACH( const interval_map::value_type & myPair, m){
    int key = myPair.first;
      i_double c = myPair.second;
      if (key != constant_idx){
	i_double r = env.getRange(key);
	c = c * r;
      }
      retVal += (c.upper() - c.lower()) * (c.upper() - c.lower()) ;
      if (debug){
	cout << " y_"<<key << ": " << (c.upper() - c.lower()) * (c.upper() - c.lower()) << endl;
      }
  }
  return retVal;
}

i_double AffineArithmeticClass::secondMoment_raw() const {
  i_double e = this-> expectation();
  i_double moment2(0.0);
  // To compute second moment, we will work 
  BOOST_FOREACH(interval_map::value_type const & p, this -> m){
    int k1 = p.first;
    i_double v1 = p.second;
    if ( k1 == constant_idx){
      moment2 = moment2 + square(v1);
      moment2 = moment2 + i_double(2.0) * v1*(e-v1);
    } else {
      i_double v = env.getMoment(k1);
      moment2 = moment2 + square(v1)*v;
      BOOST_FOREACH(interval_map::value_type const & q, this -> m){
	int k2 = q.first;
	i_double v2 = q.second;
	if (k2 != constant_idx && k2 < k1){
	  moment2 = moment2 + i_double(2.0) *  v1 * v2 * env.getProductExpectation(k1,k2);
	} 
      }

    }

  }
  
  return moment2;

}

i_double AffineArithmeticClass::variance_raw() const {
  i_double e = this-> expectation();
  i_double moment2 = this -> secondMoment_raw();
  i_double var = moment2 - square(e);
  if (var.lower() < 0.0){
    double ub = var.upper();
    var = i_double(0.0, ub);
  }
  return var;
}


void AffineArithmeticClass::linearized_multiply_assign(AffineArithmeticClass const & op2){
  // 1. Make forms for the linear part
  // 2. Collect nonlinear terms and their coefficients

  interval_map newMap; // Make a new map for our result.
  i_double c1 = this -> getConstant();
  i_double c2 = op2.getConstant();
  // Set the constant term for the new map
  newMap[constant_idx] = c1 * c2;
  // Now compute the linear terms
  BOOST_FOREACH(interval_map::value_type & p, m){
    int key = p.first;
    i_double val = p.second;
    if (key != constant_idx){
      i_double nVal = val * c2;
      newMap[key] = nVal;
    }
  }

  BOOST_FOREACH(interval_map::value_type const & q, op2.m){
    int key = q.first;
    i_double val = q.second;
    if (key != constant_idx){
      i_double nVal = val * c1;
      interval_map::const_iterator vi = newMap.find(key);
      if (vi == newMap.end()){
	newMap[key] = nVal;
      } else {
	i_double I = vi -> second;
	newMap[key] = I + nVal;
      }
    }
  }

  // Now newMap contains the linear part of the multiplication affine form.
  // Collect the nonlinear "mess"
  i_double nRange(0.0);
  i_double nExpect(0.0); // Let us not worry about the variance yet.
  set<int> depVars;
  BOOST_FOREACH(interval_map::value_type & p, m){
    interval_map::const_iterator vi, vj;
    int key1 = p.first;
    if (key1 != constant_idx){
      depVars.insert(key1);
      BOOST_FOREACH(interval_map::value_type const & q, op2.m){
	int key2 = q.first;
	if (key2 != constant_idx){
	  if (key1 != key2){
	    i_double nVal = p.second * q.second;
	    i_double pRange = env.getRange(key1) * env.getRange(key2);
	    i_double pExpect= env.getProductExpectation(key1,key2);
	    nRange = nRange + pRange * nVal;
	    nExpect = nExpect + pExpect * nVal;
	  } else {
	    i_double nVal = square(p.second);
	    i_double pRange = square(env.getRange(key1));
	    i_double pExpect = env.getProductExpectation(key1,key2);
	     nRange = nRange + pRange * nVal;
	    nExpect = nExpect + pExpect * nVal;
	  }
	  depVars.insert(key2);
	}
	
	
      }
    }
  }

  int yNew = env.makeFreshDependentVariable(nRange,depVars);
  env.setMoments(yNew, nExpect, square(nRange) );
  newMap[yNew] = i_double(1.0);
  setMap(newMap);

  
}

void AffineArithmeticClass::clusterNoiseSymbolsForChernoffHoeffding(Graph const & depGraph,  double & cSquared, i_double & variance, double & M) const {
  // 1. Compute the connected components of the dependency graph.
  // 2. For each connected component id,
  //    2.0 Compute a smaller affine form.
  //    2.1 Compute range, mean and variance for the terms
  //    2.3 Compute max deviation from mean
  int n = num_vertices(depGraph);
  std::vector<int> depComponents(n);
  int i,j;
  int nComps = getComponentsOfGraph(depGraph, depComponents);
  // Collect the set of things in each component
 
  cout << "Number of components in $\\widehat{G}$ is " << nComps << "."<<endl;
  
  std::map< int, std::set<int> > componentMap;
  
  for (i = 0; i < nComps; ++i)
    componentMap[i] = set<int>();
  
  for (j = 0; j < n; ++j){
    int k = depComponents[j];
    componentMap[k].insert(j);
  }
  // Print all the components
  if (debug){
    cout << "--- Debug message --- "<< endl;
    for (i = 0; i < nComps; ++i){
      set<int> const & s = componentMap[i];
      cout << "Component # " << i << " = { ";
      BOOST_FOREACH( int e, s){
	cout << e << ", ";
      }
      cout << "}" << endl;
    }
    cout << "--- End debug message --- "<< endl;
  }
  BOOST_ASSERT(nComps >= 1);
  cSquared = 0.0;
  variance =0.0;
  M = -1.0; // Initially, we will start this way
  
  for (i = 0; i < nComps; ++i){
    // Make a new affine form
    AffineArithmeticClass a (env);
    set<int> const & s = componentMap[i];
    BOOST_FOREACH(int e, s){
      interval_map::const_iterator vi = m.find(e);
      if (vi != m.end()){
	a.setCoefficient(vi-> first, vi -> second);
      }
    }
    // Affine form a created.
    // 1. find it's range and variance
    i_double aRange = a.range();
    cSquared += (aRange.upper() - aRange.lower()) * (aRange.upper() - aRange.lower() );
    i_double aVariance = a.variance_raw();
    variance = variance + aVariance;
    i_double e = a.expectation();
    i_double deviation = aRange - e;
    M = MAX(M, ( MAX( fabs(deviation.lower()), fabs(deviation.upper()))));
    
  }

  return;

}

ostream & operator<< (ostream & os, AffineArithmeticClass const & a){
  interval_map m = a.m;
  i_double c = a.getConstant();
  os << c ;
  BOOST_FOREACH(const interval_map::value_type & myPair, a.m){
    if (myPair.first != constant_idx) {
      os << " + "<< myPair.second << "*";
      os << a.env.varName(myPair.first);
    }
  }
  os << "  ";

  return os;

}


void AffineArithmeticClass::normalizeForm_assign() {
  i_double res(0.0);
  set<int> toErase;
  
  BOOST_FOREACH(interval_map::value_type & myPair, this -> m){
    int k = myPair.first;
    if (k != constant_idx){
      i_double v = myPair.second;
      i_double m = median(v);
      if (m.upper() != 0.0 || m.lower() != 0.0){
	this -> setCoefficient(k,m);
      } else {
	toErase.insert(k);
      }
      res = res + (v-m)*env.getRange(k);
    } 
  }
  // Erase the keys to be set to zero
  BOOST_FOREACH(int k, toErase){
    this -> m.erase(k);
  }

  i_double c = this -> getConstant();
  c = c + res;
  this -> setConstant(c);
  return;

}

void AffineArithmeticClass::performAnalysisForCMI(ostream & os) const {
  /* -- This is the main class to perform an analysis for concentration of measure inequalities.
     We will focus on two for now:
     (a) Chromatic number based bound by Jansen 2004.
     (b) Cluster dependencies and apply Chernoff Hoeffding.
     
     Possible extension to consider would be to add variance bound.
     ---*/

  set<int> relevantSymbols;
  this -> collectRelevantNoiseSymbols(relevantSymbols);
  Graph depGraph;
  this -> env.computeDependencyGraph(depGraph, relevantSymbols);

  if (debug){
    std::ofstream of;
    cout << "Dependency graph created and dumped to tmp.dot " << endl;
    of.open("tmpGraph.dot");
    outputGraphToDot(depGraph, of);
    of.close();
  }
  
  os << "\\paragraph{Chromatic number based bound} " << endl;
  
  double c0 = this -> getDenominatorBoundForConcMeasure();
  int del = findMaximumDegree(depGraph);
  i_double eI = this -> expectation();
  double eLB = eI.lower(), eUB = eI.upper();
  double xG = (del+1.0)/2.0 * c0;
  cout << "Chromatic number bound on $\\widehat{G}$ is " << del+1 <<"."<< endl; 

  
  os << " Upper tail bound ";
  os << " $ (t \\geq "<< eUB <<"): $ " << endl;
  
  os << " \\[ \\mathbb{P} ( X \\geq t ) \\leq \\mbox{exp}\\left(  \\frac{ " ;
  os << " - ( t - "<< eUB<< ")^2}{" << xG;
  os << "}\\right) \\]" << endl;

  

  os << " Lower tail bound ";
  os << " $ (t \\leq "<< eLB <<"): $ " << endl;
  
  os << " \\[ \\mathbb{P} ( X \\leq t) \\leq \\mbox{exp}\\left( \\frac{" ;
  os << " - ("<< eLB<< " -t)^2}{" << xG;
  os << "}\\right) \\]" << endl;

  applySubGaussianCMIInequalitiesOverRange(os, this -> range(), eI, xG,250);
  
  os << "\\paragraph{Clustered dependencies and Chernoff-Hoeffding}" << endl;

  double c1; i_double variance;
  double M;
  this -> clusterNoiseSymbolsForChernoffHoeffding(depGraph,c1,variance, M);
  c1 = c1/2.0;

  
  
  os << " Upper tail bound ";
  os << " $ (t \\geq "<< eUB <<"): $ " << endl;
  os << " \\[ \\mathbb{P} ( X \\geq t ) \\leq \\mbox{exp}\\left(  \\frac{ " ;
  os << " - ( t - "<< eUB<< ")^2}{" << c1;
  os << "}\\right) \\]" << endl;


  os << " Lower tail bound ";
  os << " $ (t \\leq "<< eLB <<"): $ " << endl;
  os << " \\[ \\mathbb{P} ( X \\leq t) \\leq \\mbox{exp}\\left( \\frac{" ;
  os << " - ("<< eLB<< " -t)^2}{" << c1;
  os << "}\\right) \\]" << endl;

  

  applySubGaussianCMIInequalitiesOverRange(os, this -> range(), eI, c1,250);
  
  os << "\\paragraph{Bound based on Bernstein inequality} " << endl;
  os << "We obtain $\\sigma^2$ = " << variance << endl;
  os << "$M$ = " << M << endl;

  

  os << " Upper tail bound ";
  os << " $ (t \\geq "<< eUB <<"): $ " << endl;
   os << " \\[ \\mathbb{P} ( X \\geq t ) \\leq \\mbox{exp}\\left(  \\frac{ " ;
  os << " - ( t - "<< eUB<< ")^2}{";
  os << 2*variance.upper() << " + " << 0.6666 * M << " t "; 
  os << "}\\right) \\]" << endl;


   os << " Lower tail bound ";
  os << " $ (t \\leq "<< eLB <<"): $ " << endl;
  os << " \\[ \\mathbb{P} ( X \\leq t) \\leq \\mbox{exp}\\left( \\frac{" ;
  os << " - ("<< eLB<< " -t)^2}{" ;
  os << 2*variance.upper() << " + " << 0.6666 * M << "* t "; 
  os << "}\\right) \\]" << endl;

  applyBernsteinCMIInequalityOverRange(os, this -> range(), eI, variance.upper(), M,250);
  return;
}

/*-- Overloaded Operators --*/


AffineArithmeticClass & AffineArithmeticClass::operator= (AffineArithmeticClass const & op){
  BOOST_ASSERT( this -> env.identity() == op.env.identity());
  this-> clear();
  this -> m = op.m;
  return (*this);
}


AffineArithmeticClass AffineArithmeticClass::operator- () const {
  AffineArithmeticClass ret(*this);
  ret.scale_assign(i_double(-1.0));
  return ret;
}

AffineArithmeticClass depRandomVariable(AAEnvironment & env,set<int> const & vars, i_double rng ,i_double expect ,i_double secondMoment){
  AffineArithmeticClass ret(env);
  int yNew = env.makeFreshDependentVariable(rng,vars);
  env.setMoments(yNew,expect,secondMoment);
  ret.setCoefficient(yNew,i_double(1.0));
  return ret;
}

AffineArithmeticClass randomVariable(AAEnvironment & env, i_double range, i_double mean, i_double moment2){
  AffineArithmeticClass ret(env);
  int yNew = env.makeFreshVariable(range, mean, (moment2));
  ret.setCoefficient(yNew, i_double(1.0));
  return ret;
}

AffineArithmeticClass operator+ (AffineArithmeticClass const & a, AffineArithmeticClass const & b){
  BOOST_ASSERT( a.env.identity() == b.env.identity());
  AffineArithmeticClass ret(a);
  ret.add_assign(b);
  return ret;
}


AffineArithmeticClass operator+ (AffineArithmeticClass const & a, i_double b){
  AffineArithmeticClass ret(a);
  i_double c = ret.getConstant();
  c = c + b;
  ret.setConstant (c);
  return ret;
}

AffineArithmeticClass operator+ (i_double b, AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a);
  i_double c = ret.getConstant();
  c = c + b;
  ret.setConstant (c);
  return ret;
}



AffineArithmeticClass operator- (AffineArithmeticClass const & a, AffineArithmeticClass const & b){
  BOOST_ASSERT( a.env.identity() == b.env.identity());
  AffineArithmeticClass ret(a);
  ret.sub_assign(b);
  return ret;
}


AffineArithmeticClass operator- (AffineArithmeticClass const & a, i_double b){
  AffineArithmeticClass ret(a);
  i_double c = ret.getConstant();
  c = c - b;
  ret.setConstant (c);
  return ret;
}

AffineArithmeticClass operator- (i_double b, AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a);
  ret.scale_assign(i_double(-1.0));
  i_double c = ret.getConstant();
  c = b + c;
  ret.setConstant (c);
  return ret;
}




AffineArithmeticClass operator* (AffineArithmeticClass const & a, AffineArithmeticClass const & b){
  BOOST_ASSERT( a.env.identity() == b.env.identity());
  AffineArithmeticClass ret(a);
  ret.multiply_assign(b);
  return ret;
}


AffineArithmeticClass operator* (AffineArithmeticClass const & a, i_double b){
  AffineArithmeticClass ret(a);
  ret.scale_assign(b);
  return ret;
}

AffineArithmeticClass operator* (i_double b, AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a);
  ret.scale_assign(b);
  return ret;
}


AffineArithmeticClass sine(AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a.getEnv());
   ret.sine_assign(a);
   return ret;
}



AffineArithmeticClass cosine(AffineArithmeticClass const & a){
  AffineArithmeticClass ret( a.getEnv() );
   ret.cosine_assign(a);
   return ret;
}

AffineArithmeticClass square(AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a);
  ret.square_assign();
  return ret;
}
AffineArithmeticClass tan(AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a.getEnv());
  ret.tan_assign(a);
  return ret;
}
AffineArithmeticClass inverse(AffineArithmeticClass const & a){
  AffineArithmeticClass ret(a.getEnv());
  ret.inverse_assign(a);
  return ret;
}

AffineArithmeticClass power(AffineArithmeticClass const & a, int n){
  AffineArithmeticClass ret(a.getEnv());
  ret.power_assign(a,n);
  return ret;
}


