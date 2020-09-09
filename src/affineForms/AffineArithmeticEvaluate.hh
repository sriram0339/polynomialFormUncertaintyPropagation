//
// Created by Sriram Sankaranarayanan on 6/3/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_AFFINEARITHMETICEVALUATE_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_AFFINEARITHMETICEVALUATE_HH

#include "aaEnvironment.hpp"
#include "affineArithmeticClass.hpp"
#include "ExprAST.hh"
#include "SystemDescription.hh"

using namespace PolynomialForms;

typedef std::map<int, AffineArithmeticClass> stateMap;

class AAEvaluateExpr {
private:
    AAEnvironment & env;
    stateMap const & stMap;

public:
    AAEvaluateExpr(AAEnvironment & _env, stateMap const & origMap):env(_env), stMap(origMap){};

    AffineArithmeticClass evaluate(ExprPtr what);

    stateMap evaluateUpdates(std::map<int, ExprPtr> const & updates);

};

stateMap  computeInitialMap(AAEnvironment & env, StochasticSystemPtr const & what);
stateMap  computeOneStep(AAEnvironment & env, StochasticSystemPtr const & what, stateMap const & oldMap);

void probabilityQueryEvaluateAA(ExprPtr query, AAEnvironment & env, stateMap const & what);
void expectationQueryEvaluateAA(ExprPtr query, AAEnvironment & env, stateMap const & map);


#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_AFFINEARITHMETICEVALUATE_HH
