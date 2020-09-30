//
// Created by Sriram Sankaranarayanan on 2/13/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_MODELPARSER_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_MODELPARSER_HH

#include "DistributionInfo.hh"
#include "SystemDescription.hh"
#include "ExprAST.hh"

using namespace PolynomialForms;
extern StochasticSystem * globalSystem;
extern int lineNum;

void parserMain(const char * fileName);
StateAbstractionPtr computeNSteps(int n, int maxDegree);
void computeNStepsAndSaveData(int n, int maxDegree,std::vector<std::string> & varNameVec, std::vector<MpfiWrapper> & rangeVec);
void computeNSteps_inference(int n, int maxDegree,StateAbstractionPtr & st);
StateAbstractionPtr initialize_globalSystem( int maxDegree);

void computeAffineArithmeticSteps(int n,  StochasticSystemPtr sys);
#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_MODELPARSER_HH
