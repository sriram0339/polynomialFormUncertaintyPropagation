//
// Created by Sriram Sankaranarayanan on 2/14/20.
//
#include "ModelParser.hh"
#include "affineForms/AffineArithmeticEvaluate.hh"

StochasticSystem * globalSystem;
extern "C" FILE *yyin;
extern int yyparse();

void parserMain(const char * fileName){
    globalSystem = new StochasticSystem();
    FILE * f = fopen(fileName, "r");
    if (!f) {
        std::cerr << "FATAL: failed to open file -- " << fileName << std::endl;
        exit(1);
    }

    yyin = f;
    yyparse();
    fclose(f);
    return;
}


StateAbstractionPtr computeNSteps(int n, int maxDegree){
    StateAbstractionPtr st = globalSystem -> initialize(maxDegree);
    for (int i = 0 ; i < n; ++i){
        std::cout << "After " << i << " Steps" << std::endl;
        globalSystem -> prettyPrintStateAbstraction(std::cout, st);
        globalSystem -> computeOneStep(st, maxDegree);
    }
    std::cout << "After " << n << " Steps" << std::endl;
    globalSystem -> prettyPrintStateAbstraction(std::cout, st);
    return st;
}


void computeAffineArithmeticSteps(int n,  StochasticSystemPtr sys) {
    AAEnvironment env(0);
    stateMap st = computeInitialMap(env, sys);
    for (int i = 0; i < n; ++i){
        std::cout << i << std::endl;
        st = computeOneStep(env, sys, st);
    }
    std::vector<Query> const & queries = sys -> getQueries();
    for (auto q: queries){
        std::cout << "---- Query: " << q.getID() << "---------" << std::endl;
        switch(q.getType()){
            case PolynomialForms::PROB_QUERY:
                probabilityQueryEvaluateAA(q.getExpr(), env, st);
                break;
            case PolynomialForms::EXPECT_QUERY:
                expectationQueryEvaluateAA(q.getExpr(), env, st);
                break;
        }
    }
}

void computeNStepsAndSaveData(int n, int maxDegree,std::vector<std::string> & varNameVec, std::vector<MpfiWrapper> & rangeVec){
    //initialize and clean data
    StateAbstractionPtr st = globalSystem -> initializeAndSetVarRange(maxDegree,varNameVec, rangeVec);//initialize(maxDegree);
    vector <vector <double> > ().swap(globalSystem->result_upx1);
    vector <vector <double> > ().swap(globalSystem->result_lowx1);
    //globalSystem ->setVarRange(varName, range);
    for (int i = 0 ; i < n; ++i){
        //std::cout << "After " << i << " Steps" << std::endl;
        //globalSystem -> prettyPrintStateAbstraction(std::cout, st);
        globalSystem -> computeOneStep(st, maxDegree);
        //globalSystem->evaluateQueries(st);
        globalSystem->evaluateQueriesofRange(st,1);//1:save data into result_upx1,result_lowx1

    }
//    std::cout << "After " << n << " Steps" << std::endl;
//    cout<<globalSystem -> result_upx1[globalSystem -> result_upx1.size()-1][0]<<std::endl;
//    cout<<globalSystem -> result_lowx1[globalSystem -> result_lowx1.size()-1][0]<<std::endl;
    //globalSystem -> prettyPrintStateAbstraction(std::cout, st);
    //return st;
}

void computeNSteps_inference(int n, int maxDegree,StateAbstractionPtr & st){
    for (int i = 0 ; i < n; ++i){
        globalSystem -> computeOneStep(st, maxDegree);
    }
}

StateAbstractionPtr initialize_globalSystem(int maxDegree){
    StateAbstractionPtr st = globalSystem -> initialize(maxDegree);
    return st;
}






