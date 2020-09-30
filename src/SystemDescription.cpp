//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#include "ProbabilityQueryEvaluator.hh"
#include "SystemDescription.hh"
#include <chrono>

namespace PolynomialForms{

    extern bool fourthMomentBoundCalculation;

    int StochasticSystem::addVar(std::string varName, DistributionInfoPtr dPtr) {
        int varID = numStateVars;
        numStateVars++;
        varIDs.insert(make_pair(varName, varID));
        varNames.insert(make_pair(varID, varName));
        initialDistrib.insert(make_pair(varID, dPtr));
        return  varID;
    }

    void StochasticSystem::addUpdate(int varID, ExprPtr rhs) {
        updates.insert(make_pair(varID, rhs));
    }

    int StochasticSystem::getVarIDFromName(std::string varName) {
        auto it = varIDs.find(varName);
        if (it == varIDs.end()){
            return -1; // SIGNAL FOR not found
        } else {
            return it -> second;
        }
    }

    StateAbstractionPtr StochasticSystem::initialize(int maxDegree) {
        StateAbstractionPtr st = std::make_shared<StateAbstraction>(maxDegree);
        st -> initialize(initialDistrib);
        return st;
    }

    void StochasticSystem::computeOneStep(StateAbstractionPtr st, int maxDegree) {
        std::map<int, MultivariatePoly> newStateMap;
        for (auto p: updates){
            int varID = p.first;
            ExprPtr e = p.second;
            MultivariatePoly mp2 = e -> evaluate(st);
            if (maxDegree > 0) {
                MultivariatePoly mp2_trunc = mp2.truncate(maxDegree, st->getRangeMapForNoiseSymbols());
                mp2_trunc.centerAssign(st->getRangeMapForNoiseSymbols());
                newStateMap.insert(make_pair(varID, mp2_trunc));
            } else {
                newStateMap.insert(make_pair(varID, mp2));
            }
        }
        st -> setStateMap(newStateMap);
    }

    void StochasticSystem::prettyPrintStateAbstraction(std::ostream &what, StateAbstractionPtr st) {
        std::map<int, string> name_env;

        for (int varID = 0; varID < numStateVars; ++varID){
            what << varNames[varID] << std::endl;
            MultivariatePoly const & p =  st -> getPolynomialForVar(varID);
            p.prettyPrint(what, name_env);

        }
    }

    void StochasticSystem::evaluateQueries(StateAbstractionPtr st) {
        for (auto q: queries){
            ExprPtr e = q.getExpr();
            MultivariatePoly p = e -> evaluate(st);
            std::cout << "Evaluating query " << q.getID() << std::endl;
            switch (q.getType()){
                case PROB_QUERY: {
                    ProbabilityQueryEvaluator pqe(p, st->getNoiseSymbolInfoMap());
                    pqe.separatePolynomialIntoComponents();
                    pqe.computeBestUpperTailBounds(0.0);
                };
                    break;
                case EXPECT_QUERY: {
                    MpfiWrapper r = p.expectation(st->getNoiseSymbolInfoMap());
                    std::cout << "\t RESULT: " << r << std::endl;
                    MpfiWrapper rHat = p.evaluate(st->getRangeMapForNoiseSymbols());
                    std::cout << "\t RANGE: " << rHat << std::endl;
                }
                    break;
            }

            std::cout << "Finished query " << q.getID() << std::endl;


        }
    }

//    void StochasticSystem::evaluateLogLikelihood_AtCur(std::ostream &what, StateAbstractionPtr st,int maxDegree) {
//        std::map<int, string> name_env;
//        double dy=0.2;
//
//        for (auto q: queries){
//            ExprPtr e = q.getExpr();
//            MultivariatePoly p = e -> evaluate(st);
//            p.scaleAndAddAssign(-1 , MultivariatePoly(MpfiWrapper(dy)));
//            MultivariatePoly p2=p.squarePoly();
//            if(maxDegree > 0 ) {
//                MultivariatePoly p2_trunc = p.truncate(maxDegree, st->getRangeMapForNoiseSymbols());
//                p2_trunc.centerAssign(st->getRangeMapForNoiseSymbols());
//                p2_trunc.prettyPrint(what, name_env);
//            }
//
//            cout<<" here"<<endl;
//        }
//
//    }

    void StochasticSystem::evaluateQueriesofRange(StateAbstractionPtr st, bool isSave) {
        vector<MpfiWrapper> v_rHat;
        for (auto q: queries){
            ExprPtr e = q.getExpr();
            MultivariatePoly p = e -> evaluate(st);
            //std::cout << "Evaluating query " << q.getID() << std::endl;

            MpfiWrapper rHat = p.evaluate(st->getRangeMapForNoiseSymbols());
            //std::cout << "\t RANGE: " << rHat << std::endl;
            v_rHat.push_back(rHat);
            //std::cout << "Finished query " << q.getID() << std::endl;
        }
        if(isSave){writeQueriesdata(v_rHat);}
    }

    void StochasticSystem::writeQueriesdata(const vector<MpfiWrapper> & v_rHat){
        vector<double> v_rHat_up,v_rHat_lo;
        for(auto rHat:v_rHat){
            v_rHat_up.push_back(rHat.upper());
            v_rHat_lo.push_back(rHat.lower());
        }
        result_upx1.push_back(v_rHat_up);
        result_lowx1.push_back(v_rHat_lo);

    }

    void StochasticSystem::setVarRange(std::string & varName, MpfiWrapper & range){
        int varID = varIDs[varName];
        /* Center the distribution? */
        MpfiWrapper offset = median(range);
        range = range - offset;
        initialDistrib[varID]->setRange(range);
        initialDistrib[varID]->setOffset(offset);
//        if(varName=="alpha"){
//            //cout<<dPtr->getOffset()<<endl;
//            cout<<initialDistrib[varID]->getOffset()<<endl;
//        }
    }
    StateAbstractionPtr StochasticSystem::initializeAndSetVarRange(int maxDegree, std::vector<std::string> & varNameVec, std::vector<MpfiWrapper> & rangeVec) {
        int np=varNameVec.size();
        for(int i=0;i<np;i++){
            setVarRange(varNameVec[i], rangeVec[i]);
        }
        StateAbstractionPtr st = std::make_shared<StateAbstraction>(maxDegree);
        st -> initialize(initialDistrib);
        return st;
    }




};
