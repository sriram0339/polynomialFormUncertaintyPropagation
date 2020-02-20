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


};