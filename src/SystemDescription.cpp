//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#include "SystemDescription.hh"


namespace PolynomialForms{

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

};