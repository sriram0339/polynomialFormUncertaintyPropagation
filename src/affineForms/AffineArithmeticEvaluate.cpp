//
// Created by Sriram Sankaranarayanan on 6/3/20.
//

#include "AffineArithmeticEvaluate.hh"


AffineArithmeticClass AAEvaluateExpr::evaluate(ExprPtr what){
    switch(what -> getType()){
        case VAR_TYPE: {
            std::shared_ptr<Var> whatVar = dynamic_pointer_cast<Var>(what);
            int id = whatVar -> getID();
            auto it = stMap.find(id);
            assert(it != stMap.end());
            AffineArithmeticClass retVal(env);
            //std::cout << "it -> second = " << it -> second << std::endl;
            retVal.add_assign(it -> second);
            //std::cout << "retval = " << retVal << std::endl;
            return retVal;
        }
        break;
        case CONST_TYPE: {
            std::shared_ptr<Const> whatConst = dynamic_pointer_cast<Const>(what);
            AffineArithmeticClass retVal(env);
            retVal.setConstant(whatConst -> getConst());
            return retVal;
        }
        break;
        case PLUS_TYPE:{
            std::shared_ptr<Plus> whatPlus = dynamic_pointer_cast<Plus>(what);
            std::vector<ExprPtr> const & subExprs = whatPlus -> getSubExprs();
            std::vector<MpfiWrapper> const & scaleFactors = whatPlus -> getScaleFactors();
            int numSubs = subExprs.size();
            assert(scaleFactors.size() == numSubs);
            AffineArithmeticClass retVal(env);
            for (int i = 0; i < numSubs; ++i){
                AffineArithmeticClass tmp = evaluate(subExprs[i]);
                tmp.scale_assign(scaleFactors[i]);
                retVal.add_assign(tmp);
            }
            return retVal;
        } break;
        case MULT_TYPE:{
            std::shared_ptr<Star> whatMult = dynamic_pointer_cast<Star>(what);
            std::vector<ExprPtr> const & subExprs = whatMult -> getSubExprs();
            AffineArithmeticClass retVal(env);
            retVal.setConstant(i_double(1.0));
            for (auto s: subExprs) {
                AffineArithmeticClass tmp = evaluate(s);
                //retVal.multiply_assign(tmp);
                retVal.linearized_multiply_assign(tmp);
            }
            return retVal;
        } break;
        case DIV_TYPE: {
            // Division not supported
            assert(false);
        }
        break;
        case SIN_TYPE: {
            std::shared_ptr<Trig> whatSin = dynamic_pointer_cast<Trig>(what);
            ExprPtr s = whatSin -> getSubExpr();
            AffineArithmeticClass tmp = evaluate(s);
            AffineArithmeticClass retVal(env);
            retVal.sine_assign(tmp);
            return retVal;

        } break;
        case COS_TYPE: {
            std::shared_ptr<Trig> whatCos = dynamic_pointer_cast<Trig>(what);
            ExprPtr s = whatCos -> getSubExpr();
            AffineArithmeticClass tmp = evaluate(s);
            AffineArithmeticClass retVal(env);
            retVal.cosine_assign(tmp);
            return retVal;
        }
        case DISTRIB_TYPE: {
            std::shared_ptr<Distrib> whatDistrib = dynamic_pointer_cast<Distrib>(what);
            DistributionInfoPtr dPtr = whatDistrib -> getDistributionInfo();
            MpfiWrapper rng = dPtr -> getRange();
            MpfiWrapper expect = dPtr -> getExpectation();
            MpfiWrapper moment2 = dPtr-> getMoment(2);
            AffineArithmeticClass retVal = randomVariable(env, rng, expect, moment2);
            retVal.setConstant(dPtr -> getOffset());
            return retVal;
        }
        case POW_TYPE: {
            std::shared_ptr<Pow> whatPow = dynamic_pointer_cast<Pow>(what);
            ExprPtr e = whatPow-> getSubExpr();
            int p = whatPow -> getPow();
            AffineArithmeticClass tmp = evaluate(e);
            AffineArithmeticClass retVal(env);
            retVal.power_assign(tmp, p);
            return retVal;
        } break;
        default:
            assert(false);
            break;
    }
    assert(false);
    return AffineArithmeticClass(env);
}

stateMap AAEvaluateExpr::evaluateUpdates(std::map<int, ExprPtr> const &updates) {
    stateMap newMap;
    for(auto p: updates){
        AffineArithmeticClass tmp = evaluate(p.second);
     //   std::cout << "Debug: " << tmp << std::endl;
        newMap.insert(std::make_pair(p.first, tmp));
    }
    return (newMap);
}

stateMap computeInitialMap(AAEnvironment & env, StochasticSystemPtr const & what){
    stateMap initMap;
    auto dPtr = what -> getInitialMap();
    for (auto p: dPtr){
        int varID = p.first;
        DistributionInfoPtr dPtr = p.second;
        MpfiWrapper rng = dPtr -> getRange();
        MpfiWrapper expect = dPtr -> getExpectation();
        MpfiWrapper moment2 = dPtr-> getMoment(2);
        AffineArithmeticClass retVal = randomVariable(env, rng, expect, moment2);
        retVal.setConstant(dPtr -> getOffset());
        initMap.insert(std::make_pair(varID, retVal));
    }
    return (initMap);
}

stateMap  computeOneStep(AAEnvironment & env,
        StochasticSystemPtr const & what,
        stateMap const & oldMap){

    AAEvaluateExpr aa (env, oldMap);
    auto update_map = what->getUpdates();
    return aa.evaluateUpdates(update_map);
}

void probabilityQueryEvaluateAA(ExprPtr query, AAEnvironment & env, stateMap const & what){
    AAEvaluateExpr aa(env, what);
    AffineArithmeticClass qExpr = aa.evaluate(query);
    //std::cout << qExpr << std::endl;
    std::cout << "--------------" << std::endl;
    qExpr.performAnalysisForCMI(std::cout);
    std::cout << "--------------" << std::endl;
}

void expectationQueryEvaluateAA(ExprPtr query, AAEnvironment & env, stateMap const & map){
    AAEvaluateExpr aa(env, map);
    AffineArithmeticClass qExpr = aa.evaluate(query);
    std::cout << "Expectation Range: " << qExpr.expectation() << std::endl;
}




