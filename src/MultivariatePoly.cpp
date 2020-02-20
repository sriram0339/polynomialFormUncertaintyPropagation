//
// Created by Sriram Sankaranarayanan on 2/5/20.
//

#include "MpfiWrapper.hh"
#include "MultivariatePoly.hh"
#include <cmath>

namespace PolynomialForms {


    MultivariatePoly MultivariatePoly::multiply(MultivariatePoly const &what) const {
        MultivariatePoly retPoly (what.constIntvl * constIntvl);
        // First deal with the constant term
        retPoly.scaleAndAddAssign(constIntvl, what, false);
        retPoly.scaleAndAddAssign(what.constIntvl, *this, false);
        for (auto p: terms){
            retPoly.multiplyAndAddAssign(p.second, p.first, what, false);
        }
        return retPoly;
    }

    MultivariatePoly MultivariatePoly::squarePoly() const {
        MpfiWrapper r = square(constIntvl);
        MultivariatePoly retPoly(r);
        // First deal with the constant term
        retPoly.scaleAndAddAssign(MpfiWrapper(2.0) * constIntvl, *this, false);
        for (auto const p: terms){
            retPoly.multiplyAndAddAssign(p.second, p.first, *this, false);
        }
        return retPoly;
    }

    MultivariatePoly MultivariatePoly::powPoly(int j) const {
        MultivariatePoly retPoly(MpfiWrapper(1.0));
        MultivariatePoly tmpPoly(*this);
        if (j == 0){
            return retPoly;
        }
        if (j == 1){
            return MultivariatePoly(*this);
        }
        if (j == 2) {
            return squarePoly();
        }

        while (j > 0) {
            if (j %2 == 1){
                retPoly = retPoly.multiply(tmpPoly);
            }
            j = j/2;
            if ( j > 0)
                tmpPoly = tmpPoly.squarePoly();
        }
        return retPoly;
    }


    int trig_num_terms = 5;

    trig_deriv_t nextDeriv(trig_deriv_t state) {
        switch (state){
            case PLUS_SINE:
                return PLUS_COSINE;
            case PLUS_COSINE:
                return MINUS_SINE;
            case MINUS_SINE:
                return MINUS_COSINE;
            case MINUS_COSINE:
                return PLUS_SINE;
            default:
                assert(false);
        }
    }

    MultivariatePoly MultivariatePoly::computeTrig(
                std::map<int, MpfiWrapper> const &var_env,
                trig_deriv_t derivState
            ) const {

        MpfiWrapper rng = this -> evaluate(var_env);
        // Choose the center point of this range.
        MpfiWrapper c = median(rng);
        MpfiWrapper sinc = sin(c);
        MpfiWrapper cosc = cos(c);
        MultivariatePoly retPoly;
        switch (derivState){
            case PLUS_COSINE:
                retPoly.setConst(cosc);
                break;
            case PLUS_SINE:
                retPoly.setConst(sinc);
                break;
            default:
                assert(false);
        }

        MultivariatePoly qPoly(*this);
        qPoly.addToConst(-1.0 * c);

        MultivariatePoly tmpPoly(1.0);
        MpfiWrapper fact(1.0);
        int i;
        for (i = 1; i <= trig_num_terms; ++i){
            derivState = nextDeriv(derivState); // derivative state gives us what derivative to apply
            tmpPoly = tmpPoly.multiply(qPoly); // tmpPoly = (p - c)^i
            fact = fact * MpfiWrapper((double) i);
            MpfiWrapper curDeriv;
            switch (derivState) {
                case PLUS_SINE:
                    curDeriv = sinc/fact;
                    break;
                case PLUS_COSINE:
                    curDeriv = cosc/fact;
                    break;
                case MINUS_COSINE:
                    curDeriv = -cosc/fact;
                    break;
                case MINUS_SINE:
                    curDeriv = -sinc/fact;
                    break;
            }
            // retPoly = retPoly + (f^{(i)}(c))/i! * (p-c)^i
            retPoly.scaleAndAddAssign(curDeriv, tmpPoly, true);
        }

        // Npw for the reminder term
        derivState = nextDeriv(derivState);
        tmpPoly = tmpPoly.multiply(qPoly); // tmp = (p-c)^{i+1}
        fact = fact * MpfiWrapper((double) i); // i = num_trig_terms + 1
        MpfiWrapper remRange;
        //   rng = rng - c ; // correct the range by subtracting the center point
        switch (derivState) {
            case PLUS_SINE:
                remRange = sin(rng)/fact;
                break;
            case PLUS_COSINE:
                remRange = cos(rng)/fact;
                break;
            case MINUS_COSINE:
                remRange = -cos(rng)/fact;
                break;
            case MINUS_SINE:
                remRange = -sin(rng)/fact;
                break;
        }


        //MpfiWrapper polyRange = tmpPoly.evaluate(var_env);
        retPoly.scaleAndAddAssign(remRange, tmpPoly);
        retPoly.centerAssign(var_env);
        return retPoly;

    }

    MultivariatePoly MultivariatePoly::sine(std::map<int, MpfiWrapper> const &var_env) const {
        /*
         * We will use the following idea.
         *  sine(p(w)) = sine(c + q(w)) =
         *    sin(c) + q(w) cos(c) - sin(c) q^2(w)/2 - cos(c) q^3(w)/6 + Rem(sin(c) q^4(w)/24)
         *
         */
        return computeTrig(var_env, PLUS_SINE);

    }


    MultivariatePoly MultivariatePoly::cosine(std::map<int, MpfiWrapper> const &var_env) const {
        /*
         * We will use the following idea.
         *  cos(p(w)) = cos(c + q(w)) =
         *    cos(c) - q(w) sin(c) - cos(c) q^2(w)/2 + sin(c) q^3(w)/6 + Rem(cos(c) q^4(w)/24)
         *
         */


        return computeTrig(var_env, PLUS_COSINE);
    }


    MpfiWrapper MultivariatePoly::evaluate(std::map<int, MpfiWrapper> const &var_env) const {
        MpfiWrapper reslt(constIntvl);
        for (auto const p: terms){
            reslt = reslt + p.second * p.first.evaluate(var_env);
        }
        return reslt;
    }

    MpfiWrapper MultivariatePoly::expectation(std::map<int, DistributionInfoPtr> const &env) const {
        MpfiWrapper reslt(constIntvl);
        for (auto const p: terms){
            reslt = reslt + p.second * p.first.expectation(env);
        }
        return reslt;
    }


    void MultivariatePoly::prettyPrint(ostream &out, std::map<int, string> const &name_env) const {
        out << constIntvl;
        for (auto const p: terms){
            out << " + ";
            out << p.second << "*";
            p.first.prettyPrint(out, name_env);
        }
        out << std::endl;
    }

    bool isTooSmall(MpfiWrapper s){
        return fabs(s.upper()) <= 1E-8 && fabs(s.lower()) <= 1E-8;
    }

    MultivariatePoly MultivariatePoly::truncate(int maxDegree,
            std::map<int, MpfiWrapper> const &var_env ) const {
        MpfiWrapper constPart = constIntvl;
        MultivariatePoly retPoly(constIntvl);
        for (auto p: terms){
            if (p.first.totalDegree() <= maxDegree  && !isTooSmall(p.second)){
                retPoly.setTerm(p.first, p.second);
            } else {
                MpfiWrapper intvl = p.first.evaluate(var_env);
                constPart = constPart + intvl * p.second;
            }
        }
        retPoly.setConst(constPart);
        return retPoly;
    }

    void MultivariatePoly::truncateAssign(int maxDegree,  std::map<int, MpfiWrapper> const &var_env ) {
       // MpfiWrapper constPart = constIntvl;
        for (auto it = terms.cbegin(); it != terms.cend() /* not hoisted */; /* no increment */)
            if (it -> first.totalDegree() > maxDegree || isTooSmall(it -> second)){
                MpfiWrapper intvl = it -> first.evaluate(var_env);
                constIntvl = constIntvl + intvl * it -> second;
                terms.erase(it++);
            } else {
                ++it;
            }
        }

    MultivariatePoly MultivariatePoly::rangeMultiply(const MultivariatePoly &what,
            const std::map<int, MpfiWrapper> &var_env) const {
        MultivariatePoly tmp1 = *this;
        MpfiWrapper crng1 = constIntvl - median(constIntvl);
        tmp1.setConst(median(constIntvl));
        MpfiWrapper rng1 = tmp1.evaluate(var_env);

        MultivariatePoly tmp2 = what;
        MpfiWrapper crng2 = tmp2.getConstIntvl();
        tmp2.setConst(median(crng2));
        crng2 = crng2 - median(crng2);
        MpfiWrapper rng2 = tmp2.evaluate(var_env);

        MpfiWrapper uncertaintyIntvl = rng2 * crng1 + rng1 * crng2;
        MultivariatePoly resPoly = tmp1.multiply(tmp2);
        resPoly.addToConst(uncertaintyIntvl);
        return resPoly;
    }

    void MultivariatePoly::centerAssign(std::map<int, MpfiWrapper> const & env) {
        /*
         * Make sure that all terms have coefficients are centered around their median
         * Addd the rest to the const Interval
         */

        for (auto & p : terms ){
            MpfiWrapper c = median(p.second);
            MpfiWrapper e = p.second - c;
            MpfiWrapper r = p.first.evaluate(env) * e;
            constIntvl = constIntvl + r;
            p.second = c;
        }
        return;
    }

    int num_recip_terms = 4;

    MultivariatePoly MultivariatePoly::reciprocal(std::map<int, MpfiWrapper> const & var_env) const {
        MpfiWrapper rng = this -> evaluate(var_env);
        // Choose the center point of this range.
        MpfiWrapper c = median(rng);
        if (rng.lower() <= 0.0 && rng.upper() >= 0.0){
            std::cout << "Error: taking reciprocal of a range that contains zero " << std::endl;
            assert(false);
        }
        MpfiWrapper recip = inverse(c);
        MultivariatePoly retPoly(recip);
        MultivariatePoly qPoly(*this);
        qPoly.addToConst(-1.0 * c);

        MultivariatePoly tmpPoly(1.0);
        MpfiWrapper fact(1.0);
        int i;
        for (i = 1; i <= num_recip_terms; ++i ){
            tmpPoly = tmpPoly.multiply(qPoly); // tmpPoly = (p - c)^i
            MpfiWrapper curDeriv = pow(-1.0, i) * pow(recip, (i+1));
            retPoly.scaleAndAddAssign(curDeriv, tmpPoly, true);
        }

        MpfiWrapper recipRange = inverse(rng);
        MpfiWrapper rem = pow(-1.0, i) * pow(recipRange, i+1);
        tmpPoly = tmpPoly.multiply(qPoly);
        retPoly.scaleAndAddAssign(rem, tmpPoly, true);
        retPoly.centerAssign(var_env);
        return retPoly;
    }


};