//
// Created by Sriram Sankaranarayanan on 2/5/20.
//

#include "MpfiWrapper.hh"
#include "MultivariatePoly.hh"

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

    MultivariatePoly MultivariatePoly::sine(std::map<int, MpfiWrapper> const &var_env) const {
        /*
         * We will use the following idea.
         *  sine(p(w)) = sine(c + q(w)) =
         *    sin(c) + q(w) cos(c) - sin(c) q^2(w)/2 - cos(c) q^3(w)/6 + Rem(sin(c) q^4(w)/24)
         *
         */
         MpfiWrapper rng = this -> evaluate(var_env);
         // Choose the center point of this range.
         MpfiWrapper c = median(rng);
         MpfiWrapper sinc = sin(c);
         MpfiWrapper cosc = cos(c);
         MultivariatePoly retPoly(sinc);
         MultivariatePoly tmpPoly(*this);
         tmpPoly.addToConst(-1.0 * c);
         retPoly.scaleAndAddAssign(cosc, tmpPoly);
         MultivariatePoly pows = tmpPoly.squarePoly();
         retPoly.scaleAndAddAssign(-0.5 * sinc, pows );
         MultivariatePoly pows3 = pows.multiply(*this);
         retPoly.scaleAndAddAssign(MpfiWrapper(-1.0/6.0) * cosc, pows3 );
         MultivariatePoly pows4 = pows . squarePoly();
         pows4.scaleAssign(MpfiWrapper(1.0/24.0)* sin(rng));
         MpfiWrapper rem = pows4.evaluate(var_env);
         retPoly.addToConst(rem);
         return retPoly;
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

    MultivariatePoly MultivariatePoly::cosine(std::map<int, MpfiWrapper> const &var_env) const {
        /*
         * We will use the following idea.
         *  cos(p(w)) = cos(c + q(w)) =
         *    cos(c) - q(w) sin(c) - cos(c) q^2(w)/2 + sin(c) q^3(w)/6 + Rem(cos(c) q^4(w)/24)
         *
         */

        MpfiWrapper rng = this -> evaluate(var_env);
        // Choose the center point of this range.
        MpfiWrapper c = median(rng);
        MpfiWrapper sinc = sin(c);
        MpfiWrapper cosc = cos(c);
        MultivariatePoly retPoly(cosc);
        MultivariatePoly tmpPoly(*this);
        tmpPoly.addToConst(-1.0 * c);
        retPoly.scaleAndAddAssign(MpfiWrapper(-1.0) * sinc, tmpPoly);
        MultivariatePoly pows = tmpPoly.squarePoly();
        retPoly.scaleAndAddAssign(-0.5 * cosc, pows );
        MultivariatePoly pows3 = pows.multiply(*this);
        retPoly.scaleAndAddAssign(MpfiWrapper(1.0/6.0) * sinc, pows3 );
        MultivariatePoly pows4 = pows . squarePoly();
        pows4.scaleAssign(MpfiWrapper(1.0/24.0)* cos(rng));
        MpfiWrapper rem = pows4.evaluate(var_env);
        retPoly.addToConst(rem);
        return retPoly;
    }

    MultivariatePoly MultivariatePoly::truncate(int maxDegree,
            std::map<int, MpfiWrapper> const &var_env ) const {
        MpfiWrapper constPart = constIntvl;
        MultivariatePoly retPoly(constIntvl);
        for (auto p: terms){
            if (p.first.totalDegree() <= maxDegree){
                retPoly.setTerm(p.first, p.second);
            } else {
                MpfiWrapper intvl = p.first.evaluate(var_env);
                constPart = constPart + intvl * p.second;
            }
        }
        retPoly.setConst(constPart);
        return retPoly;
    }


};