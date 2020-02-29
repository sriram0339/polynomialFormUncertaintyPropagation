//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#include "PolynomialFunctions.hh"

namespace PolynomialForms {
    MultivariatePoly computeSinXByX(MultivariatePoly const & p, std::map<int, MpfiWrapper> const &var_env, int resultDegree) {
        /*--
         * Compute power series of sine around 0.0
         *   Just make sure that the power of x is one short
         *   Note that with remainder sin(x) =
         *       x - x^3/6  - cos(rng) x^4/24
         *    sin(x)/x = 1 - x^2/6 - cos(rng) x^3/24
         */

        MultivariatePoly res(1.0);
        MultivariatePoly tmp = p.squarePoly();
        MpfiWrapper rng = p.evaluate(var_env);
        MpfiWrapper cos_rng = cos(rng);

        //res = res -tmp/6
        res.scaleAndAddAssign(MpfiWrapper(-1.0/6.0), tmp);
        //tmp3 = p^3
        MultivariatePoly tmp3 = tmp.multiply(p);
        // res = res -1/24 cos(rng)*x^3
        res.scaleAndAddAssign(rng * MpfiWrapper(-1.0/24), tmp3);
        res.centerAssign(var_env);
        res.truncateAssign(resultDegree, var_env);
        return res;
    }

    MultivariatePoly
    computeCosXMinusOneByX(MultivariatePoly const &p, std::map<int, MpfiWrapper> const &var_env, int resultDegree) {
        /*
         * cos (x) - 1= -x^2/2 + x^4/4!  -sin(rng) * x^5/5!
         * cos(x)-1/x = -x/2 + x^3/24 -sin(rng) * x^4/120
         */
        MultivariatePoly res(0.0);
        MpfiWrapper rng = p.evaluate(var_env);
        MpfiWrapper sin_rng = sin(rng);

        res.scaleAndAddAssign(MpfiWrapper(-0.5), p);
        MultivariatePoly tmp2 = p.squarePoly();
        MultivariatePoly tmp3 = p.multiply(tmp2);
        res.scaleAndAddAssign(MpfiWrapper(1.0/24), tmp3);
        MultivariatePoly tmp4 = tmp3.multiply(p);
        res.scaleAndAddAssign(MpfiWrapper(-1.0/120) * sin_rng, tmp4);
        res.centerAssign(var_env);
        res.truncateAssign(resultDegree, var_env);
        return res;
    }


};