//
// Created by Sriram Sankaranarayanan on 2/19/20.
//

#include "MultivariatePoly.hh"
#include "StateAbstraction.hh"
#include "DistributionInfo.hh"
#include "ProbabilityQueryEvaluator.hh"
#include <chrono>
#include <cmath>

using namespace PolynomialForms;

namespace PolynomialForms {
    extern bool fourthMomentBoundCalculation;
};

inline MpfiWrapper deg2rad(double ang){
    MpfiWrapper pi(3.1415);
    MpfiWrapper angRad = MpfiWrapper(ang) *pi/MpfiWrapper(180.0);
    return angRad;
}
struct ModelSimulation {
    int randomVarID = 0;
    map<int, DistributionInfoPtr> dInfo;
    map<int, MpfiWrapper> env;

    int createUniform(double a, double b) {
        int varID = randomVarID;
        DistributionInfoPtr d1 = std::make_shared<UniformDistributionInfo>(MpfiWrapper(a, b));
        dInfo.insert(std::make_pair(varID, d1));
        env.insert(std::make_pair(varID, MpfiWrapper(a, b)));
        randomVarID++;
        return varID;
    }

    int createTruncGaussian(double mean, double sd, MpfiWrapper nRange) {
        int varID = randomVarID;
        DistributionInfoPtr d1 = std::make_shared<TruncNormalDistributionInfo>
                (nRange, mean, sd);
        dInfo.insert(std::make_pair(varID, d1));
        env.insert(std::make_pair(varID, nRange));
        randomVarID++;
        return varID;

    }

    void symbolicSimulateSystem(int maxDegree, int numReps) {
        std::cout << "MaxDegree: " << maxDegree << std::endl;
        auto start = chrono::high_resolution_clock::now();
        double delta = 0.1;
        MpfiWrapper nRange(-10.0, 10.0);
        double nMean = 0.0;
        double nVar = 1.0;
        int x0  =createUniform(-0.1, 0.1);
        int theta0 = createUniform(-0.1, 0.1);
        MultivariatePoly x(1.0, x0), dxdt(0.0), theta(1.0, theta0), dthetadt(0.0);


        for (int i = 0; i < numReps; ++i) {
            cout << i << endl;
            int w1 = createTruncGaussian(nMean, nVar, nRange);
            int w2 = createTruncGaussian(nMean, nVar, nRange);
            int w3 = createTruncGaussian(nMean, nVar, nRange);
            int w4 = createTruncGaussian(nMean, nVar, nRange);
            MultivariatePoly feedbackTerm;
            feedbackTerm.scaleAndAddAssign(10.0, x);
            feedbackTerm.scaleAndAddAssign(-289.83, theta);
            feedbackTerm.scaleAndAddAssign(19.53, dxdt);
            feedbackTerm.scaleAndAddAssign(-63.25, dthetadt);
            MultivariatePoly w1Poly(0.01 * sqrt(delta), w1);
            MultivariatePoly w2Poly(0.01 * sqrt(delta), w2);
            MultivariatePoly w3Poly(0.01 * sqrt(delta), w3);
            MultivariatePoly w4Poly(0.01 * sqrt(delta), w4);
            MultivariatePoly ddx;
            ddx.scaleAndAddAssign(0.1, feedbackTerm);
            ddx.scaleAndAddAssign(0.98, theta);
            MultivariatePoly tmp0;
            tmp0.scaleAndAddAssign(-0.75, theta);
            tmp0.scaleAndAddAssign(-0.1, feedbackTerm);
            MultivariatePoly thetaSq = theta.powPoly(2);
            thetaSq = thetaSq.truncate(maxDegree, env);
            tmp0 = tmp0.multiply(theta);
            MultivariatePoly tmp1 = dthetadt.powPoly(2);
            tmp1 = tmp1.truncate(maxDegree, env);
            tmp1.scaleAssign(-0.05);
//            tmp1.linearized_multiply_assign(theta);
            tmp1 = tmp1.multiply(thetaSq);
            tmp1 = tmp1.truncate(maxDegree, env);
            //ddx = ddx + tmp0 + tmp1;
            ddx.scaleAndAddAssign(1.0, tmp0);
            ddx.scaleAndAddAssign(1.0, tmp1);
            MultivariatePoly ddtheta;
            ddtheta.scaleAndAddAssign(0.2, feedbackTerm);
            ddtheta.scaleAndAddAssign(21.56, theta);
//    ddtheta = i_double(0.2)* feedbackTerm + i_double(21.56) * theta;
            MultivariatePoly tmp2;
            tmp2.scaleAndAddAssign(-5.75, theta);
            tmp2.scaleAndAddAssign(-0.12, feedbackTerm);
//    AffineArithmeticClass tmp2(env);
            //  tmp2 = i_double(-5.75)*theta + i_double(-0.12) * feedbackTerm;
            //tmp2.linearized_multiply_assign(thetaSq);
            tmp2 = tmp2.multiply(thetaSq);
            tmp2 = tmp2.truncate(maxDegree, env);
            //AffineArithmeticClass tmp3(env);
            MultivariatePoly tmp3 = dthetadt.powPoly(2);
            tmp3 = tmp3.truncate(maxDegree, env);
//    tmp3.power_assign(dthetadt,2);
            //  tmp3.scale_assign(i_double(-0.1));
            tmp3.scaleAssign(-0.1);
            ddtheta.scaleAndAddAssign(1.0, tmp2);
            ddtheta.scaleAndAddAssign(1.0, tmp3);

            dxdt = dxdt.truncate(maxDegree, env);
            dxdt.centerAssign(env);
            x.scaleAndAddAssign(delta, dxdt);
            x.scaleAndAddAssign(1.0, w1Poly);

            theta.scaleAndAddAssign(delta, dthetadt);
            theta.scaleAndAddAssign(1.0, w2Poly);

            dxdt.scaleAndAddAssign(delta, ddx);
            dxdt.scaleAndAddAssign(1.0, w3Poly);

            dthetadt.scaleAndAddAssign(delta, ddtheta);
            dthetadt.scaleAndAddAssign(1.0, w4Poly);

            x.centerAssign(env);
            x = x.truncate(maxDegree, env);

            theta.centerAssign(env);
            theta = theta.truncate(maxDegree, env);

            dxdt.centerAssign(env);
            dxdt = dxdt.truncate(maxDegree, env);

            dthetadt.centerAssign(env);
            dthetadt = dthetadt.truncate(maxDegree, env);
//    x = x + delta * dxdt + w1;
//    theta = theta + delta * dthetadt + w2;
            //  dxdt = x + delta * ddx + w3;
//    dthetadt = dthetadt+ delta * ddtheta + w4;

        }

        auto end = chrono::high_resolution_clock::now();
        double time_taken =
                chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        time_taken *= 1e-9;

        std::cout << "E(x) = " << x.expectation(dInfo)<< std::endl;
        std::cout << "Range(x) = " << x.evaluate(env) << std::endl;

        std::cout << "E(theta) = " << theta.expectation(dInfo)<< std::endl;
        std::cout << "Range(theta) = " << theta.evaluate(env) << std::endl;

        ProbabilityQueryEvaluator pqe1(theta, dInfo);
        std::cout << "Computing tail probability P(theta <= 3.1415/6) <= ??" << std::endl;
        pqe1.separatePolynomialIntoComponents();
        pqe1.printPolyStats();
        pqe1.computeBestUpperTailBounds(3.1415/6);

        ProbabilityQueryEvaluator pqe2(x, dInfo);
        std::cout << "Computing tail probability P(x >= 2.0) <= ??" << std::endl;
        pqe2.separatePolynomialIntoComponents();
        pqe2.printPolyStats();
        pqe2.computeBestUpperTailBounds(2.0);

        auto end2 = chrono::high_resolution_clock::now();
        double time_taken2 =
                chrono::duration_cast<chrono::nanoseconds>(end2 - end).count();
        time_taken2 *= 1e-9;
        std::cout << " Time Taken: " << std::endl;
        std::cout << "Poly form calculations: " << time_taken << std::endl;
        std::cout << "Query evaluations: " << time_taken2 << std::endl;
        return;
    }

    void symbolicSimulateNonPolynomialSystem(int maxDegree, int numReps) {

        std::cout << "MaxDegree: " << maxDegree << std::endl;
        auto start = chrono::high_resolution_clock::now();
        double delta = 0.1;
        MpfiWrapper nRange(-10.0, 10.0);
        double nMean = 0.0;
        double nVar = 1.0;
        int x0 = createUniform(-0.1, 0.1);
        int theta0 = createUniform(-0.1, 0.1);
        MultivariatePoly x(1.0, x0), dxdt(0.0), theta(1.0, theta0), dthetadt(0.0);
        double L = 0.5;
        double g = 9.8;
        double mp = 1.0;
        double mc = 10.0;

        for (int i = 0; i < numReps; ++i) {
            cout << i << endl;
            int w1 = createTruncGaussian(nMean, nVar, nRange);
            int w2 = createTruncGaussian(nMean, nVar, nRange);
            int w3 = createTruncGaussian(nMean, nVar, nRange);
            int w4 = createTruncGaussian(nMean, nVar, nRange);
            MultivariatePoly feedbackTerm;
            feedbackTerm.scaleAndAddAssign(10.0, x);
            feedbackTerm.scaleAndAddAssign(-289.83, theta);
            feedbackTerm.scaleAndAddAssign(19.53, dxdt);
            feedbackTerm.scaleAndAddAssign(-63.25, dthetadt);
            MultivariatePoly w1Poly(0.01 * sqrt(delta), w1);
            MultivariatePoly w2Poly(0.01 * sqrt(delta), w2);
            MultivariatePoly w3Poly(0.01 * sqrt(delta), w3);
            MultivariatePoly w4Poly(0.01 * sqrt(delta), w4);

            MultivariatePoly tmp0 = dthetadt.powPoly(2);
            tmp0.truncateAssign(maxDegree, env);
            tmp0.scaleAssign(L);
            // tmp2 = (mc + mp sin(theta)^2)
            MultivariatePoly sinTheta = theta.sine(env);
            MultivariatePoly cosTheta = theta.cosine(env);
            MultivariatePoly tmp1 = sinTheta;
            tmp1.truncateAssign(maxDegree, env);
            MultivariatePoly tmp2 = tmp1.powPoly(2);
            tmp2.truncateAssign(maxDegree, env);
            tmp2.scaleAssign(mp);
            tmp2.addToConst(mc);
            //denTerm = 1/tmp2
            MultivariatePoly denTerm = tmp2.reciprocal(env);
            /* AffineArithmeticClass ddx(env);
               if (cheapMul){
                 AffineArithmeticClass tmp4(env);
                 tmp4 = mp * sine(theta);
                 AffineArithmeticClass tmp5(env);
                 tmp5 = tmp0 - g * cosine(theta);
                 tmp4.linearized_multiply_assign(tmp5);
                 ddx = feedbackTerm - tmp4 ;
                 ddx.linearized_multiply_assign(denTerm);
               } else {
                 ddx = (feedbackTerm - mp*sine(theta)*(tmp0 - g * cosine(theta)) ) * denTerm;
               } */
            MultivariatePoly tmp4 = tmp1;
            tmp4.scaleAssign(mp);
            MultivariatePoly tmp5 = cosTheta;
            tmp5.scaleAssign(-g);
            tmp5.scaleAndAddAssign(1.0, tmp0);
            tmp5.truncateAssign(maxDegree, env);
            tmp4 = tmp4.multiply(tmp5);
            MultivariatePoly ddx;
            ddx.scaleAndAddAssign(-1.0, tmp4);
            ddx.scaleAndAddAssign(1.0, feedbackTerm);
            ddx = ddx.multiply(denTerm);
            ddx.truncateAssign(maxDegree, env);


            /*
            AffineArithmeticClass ddtheta(env);
            if (cheapMul){
              AffineArithmeticClass tmp3(feedbackTerm);
              tmp3.linearized_multiply_assign(cosine(theta));
              tmp0.linearized_multiply_assign(sine( i_double(2.0) * theta));
              tmp0.scale_assign (i_double(0.5)* mp);
              ddtheta = tmp3 - tmp0 + (mc+mp)*g*sine(theta);
              ddtheta.linearized_multiply_assign(denTerm);
              ddtheta.scale_assign(i_double(1.0)/L);
            } else {
              ddtheta = (feedbackTerm * cosine(theta) - mp * tmp0 * i_double(0.5)*sine(i_double(2.0)*theta)  + (mc+mp)*g * sine(theta)) * denTerm * (i_double(1.0)/L);
            } */

            MultivariatePoly tmp3;
            tmp3.scaleAndAddAssign(1.0, feedbackTerm);
            tmp3 = tmp3.multiply(cosTheta);
            tmp3.truncateAssign(maxDegree, env);
            MultivariatePoly sin2Theta(theta);
            sin2Theta.scaleAssign(2.0);
            sin2Theta = sin2Theta.sine(env);

            tmp0 = tmp0.multiply(sin2Theta);
            tmp0.scaleAssign(0.5 * mp);
            tmp0.truncateAssign(maxDegree, env);
            MultivariatePoly ddtheta(tmp3);
            ddtheta.scaleAndAddAssign(-1.0, tmp0);
            ddtheta.scaleAndAddAssign((mc+mp)*g, sinTheta);
            ddtheta = ddtheta.multiply(denTerm);
            ddtheta.scaleAssign(1.0/L);



            x.scaleAndAddAssign(delta, dxdt);
            x.scaleAndAddAssign(1.0, w1Poly);

            theta.scaleAndAddAssign(delta, dthetadt);
            theta.scaleAndAddAssign(1.0, w2Poly);

            dxdt.scaleAndAddAssign(delta, ddx);
            dxdt.scaleAndAddAssign(1.0, w3Poly);

            dthetadt.scaleAndAddAssign(delta, ddtheta);
            dthetadt.scaleAndAddAssign(1.0, w4Poly);

            x.centerAssign(env);
            x.truncateAssign(maxDegree, env);

            theta.centerAssign(env);
            theta.truncateAssign(maxDegree, env);

            dxdt.centerAssign(env);
            dxdt.truncateAssign(maxDegree, env);

            dthetadt.centerAssign(env);
            dthetadt.truncateAssign(maxDegree, env);


        }
        auto end = chrono::high_resolution_clock::now();
        double time_taken =
                chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        time_taken *= 1e-9;

        std::cout << "E(x) = " << x.expectation(dInfo)<< std::endl;
        std::cout << "Range(x) = " << x.evaluate(env) << std::endl;

        std::cout << "E(theta) = " << theta.expectation(dInfo)<< std::endl;
        std::cout << "Range(theta) = " << theta.evaluate(env) << std::endl;

        ProbabilityQueryEvaluator pqe1(theta, dInfo);
        std::cout << "Computing tail probability P(theta <= 3.1415/6) <= ??" << std::endl;
        pqe1.separatePolynomialIntoComponents();
        pqe1.printPolyStats();
        pqe1.computeBestUpperTailBounds(3.1415/6);

        ProbabilityQueryEvaluator pqe2(x, dInfo);
        std::cout << "Computing tail probability P(x >= 4.0) <= ??" << std::endl;
        pqe2.separatePolynomialIntoComponents();
        pqe2.printPolyStats();
        pqe2.computeBestUpperTailBounds(4.0);

        auto end2 = chrono::high_resolution_clock::now();
        double time_taken2 =
                chrono::duration_cast<chrono::nanoseconds>(end2 - end).count();
        time_taken2 *= 1e-9;
        std::cout << " Time Taken: " << std::endl;
        std::cout << "Poly form calculations: " << time_taken << std::endl;
        std::cout << "Query evaluations: " << time_taken2 << std::endl;
        return;

    }

};

void computeCartPoleModel(int maxDegree, int numReps){
    ModelSimulation s;
    std::cout << "Running cartpole (polynomial) model for " << numReps <<" steps" << std::endl;
    s.symbolicSimulateSystem(maxDegree, numReps);
}

void computeCartPoleNonPolyModel(int maxDegree, int numReps){
    ModelSimulation s;
    std::cout << "Running cartpole (non polynomial) model for " << numReps <<" steps" << std::endl;
    s.symbolicSimulateNonPolynomialSystem(maxDegree, numReps);
}