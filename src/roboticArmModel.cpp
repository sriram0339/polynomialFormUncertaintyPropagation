//
// Created by Sriram Sankaranarayanan on 2/16/20.
//
#include "MultivariatePoly.hh"
#include "StateAbstraction.hh"
#include "DistributionInfo.hh"
#include "ProbabilityQueryEvaluator.hh"
#include <chrono>

using namespace PolynomialForms;

namespace PolynomialForms {
    extern bool fourthMomentBoundCalculation;
};

inline MpfiWrapper deg2rad(double ang){
    MpfiWrapper pi(3.1415);
    MpfiWrapper angRad = MpfiWrapper(ang) *pi/MpfiWrapper(180.0);
    return angRad;
}
struct RoboticArmModelSimulation {
    int randomVarID = 0;
    map<int, DistributionInfoPtr> dInfo;
    map<int, MpfiWrapper> env;
    int createUniform(double a, double b){
        int varID = randomVarID;
        DistributionInfoPtr d1 = std::make_shared<UniformDistributionInfo>(MpfiWrapper(a,b));
        dInfo.insert(std::make_pair(varID, d1));
        env.insert(std::make_pair(varID, MpfiWrapper(a,b)));
        randomVarID ++;
        return varID;
    }

    int createTruncGaussian(double mean, double sd, double a, double b){
        int varID = randomVarID;
        DistributionInfoPtr d1 = std::make_shared<TruncNormalDistributionInfo>
                (MpfiWrapper(a, b), mean, sd);
        dInfo.insert(std::make_pair(varID, d1));
        env.insert(std::make_pair(varID, MpfiWrapper(a,b)));
        randomVarID ++;
        return varID;

    }

    void symbolicSimulateSystem(int maxDegree, int numReps) {
        auto start = chrono::high_resolution_clock::now();
        vector<double> angles = {10, 60, 110, 160, 140, 100, 60, 20, 10, 0};
        int v0 = createTruncGaussian(0.0, 0.05, -0.5, 0.5);
        int v1 = createTruncGaussian(0.0, 0.1, -0.5, 0.5);
        MultivariatePoly x(1.0, v0);
        MultivariatePoly y(1.0, v1);

        for (int i = 0; i < numReps; ++i){
            std::cout << "i = " << i << std::endl;
            for (double ang: angles){
                int dVar = createUniform(0.98, 1.02);
                MultivariatePoly dPoly(1.0, dVar);
                int tVar = createTruncGaussian(0.0, 0.01, -0.05, 0.05);
                MultivariatePoly tPoly(1.0, tVar);
                tPoly.setConst(1.0);
                tPoly.scaleAssign(deg2rad(ang));
                MultivariatePoly tmp1 = tPoly.cosine(env);
                MultivariatePoly tmp2 = tmp1.rangeMultiply(dPoly, env);
                x.scaleAndAddAssign(1.0, tmp2);
                MultivariatePoly tmp3 = tPoly.sine(env);
                MultivariatePoly tmp4 = tmp3.rangeMultiply(dPoly, env);
                y.scaleAndAddAssign(1.0, tmp4);
                x.truncateAssign(maxDegree, env);
                y.truncateAssign(maxDegree, env);

            }

        }

        std::cout << "After " << numReps << "Cycles:" << std::endl;
        std::cout << "x = " << std::endl;
        x.prettyPrint(std::cout, std::map<int, string> ());
        std::cout << "y = " << std::endl;
        y.prettyPrint(std::cout, std::map<int, string>());
        auto end = chrono::high_resolution_clock::now();
        double time_taken =
                chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        time_taken *= 1e-9;

        std::cout << "E(x) = " << x.expectation(dInfo)<< std::endl;
        std::cout << "Range(x) = " << x.evaluate(env) << std::endl;
        ProbabilityQueryEvaluator pqe(x, dInfo);
        pqe.separatePolynomialIntoComponents();
        pqe.printPolyStats();
        std::cout << "Evaluating upper tail bounds on x " << std::endl;
        pqe.computeBestUpperTailBounds(272.0);

        pqe.computeBestUpperTailBounds(274.0);

        pqe.computeBestUpperTailBounds(275.0);

        pqe.computeBestUpperTailBounds(277.0);

        auto end2 = chrono::high_resolution_clock::now();
        double time_taken2 =
                chrono::duration_cast<chrono::nanoseconds>(end2 - end).count();
        time_taken2 *= 1e-9;
        std::cout << " Time Taken: " << std::endl;
        std::cout << "Poly form calculations: " << time_taken << std::endl;
        std::cout << "Query evaluations: " << time_taken2 << std::endl;


    }

};

void computeRoboticArmModel(int maxDegree, int numReps){
    RoboticArmModelSimulation s;
    std::cout << "Running robotic arm model for " << numReps <<" steps" << std::endl;
    s.symbolicSimulateSystem(maxDegree, numReps);
}