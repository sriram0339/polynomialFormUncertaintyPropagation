//
// Created by Sriram Sankaranarayanan on 2/17/20.
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

struct RimlessWheelMdl {
    /*-- cut and paste from roboticArmModel.cpp --*/
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

    /*-- Model actual dynamics --*/
    void symbolicSimulateSystem(int numReps) {
        auto start = chrono::high_resolution_clock::now();
        int x0 = createUniform(-0.1, 0.1); // Initialize x
        MpfiWrapper theta(30.0);
        MpfiWrapper piBy180(3.1415/180);
        MpfiWrapper cos2theta = 0.75;
        MultivariatePoly x(1.0, x0);
        MpfiWrapper gByL(10.0);
        for (int i = 0; i < numReps; ++i){
            int wVar = createTruncGaussian(8, 1.5, 0, 16);
            MultivariatePoly wPoly(1.0, wVar);
            MultivariatePoly beta1(piBy180 * theta/2.0);
            beta1.scaleAndAddAssign(piBy180, wPoly); // beta1 = pi/180 (theta/2 + w)
            MultivariatePoly beta2(piBy180 * theta/2.0);

            beta2.scaleAndAddAssign(-1.0 * piBy180, wPoly); // beta2 = pi/180 (theta/2 - w)
            MultivariatePoly cosBeta1  =beta1.cosine(env);
            MultivariatePoly cosBeta2 = beta2.cosine(env);

            cosBeta1.scaleAssign(-2.0 * gByL);
            cosBeta1.addToConst(2.0*gByL);

            cosBeta2.scaleAssign(2.0 * gByL);
            cosBeta2.addToConst(-2.0*gByL);

            x.scaleAndAddAssign(1.0, cosBeta1);
            x.scaleAssign(cos2theta);
            x.scaleAndAddAssign(1.0, cosBeta2);
        }

        auto end = chrono::high_resolution_clock::now();
        double time_taken =
                chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        time_taken *= 1e-9;

        std::cout << "E(x) = " << x.expectation(dInfo)<< std::endl;
        std::cout << "Range(x) = " << x.evaluate(env) << std::endl;

        std::cout << "Concentration of measure inequalities:" << std::endl;
        x.scaleAssign(-1.0);
        ProbabilityQueryEvaluator pqe1(x,dInfo);
        std::cout << "P(x <= 0) <= ??" << std::endl;
        pqe1.separatePolynomialIntoComponents();
        std::cout << "Chernoff Bounds: " << pqe1.computeChernoffBound() << std::endl;
        std::cout << "Chebyshev Bounds: " << pqe1.computeChebyshevBounds() << std::endl;
        if (fourthMomentBoundCalculation) {
            std::cout << "Fourth moment bounds:" << pqe1.computeFourthMomentBound() << std::endl;
        }

        auto end2 = chrono::high_resolution_clock::now();
        double time_taken2 =
                chrono::duration_cast<chrono::nanoseconds>(end2 - end).count();
        time_taken2 *= 1e-9;
        std::cout << " Time Taken: " << std::endl;
        std::cout << "Poly form calculations: " << time_taken << std::endl;
        std::cout << "Query evaluations: " << time_taken2 << std::endl;
    }


};

void computeRimlessWheel(int numReachSteps){
    RimlessWheelMdl rmw;
    std::cout << "Running rimless wheel model for " << numReachSteps <<" steps" << std::endl;
    rmw.symbolicSimulateSystem(numReachSteps);
}