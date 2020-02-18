//
// Created by Sriram Sankaranarayanan on 2/17/20.
//

#include "MultivariatePoly.hh"

using namespace PolynomialForms;

void test1(){
    // Test the polynomial manipulations
    // 1. Create a polynomial
    std::map<int, string> name_map = {{0, "x0"}, {1, "x1"}, {2, "x2"}, {3, "x3"}};
    MultivariatePoly p0(5.0, 2);
    p0.addToConst(3.0);
    std::cout << "p0 = ";
    p0.prettyPrint(std::cout, name_map);

    MultivariatePoly p1(1.0, 1);
    p1.addToConst(2.0);
    MultivariatePoly p2 = p1.squarePoly();

    std:: cout << "p1 = " ;
    p1.prettyPrint(std::cout, name_map);
    std::cout << "p2 = ";
    p2.prettyPrint(std::cout, name_map);

    MultivariatePoly p3 = p2.multiply(p0);
    std::cout << "p3 = p0 * p2 = ";
    p3.prettyPrint(std::cout, name_map);

    PowerProduct pp0;
    PowerProduct pp1;
    pp0.setPower(1, 1);
    pp1.setPower(2, 1);
    PowerProduct pp2 = pp0.multiply(pp1);
    pp2.prettyPrint(std::cout, name_map);
    std::cout << std::endl;
    p2.setTerm(pp2, MpfiWrapper(2.5));
    std::cout << "p2' = ";
    p2.prettyPrint(std::cout, name_map);

    MultivariatePoly p5 = p2.powPoly(3);
    std::cout << "p5 = p2^3 = ";
    p5.prettyPrint(std::cout, name_map);

    p5.scaleAndAddAssign(MpfiWrapper(-1.0), p2.squarePoly());
    std::cout << "p5 = p2^3  - p2^2 = ";
    p5.prettyPrint(std::cout, name_map);

    p5.scaleAssign(MpfiWrapper(0.1));
    std::cout << "0.1*p5 =  ";
    p5.prettyPrint(std::cout, name_map);

}

std::map<int, MpfiWrapper> makeSampleMap(std::map<int, MpfiWrapper> const & envMap) {
std::map<int, MpfiWrapper> retMap;
for (auto p: envMap){
double d = p.second.getRandomSample();
retMap.insert(make_pair(p.first, MpfiWrapper(d)));
}
return retMap;
}

void sineAndCosineTest(){
    /* -- Check that the sine and cosine are yielding appropriate ranges --*/
    std::map<int, string> name_map = {{0, "x0"}, {1, "x1"}, {2, "x2"}, {3, "x3"}};
    MultivariatePoly p0(5.0, 2);
    p0.addToConst(3.0);
    PowerProduct pp0, pp1, pp2;
    pp0.setPower(0, 1);
    pp0.setPower(2, 1);
    p0.setTerm(pp0, 1.0);
    pp1.setPower(1,2);
    p0.setTerm(pp1, -2.0);
    pp2.setPower(3, 1);
    p0.setTerm(pp2, -0.5);
    p0.setConst(0.1);
    std::map<int, MpfiWrapper> env_map = {{0, MpfiWrapper(-0.1, 0.1)},
                                          {1, MpfiWrapper(-0.05, 0.1)},
                                          {2, MpfiWrapper(-0.05, 0.0)},
                                          {3, MpfiWrapper(0.8, 1.2)}};
    std::cout << "p0 = ";
    p0.prettyPrint(std::cout, name_map);
    std::cout << "range of p0 = " ;
    MpfiWrapper m1 = p0.evaluate(env_map);
    std:: cout << m1 << std::endl;
    std::cout << "range of sine(p0) = " << sin(m1) << std::endl;
    std::cout << "range of cosine(p0) = " << cos(m1) << std::endl;

    MultivariatePoly p1 = p0.sine(env_map);
    std::cout << "sine(p0) = p1 = " ;
    p1.prettyPrint(std::cout, name_map);
    std::cout << "range of p1 = " << p1.evaluate(env_map) << std::endl;

    MultivariatePoly p2 = p0.cosine(env_map);
    std::cout << "cosine(p0) = p2 = " ;
    p2.prettyPrint(std::cout, name_map);
    std::cout << "range of p2 = " << p2.evaluate(env_map) << std::endl;


    // Draw Samples
    double nSamples = 10000.0;
    MpfiWrapper mean0(0.0), mean1(0.0), mean2(0.0), mean3(0.0);
    for (int i =0; i < (int) nSamples; ++i){
        std::map<int, MpfiWrapper> sampleMap = makeSampleMap(env_map);
        MpfiWrapper m0 = p0.evaluate(sampleMap);
        mean0 = mean0 + sin(m0);
        MpfiWrapper m1 = p1.evaluate(sampleMap);
        mean1 = mean1  + m1;
        MpfiWrapper m2 = p2.evaluate(sampleMap);
        mean2 = mean2  +m2;
        mean3 = mean3  +cos(m0);

    }

    std::cout << " sine of samples: " << mean0/nSamples << std::endl;
    std::cout << " expectation of sine approximation" << mean1/nSamples << std::endl;
    std::cout << " cosine of samples: " << mean3/nSamples << std::endl;
    std::cout << " expectation of cosine approximation" << mean2/nSamples << std::endl;



}
