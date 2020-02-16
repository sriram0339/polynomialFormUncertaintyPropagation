//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_PROBABILITYQUERYEVALUATOR_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_PROBABILITYQUERYEVALUATOR_HH

#include "MultivariatePoly.hh"
#include "DistributionInfo.hh"

namespace PolynomialForms {

    class ProbabilityQueryEvaluator {
    protected:
        MultivariatePoly mp;
        std::map<int, DistributionInfoPtr > distributionInfo;
        std::vector<MultivariatePoly> splitComponents;

        double t;
    public:
        ProbabilityQueryEvaluator(MultivariatePoly const & mp_, std::map<int, DistributionInfoPtr > const & dInfo_):
        mp(mp_), distributionInfo(dInfo_){
            MpfiWrapper cTerm = mp.constantIntvl();
            mp.setConst(0);
            t = - cTerm.lower(); // p + [l,u] >= 0 only if p >= -l
        }

        void separatePolynomialIntoComponents();


    };

};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_PROBABILITYQUERYEVALUATOR_HH
