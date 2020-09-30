//
// Created by Sriram Sankaranarayanan on 2/10/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_DISTRIBUTIONINFO_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_DISTRIBUTIONINFO_HH

#include <string>
#include <map>
#include <vector>
#include "MpfiWrapper.hh"
#include "truncatedNormalCode.hh"
#include <memory>

namespace PolynomialForms {
    extern bool doCenteringOfRVs ;

    class VariableDistributionInfo {
    protected:

        MpfiWrapper range;
        MpfiWrapper offset;
    public:
        VariableDistributionInfo( MpfiWrapper range_, MpfiWrapper offset_):  range(range_), offset(offset_){};
        virtual MpfiWrapper getRange() const { return range; };
        virtual MpfiWrapper getOffset() const {return offset;};
        virtual MpfiWrapper getExpectation() const = 0;
        virtual MpfiWrapper getMoment(int j) const = 0;
        virtual ~VariableDistributionInfo() = default;
        virtual void setRange(MpfiWrapper & range_) { range.set(range_) ; };
        virtual void setOffset(MpfiWrapper & offset_) { offset.set(offset_) ; }
    };
    typedef std::shared_ptr<VariableDistributionInfo> DistributionInfoPtr;


    class TruncNormalDistributionInfo: public VariableDistributionInfo {
    protected:
        MpfiWrapper mean;
        MpfiWrapper sdev;
    public:
        TruncNormalDistributionInfo( MpfiWrapper range_, MpfiWrapper mean_, MpfiWrapper sdev_):
        VariableDistributionInfo(range_, 0.0),
        mean(mean_),
        sdev(sdev_) {
            if (doCenteringOfRVs) {
                /* Center the distribution? */
                offset = mean;
                mean = MpfiWrapper(0.0);
                range = range - offset;
            }

        };

        virtual MpfiWrapper getExpectation() const { return mean; };

        virtual MpfiWrapper getMoment(int j) const {
            if (j == 1) { return getExpectation(); }
            double a = range.lower();
            double b = range.upper();
            double mu = median(mean);
            double sigma = median(sdev);
            return MpfiWrapper(truncated_normal_ab_moment (j, mu, sigma, a, b));
        }

    };


    class UniformDistributionInfo: public VariableDistributionInfo {
    protected:
        MpfiWrapper mean;
    public:
        UniformDistributionInfo(MpfiWrapper range_):
                VariableDistributionInfo(range_, 0.0),
                mean(median(range_))
                {
                if (doCenteringOfRVs) {
                    /* Center the distribution? */
                    offset = mean;
                    mean = MpfiWrapper(0.0);
                    range = range - offset;
                }
                };

        virtual MpfiWrapper getExpectation() const { return mean; };

        virtual MpfiWrapper getMoment(int j) const {
            if (j == 1) { return getExpectation(); }
            MpfiWrapper den(1.0 + (double) j);

            MpfiWrapper a = range.upper();
            MpfiWrapper b = range.lower();
            den = den * (b-a);
            MpfiWrapper v = pow(b, j+1) - pow(a, j+1);
            MpfiWrapper retVal = v/den;
            return retVal;
        }



    };


};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_DISTRIBUTIONINFO_HH
