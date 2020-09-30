#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <MpfiWrapper.hh>
#include <MultivariatePoly.hh>
#include "ModelParser.hh"
#include <MpfiWrapper.hh>

using namespace std;
typedef vector<double> (*FnPtr)(vector<double>,vector<double>,double);
#ifndef inference
#define inference
void mono_next( int m, vector <int> & x );
double CalculateMean(vector<double> & value);
void getsubcell (const int m,const int k_order,vector<vector<int> > & v);
void cell_generator (const int n,const int kmax, vector<vector<int> > & table );
double log_mvnpdf(const vector<double> & mu,const vector<double> & x,const vector <vector<double>>& inv_cov);
void computePosterior(const int numInt,const int sdim,const string model,const double dt,const int t,const vector<double> & cur_x,const vector<double> & next_x,const vector <vector<double>>& inv_cov,const vector <vector<double>> theta_list,map<string, FnPtr> & modelMap,vector<double> & logPrior,vector<double> & logPosterior);
//void UpAndLowGaussian_flowstar(const int & j,const double & dy, const vector <double> & cov,const double &  maxdx,const double &  mindx,double & sum1_lo,double & sum1_up);
//void CellTarget_pdf_flowstar(const vector <double> & X,const vector <Interval> & x0,const double tf,const double delta_t,const std::string & y_str,const Variables & stateVars,const Polynomial_ODE & ode,Continuous_Reachability_Setting & crs, const vector<vector <double > > & obs_y, const vector <double> & cov,const vector <double> & delta_p,const double tol_uplo,const vector <int> & indSamp, vector <double> & logp);
//void mostLikelyCell_flowstar(const vector <Interval> & x0_initial,const double tf, const double delta_t,const std::string & y_str,const Variables & stateVars,const Polynomial_ODE & ode, Continuous_Reachability_Setting & crs,const vector <vector <double>> & obs_y,const vector <double> & cov,const vector <double> & maxtheta,const vector <double> & mintheta,const int numInt,const double tol_uplo,const double tor, const vector <int> & indSamp, vector <vector <double>> & theta_list, vector <double> & delta_p,vector <vector <double>> & thetaProb_list, vector <vector <double>> & res_deltap_list,vector <vector <double>> & res_theta_list,vector <vector <double>> & res_thetaProb_list);
namespace PolynomialForms {
    class bayesianPolyForm {
    protected:

        double tf;
        double delta_t;
        double Tdel;
        vector <vector <double>> obs_y;
        vector <double> cov;
        int numInt;
        vector <int> indSamp;
        int maxDegree;
        vector<double> tol_uplo_lay;
        vector<std::string> varNameVec;
        //void initialize();

    public:

        bayesianPolyForm(double tf_,double delta_t_,double Tdel_, vector <vector <double>> obs_y_,
                         vector <double> cov_,int numInt_,vector <int> indSamp_,int maxDegree_, vector<double> tol_uplo_lay_,vector<std::string> varNameVec_) :
                tf(tf_),
                delta_t(delta_t_),
                Tdel(Tdel_),
                obs_y(obs_y_),
                cov(cov_),
                numInt(numInt_),
                indSamp(indSamp_),
                maxDegree(maxDegree_),
                tol_uplo_lay(tol_uplo_lay_),
                varNameVec(varNameVec_){};


        void computeBayesian(const vector <double> & maxtheta,const vector <double> & mintheta);
        void UpAndLowGaussian_polyform(const int & j,const double & dy,const double &  maxdx,const double &  mindx,double & sum1_lo,double & sum1_up);
        void CellTarget_pdf_polyform(const vector <double> & X, const vector <double> & delta_p,const double tol_uplo, vector <double> & logp);
        void CellTarget_pdf_polyform_refine(const vector <double> & X, const vector <double> & delta_p,const double tol_uplo, vector <double> & logp);
        void allCell_polyform( const vector <double> & maxtheta,const vector <double> & mintheta, const double tol_uplo, vector <vector <double>> & theta_list, vector <double> & delta_p,vector <vector <double>> & thetaProb_list, vector <vector <double>> & res_deltap_list,vector <vector <double>> & res_theta_list,vector <vector <double>> & res_thetaProb_list);

        //void computeNStepForm(bool debug=false);
        //void printExpectationsAndRanges();
        //void printSpreadsheetRow(ostream & outFileHandle);
        void evaluateLogLikelihood_AtCur( StateAbstractionPtr st,int maxDegree, double dy, MultivariatePoly & logp);
        void computeLikelihoodPolyFormBound(int maxDegree,const vector <double> & maxtheta,const vector <double> & mintheta);
        StateAbstractionPtr fullSpace_pdf_polyform( MultivariatePoly & logp);
        void computePolyIntegralAndSaveResults(const vector <double> & maxtheta,const vector <double> & mintheta, const MultivariatePoly & logp,StateAbstractionPtr & st);
            std::map<int, MpfiWrapper> ConvertToRangeMapForNoiseSymbols(std::map<int, DistributionInfoPtr> & Distrib);
        void computePolyIntegralAndSaveResults_Taylor(const vector <double> & maxtheta,const vector <double> & mintheta, const MultivariatePoly & logp, StateAbstractionPtr & st);
        void test(int n);
        void normalizationOfInterval_discrete( vector <vector <double>> const & res_thetaProb_list, vector <vector <double>> & res_thetaProb_listNorm );

    };
}
#endif
