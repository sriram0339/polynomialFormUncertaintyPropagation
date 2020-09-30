//
// Created by Sriram Sankaranarayanan on 2/3/20.
//

#include "MultivariatePoly.hh"
#include "ModelParser.hh"
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream
#include "inference.h"

namespace PolynomialForms{
    bool debug = false;
    bool fourthMomentBoundCalculation = false;
    bool doCenteringOfRVs =true;
};

using namespace PolynomialForms;

extern void computeRoboticArmModel(int maxDegree, int numReachSteps);
extern void computeRimlessWheel(int numReachSteps, double meanParam);
extern void computeCartPoleModel(int maxDegree, int numReps);
void computeCartPoleNonPolyModel(int maxDegree, int numReps);
void printHelpMessage(const char * progName){
    std::cout << "Usage: " << progName << " [options] [optional file name to parse]" << std::endl;
    std::cout << "Options: "<<std::endl;
    std::cout << "\t -n <number of time steps> -- Number of steps to run " << std::endl;
    std::cout << "\t -d <max poly form degree> -- Set max polynomial form degree " << std::endl;
    std::cout << "\t -r                        -- Run robotic Arm Ex. (do not provide a filename to parse)" << std::endl;
    std::cout << "\t -w <meanParameterValue>   -- Run rimless Wheel Ex. (do not provide a filename to parse)" << std::endl;
    std::cout << "\t -a                        -- Run affine arithmetic " << std::endl;
}
//const string modelDataFile = "./data/obsy.txt";
void bayesianInference_neuronmodel(int numReachSteps,int maxDegree){
    vector<string> varNameVec;
    string varName1 = "alpha";
    string varName2 = "beta";
    varNameVec.push_back(varName1);
    varNameVec.push_back(varName2);
    double delta_t = 0.02;
    double tf=delta_t*numReachSteps;//1.0;//final time
    double Tdel=0.2; //observation time step
    int n_samp=tf/Tdel;
    std::vector <int> indSamp;
    for (int i=1;i<=n_samp;i++){
        int ind_temp=i*(Tdel/delta_t);
        cout<<ind_temp<<endl;
        indSamp.push_back(ind_temp);
    }



//    double x10 = -1.0;
//    double x20 = 1.0;
//    std::vector<double> x0_initial(2);
//    x0_initial[0]=x10;
//    x0_initial[1]=x20;

    //define the range of parameter's space
    std::vector <double> maxtheta,mintheta;
    std::map<int, DistributionInfoPtr> Distrib = globalSystem->getInitialMap();
    std::map<string, int> varIDs = globalSystem->getvarIDs();
    for(int i=0;i<varNameVec.size();i++){
        int varID = varIDs[varNameVec[i]];
        MpfiWrapper range = Distrib[varID]->getRange();
        MpfiWrapper offset = Distrib[varID]->getOffset();
        range = range + offset;
        maxtheta.push_back(range.upper());
        mintheta.push_back(range.lower());
    }



    int numInt=20;
    std::vector <double> tol_uplo_lay;//the tolerance of |Upper-Lower| for refining in each layer [first layer, second layer ...]
    tol_uplo_lay.push_back(0.3);//set 1: no refine

    double sigma2=0.01;                   // noise variance
    vector<double> cov;
    cov.push_back(sigma2);
// add noise on x1 state
    std::vector<vector <double>> obs_y;
    ////read observed data
    const string modelDataFile = "./data/obsy_neuron.txt";
    std::ifstream infile(modelDataFile);
    std::string line;
    while (std::getline(infile, line))
    {
        vector <double> row;
        double a;
        std::istringstream iss(line);
        while((iss >> a)){// Read in from line stream
            row.push_back(a);
        }
        obs_y.push_back(row);
    }


    auto start = chrono::high_resolution_clock::now();
    bayesianPolyForm inf_dynamic(tf,delta_t,Tdel, obs_y, cov, numInt, indSamp,maxDegree, tol_uplo_lay,varNameVec);
    //inf_dynamic.computeBayesian(maxtheta,mintheta);
    //MultivariatePoly logp(MpfiWrapper(0.0));
    //inf_dynamic.test(numReachSteps);

    inf_dynamic.computeLikelihoodPolyFormBound( maxDegree, maxtheta, mintheta);

    //computeNSteps(numReachSteps, maxDegree);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    std::cout << "Time Taken: " << time_taken << std::endl;
}

void bayesianInference_ebola(int numReachSteps,int maxDegree){
    vector<string> varNameVec;
    string varName1 = "p1";
    string varName2 = "p2";
    string varName3 = "p3";
    varNameVec.push_back(varName1);
    varNameVec.push_back(varName2);
    varNameVec.push_back(varName3);
    double delta_t = 0.5;
    double tf=delta_t*numReachSteps;//1.0;//final time
    double Tdel=0.5; //observation time step
    int n_samp=tf/Tdel;
    std::vector <int> indSamp;
    for (int i=1;i<=n_samp;i++){
        int ind_temp=i*(Tdel/delta_t);
        cout<<ind_temp<<endl;
        indSamp.push_back(ind_temp);
    }




//    std::vector<double> x0_initial(5);
//    x0_initial[0]=0.8;
//    x0_initial[1]=0.3;
//    x0_initial[2]=0.2;
//    x0_initial[3]=0.0;
//    x0_initial[4]=0.02;

    //define the range of parameter's space
    std::vector <double> maxtheta,mintheta;
    std::map<int, DistributionInfoPtr> Distrib = globalSystem->getInitialMap();
    std::map<string, int> varIDs = globalSystem->getvarIDs();
    for(int i=0;i<varNameVec.size();i++){
        int varID = varIDs[varNameVec[i]];
        MpfiWrapper range = Distrib[varID]->getRange();
        MpfiWrapper offset = Distrib[varID]->getOffset();
        range = range + offset;
        maxtheta.push_back(range.upper());
        mintheta.push_back(range.lower());
    }


    int numInt=30;//20;
    std::vector <double> tol_uplo_lay;//the tolerance of |Upper-Lower| for refining in each layer [first layer, second layer ...]
    tol_uplo_lay.push_back(0.3);//set 1: no refine

    double sigma2=0.000025;                   // noise variance
    vector<double> cov;
    cov.push_back(sigma2);
// add noise on x1 state
    std::vector<vector <double>> obs_y;
    ////read observed data
    const string modelDataFile = "./data/obsy_ebola_i.txt";
    //const string modelDataFile = "./data/obsy_ebola_s.txt";
    std::ifstream infile(modelDataFile);
    std::string line;
    while (std::getline(infile, line))
    {
        vector <double> row;
        double a;
        std::istringstream iss(line);
        while((iss >> a)){// Read in from line stream
            row.push_back(a);
        }
        obs_y.push_back(row);
    }


    auto start = chrono::high_resolution_clock::now();
    bayesianPolyForm inf_dynamic(tf,delta_t,Tdel, obs_y, cov, numInt, indSamp,maxDegree, tol_uplo_lay,varNameVec);
    //inf_dynamic.computeBayesian(maxtheta,mintheta);
    //MultivariatePoly logp(MpfiWrapper(0.0));
    //inf_dynamic.test(numReachSteps);
    inf_dynamic.computeLikelihoodPolyFormBound( maxDegree, maxtheta, mintheta);

    //computeNSteps(numReachSteps, maxDegree);
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    std::cout << "Time Taken: " << time_taken << std::endl;
}

int main(int argc, char * argv[]){
    int c;
    int numReachSteps = 15;
    int maxDegree = 4;
    bool affineArithmeticDo = false;
    opterr = 0;

    while ((c = getopt (argc, argv, "n:d:hacetrw:4")) != -1)
        switch (c)
        {
            case 'a':
                affineArithmeticDo = true;
                break;
            case 'n':
                numReachSteps = atoi(optarg);
                break;
            case 'd':
                maxDegree = atoi(optarg);
                break;
            case 'h':
                printHelpMessage(argv[0]);
                exit(1);
                break;
            case 'r':
                computeRoboticArmModel(maxDegree, numReachSteps);
                exit(1);
                break;
            case 'c':
                computeCartPoleModel(maxDegree, numReachSteps);
                exit(1);
            case 'e':
                computeCartPoleNonPolyModel(maxDegree, numReachSteps);
                exit(1);
            case 't':
                doCenteringOfRVs = false;
                break;
            case 'w': {
                double meanParam = atof(optarg);
                std::cout << "Rimless Wheel mean Param: " << meanParam << std::endl;
                computeRimlessWheel(numReachSteps, meanParam);
                exit(1);
            }
                break;
            case '4':
                fourthMomentBoundCalculation = true;
                std::cerr << "WARNING: You are turning on 4th moment bounds -- expensive!" << std::endl;
                break;
            case '?':
                if (optopt == 'n' || optopt == 's')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",
                             optopt);
                return 1;
            default:
                abort ();
        }

    if (optind < argc){

        const char * fileName = argv[optind];
        std::cout << "Parsing file name: " << fileName << std::endl;
        parserMain(fileName);

        if (!affineArithmeticDo) {
            bayesianInference_neuronmodel(numReachSteps, maxDegree);
            //bayesianInference_ebola(numReachSteps, maxDegree);
//            auto start = chrono::high_resolution_clock::now();

//            vector<MpfiWrapper> rangeVec;
//            MpfiWrapper range;
//            string varName1 = "alpha";
//            range.set(0.295, 0.305);
//            rangeVec.push_back(range);
//            string varName2 = "beta";
//            range.set(0.147, 0.153);
//            rangeVec.push_back(range);
//            vector<string> varNameVec;
//            varNameVec.push_back(varName1);
//            varNameVec.push_back(varName2);
//
//            StateAbstractionPtr st = computeNStepsAndSaveData(numReachSteps, maxDegree,varNameVec ,rangeVec);
//            //range.set(0.325, 0.335);


//            auto end = chrono::high_resolution_clock::now();
//            double time_taken =
//                    chrono::duration_cast<chrono::nanoseconds>(end - start).count();
//            time_taken *= 1e-9;
////            std::cout << "Evaluating Queries" << std::endl;
////            globalSystem->evaluateQueries(st);
////            auto end2 = chrono::high_resolution_clock::now();
////            double time_taken2 =
////                    chrono::duration_cast<chrono::nanoseconds>(end2 - end).count();
////            time_taken2 *= 1e-9;
//            std::cout << " Time Taken: " << std::endl;
//            std::cout << "Bayesian calculations: " << time_taken << std::endl;
////            std::cout << "Poly form calculations: " << time_taken << std::endl;
////            std::cout << "Query evaluations: " << time_taken2 << std::endl;
        } else {
            auto start = chrono::high_resolution_clock::now();
            computeAffineArithmeticSteps(numReachSteps, std::make_shared<StochasticSystem>(*globalSystem));
            auto end = chrono::high_resolution_clock::now();
            double time_taken = chrono::duration_cast<chrono::nanoseconds>(end-start).count();
            time_taken *= 1e-9;
            std::cout << " Time Taken:" << time_taken << std::endl;
        }
    } else {
        printHelpMessage(argv[0]);
        std::cerr << "No file provided" << std::endl;
    }
    return 1;


    //computeRoboticArmModel();
    //computeRimlessWheel();

    //sineAndCosineTest();
    return 1;
}
