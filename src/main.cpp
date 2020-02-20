//
// Created by Sriram Sankaranarayanan on 2/3/20.
//

#include "MultivariatePoly.hh"
#include "ModelParser.hh"
#include <cstdlib>
#include <unistd.h>
#include <cstdio>

namespace PolynomialForms{
    bool debug = false;
    bool fourthMomentBoundCalculation = false;
};

using namespace PolynomialForms;

extern void computeRoboticArmModel(int numReachSteps);
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
}

int main(int argc, char * argv[]){
    int c;
    int numReachSteps = 15;
    int maxDegree = 4;
    opterr = 0;

    while ((c = getopt (argc, argv, "n:d:hcerw:4")) != -1)
        switch (c)
        {
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
                computeRoboticArmModel(numReachSteps);
                exit(1);
                break;
            case 'c':
                computeCartPoleModel(maxDegree, numReachSteps);
                exit(1);
            case 'e':
                computeCartPoleNonPolyModel(maxDegree, numReachSteps);
                exit(1);
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


        auto start = chrono::high_resolution_clock::now();
        StateAbstractionPtr st = computeNSteps(numReachSteps, maxDegree);
        auto end = chrono::high_resolution_clock::now();
        double time_taken =
                chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        time_taken *= 1e-9;
        std::cout << "Evaluating Queries"<< std::endl;
        globalSystem -> evaluateQueries(st);
        auto end2 = chrono::high_resolution_clock::now();
        double time_taken2 =
                chrono::duration_cast<chrono::nanoseconds>(end2 - end).count();
        time_taken2 *= 1e-9;
        std::cout << " Time Taken: " << std::endl;
        std::cout << "Poly form calculations: " << time_taken << std::endl;
        std::cout << "Query evaluations: " << time_taken2 << std::endl;
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