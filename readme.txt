Project:
Bayesian Parameter Estimation for Nonlinear Dynamical Systems using Polynomial Expansions

Prerequisites:
MPFI library

Getting Started:
It uses cmake.  To compile:
- cmake ./CMakeLists.txt
- make

Running the examples:
You can try a command like 
./polynomialFormUncertaintyPropagation -n 250 -d 6 test/neuronmodel.sys

(The commands for other benchmarks in this paper:)
./polynomialFormUncertaintyPropagation -n 20 -d 6 test/p53-model.sys
./polynomialFormUncertaintyPropagation -n 30 -d 10 test/ebola-with-params.sys
./polynomialFormUncertaintyPropagation -n 5 -d 4 test/honeybee-with-params.sys
./polynomialFormUncertaintyPropagation -n 5 -d 4 test/LaubLoomis-with-params.sys

Output files of the polynomial for the log-likelihood: logPoly_trunc.txt and logPoly.txt

Next, to perform optimization for Variational Inference:

1. Mean-field univariate Gaussian:
run VI_meanfield_univariateGaussian.m

2. Multivariate Gaussian:
run VI_multivariateGaussian.m

3. Mean-field histogram approximation
run VI_meanfield_hist.m 


