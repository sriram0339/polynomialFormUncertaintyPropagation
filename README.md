This package contains code and appendix.

## Appendix to the paper

Containing proofs and benchmark details.

See file `appendix.PDF`

## Dependencies:

The code uses the CMake compilation system. It has been tested on a MAC OSX
machine.



To compile, you will need to install the libraries
- C++ BOOST library
- GNU MPFR library for multi precision floating point
- MPFI library
- glpk library

You will also need parser generator and lexical analyzers
- Flex
- Bison 

These are installed from packages on homebrew (MAC OSX) or
apt (Linux/Ubuntu). 

~~~
$ apt-get install -y cmake libmpfr-dev libmpfi-dev flex bison libglpk
~~~

## Compile:

Type

~~~
$ cmake ./CMakeLists.txt
$ make
~~~

If your system has libraries installed at non-standard locations, you will
need to edit the files in the directory `cmake-includes`. Specifically the files FindMPFR.cmake and FindMPFI.cmake specify where to search for these libraries.


## Running to reproduce table entries

First go to the test directory

~~~
$ cd ./tests
~~~

### Rimless Wheel Model (run using polynomial form API functions)

~~~
../polynomialFormUncertaintyPropagation -n 2000 -t -d 2 -w 4
~~~

### 2D Robot model (run using polynomial form API functions)

~~~
../polynomialFormUncertaintyPropagation  -n 100 -d 2 -r
~~~

### Cartpole model (run using polynomial form API functions)


~~~
../polynomialFormUncertaintyPropagation  -n 8 -d 2 -c 
~~~

### Ebola model (run from model file)

~~~
$ cat ebola-model.sys
$ ../polynomialFormUncertaintyPropagation -d 6 -n 25 -4 ebola-model.sys 
~~~

### Honeybee model (run from model file)

~~~
$ cat honeybee.sys
$ ../polynomialFormUncertaintyPropagation -n 25 -d 4 -4 honeybee.sys 
~~~

### Coupled Vanderpol oscillator

~~~
$ cat coupledVanderpol-3.sys
$ ../polynomialFormUncertaintyPropagation -n 15 -d 2  coupledVanderPolOscillators-3.sys 
~~~

### Laub Loomis Model

~~~
$ cat laub-loomis.sys
$ ../polynomialFormUncertaintyPropagation -n 25 -d 4 laub-loomis.sys 
~~~

### Lattice 10x1 Model

~~~
$ cat 1d-lattice-10.sys
$ ../polynomialFormUncertaintyPropagation -n 10 -d 4 1d-lattice-10.sys 
~~~


# Reproducing Affine Arithmetic Results

NOTE: Code requested and provided as is by Bouissou et al. Used as a library
inside our tool. Rimless wheel, 2D robot and Cartpole model results reported
directly from Bouissou et al.


### Ebola model (run from model file)

~~~
$ cat ebola-model.sys
$ ../polynomialFormUncertaintyPropagation  -n 25 -a ebola-model.sys 
~~~

### Honeybee model (run from model file)

~~~
$ cat honeybee.sys
$ ../polynomialFormUncertaintyPropagation -n 25 -a honeybee.sys 
~~~

### Coupled Vanderpol oscillator

~~~
$ cat coupledVanderpol-3.sys
$ ../polynomialFormUncertaintyPropagation -n 15 -a coupledVanderPolOscillators-3.sys 
~~~

### Laub Loomis Model

~~~
$ cat laub-loomis.sys
$ ../polynomialFormUncertaintyPropagation -n 25 -a laub-loomis.sys 
~~~

### Lattice 10x1 Model

~~~
$ cat 1d-lattice-10.sys
$ ../polynomialFormUncertaintyPropagation -n 10 -a  1d-lattice-10.sys 
~~~

# Running Simulations

~~~
cd ../python-simulation
~~~

Each python file corresponds to the model of the same name

~~~
python3 rimlessWheel.py
python3 roboticArm.py
python3 cartpole.py
python3 ebola-model.py
python3 honeybee.py
python3 coupledVanderPolOscillators-3.py
python3 laub-loomis.py
python3 1d-lattice-10.py
~~~



