cmake_minimum_required(VERSION 3.10)
project(polynomialFormUncertaintyPropagation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-includes)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g  -Wall  -O3 -Dregister= ")
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
include_directories(.)
include_directories(src)
include_directories(/usr/local/include)
link_directories(/usr/local/lib /usr/lib /usr/local/opt/llvm@9/lib)

include(FindMPFI)
if (NOT MPFI_FOUND)
    message(FATAL_ERROR "ERROR: MPFI not found")
endif(NOT MPFI_FOUND)

if (NOT MPFR_FOUND)
    message(FATAL_ERROR "ERROR: MPFR not found.")
endif(NOT MPFR_FOUND)

if (NOT glpk_FOUND)
    message(FATAL_ERROR "ERROR: glpk not found.")
endif(NOT glpk_FOUND)

find_package(MPFI REQUIRED)

find_package(BISON)
find_package(FLEX)

BISON_TARGET(MyParser src/modelParser.yacc ${CMAKE_SOURCE_DIR}/modelParser.yy.cc)
FLEX_TARGET(MyScanner src/modelLexer.lex ${CMAKE_SOURCE_DIR}/modelLexer.cc)

ADD_FLEX_BISON_DEPENDENCY(MyScanner MyParser)


add_executable(polynomialFormUncertaintyPropagation src/main.cpp  src/MultivariatePoly.hh
               src/MpfiWrapper.hh src/MpfiWrapper.cpp src/MultivariatePoly.cpp
               src/DistributionInfo.hh src/truncatedNormalCode.cpp src/truncatedNormalCode.hh
               ${BISON_MyParser_OUTPUTS}  ${FLEX_MyScanner_OUTPUTS}
               src/ExprAST.hh src/ExprAST.cpp src/SystemDescription.hh
               src/ModelParser.hh src/SystemDescription.cpp
               src/ModelParser.cpp src/Queries.hh src/StateAbstraction.hh
               src/affineForms/utilities.hpp src/affineForms/utilities.cpp
               src/affineForms/affineArithmeticClass.hpp src/affineForms/aaEnvironment.hpp
               src/affineForms/affineArithmeticClass.cpp src/affineForms/aaEnvironment.cpp
               src/ProbabilityQueryEvaluator.hh src/ProbabilityQueryEvaluator.cpp
               src/roboticArmModel.cpp src/rimlessWheelModel.cpp src/simpleTests.cpp
               src/cartPoleBalanceModel.cpp src/PolynomialFunctions.hh src/PolynomialFunctions.cpp src/affineForms/AffineArithmeticEvaluate.hh src/affineForms/AffineArithmeticEvaluate.cpp)

add_library(polynomialForm  src/MultivariatePoly.hh
        src/MpfiWrapper.hh src/MpfiWrapper.cpp src/MultivariatePoly.cpp
        src/DistributionInfo.hh src/truncatedNormalCode.cpp src/truncatedNormalCode.hh
        src/ExprAST.hh src/ExprAST.cpp src/SystemDescription.hh src/ModelParser.hh
        src/SystemDescription.cpp src/ModelParser.cpp src/Queries.hh
        src/StateAbstraction.hh src/ProbabilityQueryEvaluator.hh src/ProbabilityQueryEvaluator.cpp
        src/PolynomialFunctions.hh src/PolynomialFunctions.cpp src/libsetup.cpp
        )

install (TARGETS polynomialForm DESTINATION lib)
install (FILES src/MultivariatePoly.hh src/MpfiWrapper.hh src/DistributionInfo.hh
        src/truncatedNormalCode.hh src/ExprAST.hh src/SystemDescription.hh src/ModelParser.hh
        src/PolynomialFunctions.hh
        src/Queries.hh src/StateAbstraction.hh src/ProbabilityQueryEvaluator.hh DESTINATION include)

target_link_libraries(polynomialFormUncertaintyPropagation
        ${MPFR_LIBRARIES}
        ${MPFI_LIBRARIES}
        ${glpk_LIBRARIES})
