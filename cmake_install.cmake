# Install script for directory: /Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/libpolynomialForm.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libpolynomialForm.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libpolynomialForm.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libpolynomialForm.a")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/MultivariatePoly.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/MpfiWrapper.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/DistributionInfo.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/truncatedNormalCode.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/ExprAST.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/SystemDescription.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/ModelParser.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/PolynomialFunctions.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/Queries.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/StateAbstraction.hh"
    "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/src/ProbabilityQueryEvaluator.hh"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/CU-mac/Library/Mobile Documents/com~apple~CloudDocs/workspaceC/bayesianPolyForm/polynomialFormPropagationBayesian/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
