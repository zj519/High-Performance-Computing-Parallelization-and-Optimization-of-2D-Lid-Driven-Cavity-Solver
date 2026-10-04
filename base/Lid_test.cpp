#include "LidDrivenCavity.h"
#include "SolverCG.h"
#include <cmath>
#include <iostream>
#include "cblas.h"
using namespace std;

#define BOOST_TEST_MODULE PopulateTestSuite
#include <boost/test/included/unit_test.hpp>
#define IDX(I,J) ((J)*lid.Nx + (I))

//test the correct initialization of the LidDrivencavity class
BOOST_AUTO_TEST_SUITE(LidDrivenCavityTest)

BOOST_AUTO_TEST_CASE(Initialization){
    double tol=1e-6;
    
    LidDrivenCavity lid;
    lid.SetDomainSize(1.0,1.0);
    lid.SetGridSize(11,11);
    lid.SetTimeStep(0.01);
    lid.SetFinalTime(1.0);
    lid.SetReynoldsNumber(100.0);
    lid.Initialise();
    
    BOOST_CHECK_EQUAL(lid.Lx, 1.0);
    BOOST_CHECK_EQUAL(lid.Ly, 1.0);
    BOOST_CHECK_EQUAL(lid.Nx, 11);
    BOOST_CHECK_EQUAL(lid.Ny, 11);
    BOOST_CHECK_EQUAL(lid.dt, 0.01);
    BOOST_CHECK_EQUAL(lid.T, 1.0);
    BOOST_CHECK_EQUAL(lid.Re, 100.0);
    BOOST_CHECK_CLOSE(lid.nu, 0.01,tol);
    BOOST_CHECK_CLOSE(lid.dx, 0.1,tol);  
    BOOST_CHECK_CLOSE(lid.dx, 0.1,tol);  
    BOOST_CHECK_EQUAL(lid.Npts, 121);
    
    BOOST_CHECK(lid.v != nullptr);
    BOOST_CHECK(lid.s != nullptr);
    BOOST_CHECK(lid.tmp != nullptr);
    BOOST_CHECK(lid.cg != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()




