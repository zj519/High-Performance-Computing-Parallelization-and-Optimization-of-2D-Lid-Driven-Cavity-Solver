#include <iostream>
#include "SolverCG.h"
#include <cblas.h>
using namespace std;

#define BOOST_TEST_MODULE PopularTestSuite
#include <boost/test/included/unit_test.hpp>

//define the 3x3 matrix and expected results
int testNx=3;
int testNy=3;
int testNo=testNx*testNy;
double testDx=1.0;
double testDy=1.0;
double tol=1e-6;
//input array
double testIn[9]={1,2,3,4,5,6,7,8,9};
//expected output
double expectedApplyOperator[9]={0,0,0,0,0,0,0,0,0};
double expectedPrecondition[9]={1,2,3,4,1.25,6,7,8,9};
double expectedImposeBC[9]={0,0,0,0,5,0,0,0,0};
double testOut[9]={0};
    
//test the correct initialization of the SolverCG class
BOOST_AUTO_TEST_CASE(testConstructor){
    SolverCG solver(testNx,testNy,testDx,testDy);
    
    BOOST_CHECK_EQUAL(solver.Nx, testNx);
    BOOST_CHECK_EQUAL(solver.Ny, testNy);
    BOOST_CHECK_EQUAL(solver.dx, testDx);
    BOOST_CHECK_EQUAL(solver.dy, testDy);
}


//test the ApplyOperator function outputs expected value 
BOOST_AUTO_TEST_CASE(testApplyOperator){
    SolverCG solver(testNx,testNy,testDx,testDy);
    
    solver.ApplyOperator(testIn, testOut);
    
    //check each array member aginst the expected
    for (int i=0;i<testNo;i++){
        BOOST_CHECK_CLOSE(testOut[i],expectedApplyOperator[i], tol);
        
    }
}

//test whether Precondition function outputs expected value
BOOST_AUTO_TEST_CASE(testPrecondition){
    SolverCG solver(testNx,testNy,testDx,testDy);
    
    solver.Precondition(testIn, testOut);
    
    for (int i=0;i<testNo; i++){
        BOOST_CHECK_CLOSE(testOut[i],expectedPrecondition[i], tol);
    }
}

//test whether ImposeBC function outputs expected value
BOOST_AUTO_TEST_CASE(testImposeBC){
    SolverCG solver(testNx,testNy,testDx,testDy);
    solver.ImposeBC(testIn);
    
    for (int i=0;i<testNo; i++){
        BOOST_CHECK_CLOSE(testIn[i],expectedImposeBC[i], tol);
    }
}






