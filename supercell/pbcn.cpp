#include "pbcn.hpp"
#include <iostream>

using namespace std;
using namespace Eigen;

PbcnSupercell::PbcnSupercell(void){}

void PbcnSupercell::calcUnitCellFractionalPositions(void){
    double x,y,z;
    x = atomicPos.x();
    y = atomicPos.y();
    z = atomicPos.z();
    unitCellFractionalPositions.push_back( Vector3d(x,y,z) );
    unitCellFractionalPositions.push_back( Vector3d(-x, y, -z+0.5) );
    unitCellFractionalPositions.push_back( Vector3d(-x+0.5, y+0.5, z) );
    unitCellFractionalPositions.push_back( Vector3d(x+0.5, y+0.5, -z+0.5) );
    unitCellFractionalPositions.push_back( Vector3d(-x+0.5, -y+0.5, z+0.5) );
    unitCellFractionalPositions.push_back( Vector3d(x+0.5, -y+0.5, -z) );
    unitCellFractionalPositions.push_back( Vector3d(x, -y, z+0.5) );
    unitCellFractionalPositions.push_back( Vector3d(-x, -y, -z) );
}

void PbcnSupercell::calcSpaceGroupSymmetryMatrices(void){
    SMatrixXi alphaB, betaC, gammaN;
    SMatrixXi identitySMatrix(getNumPositions(), getNumPositions());
    identitySMatrix.setIdentity();
    alphaB = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0.5,0), Vector3d(0.25,0,0)) );
    betaC = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0,0.5), Vector3d(0,0.5,0)) );
    gammaN = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0.5,0.5,0), Vector3d(0,0,0.25))) ;
    spaceGroupSymmetryMatrices.push_back( identitySMatrix );
    spaceGroupSymmetryMatrices.push_back( alphaB );
    spaceGroupSymmetryMatrices.push_back( betaC );
    spaceGroupSymmetryMatrices.push_back( gammaN );
    spaceGroupSymmetryMatrices.push_back( alphaB * betaC );
    spaceGroupSymmetryMatrices.push_back( betaC * gammaN );
    spaceGroupSymmetryMatrices.push_back( gammaN * alphaB );
    spaceGroupSymmetryMatrices.push_back( alphaB * betaC * gammaN );
}
