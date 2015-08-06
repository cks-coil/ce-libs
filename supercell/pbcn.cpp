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

void PbcnSupercell::calcSpaceGroupSymOpMatrices(void){
    SMatrixXi alphaB, betaC, gammaN;
    SMatrixXi identitySMatrix(getNumPositions(), getNumPositions());
    identitySMatrix.setIdentity();
    alphaB = getSymOpMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0.5,0), Vector3d(0.25,0,0)) );
    betaC = getSymOpMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0,0.5), Vector3d(0,0.5,0)) );
    gammaN = getSymOpMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0.5,0.5,0), Vector3d(0,0,0.25))) ;
    spaceGroupSymOpMatrices.push_back( identitySMatrix );
    spaceGroupSymOpMatrices.push_back( alphaB );
    spaceGroupSymOpMatrices.push_back( betaC );
    spaceGroupSymOpMatrices.push_back( gammaN );
    spaceGroupSymOpMatrices.push_back( alphaB * betaC );
    spaceGroupSymOpMatrices.push_back( betaC * gammaN );
    spaceGroupSymOpMatrices.push_back( gammaN * alphaB );
    spaceGroupSymOpMatrices.push_back( alphaB * betaC * gammaN );
}
