#include "pbcn.hpp"

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

void PbcnSupercell::calcSymmetryMatrices(void){
    int maxN = getCellSize().maxCoeff() * 2;
    MatrixXi alphaBMatrix, betaCMatrix, gammaNMatrix;
    alphaBMatrix = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0.5,0), Vector3d(0.25,0,0)) );
    betaCMatrix = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0,0.5), Vector3d(0,0.5,0)) );
    gammaNMatrix = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0.5,0.5,0), Vector3d(0,0,0.25))) ;

    for( auto alpha : getPoweredMatrices(alphaBMatrix, maxN) ){
        for( auto beta : getPoweredMatrices(betaCMatrix, maxN) ){
            for( auto gamma : getPoweredMatrices(gammaNMatrix, maxN) ){
                symmetriMatrices.push_back( alpha * beta * gamma );
            }
        }
    }
}
