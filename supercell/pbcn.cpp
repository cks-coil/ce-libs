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

void PbcnSupercell::calcSymmetryMatrices(void){}
