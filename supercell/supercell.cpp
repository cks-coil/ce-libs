#include "supercell.hpp"
#include <algorithm>
#include <math.h>
#include <iostream>

#define ALLOWABLE_ERROR pow(10.0, -7)

using namespace Eigen;
using namespace std;

Supercell::Supercell(void){
    cellSize = Vector3i::Zero();
    atomicPos = Vector3d::Zero();
    crystalAxisMatrix = Matrix3d::Zero();
    orthogonalPositions.push_back( Vector3d::Zero() );
    fractionalPositions.push_back( Vector3d::Zero() );
    symmetriMatrices.push_back( MatrixXi::Zero(1,1) );
}

Supercell::~Supercell(void){
    fractionalPositions.clear();
    orthogonalPositions.clear();
    symmetriMatrices.clear();
}

void Supercell::setCellSize(Vector3i cellSize){
    this->cellSize = cellSize;
    updateVariables();
}

void Supercell::setAtomicPos(Vector3d atomicPos){
    this->atomicPos = atomicPos;
    updateVariables();
}

void Supercell::setCrystalAxis(Vector3d a, Vector3d b, Vector3d c){
    crystalAxisMatrix << a,b,c;
    updateVariables();
}


int Supercell::getNumAtoms(void){
    return fractionalPositions.size();
}

int Supercell::getNumUnitCellAtoms(void){
    return fractionalPositions.size() / ( cellSize(0) * cellSize(1) * cellSize(2) );
}

Vector3i Supercell::getCellSize(void){
    return cellSize;
}

Matrix3d Supercell::getCrystalAxisMatrix(void){
    return crystalAxisMatrix;
}

Vector3d Supercell::getFractionalPos(int index){
    return fractionalPositions[index];
}

Vector3d Supercell::getOrthogonalPos(int index){
    return orthogonalPositions[index];
}

vector<Vector3d> Supercell::getFractionalPositions(void){
    return fractionalPositions;
}
vector<Vector3d> Supercell::getOrthogonalPositions(void){
    return orthogonalPositions;
}
vector<MatrixXi> Supercell::getSymmetryMatrices(void){
    return symmetriMatrices;
}

void Supercell::updateVariables(void){
    if( cellSize == Vector3i::Zero() ) return;
    if( crystalAxisMatrix == Matrix3d::Zero() ) return;
    
    fractionalPositions.clear();
    calcFractionalPositions();
    periodicBoundaryCondition( fractionalPositions );

    orthogonalPositions.clear();
    calcOrthogonalPositions();

    symmetriMatrices.clear();
    calcSymmetryMatrices();
}

void Supercell::calcFractionalPositions(void){}
void Supercell::calcSymmetryMatrices(void){}


void Supercell::calcOrthogonalPositions(void){
    for(auto pos : fractionalPositions){
        orthogonalPositions.push_back( crystalAxisMatrix * pos ) ;
    }
}

void Supercell::periodicBoundaryCondition(vector<Eigen::Vector3d> &positions){
    for(auto &pos : positions){
        for(int i=0; i<3; i++){
            while(pos(i) >= (double)cellSize(i) ) pos(i) -= (double)cellSize(i);
            while(pos(i) < 0.0) pos(i) += (double)cellSize(i);
        }
    }
}
vector<Vector3d> Supercell::glideReflection(vector<Vector3d> positions, Vector3d transVector, Vector3d reflectionPos){
    vector<Vector3d> newPositions;
    Matrix3d reflectionMatrix = Matrix3d::Zero();
    for(int i=0;i<3;i++){
        if(reflectionPos(i) != 0) reflectionMatrix(i,i) = -1;
        else reflectionMatrix(i,i)=1;
    }
    for(auto pos: positions){
        Vector3d newPos;
        newPos = reflectionMatrix * pos + 2*reflectionPos + transVector;
        newPositions.push_back(newPos);
    }
    periodicBoundaryCondition(newPositions);
    return newPositions;
}

MatrixXi Supercell::calcSymmetryMatrix(vector<Vector3d> arr1, vector<Vector3d> arr2){
    MatrixXi symmetriMatrix = MatrixXi::Zero(getNumAtoms(), getNumAtoms());
    for(int i=0; i<(int)arr1.size(); i++){
        Vector3d v1 = arr1[i];
        int j = find_if(arr2.begin(), arr2.end(), [v1](const Vector3d &v2){ return (v1-v2).norm() < ALLOWABLE_ERROR; }) - arr2.begin();
        if ( j >= getNumAtoms() ){
            cerr << "ERROR: Faild to Find Pair (calcSymmetryMatrix): " << v1.transpose() << endl;
            for(auto v2: arr2) cout << v2.transpose() << endl;
            exit(1);
        }
        symmetriMatrix(j ,i) = 1;
    }
    return symmetriMatrix;
}
