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
    unitCellFractionalPositions.push_back( Vector3d::Zero() );
}

Supercell::~Supercell(void){
    fractionalPositions.clear();
    orthogonalPositions.clear();
    symmetriMatrices.clear();
    unitCellFractionalPositions.clear();
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


int Supercell::getNumPositions(void){
    return fractionalPositions.size();
}

int Supercell::getNumUnitCellPositions(void){
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
    
    calcFractionalPositions();
    calcOrthogonalPositions();

    symmetriMatrices.clear();
    calcSymmetryMatrices();
}

void Supercell::calcFractionalPositions(void){
    unitCellFractionalPositions.clear();
    calcUnitCellFractionalPositions();
    periodicBoundaryCondition(unitCellFractionalPositions, Vector3i::Ones());

    fractionalPositions.clear();
    for(int a=0; a<cellSize(0); a++){
        for(int b=0; b<cellSize(1); b++){
            for(int c=0; c<cellSize(2); c++){
                for(auto pos: unitCellFractionalPositions){
                    fractionalPositions.push_back( pos + Vector3d(a,b,c) );
                }
            }
        }
    }
}

void Supercell::calcUnitCellFractionalPositions(void){}
void Supercell::calcSymmetryMatrices(void){}


void Supercell::calcOrthogonalPositions(void){
    orthogonalPositions.clear();
    for(auto pos : fractionalPositions){
        orthogonalPositions.push_back( crystalAxisMatrix * pos ) ;
    }
}

void Supercell::periodicBoundaryCondition(vector<Vector3d> &positions, Vector3i boundary){
    for(auto &pos : positions){
        for(int i=0; i<3; i++){
            while(pos(i) >= (double)boundary(i) ) pos(i) -= (double)boundary(i);
            while(pos(i) < 0.0) pos(i) += (double)boundary(i);
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
    periodicBoundaryCondition(newPositions, cellSize);
    return newPositions;
}

MatrixXi Supercell::getSymmetryMatrix(vector<Vector3d> arr1, vector<Vector3d> arr2){
    MatrixXi symmetriMatrix = MatrixXi::Zero(getNumPositions(), getNumPositions());
    for(int i=0; i<(int)arr1.size(); i++){
        Vector3d v1 = arr1[i];
        int j = find_if(arr2.begin(), arr2.end(), [v1](const Vector3d &v2){ return (v1-v2).norm() < ALLOWABLE_ERROR; }) - arr2.begin();
        if ( j >= getNumPositions() ){
            cerr << "ERROR: Faild to Find Pair (getSymmetryMatrix): " << v1.transpose() << endl;
            for(auto v2: arr2) cout << v2.transpose() << endl;
            exit(1);
        }
        symmetriMatrix(j ,i) = 1;
    }
    return symmetriMatrix;
}

vector<MatrixXi> Supercell::getPoweredMatrices(MatrixXi matrix, int maxN){
    vector<MatrixXi> matrices;
    matrices.push_back( MatrixXi::Identity(getNumPositions(), getNumPositions()) );
    for(int i=1; i<=maxN; i++){
        MatrixXi tmpMatrix = matrices.back() * matrix;
        if( tmpMatrix == matrices.front() ){
            return matrices;
        }
        matrices.push_back(tmpMatrix);
    }
    cerr << "ERROR: Faild to Find A^N==E (N<" << maxN << ") (getPowerdMatrices)" << endl;
    cerr << matrix << endl;
    exit(1);
}
