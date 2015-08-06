#include "supercell.hpp"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>

#define ALLOWABLE_ERROR pow(10.0, -7)

using namespace Eigen;
using namespace std;


Supercell::Supercell(void){
    cellSize = Vector3i::Zero();
    atomicPos = Vector3d::Zero();
    crystalAxisMatrix = Matrix3d::Zero();
}

Supercell::~Supercell(void){
    fractionalPositions.clear();
    orthogonalPositions.clear();
    symmetryMatrices.clear();
    unitCellFractionalPositions.clear();
    spaceGroupSymmetryMatrices.clear();
}

void Supercell::setCellSize(Vector3i cellSize){
    this->cellSize = cellSize;
}

void Supercell::setAtomicPos(Vector3d atomicPos){
    this->atomicPos = atomicPos;
}

void Supercell::setCrystalAxis(Vector3d a, Vector3d b, Vector3d c){
    crystalAxisMatrix << a,b,c;
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

Vector3d Supercell::getFractionalPos(int supercellIndex){
    return fractionalPositions[supercellIndex];
}
Vector3d Supercell::getFractionalPos(int unitCellIndex, Vector3i cell){
    return getFractionalPos( getSupercellIndex(unitCellIndex, cell) );
}

Vector3d Supercell::getOrthogonalPos(int supercellIndex){
    return orthogonalPositions[supercellIndex];
}
Vector3d Supercell::getOrthogonalPos(int unitCellIndex, Vector3i cell){
    return getOrthogonalPos( getSupercellIndex(unitCellIndex, cell) );
}

vector<Vector3d> Supercell::getFractionalPositions(void){
    return fractionalPositions;
}
vector<Vector3d> Supercell::getOrthogonalPositions(void){
    return orthogonalPositions;
}
vector<SMatrixXi> Supercell::getSymmetryMatrices(void){
    return symmetryMatrices;
}

int Supercell::getSupercellIndex(int unitCellIndex, Vector3i cell){
    int supercellIndex=0;
    supercellIndex += cell(0) * cellSize(1) * cellSize(2) * getNumUnitCellPositions();
    supercellIndex += cell(1) * cellSize(2) * getNumUnitCellPositions();
    supercellIndex += cell(2) * getNumUnitCellPositions();
    supercellIndex += unitCellIndex;
    return supercellIndex;
}


void Supercell::calcPositions(void){
    if( cellSize == Vector3i::Zero() ) return;
    if( crystalAxisMatrix == Matrix3d::Zero() ) return;
    
    calcFractionalPositions();
    calcOrthogonalPositions();
}

void Supercell::calcSymmetryMatrices(void){
    if( fractionalPositions.empty() ) return;

    spaceGroupSymmetryMatrices.clear();
    calcSpaceGroupSymmetryMatrices();

    symmetryMatrices.clear();
    for(auto slide : getSlideMatrices() ){
        for(auto spaceGroup: spaceGroupSymmetryMatrices ){
            symmetryMatrices.push_back(slide * spaceGroup);
        }
    }
}

void Supercell::checkSymmetryMatrices(void){
    if( (int)symmetryMatrices.size() != getNumPositions() ){
        cerr << "ERROR: Number of Symmetry Matrices (" << symmetryMatrices.size() << ") is NOT Equal to Number of Positions (" << getNumPositions() << ")" << endl;
        exit(1);
    }

    for(auto m: symmetryMatrices){
        if( m * VectorXi::Ones(getNumPositions()) != VectorXi::Ones(getNumPositions()) ){
            cerr << "ERROR: Wrong Symmetriy Matrix" << endl;
            cerr << m << endl;
            exit(1);
        }
    }
    
    SMatrixXi zeroMatrix(getNumPositions(), getNumPositions());
    zeroMatrix.setZero();
    MatrixXi sum = accumulate(symmetryMatrices.begin(), symmetryMatrices.end(), zeroMatrix);
    if(sum != MatrixXi::Ones(getNumPositions(), getNumPositions())){
        cerr << "ERROR: Wrong Symmetriy Matrices" << endl;
        exit(1);
    }

    cerr << "PASS: Symmetry Matrices Test" << endl;
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

void Supercell::calcOrthogonalPositions(void){
    orthogonalPositions.clear();
    for(auto pos : fractionalPositions){
        orthogonalPositions.push_back( crystalAxisMatrix * pos ) ;
    }
}

void Supercell::calcUnitCellFractionalPositions(void){}
void Supercell::calcSpaceGroupSymmetryMatrices(void){}

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

SMatrixXi Supercell::getSymmetryMatrix(vector<Vector3d> arr1, vector<Vector3d> arr2){
    SMatrixXi symmetryMatrix(getNumPositions(), getNumPositions());
    for(int j=0; j<(int)arr1.size(); j++){
        Vector3d v1 = arr1[j];
        int i = find_if(arr2.begin(), arr2.end(), [v1](const Vector3d &v2){ return (v1-v2).norm() < ALLOWABLE_ERROR; }) - arr2.begin();
        if ( i >= getNumPositions() ){
            cerr << "ERROR: Faild to Find Pair (getSymmetryMatrix): " << v1.transpose() << endl;
            for(auto v2: arr2) cerr << v2.transpose() << endl;
            exit(1);
        }
        symmetryMatrix.insert(i ,j) = 1;
    }
    return symmetryMatrix;
}

vector<SMatrixXi> Supercell::getPoweredMatrices(SMatrixXi matrix, int maxN){
    vector<SMatrixXi> matrices;
    SMatrixXi identityMatrix(getNumPositions(), getNumPositions());
    identityMatrix.setIdentity();
    matrices.push_back( identityMatrix );
    for(int i=1; i<=maxN; i++){
        SMatrixXi tmpMatrix = matrices.back() * matrix;
        if( tmpMatrix == identityMatrix ){
            return matrices;
        }
        matrices.push_back(tmpMatrix);
    }
    cerr << "ERROR: Faild to Find A^N==E (N<=" << maxN << ") (getPowerdMatrices)" << endl;
    cerr << matrix << endl;
    exit(1);
}

vector<SMatrixXi> Supercell::getSlideMatrices(void){
    vector<SMatrixXi> slideMatrices;
    SMatrixXi aSlideMatrix, bSlideMatrix, cSlideMatrix;
    aSlideMatrix = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(1,0,0), Vector3d(0,0,0)) );
    bSlideMatrix = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,1,0), Vector3d(0,0,0)) );
    cSlideMatrix = getSymmetryMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0,1), Vector3d(0,0,0)) );
    for(auto a: getPoweredMatrices(aSlideMatrix, getCellSize().maxCoeff()) ){
        for(auto b: getPoweredMatrices(bSlideMatrix, getCellSize().maxCoeff()) ){
            for(auto c: getPoweredMatrices(cSlideMatrix, getCellSize().maxCoeff()) ){
                slideMatrices.push_back(a*b*c);
            }
        }
    }
    return slideMatrices;
}

bool operator == (const SMatrixXi &obj1, const SMatrixXi &obj2){
    Eigen::SMatrixXi tmpMatrix = obj1- obj2;
    return tmpMatrix.norm() == 0;
}
