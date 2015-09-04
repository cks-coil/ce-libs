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
    numUnitCellPositions = 0;
}

Supercell::~Supercell(void){
    fractionalPositions.clear();
    orthogonalPositions.clear();
    symOpMatrices.clear();
    unitCellFractionalPositions.clear();
    spaceGroupSymOpMatrices.clear();
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

int Supercell::getNumPositions(void) const{
    return getNumUnitCellPositions() * cellSize.prod();
}

int Supercell::getNumUnitCellPositions(void) const{
    if( numUnitCellPositions ) return numUnitCellPositions;
    else return fractionalPositions.size() / cellSize.prod();
}

Vector3i Supercell::getCellSize(void) const{
    return cellSize;
}

Matrix3d Supercell::getCrystalAxisMatrix(void) const{
    return crystalAxisMatrix;
}

Vector3d Supercell::getFractionalPos(int supercellIndex) const{
    return fractionalPositions[supercellIndex];
}
Vector3d Supercell::getFractionalPos(int unitCellIndex, Vector3i cellPos) const{
    return getFractionalPos( getSupercellIndex(unitCellIndex, cellPos) );
}

Vector3d Supercell::getOrthogonalPos(int supercellIndex) const{
    return orthogonalPositions[supercellIndex];
}
Vector3d Supercell::getOrthogonalPos(int unitCellIndex, Vector3i cellPos) const{
    return getOrthogonalPos( getSupercellIndex(unitCellIndex, cellPos) );
}

vector<Vector3d> Supercell::getFractionalPositions(void) const{
    return fractionalPositions;
}
vector<Vector3d> Supercell::getOrthogonalPositions(void) const{
    return orthogonalPositions;
}
vector<SMatrixXi> Supercell::getSymOpMatrices(void) const{
    return symOpMatrices;
}

int Supercell::getSupercellIndex(int unitCellIndex, Vector3i cellPos) const{
    int supercellIndex=0;
    supercellIndex += cellPos(0) * cellSize(1) * cellSize(2) * getNumUnitCellPositions();
    supercellIndex += cellPos(1) * cellSize(2) * getNumUnitCellPositions();
    supercellIndex += cellPos(2) * getNumUnitCellPositions();
    supercellIndex += unitCellIndex;
    return supercellIndex;
}

int Supercell::getUnitCellIndex(int supercellIndex) const{
    return supercellIndex % getNumUnitCellPositions();
}

Vector3i Supercell::getCellPos(int supercellIndex) const{
    Vector3i cellPos;
    for(int i=0;i<3;i++){
        int tmp = getNumUnitCellPositions();
        for(int j=i+1;j<3;j++) tmp*=cellSize(j);
        cellPos(i) = supercellIndex / tmp;
        supercellIndex %= tmp;
    }
    return cellPos;
}

void Supercell::calcPositions(void){
    if( cellSize == Vector3i::Zero() ) return;
    if( crystalAxisMatrix == Matrix3d::Zero() ) return;
    
    calcFractionalPositions();
    calcOrthogonalPositions();
}

void Supercell::calcSymOpMatrices(void){
    if( fractionalPositions.empty() ) return;

    spaceGroupSymOpMatrices.clear();
    calcSpaceGroupSymOpMatrices();

    symOpMatrices.clear();
    for(auto slide : getSlideMatrices() ){
        for(auto spaceGroup: spaceGroupSymOpMatrices ){
            symOpMatrices.push_back(slide * spaceGroup);
        }
    }
}

void Supercell::checkSymOpMatrices(void){
    if( (int)symOpMatrices.size() != getNumPositions() ){
        cerr << "ERROR: Number of SymOp Matrices (" << symOpMatrices.size() << ") is NOT Equal to Number of Positions (" << getNumPositions() << ")" << endl;
        exit(1);
    }

    for(auto m: symOpMatrices){
        if( m * VectorXi::Ones(getNumPositions()) != VectorXi::Ones(getNumPositions()) ){
            cerr << "ERROR: Wrong Symmetriy Matrix" << endl;
            cerr << m << endl;
            exit(1);
        }
    }
    
    MatrixXi sum = MatrixXi::Zero(getNumPositions(), getNumPositions());
    for(auto matrix:symOpMatrices){
        /*
          summation
          eigen3's defalut MatrixXi += SMatrixXi is too slow
         */
        for (int k=0; k < matrix.outerSize(); ++k){
            for (SparseMatrix<int>::InnerIterator it(matrix,k); it; ++it){
                sum(it.row(), it.col()) += it.value();
            }
        }
    }
    if(sum != MatrixXi::Ones(getNumPositions(), getNumPositions())){
        cerr << "ERROR: Wrong Symmetriy Matrices" << endl;
        exit(1);
    }

    cerr << "PASS: SymOp Matrices Test" << endl;
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
void Supercell::calcSpaceGroupSymOpMatrices(void){}

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

SMatrixXi Supercell::getSymOpMatrix(vector<Vector3d> arr1, vector<Vector3d> arr2){
    SMatrixXi symOpMatrix(getNumPositions(), getNumPositions());
    for(int j=0; j<(int)arr1.size(); j++){
        Vector3d v1 = arr1[j];
        int i = find_if(arr2.begin(), arr2.end(), [v1](const Vector3d &v2){ return (v1-v2).norm() < ALLOWABLE_ERROR; }) - arr2.begin();
        if ( i >= getNumPositions() ){
            cerr << "ERROR: Faild to Find Pair (getSymOpMatrix): " << v1.transpose() << endl;
            for(auto v2: arr2) cerr << v2.transpose() << endl;
            exit(1);
        }
        symOpMatrix.insert(i ,j) = 1;
    }
    return symOpMatrix;
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
    aSlideMatrix = getSymOpMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(1,0,0), Vector3d(0,0,0)) );
    bSlideMatrix = getSymOpMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,1,0), Vector3d(0,0,0)) );
    cSlideMatrix = getSymOpMatrix( fractionalPositions, glideReflection(fractionalPositions, Vector3d(0,0,1), Vector3d(0,0,0)) );
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
