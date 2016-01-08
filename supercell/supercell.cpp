#include <Eigen/Dense>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <numeric>
#include "supercell.hpp"

#define ALLOWABLE_ERROR pow(10.0, -7)

using namespace Eigen;
using namespace std;


Supercell::Supercell(void){
    cellSize = Vector3i::Zero();
    atomicPos = Vector3d::Zero();
    crystalAxisMatrix = Matrix3d::Zero();
    numPositions = 0;
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
    numPositions = cellSize.prod() * numUnitCellPositions;
}

void Supercell::setAtomicPos(Vector3d atomicPos){
    this->atomicPos = atomicPos;
}

void Supercell::setCrystalAxis(Vector3d a, Vector3d b, Vector3d c){
    crystalAxisMatrix << a,b,c;
}

int Supercell::getNumPositions(void) const{
    return numPositions;
}

int Supercell::getNumUnitCellPositions(void) const{
    return numUnitCellPositions;
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

const vector<Vector3d>& Supercell::getFractionalPositions(void) const{
    return fractionalPositions;
}
const vector<Vector3d>& Supercell::getOrthogonalPositions(void) const{
    return orthogonalPositions;
}
const vector<SMatrixXi>& Supercell::getSymOpMatrices(void) const{
    return symOpMatrices;
}
const vector<SMatrixXi>& Supercell::getSpaceGroupSymOpMatrices(void) const{
    return spaceGroupSymOpMatrices;
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

double Supercell::getVolume(void) const{
    return crystalAxisMatrix.col(0).dot( crystalAxisMatrix.col(1).cross(crystalAxisMatrix.col(2)) ) * cellSize.prod();;
}

void Supercell::calcPositions(void){
    if( cellSize == Vector3i::Zero() ) return;
    if( crystalAxisMatrix == Matrix3d::Zero() ) return;
    
    calcFractionalPositions();
    calcOrthogonalPositions();

    numPositions = fractionalPositions.size();
    numUnitCellPositions = unitCellFractionalPositions.size();
}

void Supercell::calcSymOpMatrices(void){
    if( fractionalPositions.empty() ) return;

    spaceGroupSymOpMatrices.clear();
    calcSpaceGroupSymOpMatrices();

    symOpMatrices.clear();
    symOpMatrices.reserve(numPositions);
    for(auto slide : getSlideMatrices() ){
        for(auto spaceGroup: spaceGroupSymOpMatrices ){
            symOpMatrices.push_back( (slide * spaceGroup).pruned() );
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


VectorXi getConvertedConfiguration(VectorXi configuration, const Supercell &source, const Supercell &dest){
    VectorXi convertedConfiguration(dest.getNumPositions());
    Vector3i scale = (dest.getCellSize().array() / source.getCellSize().array() ).matrix();
    if( (source.getCellSize().array() * scale.array()).matrix() != dest.getCellSize()  ){
        cerr << "Error: Can NOT Convert Configuration" << endl;
        exit(1);
    }

    for(int i=0; i<scale(0); i++){
        for(int j=0; j<scale(1); j++){
            for(int k=0; k<scale(2); k++){
                Vector3i cellSlide;
                cellSlide = (Vector3i(i,j,k).array() * source.getCellSize().array()).matrix();
                for(int row=0; row<source.getNumPositions(); row++){
                    int unitCellIndex = source.getUnitCellIndex(row);
                    Vector3i cellPos = source.getCellPos(row) + cellSlide;
                    convertedConfiguration( dest.getSupercellIndex(unitCellIndex, cellPos ) ) = configuration(row);
                }
            }
        }
    }
    return convertedConfiguration;
}

VectorXi getUnitCellConfiguration(VectorXi configuration, const Supercell &supercell, Vector3i cellPos){
    VectorXi unitCellConf(supercell.getNumUnitCellPositions());
    for(int i=0;i<supercell.getNumUnitCellPositions();i++){
        int supercellIndex = supercell.getSupercellIndex(i,cellPos);
        unitCellConf(i) = configuration(supercellIndex);
    }
    return unitCellConf;
}

string getConfCombinedStr(VectorXi configuration){
    string str = "";
    for(int i=0; i<configuration.rows(); i++) str+= to_string(configuration(i));
    return str;
}


string getConfSplitedStr(VectorXi configuration, const Supercell &supercell){
    string str = "";
    Vector3i cellSize = supercell.getCellSize();
    str += to_string(cellSize(0)) + "," + to_string(cellSize(1)) + "," + to_string(cellSize(2));
    for(int i=0; i<configuration.rows(); i++){
        if(i%8==0){
            Vector3i cellPos = supercell.getCellPos(i);
            str  += " " + to_string(cellPos(0)) + "," + to_string(cellPos(1)) + "," + to_string(cellPos(2)) + "-";
        }
        str+= to_string(configuration(i));
    }
    return str;
}
