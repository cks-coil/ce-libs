#include "clusters.hpp"
#include <iostream>
#include <algorithm>
#include <set>

using namespace std;
using namespace Eigen;

bool operator< (const SVectorXi &obj1, const SVectorXi &obj2){
    SVectorXi diff = obj1- obj2;
    diff.prune(0);
    if ( diff.nonZeros()==0 ) return false;
    SVectorXi::InnerIterator it(diff);
    return (it.value() < 0);
}

Clusters::Clusters(void){
    tgt = nullptr;
    maxDistance = 0;
    maxNum = 0;
}

Clusters::~Clusters(void){
    uniqueClusters.clear();
}

void Clusters::setSupercell(Supercell *tgt){
    if(  ( tgt->getCellSize().prod() ) %2 == 0  ){
        cerr << "ERROR: CellSize of Supercell is NOT Odd ( " << tgt->getCellSize().transpose() << " )" << endl;
        exit(1);
    } 
    this->tgt = tgt;
}

void Clusters::setMaxDistance(double maxDistance){
    this->maxDistance = maxDistance;
}

void Clusters::setMaxNum(int maxNum){
    this->maxNum = maxNum;
}

vector<SVectorXi> Clusters::getUniqueClusters(void){
    return uniqueClusters;
}
int Clusters::getNumUniqueClusters(void){
    return uniqueClusters.size();
}

void Clusters::findUniqueClusters(void){
    if(tgt == nullptr || maxNum == 0 || maxDistance == 0) return;
    SVectorXi emptyCluster(tgt->getNumPositions());
    SVectorXi singleCluster(tgt->getNumPositions());
    vector<SVectorXi> preClusters;
    emptyCluster.setZero();
    singleCluster.setZero();
    singleCluster.insert( tgt->getSupercellIndex(0, tgt->getCellSize()/2) ) = 1;
    uniqueClusters.push_back(emptyCluster);
    uniqueClusters.push_back(singleCluster);
    preClusters.push_back(singleCluster);
    for(int n=2; n<= maxNum; n++){
        preClusters = getNextUniqueClusters(preClusters);
        for(auto cluster : preClusters) uniqueClusters.push_back(cluster);
    }
}


bool Clusters::isRange(SVectorXi cluster){
    cluster.prune(0);
    for(SVectorXi::InnerIterator it1(cluster); it1; ++it1){
        for(SVectorXi::InnerIterator it2=it1; it2; ++it2){
            Vector3d diffVector = tgt->getOrthogonalPos(it1.index()) - tgt->getOrthogonalPos(it2.index());
            if ( diffVector.norm() > maxDistance ) return false;
        }
    }
    return true;
}

bool Clusters::isUnique(SVectorXi cluster, const vector<SVectorXi> &refs){
    auto comp = [](const SVectorXi &obj1, const SVectorXi &obj2){return obj1<obj2;};
    set<SVectorXi, decltype(comp)> clusterSet(comp);
    for(auto matrix: tgt->getSymOpMatrices()){
        SVectorXi tmpCluster = matrix*cluster;
        if( tmpCluster.coeff( tgt->getSupercellIndex(0,tgt->getCellSize()/2) ) == 1 ){
            clusterSet.insert( tmpCluster );
        }
    }
    for(auto ref: refs){
        if( clusterSet.find(ref) != clusterSet.end() ) return false;
    }
    return true;
}


vector<SVectorXi> Clusters::getNextUniqueClusters(const vector<SVectorXi> &seeds){
    vector<SVectorXi> nextUniqueClusters;
    for(auto seed: seeds){
        for(int i=0; i<tgt->getNumPositions(); i++){
            SVectorXi tmpCluster = seed;;
            if( seed.coeff(i)!=0 ) continue;
            tmpCluster.coeffRef(i) = 1;
            if( isRange(tmpCluster) && isUnique(tmpCluster, nextUniqueClusters)) nextUniqueClusters.push_back(tmpCluster);
        }
    }
    return nextUniqueClusters;
}

