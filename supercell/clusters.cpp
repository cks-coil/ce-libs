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

Clusters::Clusters(void){}

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

void Clusters::findClusters(void){
    SVectorXi emptyCluster(tgt->getNumPositions());
    SVectorXi singleCluster(tgt->getNumPositions());
    vector<SVectorXi> preClusters;
    emptyCluster.setZero();
    singleCluster.setZero();
    singleCluster.insert( tgt->getSupercellIndex(0, tgt->getCellSize()/2) ) = 1;
    uniqueClusters.push_back(emptyCluster);
    uniqueClusters.push_back(singleCluster);
    preClusters.push_back(singleCluster);
    getNextClusters(preClusters);
}

bool Clusters::isRange(SVectorXi cluster){
    cluster.prune(0);
    for(SVectorXi::InnerIterator it1(cluster); it1; ++it1){
        for(SVectorXi::InnerIterator it2=it1; it2; ++it2){
            Vector3d diffVector;
            if ( diffVector.norm() > maxDistance ) return false;
        }
    }
    return true;
}
bool Clusters::isUnique(SVectorXi cluster, const vector<SVectorXi> &refs){
    vector<SVectorXi> clusters;
    for(auto matrix: tgt->getSymmetryMatrices()){
        SVectorXi tmpCluster = matrix*cluster;
        if( tmpCluster.coeff( tgt->getSupercellIndex(0,tgt->getCellSize()/2) ) == 1 ){
            clusters.push_back( tmpCluster );
        }
    }
    sort(clusters.begin(), clusters.end(),
         [](const SVectorXi &obj1, const SVectorXi &obj2){return obj1<obj2;} );
    return true;
}


vector<SVectorXi> Clusters::getNextClusters(const vector<SVectorXi> &seeds){
    vector<SVectorXi> nextClusters;
    for(auto seed: seeds){
        for(int i=0; i<tgt->getNumPositions(); i++){
            SVectorXi tmpCluster = seed;;
            tmpCluster.coeffRef(tgt->getNumPositions()-1)=1;
            if( seed.coeff(i)!=0 ) continue;
            tmpCluster.coeffRef(i) = 1;
            if( isRange(tmpCluster) && isUnique(tmpCluster, nextClusters)) nextClusters.push_back(tmpCluster);
        }
    }
    return nextClusters;
}

