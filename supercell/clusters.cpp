#include <iostream>
#include "clusters.hpp"

using namespace std;
using namespace Eigen;

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

void findClusters(void){
}
