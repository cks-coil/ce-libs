#include "cluster-expansion.hpp"
#include <set>

using namespace std;
using namespace Eigen;

ClusterExpansion::ClusterExpansion(void){}

ClusterExpansion::~ClusterExpansion(void){
    uniqueClusters.clear();
    expandedClusters.clear();
}

void ClusterExpansion::setSupercell(const Supercell *supercell){
    this->supercell=supercell;
}

void ClusterExpansion::setUniqueClusters(vector<Eigen::SVectorXi> uniqueClusters){
    this->uniqueClusters = uniqueClusters;
}

void ClusterExpansion::expandClusters(void){
    expandedClusters.clear();
    for(auto uniqueCluster : uniqueClusters){
        auto less = [](const SVectorXi &obj1, const SVectorXi &obj2){return obj1<obj2;};
        auto equal = [](const SVectorXi &obj1, const SVectorXi &obj2){return obj1==obj2;};
        vector<SVectorXi> clusters;
        for(auto symOpMatrix : supercell->getSymOpMatrices()){
            clusters.push_back( symOpMatrix * uniqueCluster );
        }
        sort(clusters.begin(),clusters.end(), less);
        clusters.erase(unique(clusters.begin(),clusters.end(), equal),clusters.end());
        expandedClusters.push_back(clusters);
    }
}
