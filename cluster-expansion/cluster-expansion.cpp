#include "cluster-expansion.hpp"
#include <set>
#include <iostream>

using namespace std;
using namespace Eigen;

ClusterExpansion::ClusterExpansion(void){
    supercell = nullptr;
}

ClusterExpansion::~ClusterExpansion(void){
    effectiveClusters.clear();
    expandedClusters.clear();
}

void ClusterExpansion::setSupercell(const Supercell *supercell){
    this->supercell=supercell;
}

void ClusterExpansion::setEffectiveClusters(vector<Eigen::SVectorXi> effectiveClusters){
    this->effectiveClusters = effectiveClusters;
}

void ClusterExpansion::setEffectiveClusterInteractions(VectorXd effectiveClusterInteractions){
    this->effectiveClusterInteractions = effectiveClusterInteractions;
}

void ClusterExpansion::expandClusters(void){
    expandedClusters.clear();
    for(auto effectiveCluster : effectiveClusters){
        auto less = [](const SVectorXi &obj1, const SVectorXi &obj2){return obj1<obj2;};
        auto equal = [](const SVectorXi &obj1, const SVectorXi &obj2){return obj1==obj2;};
        vector<SVectorXi> clusters;
        for(auto symOpMatrix : supercell->getSymOpMatrices()){
            clusters.push_back( symOpMatrix * effectiveCluster );
        }
        sort(clusters.begin(),clusters.end(), less);
        clusters.erase(unique(clusters.begin(),clusters.end(), equal),clusters.end());
        expandedClusters.push_back(clusters);
    }
}

int ClusterExpansion::getNumEffectiveClusters(void) const{
    return effectiveClusters.size();
}

double ClusterExpansion::getEnergy(VectorXi configuration) const{
    return effectiveClusterInteractions.transpose() * getClusterCountVector(configuration).cast<double>();
}

VectorXi ClusterExpansion::getClusterCountVector(VectorXi configuration) const{
    VectorXi clusterCountVector = VectorXi::Zero(getNumEffectiveClusters());
    configuration -= VectorXi::Ones(supercell->getNumPositions());
    for(int i=0; i<getNumEffectiveClusters(); i++){
        for(auto cluster: expandedClusters[i]){
            if( (int)( - cluster.dot(configuration)) % 2 == 0 ) clusterCountVector(i) += 1;
            else clusterCountVector(i) -= 1;
        }
    }
    return clusterCountVector;
}



VectorXi getConvertedConfiguration(const VectorXi configuration, const Supercell &source, const Supercell &dest){
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

