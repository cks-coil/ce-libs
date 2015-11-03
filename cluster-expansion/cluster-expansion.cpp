#include "cluster-expansion.hpp"
#include <set>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace Eigen;

ClusterExpansion::ClusterExpansion(void){
    supercell = nullptr;
}

ClusterExpansion::~ClusterExpansion(void){
    effectiveClusters.clear();
    expandedClusters.clear();
    mappedClusters.clear();
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
            SVectorXi configuration = symOpMatrix * effectiveCluster;
            configuration.prune(0);
            clusters.push_back( configuration );
        }
        sort(clusters.begin(),clusters.end(), less);
        clusters.erase(unique(clusters.begin(),clusters.end(), equal),clusters.end());
        expandedClusters.push_back(clusters);
    }
}

void ClusterExpansion::mapClusters(void){
    mappedClusters.clear();
    mappedClusters.resize(supercell->getNumPositions());
    for(int clusterIndex=0; clusterIndex<getNumEffectiveClusters(); clusterIndex++){
        for(auto cluster: expandedClusters[clusterIndex]){
            for(SVectorXi::InnerIterator it(cluster); it; ++it){
                mappedClusters[it.index()].push_back(make_pair(clusterIndex, cluster));
            }
        }
    }
}

int ClusterExpansion::getNumEffectiveClusters(void) const{
    return effectiveClusters.size();
}

double ClusterExpansion::getEnergy(VectorXi configuration) const{
    return effectiveClusterInteractions.transpose() * getClusterCountVector(configuration).cast<double>();
}

pair<double,VectorXi> ClusterExpansion::getEnergyInitial(VectorXi configuration) const{
    VectorXi clusterCountVector = getClusterCountVector(configuration);
    double energy = effectiveClusterInteractions.transpose() * clusterCountVector.cast<double>();
    return make_pair(energy, clusterCountVector);
}

pair<double,VectorXi> ClusterExpansion::getEnergyDifferential(const VectorXi &configuration, const VectorXi &oldConfiguration, const VectorXi &oldClusterCountVector, const vector<pair<int,int>> &changes) const{
    VectorXi clusterCountVector = getClusterCountVectorDifferential(configuration, oldConfiguration, oldClusterCountVector, changes);
    double energy = effectiveClusterInteractions.transpose() * clusterCountVector.cast<double>();
    return make_pair(energy, clusterCountVector);
}

VectorXi ClusterExpansion::getClusterCountVector(VectorXi configuration) const{
    VectorXi clusterCountVector = VectorXi::Zero(getNumEffectiveClusters());
    configuration -= VectorXi::Ones(supercell->getNumPositions());
    for(int i=0; i<getNumEffectiveClusters(); i++){
        for(auto cluster: expandedClusters[i]){
            // for empty cluster
            if(cluster.nonZeros()==0){
                clusterCountVector(i) = supercell->getNumPositions();
                break;
            }

            if( (int)( - cluster.dot(configuration)) % 2 == 0 ) clusterCountVector(i) += 1;
            else clusterCountVector(i) -= 1;
        }
    }
    return clusterCountVector;
}

const Supercell *ClusterExpansion::getSupercell(void) const{
    return supercell;
}

VectorXd ClusterExpansion::getEffectiveClusterInteractions(void) const{
    return effectiveClusterInteractions;
}

void ClusterExpansion::output(ostream &out) const{
    for(int i=0; i<getNumEffectiveClusters(); i++){
        SVectorXi cluster = effectiveClusters[i];
        double eci = effectiveClusterInteractions.coeff(i);
        out << cluster.nonZeros();
        for(SVectorXi::InnerIterator it(cluster); it; ++it){
            int supercellIndex = it.index();
            int unitcellIndex = supercell->getUnitCellIndex(supercellIndex);
            Vector3i cellPos = supercell->getCellPos(supercellIndex);
            out << " " << cellPos(0) << "," << cellPos(1) << "," << cellPos(2) << "-" << unitcellIndex;
        }
        out << " " << scientific << setprecision(15) << eci << " #ECI" << endl;
    }
}

VectorXi ClusterExpansion::getClusterCountVectorDifferential(const VectorXi &configuration, const VectorXi &oldConfiguration, const VectorXi &oldClusterCountVector, const vector<pair<int,int>> &changes) const{
    int numPositions = supercell->getNumPositions();
    VectorXi clusterCountVectorDiff = VectorXi::Zero(getNumEffectiveClusters());
    SVectorXi diffConf(numPositions);
    for(auto change: changes) diffConf.coeffRef(change.first)=change.second;

    for(auto change: changes){
        int posIndex = change.first;
        for( auto cluster: mappedClusters[posIndex] ){
            if( cluster.second.dot(diffConf)%2==0 ) continue;
            
            int firstIndex;
            SVectorXi::InnerIterator it1(cluster.second);
            while( diffConf.coeff(firstIndex=it1.index()) == 0 ) ++it1;
            if(firstIndex != posIndex) continue;

            int prod = 1;
            for(SVectorXi::InnerIterator it2(cluster.second); it2; ++it2) if(configuration[it2.index()] ==0 ) prod *= -1;
            if( prod == 1 ) clusterCountVectorDiff(cluster.first) += 2;
            else clusterCountVectorDiff(cluster.first) -= 2;
        }
    }
    return clusterCountVectorDiff + oldClusterCountVector;
}

ostream &operator<<(std::ostream &out, const ClusterExpansion &tgt){
    tgt.output(out);
    return out;
}
