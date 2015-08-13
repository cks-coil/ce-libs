#ifndef __INCLUDED_CLUSTER_EXPANSION_HPP__
#define __INCLUDED_CLUSTER_EXPANSION_HPP__

#include "../supercell/supercell.hpp"
#include "../supercell/clusters.hpp"

class ClusterExpansion{
public:
    ClusterExpansion(void);
    ~ClusterExpansion(void);
    void setSupercell(const Supercell *supercell);
    void setEffectiveClusters(std::vector<Eigen::SVectorXi> effectiveClusters);
    void expandClusters(void);
    int getNumEffectiveClusters(void) const;
    Eigen::VectorXi getClusterCountVector(Eigen::VectorXi configuration) const;
private:
    std::vector<Eigen::SVectorXi> effectiveClusters;
    std::vector<std::vector<Eigen::SVectorXi>>  expandedClusters;
    const Supercell *supercell;
};

#endif
