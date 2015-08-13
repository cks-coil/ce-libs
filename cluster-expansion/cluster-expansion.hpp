#ifndef __INCLUDED_CLUSTER_EXPANSION_HPP__
#define __INCLUDED_CLUSTER_EXPANSION_HPP__

#include "../supercell/supercell.hpp"
#include "../supercell/clusters.hpp"

class ClusterExpansion{
public:
    ClusterExpansion(void);
    ~ClusterExpansion(void);
    void setSupercell(const Supercell *supercell);
    void setUniqueClusters(std::vector<Eigen::SVectorXi> uniqueClusters);
    void expandClusters(void);
private:
    std::vector<Eigen::SVectorXi> uniqueClusters;
    std::vector<std::vector<Eigen::SVectorXi>>  expandedClusters;
    const Supercell *supercell;
};

#endif
