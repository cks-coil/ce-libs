#ifndef __INCLUDED_CLUSTER_EXPANSION_HPP__
#define __INCLUDED_CLUSTER_EXPANSION_HPP__

#include <ostream>
#include "../supercell/supercell.hpp"
#include "../supercell/clusters.hpp"

class ClusterExpansion{
public:
    ClusterExpansion(void);
    ~ClusterExpansion(void);
    void setSupercell(const Supercell *supercell);
    void setEffectiveClusters(std::vector<Eigen::SVectorXi> effectiveClusters);
    void setEffectiveClusterInteractions(Eigen::VectorXd effectiveClusterInteractions);
    void expandClusters(void);
    int getNumEffectiveClusters(void) const;
    double getEnergy(Eigen::VectorXi configuration) const;
    Eigen::VectorXi getClusterCountVector(Eigen::VectorXi configuration) const;
    const Supercell *getSupercell(void) const;
    void output(std::ostream &out) const;
private:
    const Supercell *supercell;
    std::vector<Eigen::SVectorXi> effectiveClusters;
    Eigen::VectorXd effectiveClusterInteractions;
    std::vector<std::vector<Eigen::SVectorXi>>  expandedClusters;
};

std::ostream &operator<<(std::ostream &out, const ClusterExpansion &tgt);

#endif
