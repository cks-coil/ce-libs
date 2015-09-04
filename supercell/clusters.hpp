#ifndef __INCLUDED_CLUSTERS_HPP__
#define __INCLUDED_CLUSTERS_HPP__

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "eigen-extension.hpp"
#include "supercell.hpp"

class Clusters{
public:
    Clusters(void);
    ~Clusters(void);
    void setSupercell(const Supercell *tgt);
    void setMaxDistance(double maxDistance);
    void setMaxNum(int maxNum);
    std::vector<Eigen::SVectorXi> getUniqueClusters(void) const;
    int getNumUniqueClusters(void) const;
    void findUniqueClusters(void);
private:
    bool isRange(Eigen::SVectorXi cluster);
    bool isUnique(Eigen::SVectorXi cluster, const std::vector<Eigen::SVectorXi> &refs);
    std::vector<Eigen::SVectorXi> getNextUniqueClusters(const std::vector<Eigen::SVectorXi> &seeds);
    std::vector<Eigen::SVectorXi> uniqueClusters;
    const Supercell *tgt;
    int maxNum;
    int maxDistance;
};

Eigen::SVectorXi getConvertedCluster(const Eigen::SVectorXi cluster, const Supercell &source, const Supercell &dest);

#endif
