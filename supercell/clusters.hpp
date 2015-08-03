#ifndef __INCLUDED_CLUSTERS_HPP__
#define __INCLUDED_CLUSTERS_HPP__

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "supercell.hpp"

namespace Eigen{
    typedef SparseVector<int> SVectorXi;
}

class Clusters{
public:
    Clusters(void);
    void setSupercell(Supercell *tgt);
    void setMaxDistance(double maxDistance);
    void setMaxNum(int maxNum);
    void findClusters(void);
private:
    std::vector<Eigen::SVectorXi> uniqueClusters;
    std::vector<Eigen::SVectorXi> allClusters;
    Supercell *tgt;
    int maxNum;
    int maxDistance;
};

#endif
