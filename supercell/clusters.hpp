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
    bool isRange(Eigen::SVectorXi cluster);
    bool isUnique(Eigen::SVectorXi cluster, const std::vector<Eigen::SVectorXi> &refs);
    std::vector<Eigen::SVectorXi> getNextClusters(const std::vector<Eigen::SVectorXi> &seeds);
    std::vector<Eigen::SVectorXi> uniqueClusters;
    std::vector<Eigen::SVectorXi> allClusters;
    Supercell *tgt;
    int maxNum;
    int maxDistance;
};

inline bool operator < (const Eigen::SparseVector<int> &obj1, const  Eigen::SparseVector<int> &obj2);

#endif
