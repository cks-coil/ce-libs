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
    ~Clusters(void);
    void setSupercell(Supercell *tgt);
    void setMaxDistance(double maxDistance);
    void setMaxNum(int maxNum);
    std::vector<Eigen::SVectorXi> getUniqueClusters(void);
    int getNumUniqueClusters(void);
    void findUniqueClusters(void);
private:
    bool isRange(Eigen::SVectorXi cluster);
    bool isUnique(Eigen::SVectorXi cluster, const std::vector<Eigen::SVectorXi> &refs);
    std::vector<Eigen::SVectorXi> getNextUniqueClusters(const std::vector<Eigen::SVectorXi> &seeds);
    std::vector<Eigen::SVectorXi> uniqueClusters;
    Supercell *tgt;
    int maxNum;
    int maxDistance;
};

inline bool operator < (const Eigen::SparseVector<int> &obj1, const  Eigen::SparseVector<int> &obj2);

Eigen::SVectorXi getConvertedCluster(const Eigen::SVectorXi cluster, const Supercell &source, const Supercell &dest);

#endif
