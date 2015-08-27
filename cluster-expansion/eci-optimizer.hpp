#ifndef __INCLUDED_ECI_OPTIMIZER_HPP__
#define __INCLUDED_ECI_OPTIMIZER_HPP__

#include <vector>
#include "cluster-expansion.hpp"

class ECIOptimizer{
public:
    ECIOptimizer(void);
    ~ECIOptimizer(void);
    void setTarget(ClusterExpansion *target);
    void setSample(Eigen::SVectorXi configuration, double energy);
    void optimizeECI(void);
    double getLOOCVScore(void);
private:
    ClusterExpansion *tgt;
    std::vector<std::pair<Eigen::VectorXi, double>> samples;
};

#endif
