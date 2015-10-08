#ifndef __INCLUDED_ECI_OPTIMIZER_HPP__
#define __INCLUDED_ECI_OPTIMIZER_HPP__

#include <vector>
#include <ostream>
#include "cluster-expansion.hpp"

class ECIOptimizer{
public:
    ECIOptimizer(void);
    ~ECIOptimizer(void);
    void setTarget(ClusterExpansion *target);
    void addSample(Eigen::VectorXi configuration, double energy);
    void optimizeECI(void);
    double getLOOCVScore(void);
    void output(std::ostream &out) const;
private:
    ClusterExpansion *tgt;
    std::vector<std::pair<Eigen::VectorXi, double>> samples;
};

std::ostream &operator<<(std::ostream &out, const ECIOptimizer &tgt);

#endif
