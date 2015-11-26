#ifndef __INCLUDED_ECI_OPTIMIZER_HPP__
#define __INCLUDED_ECI_OPTIMIZER_HPP__

#include <vector>
#include <ostream>
#include "cluster-expansion.hpp"
#include "eigen-extension.hpp"

class ECIOptimizer{
public:
    ECIOptimizer(void);
    ~ECIOptimizer(void);
    void setTarget(ClusterExpansion *target);
    void addSample(Eigen::VectorXi configuration, double energy);
    void optimizeECI(Eigen::SVectorXi flag);
    double getLOOCVScore(Eigen::SVectorXi flag) const;
    void output(std::ostream &out) const;
private:
    ClusterExpansion *tgt;
    std::vector< std::pair<Eigen::VectorXi, double> > samples;
    std::vector< std::pair<Eigen::VectorXi, double> > getCurrentSamples(const Eigen::SVectorXi &flag) const;
};

std::ostream &operator<<(std::ostream &out, const ECIOptimizer &tgt);

#endif
