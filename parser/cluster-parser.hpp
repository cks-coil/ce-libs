#ifndef __INCLUDED_CLUSTER_PARSER_HPP__
#define __INCLUDED_CLUSTER_PARSER_HPP__

#include <string>
#include <vector>
#include "parser.hpp"
#include "supercell.hpp"
#include "eigen-extension.hpp"

class ClusterParser: public Parser{
public:
    ClusterParser(void);
    void setSupercell(const Supercell *supercell);
    std::vector<Eigen::SVectorXi> getClusters(void);
    Eigen::VectorXd getInteractions(void);
private:
    void parseLine(std::vector<std::string> strs);
    void resetResults(void);
    const Supercell *supercell;
    std::vector<Eigen::SVectorXi> clusters;
    std::vector<double> interactions;
};

#endif

