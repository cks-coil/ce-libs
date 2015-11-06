#ifndef __INCLUDED_CONF_ANALYZER_HPP__
#define __INCLUDED_CONF_ANALYZER_HPP__

#include <vector>
#include <ostream>
#include "eigen-extension.hpp"
#include "supercell.hpp"

class ConfAnalyzer{
public:
    ConfAnalyzer(void);
    ~ConfAnalyzer(void);
    void setUnitCell(const Supercell *unitCell);
    void setSupercell(const Supercell *supercell);
    void addUnitCellConf(Eigen::VectorXi unitCellConf);
    int getNumUnitCellConf(void) const;
    void expandUnitCellConfigurations(void);
    Eigen::VectorXi getConfCountVector(const Eigen::VectorXi &configuration) const;
    void getConfCountVectorDifferential(Eigen::VectorXi *confCountVector, const Eigen::VectorXi &configuration, const Eigen::VectorXi &oldConfiguration, const std::vector<std::pair<int,int>> &changes) const;
    void output(std::ostream &out) const;
private:
    int getUnitCellConfIndex(const Eigen::VectorXi &unitCellConf) const;
    Eigen::VectorXi getUnitCellConf(const Eigen::VectorXi &configuration, int pos) const;
    std::vector<std::vector<std::pair<int,std::vector<Eigen::VectorXi>>>> expandedUnitCellConfigurations;
    std::vector<Eigen::VectorXi> unitCellConfigurations;
    const Supercell *unitCell;
    const Supercell *supercell;
    int numUnitCellConf;
    int numUnitCellPos;
};

std::ostream &operator<<(std::ostream &out, const ConfAnalyzer &tgt);

#endif
