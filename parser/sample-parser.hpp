#ifndef __INCLUDED_SAMPLE_PARSER_HPP__
#define __INCLUDED_SAMPLE_PARSER_HPP__

#include <string>
#include <vector>
#include "supercell.hpp"
#include "eigen-extension.hpp"

class SampleParser{
public:
    SampleParser(void);
    void setFileName(std::string fileName);
    void setTargetSupercell(const Supercell *targetSupercell);
    void setParseSupercell(Supercell *parseSupercell);
    void parse(void);
    std::vector< std::pair<Eigen::VectorXi, double> > getSamples(void);
private:
    std::string fileName;
    const Supercell *targetSupercell;
    Supercell *parseSupercell;
    std::vector< std::pair<Eigen::VectorXi, double> > samples;
};

#endif
