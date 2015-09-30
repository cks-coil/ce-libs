#ifndef __INCLUDED_SAMPLE_PARSER_HPP__
#define __INCLUDED_SAMPLE_PARSER_HPP__

#include <string>
#include <vector>
#include "parser.hpp"
#include "supercell.hpp"
#include "eigen-extension.hpp"

class SampleParser : public Parser{
public:
    SampleParser(void);
    void setTargetSupercell(const Supercell *targetSupercell);
    void setParseSupercell(Supercell *parseSupercell);
    std::vector< std::pair<Eigen::VectorXi, double> > getSamples(void);
private:
    void parseLine(std::vector<std::string> strs);
    void resetResults(void);
    const Supercell *targetSupercell;
    Supercell *parseSupercell;
    std::vector< std::pair<Eigen::VectorXi, double> > samples;
};

#endif
