#include <fstream>
#include <iostream>
#include <exception>
#include <boost/algorithm/string.hpp>
#include "sample-parser.hpp"

using namespace std;
using namespace Eigen;

SampleParser::SampleParser(void){
    targetSupercell = nullptr;
    parseSupercell = nullptr;
}

void SampleParser::setTargetSupercell(const Supercell *targetSupercell){ this->targetSupercell = targetSupercell; }
void SampleParser::setParseSupercell(Supercell *parseSupercell){ this->parseSupercell = parseSupercell; }
vector< pair<VectorXi, double> > SampleParser::getSamples(void){ return samples; }


void SampleParser::parseLine(vector<string> strs){
    if(strs.size() < 2) throw runtime_error("SampleParser::parseLine");

    vector<string> confStrs;
    copy(strs.begin(), strs.end()-1, back_inserter(confStrs) );
    VectorXi configuration = parseConfiguration(confStrs, parseSupercell);
    double energy = stod(strs.back());
    VectorXi convertedConfiguration = getConvertedConfiguration(configuration, *parseSupercell, *targetSupercell);
    double convertedEnergy = energy * targetSupercell->getCellSize().prod() / parseSupercell->getCellSize().prod();
    samples.push_back( make_pair(convertedConfiguration, convertedEnergy) );
}



void SampleParser::resetResults(void){
    samples.clear();
}
