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
    Vector3i cellSize = parseCellPos(strs.front());
    if(strs.size() != (cellSize.prod()+2) ) throw runtime_error("SampleParser::parseLine");
    parseSupercell->setCellSize( cellSize );

    double energy = stod(strs.back());
    VectorXi configuration(parseSupercell->getNumPositions());

    for(int i=1; i<(cellSize.prod()+1); i++){
        VectorXi unitCellConf;
        Vector3i cellPos;
        vector<string> v;
        boost::algorithm::split(v, strs[i], boost::algorithm::is_any_of("-"), boost::algorithm::token_compress_on);
        if(v.size() != 2) throw runtime_error("SampleParser::pareLine");
        
        cellPos = parseCellPos(v[0]);
        unitCellConf = parseUnitCellConf(v[1], parseSupercell->getNumUnitCellPositions());
        for(int unitCellIndex=0; unitCellIndex<parseSupercell->getNumUnitCellPositions(); unitCellIndex++){
            int supercellIndex = parseSupercell->getSupercellIndex(unitCellIndex, cellPos);
            configuration(supercellIndex) = unitCellConf(unitCellIndex);
        }
    }

    VectorXi convertedConfiguration = getConvertedConfiguration(configuration, *parseSupercell, *targetSupercell);
    double convertedEnergy = energy * targetSupercell->getCellSize().prod() / parseSupercell->getCellSize().prod();
    samples.push_back( make_pair(convertedConfiguration, convertedEnergy) );
}



void SampleParser::resetResults(void){
    samples.clear();
}
