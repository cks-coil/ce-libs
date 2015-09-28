#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "sample-parser.hpp"

using namespace std;
using namespace Eigen;

SampleParser::SampleParser(void){
    fileName = "";
    targetSupercell = nullptr;
    parseSupercell = nullptr;
}

void SampleParser::setFileName(string fileName){ this->fileName = fileName; }
void SampleParser::setTargetSupercell(const Supercell *targetSupercell){ this->targetSupercell = targetSupercell; }
void SampleParser::setParseSupercell(Supercell *parseSupercell){ this->parseSupercell = parseSupercell; }

void SampleParser::parse(void){
    ifstream ifs(fileName.c_str());
    string buf;
    vector <string> v;

    while( getline(ifs,buf) ){
        boost::algorithm::split(v, buf, boost::algorithm::is_space(), boost::algorithm::token_compress_on);
        parseSupercell->setCellSize( Vector3i(1,1,1) );
        VectorXi configuration(parseSupercell->getNumUnitCellPositions());

        for(int i=0; i<(int)v.front().size(); i++){
            switch(v.front()[i]){
            case '0':
                configuration(i) = 0;
                break;
            case '1':
                configuration(i) = 1;
                break;
            default:
                cerr << "ERROR: Can NOT Parse Sample: " <<  v.front() << endl;
                exit(1);
            }
        }
        double energy;
        char *endptr;
        energy = strtod(v.back().c_str(), &endptr);
        if( (errno !=0) || (*endptr != '\0') ){ // 変換失敗
            cerr << "ERROR: Can NOT Parse Sample: " <<  v.back() << endl;
            exit(1);
        }
        
        VectorXi convertedConfiguration = getConvertedConfiguration(configuration, *parseSupercell, *targetSupercell);
        double convertedEnergy = energy * targetSupercell->getCellSize().prod() / parseSupercell->getCellSize().prod();
        samples.push_back( make_pair(convertedConfiguration, convertedEnergy) );
    }
}

vector< pair<VectorXi, double> > SampleParser::getSamples(void){
    return samples;
}
