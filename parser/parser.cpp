#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "parser.hpp"
#include <exception>

using namespace std;
using namespace Eigen;

Parser::Parser(void){
    fileName = "";
}

void Parser::setFileName(string fileName){ this->fileName = fileName; }

void Parser::parse(void){
    string buf;
    vector<string> v;
    ifstream ifs(fileName.c_str());
    if( !ifs ) cerr << "ERROR: Could NOT Open File: " << fileName << endl;

    resetResults();

    while( getline(ifs,buf) ){
        boost::algorithm::split(v, buf, boost::algorithm::is_space(), boost::algorithm::token_compress_on);

        if( v.size() == 0 || v.front().front() == '#') continue;
        for( auto it = v.begin(); it != v.end(); ++it){
            if((*it).front() == '#'){
                v.erase(it, v.end());
                break;
            }
        }

        try{  parseLine(v); } 
        catch(exception& e){
            cerr << "ERROR: Could NOT Parse Line ( " << e.what() << " ) : " << buf << endl;
            exit(1);
        }
        catch(...){
            cerr << "ERROR: Could NOT Parse Line : " << buf << endl;
            exit(1);
        }
    }
}

void Parser::resetResults(void){}
void Parser::parseLine(vector<string> strs){}

Vector3i Parser::parseCellPos(string str){
    Vector3i cellPos;
    vector<string> v;
    boost::algorithm::split(v, str, boost::algorithm::is_any_of(","), boost::algorithm::token_compress_on);
    if(v.size() != 3) throw runtime_error("Parser::parseCellPos");
    for(int i=0; i<3; i++) cellPos(i) = stoi(v[i]);
    return cellPos;
}

VectorXi Parser::parseUnitCellConf(string str, int numPos){
    VectorXi conf(numPos);
    if(str.size() != numPos) throw runtime_error("Parser::parseUnitCellConf");
    for(int i=0; i<numPos; i++){
        switch(str[i]){
        case '0':
            conf(i) = 0;
            break;
        case '1':
            conf(i) = 1;
            break;
        default:
            throw runtime_error("Parser::parseUnitCellConf");
        }
    }
    return conf;
}

VectorXi Parser::parseConfiguration(vector<string> strs, Supercell *parseSupercell){
    Vector3i cellSize = parseCellPos(strs.front());
    if(strs.size() != (cellSize.prod()+1) ) throw runtime_error("Parser::parseUnitCellConf");
    parseSupercell->setCellSize( cellSize );

    VectorXi configuration(parseSupercell->getNumPositions());

    for(int i=1; i<(cellSize.prod()+1); i++){
        VectorXi unitCellConf;
        Vector3i cellPos;
        vector<string> v;
        boost::algorithm::split(v, strs[i], boost::algorithm::is_any_of("-"), boost::algorithm::token_compress_on);
        if(v.size() != 2) throw runtime_error("Parser::parseUnitCellConf");

        cellPos = parseCellPos(v[0]);
        unitCellConf = parseUnitCellConf(v[1], parseSupercell->getNumUnitCellPositions());
        for(int unitCellIndex=0; unitCellIndex<parseSupercell->getNumUnitCellPositions(); unitCellIndex++){
            int supercellIndex = parseSupercell->getSupercellIndex(unitCellIndex, cellPos);
            configuration(supercellIndex) = unitCellConf(unitCellIndex);
        }
    }

    return configuration;
}
