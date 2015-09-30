#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include "cluster-parser.hpp"

using namespace std;
using namespace Eigen;

ClusterParser::ClusterParser(void){ supercell = nullptr; }

void ClusterParser::setSupercell(const Supercell *supercell){ this->supercell = supercell; }

vector<SVectorXi> ClusterParser::getClusters(void){ return clusters; }

VectorXd ClusterParser::getInteractions(void){
    VectorXd ECI(interactions.size());
    for(int i=0; i<interactions.size(); i++){
        ECI(i) = interactions[i];
    }
    return ECI;
}

void ClusterParser::parseLine(vector<string> strs){
    SVectorXi cluster(supercell->getNumPositions());
    int num = stoi(strs.front());
    if( num != (strs.size()-2) ) throw runtime_error("ClusterParser::parseLine");
    for(int i=0; i<num; i++){
        vector<string> v;
        Vector3i cellPos;
        int unitCellIndex;
        int supercellIndex;
        boost::algorithm::split(v, (strs[i+1]), boost::algorithm::is_any_of("-"), boost::algorithm::token_compress_on);
        if( v.size() != 2 ) throw runtime_error("ClusterParser::parseLine");
        cellPos = parseCellPos(v.front());
        unitCellIndex = stoi(v.back());
        supercellIndex = supercell->getSupercellIndex(unitCellIndex, cellPos);
        cluster.coeffRef(supercellIndex) = 1;
    }
    interactions.push_back( stod(strs.back()) );
    clusters.push_back(cluster);
}

void ClusterParser::resetResults(void){
    clusters.clear();
    interactions.clear();
}
