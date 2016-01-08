#include "conf-analyzer.hpp"

using namespace std;
using namespace Eigen;

ConfAnalyzer::ConfAnalyzer(void){}
ConfAnalyzer::~ConfAnalyzer(void){}

void ConfAnalyzer::setUnitCell(const Supercell *unitCell){this->unitCell = unitCell;}
void ConfAnalyzer::setSupercell(const Supercell *supercell){this->supercell = supercell;}
void ConfAnalyzer::addUnitCellConf(VectorXi unitCellConf){ unitCellConfigurations.push_back(unitCellConf);}

int ConfAnalyzer::getNumUnitCellConf(void) const{ return numUnitCellConf; }

void ConfAnalyzer::expandUnitCellConfigurations(void){
    expandedUnitCellConfigurations.clear();
    expandedUnitCellConfigurations.resize(unitCell->getNumPositions()+1);
    auto less = [](const VectorXi &obj1, const VectorXi &obj2){return obj1<obj2;};
    auto equal = [](const VectorXi &obj1, const VectorXi &obj2){return obj1==obj2;};
    int index=0;
    for(auto unitCellConf: unitCellConfigurations){
        vector<VectorXi> expanded;
        int sum = unitCellConf.sum();
        for(auto m: unitCell->getSymOpMatrices()){
            expanded.push_back( m * unitCellConf);
        }
        sort(expanded.begin(),expanded.end(),less);
        expanded.erase(unique(expanded.begin(),expanded.end(), equal),expanded.end());
        expandedUnitCellConfigurations[sum].push_back(make_pair(index, expanded));
        index++;
    }
    numUnitCellConf = unitCellConfigurations.size();
    numUnitCellPos = unitCell->getNumUnitCellPositions();
}

VectorXi ConfAnalyzer::getConfCountVector(const VectorXi &configuration) const{
    VectorXi confCountVector = VectorXi::Zero(numUnitCellConf);
    for(int i=0; i< supercell->getCellSize().prod();i++){
        int index = getUnitCellConfIndex(getUnitCellConf(configuration, i*numUnitCellPos));
        if(index != -1) confCountVector[index]++;
    }
    return confCountVector;
}

void ConfAnalyzer::getConfCountVectorDifferential(VectorXi *confCountVector, const VectorXi &configuration, const VectorXi &oldConfiguration, const vector<pair<int, int> > &changes) const{
    for(auto const &change: changes){
        int index1 = getUnitCellConfIndex(getUnitCellConf(configuration, change.first));
        int index2 = getUnitCellConfIndex(getUnitCellConf(oldConfiguration, change.first));
        if(index1 != -1) (*confCountVector)[index1]++;
        if(index2 != -1) (*confCountVector)[index2]--;
    }
}

VectorXi ConfAnalyzer::getMaxParcentageConfCountVector(const VectorXi &configuration, const vector<VectorXi> &unitCellConfs) const{
    vector<int> unitCellConfIndexes;
    for(auto conf:unitCellConfs){ unitCellConfIndexes.push_back(getUnitCellConfIndex(conf)); }

    int bestSum = -1;
    VectorXi maxParcentageConfCountVector;
    for(auto const &m: supercell->getSpaceGroupSymOpMatrices()){
        VectorXi converted = m*configuration;
        VectorXi confCountVector = getConfCountVector(converted);
        int sum=0;
        for(auto index: unitCellConfIndexes) sum += confCountVector[index];
        if(sum > bestSum){
            bestSum = sum;
            maxParcentageConfCountVector = confCountVector;
        }
    }
    return maxParcentageConfCountVector;
}

void ConfAnalyzer::output(ostream &out) const{
    string str;
    int index=0;
    for(auto conf:unitCellConfigurations){
        str +=  getConfCombinedStr(conf) + "(" + to_string(index) + ") ";
        index ++;
    }
    str.erase(--str.end());
    out << str;
}

int ConfAnalyzer::getUnitCellConfIndex(const VectorXi &unitCellConf) const{
    int sum = unitCellConf.sum();
    for(auto const &expandedConf: expandedUnitCellConfigurations[sum]){
        for(auto const &conf:expandedConf.second){
            if( conf == unitCellConf ) return expandedConf.first;
        }
    }
    return -1;
}

VectorXi ConfAnalyzer::getUnitCellConf(const Eigen::VectorXi &configuration, int pos) const{
    VectorXi unitCellConf(numUnitCellPos);
    int startPos = pos - pos%numUnitCellPos;
    for(int i=0; i<numUnitCellPos; i++){
        unitCellConf[i] = configuration[i+startPos];
    }
    return unitCellConf;
}

ostream &operator<<(ostream &out, const ConfAnalyzer &tgt){
    tgt.output(out);
    return out;
}
