#include "eci-optimizer.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <limits>

using namespace std;
using namespace Eigen;

ECIOptimizer::ECIOptimizer(void){
    tgt = nullptr;
    samples.clear();
}

ECIOptimizer::~ECIOptimizer(void){
    samples.clear();
}

void ECIOptimizer::setTarget(ClusterExpansion *target){
    tgt = target;
}

void ECIOptimizer::addSample(VectorXi configuration, double energy){
    auto sample = make_pair(tgt->getClusterCountVector(configuration), energy);
    samples.push_back(sample);
}

void ECIOptimizer::optimizeECI(SVectorXi flag){
    flag.prune(0);
    vector< pair<VectorXi, double> > currentSamples = getCurrentSamples(flag);

    MatrixXd X(currentSamples.size(), flag.nonZeros());
    VectorXd Y(currentSamples.size());
    for(int i=0; i<currentSamples.size(); i++){
        X.row(i) = currentSamples[i].first.cast<double>();
        Y(i) = currentSamples[i].second;
    }
    VectorXd eci = X.jacobiSvd(ComputeThinU | ComputeThinV).solve(Y);
    VectorXd allECI = VectorXd::Zero(tgt->getNumEffectiveClusters());
    int j=0;
    for(SVectorXi::InnerIterator it(flag); it; ++it){
        allECI(it.index()) = eci(j);
        j++;
    }
    tgt->setEffectiveClusterInteractions(allECI);
}

double ECIOptimizer::getLOOCVScore(SVectorXi flag) const{
    flag.prune(0);
    vector< pair<VectorXi, double> > currentSamples = getCurrentSamples(flag);
    int num = flag.nonZeros();

    MatrixXd allX(currentSamples.size(), num);
    VectorXd allY(currentSamples.size());
    for(int i=0; i<currentSamples.size(); i++){
        allX.row(i) = currentSamples[i].first.cast<double>();
        allY(i) = currentSamples[i].second;
    }

    JacobiSVD<MatrixXd> svd(allX, ComputeThinU | ComputeThinV);
    int rankDiff = svd.rank() - num;
    if( rankDiff != 0 ) return rankDiff;

    MatrixXd hat = allX * (allX.transpose()*allX).inverse() * allX.transpose();
    VectorXd eci = svd.solve(allY);

    double score=0;

    for(int i=0; i<currentSamples.size(); i++){
        VectorXd testX = currentSamples[i].first.cast<double>();
        double testY = currentSamples[i].second;
        double predictY = eci.dot(testX);
        score += pow( (testY - predictY) / (1-hat(i,i)), 2);
    }
    if( std::isnan(score) || std::isinf(score) ) return -1;
    return score / (double)currentSamples.size();
}

vector< pair<VectorXi, double> > ECIOptimizer::getCurrentSamples(const SVectorXi &flag) const{
    vector< pair<VectorXi, double> > currentSamples;
    for(auto sample: samples){
        VectorXi tmp(flag.nonZeros());
        int j=0;
        for(SVectorXi::InnerIterator it(flag); it; ++it){
            tmp(j)= sample.first(it.index());
            j++;
        }
        currentSamples.push_back( make_pair(tmp, sample.second));
    }
    return currentSamples;
}
    
void ECIOptimizer::output(ostream &out) const{
    for(auto sample:samples){
        double cellNum = tgt->getSupercell()->getCellSize().prod();
        double energy = tgt->getEffectiveClusterInteractions().transpose() * sample.first.cast<double>();
        double diff = energy - sample.second;
        out << sample.second/cellNum << " " << energy/cellNum << " " << diff/cellNum << endl;
    }
}

ostream &operator<<(std::ostream &out, const ECIOptimizer &tgt){
    tgt.output(out);
    return out;
}
