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

void ECIOptimizer::optimizeECI(void){
    MatrixXd X(samples.size(), tgt->getNumEffectiveClusters());
    VectorXd Y(samples.size());
    VectorXd ECI;
    for(int i=0; i<samples.size(); i++){
        X.row(i) = samples[i].first.cast<double>();
        Y(i) = samples[i].second;
    }
    ECI = X.jacobiSvd(ComputeThinU | ComputeThinV).solve(Y);
    tgt->setEffectiveClusterInteractions(ECI);
}

double ECIOptimizer::getLOOCVScore(void){
    double score=0;

    for(auto testSample: samples){
        VectorXd testX = testSample.first.cast<double>();
        double testY = testSample.second;
        MatrixXd trainingX(samples.size()-1, tgt->getNumEffectiveClusters());
        VectorXd trainingY(samples.size()-1);
        VectorXd ECI;
        int i=0;
        for(auto trainingSample: samples){
            if( testSample == trainingSample ) continue;
            trainingX.row(i) = trainingSample.first.cast<double>();
            trainingY(i) = trainingSample.second;
            i++;
        }
        ECI = trainingX.jacobiSvd(ComputeThinU | ComputeThinV).solve(trainingY);
        score += pow( (ECI.transpose() * testX - testY), 2);
        if( std::isnan(score) ){ return -1.0;}
    }
    return score / (double)samples.size();
}

void ECIOptimizer::output(ostream &out) const{
    for(auto sample:samples){
        double energy = tgt->getEffectiveClusterInteractions().transpose() * sample.first.cast<double>();
        double diff = energy - sample.second;
        out << sample.second << " " << energy << " " << diff << endl;
    }
}

ostream &operator<<(std::ostream &out, const ECIOptimizer &tgt){
    tgt.output(out);
    return out;
}
