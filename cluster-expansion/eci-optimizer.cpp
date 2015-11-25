#include "eci-optimizer.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/LU>
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
    MatrixXd allX(samples.size(), tgt->getNumEffectiveClusters());
    VectorXd allY(samples.size());
    for(int i=0; i<samples.size(); i++){
        allX.row(i) = samples[i].first.cast<double>();
        allY(i) = samples[i].second;
    }

    JacobiSVD<MatrixXd> svd(allX, ComputeThinU | ComputeThinV);
    int rankDiff = svd.rank() - tgt->getNumEffectiveClusters();
    if( rankDiff != 0 ) return 1.0/(double)rankDiff;

    MatrixXd hat = allX * (allX.transpose()*allX).inverse() * allX.transpose();
    VectorXd eci = svd.solve(allY);

    double score=0;

    for(int i=0; i<samples.size(); i++){
        VectorXd testX = samples[i].first.cast<double>();
        double testY = samples[i].second;
        double predictY = eci.dot(testX);
        score += pow( (testY - predictY) / (1-hat(i,i)), 2);
    }
    return score / (double)samples.size();
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
