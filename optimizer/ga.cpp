#include <stdlib.h>
#include "ga.hpp"

using namespace std;
using namespace Eigen;

GA::GA(void){
    std::random_device raddev;
    engine = new std::mt19937(randdev());
    generation = 0;
    numTotalGene = 0;
    numEffectiveGene = 0;
    numParents = 0;
    numChildren = 0;
    numElites = 0;
    generationLimit = 0;
    evalLimit = 0.0;
    mutationP = 0.0;
    evalFunc = nullptr;
    parents.clear();
    children.clear();
}

GA::~GA(void){
    delete engine;
    parents.clear();
    children.clear();
}

void GA::setNumTotalGene(int numTotalGene){ this->numTotalGene = numTotalGene; }
void GA::setNumEffectiveGene(int numEffectiveGene){ this->numEffectiveGene = numEffectiveGene; }
void GA::setNumParents(int numParents){ this->numParents = numParents; }
void GA::setNumChildren(int numChildren){ this->numChildren = numChildren; }
void GA::setNumElites(int numElites){ this->numElites = numElites; }
void GA::setGenerationLimit(int generationLimit){ this->generationLimit = generationLimit; }
void GA::setEvalLimit(double evalLimit){ this->evalLimit = evalLimit; }
void GA::setMutationP(double mutationP){ this->mutationP = mutationP; }
void GA::setEvalFunc(double (*evalFunc)(Eigen::SVectorXi chromosome)){ this->evalFunc = evalFunc; }

void GA::run(void){}

void GA::getResults(int *generation, double *eval, SVectorXi *chromosome) const{
    *generation = this->generation;
    *eval = parents[0].second;
    *gene = parents[0].first;
}


void GA::initializeParents(void){
}
