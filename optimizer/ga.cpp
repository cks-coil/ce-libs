#include <set>
#include <stdlib.h>
#include "ga.hpp"

using namespace std;
using namespace Eigen;

GA::GA(void){
    std::random_device randdev;
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
    elites.clear();
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

void GA::run(void){
    generation = 0;
    changes.clear();
    initializeParents();
    sortParents();
    changes.push_back(parents[0]);
    while(generation < generationLimit && parents[0].second < evalLimit ){
        selectElite();
        generateChildren();
        mutateChildren();
        evaluateChildren();
        selectParents();
        insertElite();
        sortParents();
        changes.push_back(parents[0]);
        generation++;
    }
}

void GA::getResults(int *generation, double *eval, SVectorXi *chromosome) const{
    *generation = this->generation;
    *eval = parents[0].second;
    *chromosome = parents[0].first;
}


void GA::initializeParents(void){
    parents.clear();
    uniform_int_distribution<int> rnd(0,numTotalGene-1);

    for(int i=0; i<numParents; i++){
        SVectorXi chromosome(numTotalGene);
        for(int j=0; j<numEffectiveGene; j++){
            int gene;
            while( chromosome.coeff(gene=rnd(*engine)) != 0 );
            chromosome.coeffRef(gene) = 1;
        }
        parents.push_back( make_pair(chromosome, evalFunc(chromosome)) );
    }
}

void GA::sortParents(void){
    auto less = [](const pair<SVectorXi,double> &obj1, const pair<SVectorXi,double> &obj2){ return obj1.second < obj2.second; };
    sort(parents.begin(),parents.end(),less);
}

void GA::selectElite(void){
    elites.clear();
    for(int i=0; i<numElites; i++){
        elites.push_back(parents[i]);
    }
}

// (modified) uniform crossover
void GA::generateChildren(void){
    children.clear();
    uniform_int_distribution<int> rnd(0,numChildren - numElites-1);

    while(children.size() < numChildren - numElites){
        SVectorXi parent1, parent2;
        parent1 = parents[rnd(*engine)].first;
        while( (parent2=parents[rnd(*engine)].first) == parent1 );

        SVectorXi child1, child2;
        SVectorXi sum, common, unique;
        sum = parent1 + parent2;
        common = (sum/2);
        unique = sum - common*2;
        common.prune(0);
        unique.prune(0);
        child1 = common;

        set<int> enableIndexSet;
        uniform_int_distribution<int> rnd(0,unique.nonZeros()/2);
        while(enableIndexSet.size() < unique.nonZeros()/2){
            int enableIndex = rnd(*engine);
            if( unique.coeff( enableIndex ) == 1 ) enableIndexSet.insert( enableIndex );
        }
        for(auto enableIndex: enableIndexSet) child1.coeffRef(enableIndex) = 1;
        child2 = sum - child1;
        child2.prune(0);
        
        children.push_back( make_pair(child1, 0) );
        children.push_back( make_pair(child2, 0) );
    }
}

// (modified) swap
void GA::mutateChildren(void){
    uniform_real_distribution<double> judgeRnd(0,1);
    uniform_int_distribution<int> selectRnd(0,numChildren - numElites-1);

    for(auto &child: children){
        if(judgeRnd(*engine) > mutationP ) continue;
        int disableIndex, enableIndex;
        while( child.first.coeff( disableIndex=selectRnd(*engine) ) != 1 );
        while( child.first.coeff( enableIndex=selectRnd(*engine) ) != 0 );
        child.first.coeffRef(disableIndex) = 0;
        child.first.coeffRef(enableIndex) = 1;
    }
}

void GA::selectParents(void){
    vector<double> ruletteTable;
    double sum = 0;
    parents.clear();

    for(auto child: children){
        sum += child.second;
        ruletteTable.push_back(sum);
    }
    
    uniform_real_distribution<double> rnd(0,sum);
    while(parents.size() < numParents-numElites ){
        double selector = rnd(*engine);
        for(int i=0; i<ruletteTable.size(); i++){
            if( ruletteTable[i] >= selector ) parents.push_back(children[i]);
        }
    }
}

void GA::evaluateChildren(void){
    for(auto &child: children){
        child.second = evalFunc(child.first);
    }
}

void GA::insertElite(void){
    for( auto elite: elites) parents.push_back(elite);
}
