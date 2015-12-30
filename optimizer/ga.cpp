#include <set>
#include <stdlib.h>
#include "ga.hpp"
#include <iostream>

using namespace std;
using namespace Eigen;

#define ALLOWABLE_ERROR pow(10.0, -7)

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
void GA::setEvalFunc(std::function<double(Eigen::SVectorXi)> evalFunc){ this->evalFunc = evalFunc; }

void GA::run(void){
    generation = 0;
    changes.clear();
    initializeParents();
    sortParents();
    changes.push_back(parents.front());
    cerr << generation << " "<< parents.front().second << endl;
    while(generation < generationLimit && parents.front().second < evalLimit ){
        generation++;
        selectElite();
        generateChildren();
        mutateChildren();
        evaluateChildren();
        selectParents();
        insertElite();
        sortParents();
        changes.push_back(parents.front());
        cerr << generation << " "<< parents.front().second << endl;
    }
}

vector<gaResult> GA::getResults(void) const{
    vector<gaResult> results;
    double eval = parents.front().second;
    for(auto parent:parents){
        if( eval-parent.second > ALLOWABLE_ERROR ) break;
        gaResult result;
        result.generation = generation;
        result.eval = parent.second;
        result.chromosome = parent.first;
        results.push_back(result);
    }
    auto less = [](const gaResult &obj1, const gaResult &obj2){ return obj1.chromosome < obj2.chromosome; };
    auto equal = [](const gaResult &obj1, const gaResult &obj2){ return obj1.chromosome == obj2.chromosome; };
    sort(results.begin(),results.end(),less);
    results.erase(unique(results.begin(), results.end(),equal), results.end());
    return results;
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
    auto greater = [](const pair<SVectorXi,double> &obj1, const pair<SVectorXi,double> &obj2){ return obj1.second > obj2.second; };
    sort(parents.begin(),parents.end(),greater);
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
    children.reserve(numChildren);
    children.resize(numChildren-numElites);
    uniform_int_distribution<int> parentRnd(0,numParents-1);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<(numChildren-numElites+1)/2; i++){
        SVectorXi parent1, parent2;
        parent1 = parents[parentRnd(*engine)].first;
        while( (parent2=parents[parentRnd(*engine)].first) == parent1 );

        SVectorXi child1, child2;
        SVectorXi sum, common, unique;
        sum = parent1 + parent2;
        common = (sum/2);
        unique = sum - common*2;
        common.prune(0);
        unique.prune(0);
        child1 = common;

        set<int> enableGeneSet;
        uniform_int_distribution<int> geneRnd(0,numTotalGene-1);
        while(enableGeneSet.size() < unique.nonZeros()/2){
            int enableGene = geneRnd(*engine);
            if( unique.coeff( enableGene ) == 1 ) enableGeneSet.insert( enableGene );
        }
        for(auto enableGene: enableGeneSet) child1.coeffRef(enableGene) = 1;
        child2 = sum - child1;
        child2.prune(0);
        
        children[i*2] = make_pair(child1, 0);
        if( i*2+1 < numChildren-numElites ) children[i*2+1] = make_pair(child2, 0);
    }
}

// (modified) swap
void GA::mutateChildren(void){
    uniform_real_distribution<double> judgeRnd(0,1);
    uniform_int_distribution<int> selectRnd(0,numTotalGene-1);

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<children.size();i++){
        auto &child = children[i];
        if(judgeRnd(*engine) < mutationP ) continue;
        int disableGene, enableGene;
        while( child.first.coeff( disableGene=selectRnd(*engine) ) != 1 );
        while( child.first.coeff( enableGene=selectRnd(*engine) ) != 0 );
        child.first.coeffRef(disableGene) = 0;
        child.first.coeffRef(enableGene) = 1;
        child.first.prune(0);
    }
}

void GA::selectParents(void){
    vector<double> ruletteTable;
    parents.clear();
    parents.reserve(numParents);
    parents.resize(numParents-numElites);

    double sum = 0;
    for(auto child: children){
        sum += child.second;
        ruletteTable.push_back(sum);
    }
    
    uniform_real_distribution<double> rnd(0,sum);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0; i<numParents-numElites; i++){
        double selector = rnd(*engine);
        vector<double>::iterator it;
        while( (  it = lower_bound(ruletteTable.begin(), ruletteTable.end(), selector ) ) == ruletteTable.end() );
        int index = it-ruletteTable.begin();
        parents[i] = children[index];
    }
}

void GA::evaluateChildren(void){
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0;i<children.size();i++){
        children[i].second = evalFunc(children[i].first);
    }
}

void GA::insertElite(void){
    for( auto elite: elites) parents.push_back(elite);
}
