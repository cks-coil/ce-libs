#ifndef __INCLUDED_GA_HPP__
#define __INCLUDED_GA_HPP__

#include <vector>
#include <random>
#include <functional>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "eigen-extension.hpp"

class GA{
public:
    GA(void);
    ~GA(void);
    void setNumTotalGene(int numTotalGene);
    void setNumEffectiveGene(int numEffectiveGene);
    void setNumParents(int numParents);
    void setNumChildren(int numChildren);
    void setNumElites(int numElites);
    void setGenerationLimit(int generationLimit);
    void setEvalLimit(double evalLimit);
    void setMutationP(double mutationP);
    void setEvalFunc(std::function<double(Eigen::SVectorXi)> evalFunc);
    void run(void);
    void getResults(int *generation, double *eval, Eigen::SVectorXi *chromosome) const;
private:
    std::mt19937 *engine;
    void initializeParents(void);
    void sortParents(void);
    void selectElite(void);
    void generateChildren(void);
    void mutateChildren(void);
    void evaluateChildren(void);
    void selectParents(void);
    void insertElite(void);
    int generation;
    int numTotalGene;
    int numEffectiveGene;
    int numParents;
    int numChildren;
    int numElites;
    int generationLimit;
    double evalLimit;
    double mutationP;
    std::function<double(Eigen::SVectorXi)> evalFunc;
    std::vector< std::pair<Eigen::SVectorXi, double> > parents;
    std::vector< std::pair<Eigen::SVectorXi, double> > children;
    std::vector< std::pair<Eigen::SVectorXi, double> > elites;
    std::vector< std::pair<Eigen::SVectorXi, double> > changes;
};

#endif
