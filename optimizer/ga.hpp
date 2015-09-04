#ifndef __INCLUDED_GA_HPP__
#define __INCLUDED_GA_HPP__

#include <vector>
#include <random>
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
    void setEvalFunc(double (*evalFunc)(Eigen::SVectorXi chromosome));
    void setFixedGene(Eigen::SVectorXi fixedGene);
    void run(void);
    void getResults(int *generation, double *eval, Eigen::SVectorXi *chromosome) const;
private:
    std::mt19937 *engine;
    void initializeParents(void);
    void evaluate(void);
    void selecte(void);
    void mutate(void);
    Eigen::SVectorXi getDecompressedChromosome(Eigen:: SVectorXi chromosome) const;
    int generation;
    int numTotalGene;;
    int numEffectiveGene;
    int numParents;
    int numChildren;
    int numElites;
    int generationLimit;
    double evalLimit;
    double mutationP;
    double (*evalFunc)(Eigen::SVectorXi chromosome);
    Eigen::SVectorXi fixedGene;
    std::vector< std::pair<Eigen::SVectorXi, double> > parents;
    std::vector< std::pair<Eigen::SVectorXi, double> > children;
};

#endif
