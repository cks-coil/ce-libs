#ifndef __INCLUDED_SUPERCELL_HPP__
#define __INCLUDED_SUPERCELL_HPP__

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "eigen-extension.hpp"

class Supercell{
public:
    Supercell(void);
    ~Supercell(void);
    void setCellSize(Eigen::Vector3i cellSize);
    void setAtomicPos(Eigen::Vector3d atomicPos);
    void setCrystalAxis(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c);
    int getNumPositions(void) const;
    int getNumUnitCellPositions(void) const;
    Eigen::Vector3i getCellSize(void) const;
    Eigen::Matrix3d getCrystalAxisMatrix(void) const;
    Eigen::Vector3d getFractionalPos(int supercellIndex) const;
    Eigen::Vector3d getFractionalPos(int unitCellIndex, Eigen::Vector3i cellPos) const;
    Eigen::Vector3d getOrthogonalPos(int supercellIndex) const;
    Eigen::Vector3d getOrthogonalPos(int unitCellIndex, Eigen::Vector3i cellPos) const;
    const std::vector<Eigen::Vector3d>& getFractionalPositions(void) const;
    const std::vector<Eigen::Vector3d>& getOrthogonalPositions(void) const;
    const std::vector<Eigen::SMatrixXi>& getSymOpMatrices(void) const;
    const std::vector<Eigen::SMatrixXi>& getSpaceGroupSymOpMatrices(void) const;
    int getSupercellIndex(int unitCellIndex, Eigen::Vector3i cellPos) const;
    int getUnitCellIndex(int supercellIndex) const;
    Eigen::Vector3i getCellPos(int supercellIndex) const;
    double getVolume(void) const;
    void calcPositions(void);
    void calcSymOpMatrices(void);
    void checkSymOpMatrices(void);
protected:
    int numPositions;
    int numUnitCellPositions;
    void calcFractionalPositions(void);
    void calcOrthogonalPositions(void);
    virtual void calcUnitCellFractionalPositions(void);
    virtual void calcSpaceGroupSymOpMatrices(void);
    void periodicBoundaryCondition(std::vector<Eigen::Vector3d> &positions, Eigen::Vector3i cellSize);
    std::vector<Eigen::Vector3d> glideReflection(std::vector<Eigen::Vector3d> positions, Eigen::Vector3d transVector, Eigen::Vector3d reflectionPos);
    Eigen::SMatrixXi getSymOpMatrix(std::vector<Eigen::Vector3d> arr1, std::vector<Eigen::Vector3d> arr2);
    std::vector<Eigen::SMatrixXi> getPoweredMatrices(Eigen::SMatrixXi matrix, int maxN);
    std::vector<Eigen::SMatrixXi> getSlideMatrices(void);
    Eigen::Vector3i cellSize;
    Eigen::Vector3d atomicPos;
    Eigen::Matrix3d crystalAxisMatrix;
    std::vector<Eigen::Vector3d> fractionalPositions;
    std::vector<Eigen::Vector3d> orthogonalPositions;
    std::vector<Eigen::SMatrixXi> symOpMatrices;
    std::vector<Eigen::Vector3d> unitCellFractionalPositions;
    std::vector<Eigen::SMatrixXi> spaceGroupSymOpMatrices;
};

Eigen::VectorXi getConvertedConfiguration(Eigen::VectorXi configuration, const Supercell &source, const Supercell &dest);
Eigen::VectorXi getUnitCellConfiguration(Eigen::VectorXi configuration, const Supercell &supercell, Eigen::Vector3i cellPos);


std::string getConfCombinedStr(Eigen::VectorXi configuration);
std::string getConfSplitedStr(Eigen::VectorXi configuration, const Supercell &supercell);


#endif
