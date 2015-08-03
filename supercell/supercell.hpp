#ifndef __INCLUDED_SUPERCELL_HPP__
#define __INCLUDED_SUPERCELL_HPP__

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Eigen{
    typedef SparseMatrix<int> SMatrixXi;
}

class Supercell{
public:
    Supercell(void);
    ~Supercell(void);
    void setCellSize(Eigen::Vector3i cellSize);
    void setAtomicPos(Eigen::Vector3d atomicPos);
    void setCrystalAxis(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c);
    int getNumPositions(void);
    int getNumUnitCellPositions(void);
    Eigen::Vector3i getCellSize(void);
    Eigen::Matrix3d getCrystalAxisMatrix(void);
    Eigen::Vector3d getFractionalPos(int index);
    Eigen::Vector3d getOrthogonalPos(int index);
    std::vector<Eigen::Vector3d> getFractionalPositions(void);
    std::vector<Eigen::Vector3d> getOrthogonalPositions(void);
    std::vector<Eigen::SMatrixXi> getSymmetryMatrices(void);
protected:
    void updateVariables(void);
    void calcFractionalPositions(void);
    void calcSymmetryMatrices(void);
    virtual void calcUnitCellFractionalPositions(void);
    virtual void calcSpaceGroupSymmetryMatrices(void);
    void checkSymmetryMatrices(void);
    void calcOrthogonalPositions(void);
    void periodicBoundaryCondition(std::vector<Eigen::Vector3d> &positions, Eigen::Vector3i cellSize);
    std::vector<Eigen::Vector3d> glideReflection(std::vector<Eigen::Vector3d> positions, Eigen::Vector3d transVector, Eigen::Vector3d reflectionPos);
    Eigen::SMatrixXi getSymmetryMatrix(std::vector<Eigen::Vector3d> arr1, std::vector<Eigen::Vector3d> arr2);
    std::vector<Eigen::SMatrixXi> getPoweredMatrices(Eigen::SMatrixXi matrix, int maxN);
    std::vector<Eigen::SMatrixXi> getSlideMatrices(void);
    Eigen::Vector3i cellSize;
    Eigen::Vector3d atomicPos;
    Eigen::Matrix3d crystalAxisMatrix;
    std::vector<Eigen::Vector3d> fractionalPositions;
    std::vector<Eigen::Vector3d> orthogonalPositions;
    std::vector<Eigen::SMatrixXi> symmetryMatrices;
    std::vector<Eigen::Vector3d> unitCellFractionalPositions;
    std::vector<Eigen::SMatrixXi> spaceGroupSymmetryMatrices;
};


bool operator == (const Eigen::SMatrixXi &obj1, const Eigen::SMatrixXi &obj2);

#endif
