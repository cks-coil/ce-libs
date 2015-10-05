#ifndef __INCLUDED_EIGEN_EXTENSION_HPP__
#define __INCLUDED_EIGEN_EXTENSION_HPP__

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace Eigen{
    typedef SparseMatrix<int> SMatrixXi;
    typedef SparseVector<int> SVectorXi;
}

bool operator< (const Eigen::VectorXi &obj1, const  Eigen::VectorXi &obj2);

bool operator== (const Eigen::SMatrixXi &obj1, const Eigen::SMatrixXi &obj2);
bool operator< (const Eigen::SVectorXi &obj1, const  Eigen::SVectorXi &obj2);


#endif
