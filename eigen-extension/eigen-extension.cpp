#include "eigen-extension.hpp"

using namespace Eigen;

bool operator == (const SMatrixXi &obj1, const SMatrixXi &obj2){
    Eigen::SMatrixXi tmpMatrix = obj1- obj2;
    return tmpMatrix.norm() == 0;
}

bool operator< (const SVectorXi &obj1, const SVectorXi &obj2){
    SVectorXi diff = obj1- obj2;
    diff.prune(0);
    if ( diff.nonZeros()==0 ) return false;
    SVectorXi::InnerIterator it(diff);
    return (it.value() < 0);
}
