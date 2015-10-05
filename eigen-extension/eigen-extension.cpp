#include "eigen-extension.hpp"

using namespace Eigen;

bool operator == (const SMatrixXi &obj1, const SMatrixXi &obj2){
    Eigen::SMatrixXi tmpMatrix = obj1- obj2;
    return tmpMatrix.norm() == 0;
}

bool operator< (const SVectorXi &obj1, const SVectorXi &obj2){
    SVectorXi::InnerIterator it1(obj1), it2(obj2);
    while(it1 && it2){
        if(it1.index() != it2.index()) return (it1.index() > it2.index());
        if(it1.value() != it2.value()) return (it1.value() < it2.value());
        ++it1;
        ++it2;
    }
    return false;
}
