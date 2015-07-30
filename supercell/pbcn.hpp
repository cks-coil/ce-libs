#ifndef __INCLUDED_PBCN_HPP__
#define __INCLUDED_PBCN_HPP__

#include "supercell.hpp"

class PbcnSupercell : public Supercell{
public:
    PbcnSupercell(void);
protected:
    void calcUnitCellFractionalPositions(void);
    void calcSpaceGroupSymmetryMatrices(void);
};

#endif
