#ifndef __INCLUDED_PARSER_HPP__
#define __INCLUDED_PARSER_HPP__

#include <string>
#include <vector>
#include <Eigen/Core>


class Parser{
public:
    Parser(void);
    void setFileName(std::string fileName);
    void parse(void);
protected:
    virtual void resetResults(void);
    virtual void parseLine(std::vector<std::string> strs);
    Eigen::Vector3i parseCellPos(std::string str);
    Eigen::VectorXi parseUnitCellConf(std::string str, int numPos);
private:
    std::string fileName;
};

#endif
