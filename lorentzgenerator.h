#ifndef LORENTZGENERATOR_H
#define LORENTZGENERATOR_H
#include <armadillo>
#include <vector>
using namespace arma;

class LorentzGenerator
{
private:
    std::vector<std::vector<cx_mat44> > m_generators;
    static cx_mat44& getGenerator(int mu, int nu);
    mat44 m_metricTensor;

    LorentzGenerator();
    LorentzGenerator(LorentzGenerator const&) = delete;
    void operator=(LorentzGenerator const&)  = delete;

public:
    static LorentzGenerator& getInstance() {
        static LorentzGenerator instance;
        return instance;
    }

    static mat44 generateElement(const mat44 &omega);
};

#endif // LORENTZGENERATOR_H
