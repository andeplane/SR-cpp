#include "lorentzgenerator.h"
#include <iostream>
using namespace std;
cx_mat44 &LorentzGenerator::getGenerator(int mu, int nu)
{
    return LorentzGenerator::getInstance().m_generators[mu][nu];
}

LorentzGenerator::LorentzGenerator()
{
    m_metricTensor.zeros();
    m_metricTensor(0,0) = 1;
    m_metricTensor(1,1) = -1;
    m_metricTensor(2,2) = -1;
    m_metricTensor(3,3) = -1;

    m_generators.resize(4);
    for(auto &element : m_generators) {
        element.resize(4);
    }

    for(int mu=0; mu<4; mu++) {
        for(int nu=0; nu<4; nu++) {
            cx_mat44 &matrix = m_generators[mu][nu];
            matrix.zeros();
            bool isBoostGenerator = mu==0 || nu == 0;

            for(int a=0; a<4; a++) {
                for(int b=0; b<4; b++) {
                    cx_double value = cx_double(0,1);
                    if(a==b) continue; // These matrices have diagonal zero

                    if(isBoostGenerator) {
                        // Boosts have symmetric generators
                        value *= ((a==mu&&b==nu) + (a==nu&&b==mu));
                    } else {
                        // Rotations have antisymmetric generators
                        value *= ((a==mu&&b==nu) - (a==nu&&b==mu));
                    }

                    matrix(a,b) = value;
                }
            }
        }
    }
}

mat44 LorentzGenerator::generateElement(const mat44 &omega)
{
    cx_mat44 exponentialArgument;
    exponentialArgument.zeros();

    // We will compute element = exp(i/2*omega*J). First we compute omega*J
    for(int mu=0; mu<4; mu++) {
        for(int nu=0; nu<4; nu++) {
            cx_mat44 &generator = getGenerator(mu,nu);
            exponentialArgument += generator*omega(mu,nu);
        }
    }

    // Now multiply by i/2
    exponentialArgument *= 0.5*cx_double(0,1);

    // Then take the exponential
    cx_mat44 complexElement = expmat(exponentialArgument);

    return real(complexElement); // This matrix should be real
}
