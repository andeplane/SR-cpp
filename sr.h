#ifndef SR_H
#define SR_H
#include <vector>
#include <iostream>
#include <armadillo>
#include "lorentzgenerator.h"

class SR
{
private:
    SR();
    SR(SR const&) = delete;
    void operator=(SR const&)  = delete;
public:
    static SR& getInstance() {
        static SR instance;
        return instance;
    }

    class vec {
    public:
        arma::vec4 m_vector;
        vec() { }
        vec(double x0, double x1, double x2, double x3) { m_vector(0) = x0; m_vector(1) = x1; m_vector(2) = x2; m_vector(3) = x3; }
        double length();
        void applyOperator(const mat44 &op);
        friend std::ostream& operator<<(std::ostream &stream, const vec vector);
    private:
    };
    static void boost(vec &fourVector, vec3 direction, double beta);
    static void rotate(vec &fourVector, vec3 axis, double angle);
    static void boostAndRotate(vec &fourVector, vec3 boostDirection, double beta, vec3 rotationAxis, double angle);
};

#endif // SR_H
