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
    private:
        arma::vec4 m_vector;
    public:
        vec() { m_vector.zeros(); }
        vec(arma::vec4 v) { m_vector = v; }
        vec(double t, arma::vec3 r) { m_vector[0] = t; m_vector[1] = r(0); m_vector[2] = r(1); m_vector[3] = r(2);}
        vec(double x0, double x1, double x2, double x3) { m_vector(0) = x0; m_vector(1) = x1; m_vector(2) = x2; m_vector(3) = x3; }
        double lengthSquared();
        void applyOperator(const arma::mat44 &op);
        vec operator +(vec rhs);
        vec operator +=(vec rhs);
        vec operator -=(vec rhs);
        vec operator -(vec rhs);
        vec operator*(double scalar);
        vec operator+(double scalar);
        vec operator-(double scalar);
        vec operator/(double scalar);
        vec operator*=(double scalar);
        vec operator+=(double scalar);
        vec operator-=(double scalar);
        vec operator/=(double scalar);
        double operator *(vec rhs);
        double operator[](int index) const { return m_vector(index); }
        double operator()(int index) const { return m_vector(index); }

        void add(vec rhs);
        void addAndMultiply(vec rhs, double scalar);
        double t() { return m_vector[0]; }
        double x() { return m_vector[1]; }
        double y() { return m_vector[2]; }
        double z() { return m_vector[3]; }
        double E() { return m_vector[0]; }
        double px() { return m_vector[1]; }
        double py() { return m_vector[2]; }
        double pz() { return m_vector[3]; }
        arma::vec3 p() { arma::vec3 v; v(0) = m_vector[1]; v(1) = m_vector[2]; v(2) = m_vector[3]; return v; }
        arma::vec3 r() { return p(); }

        friend std::ostream& operator<<(std::ostream &stream, const vec vector);
    };
    static void boost(vec &fourVector, arma::vec3 direction, double beta);
    static void rotate(vec &fourVector, arma::vec3 axis, double angle);
    static void boostAndRotate(vec &fourVector, vec3 boostDirection, double beta, vec3 rotationAxis, double angle);
};

#endif // SR_H
