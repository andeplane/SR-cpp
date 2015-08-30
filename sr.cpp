#include "sr.h"

SR::SR()
{

}

void SR::boost(SR::vec &fourVector, vec3 direction, double beta)
{
    direction = normalise(direction);
    mat44 omega;
    omega.zeros();
    double gamma = 1.0/sqrt(1 - beta*beta);
    double chi = acosh(gamma);

    for(int mu=1; mu<=3; mu++) {
        // direction is a 3-vector, so subtract 1 to go from 4-vector indices to 3-vector index
        omega(0,mu) = direction(mu-1)*chi;
        omega(mu,0) = direction(mu-1)*chi;
    }

    mat44 boostOperator = LorentzGenerator::generateElement(omega);
    fourVector.applyOperator(boostOperator);
}

void SR::rotate(SR::vec &fourVector, vec3 axis, double angle)
{
    axis = normalise(axis);
    mat44 omega;
    omega.zeros();

    for(int mu=1; mu<=3; mu++) {
        for(int nu=1; nu<=3; nu++) {
            if(mu==nu) continue;
            // Sum of all indices is 1+2+3 = 6.
            // So the third index is then 6 - mu - nu = 6 - (mu+nu)
            int axisOfRotation = 6 - (mu+nu);
            omega(mu,nu) = angle*axis(axisOfRotation-1);
            omega(nu,mu) = -angle*axis(axisOfRotation-1);
        }
    }

    mat44 rotationOperator = LorentzGenerator::generateElement(omega);
    fourVector.applyOperator(rotationOperator);
}

void SR::boostAndRotate(SR::vec &fourVector, vec3 boostDirection, double beta, vec3 rotationAxis, double angle)
{
    boostDirection = normalise(boostDirection);
    rotationAxis = normalise(rotationAxis);
    mat44 omega;
    omega.zeros();

    double gamma = 1.0/sqrt(1 - beta*beta);
    double chi = acosh(gamma);

    for(int mu=1; mu<=3; mu++) {
        // direction is a 3-vector, so subtract 1 to go from 4-vector indices to 3-vector index
        omega(0,mu) = boostDirection(mu-1)*chi;
        omega(mu,0) = boostDirection(mu-1)*chi;

        for(int nu=1; nu<=3; nu++) {
            if(mu==nu) continue;
            // Sum of all indices is 1+2+3 = 6.
            // So the third index is then 6 - mu - nu = 6 - (mu+nu)
            int axisOfRotation = 6 - (mu+nu);
            omega(mu,nu) = angle*rotationAxis(axisOfRotation);
            omega(nu,mu) = -angle*rotationAxis(axisOfRotation);
        }
    }

    mat44 rotationOperator = LorentzGenerator::generateElement(omega);
    fourVector.applyOperator(rotationOperator);

}

std::ostream& operator<<(std::ostream &stream, const SR::vec vector) {
    return stream << vector.m_vector;
}

double SR::vec::length()
{
    double len = m_vector(0)*m_vector(0);
    for(int i=1; i<=3; i++) {
        len -= m_vector(i)*m_vector(i);
    }

    return len;
}

void SR::vec::applyOperator(const mat44 &op)
{
    m_vector = op*m_vector;
}
