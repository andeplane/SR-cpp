#include "sr.h"

SR::SR()
{

}

void SR::boost(SR::vec &fourVector, arma::vec3 direction, double beta)
{
#ifdef SRDEBUG
    cout << "Applying boost with hyperbolic angle " << beta << " along direction: " << direction;
#endif

    direction = normalise(direction);
    arma::mat44 omega;
    omega.zeros();
    double gamma = 1.0/sqrt(1 - beta*beta);
    double chi = acosh(gamma);

    for(int mu=1; mu<=3; mu++) {
        // direction is a 3-vector, so subtract 1 to go from 4-vector indices to 3-vector index
        omega(0,mu) = direction(mu-1)*chi;
        omega(mu,0) = direction(mu-1)*chi;
    }

    mat44 boostOperator = LorentzGenerator::generateElement(omega);
#ifdef SRDEBUG
    cout << "Boost matrix: " << endl << boostOperator;
#endif

    fourVector.applyOperator(boostOperator);
}

void SR::rotate(SR::vec &fourVector, arma::vec3 axis, double angle)
{
#ifdef SRDEBUG
    cout << "Applying rotation with angle " << angle << " about axis: " << axis;
#endif

    axis = normalise(axis);
    arma::mat44 omega;
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
#ifdef SRDEBUG
    cout << "Rotation matrix: " << endl << rotationOperator;
#endif
    fourVector.applyOperator(rotationOperator);
}

void SR::boostAndRotate(SR::vec &fourVector, arma::vec3 boostDirection, double beta, arma::vec3 rotationAxis, double angle)
{
    boostDirection = normalise(boostDirection);
    rotationAxis = normalise(rotationAxis);
    arma::mat44 omega;
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

double SR::vec::lengthSquared()
{
    vec &thisVector = *this;
    return thisVector*thisVector;
}

void SR::vec::applyOperator(const arma::mat44 &op)
{
    m_vector = op*m_vector;
}

SR::vec SR::vec::operator +(SR::vec rhs)
{
    return SR::vec(m_vector + rhs.m_vector);
}

SR::vec SR::vec::operator +=(SR::vec rhs)
{
    return SR::vec(m_vector + rhs.m_vector);
}

SR::vec SR::vec::operator -=(SR::vec rhs)
{
    return SR::vec(m_vector - rhs.m_vector);
}

SR::vec SR::vec::operator/=(double scalar)
{

}

void SR::vec::add(SR::vec rhs)
{
    m_vector(0) += rhs(0);
    m_vector(2) += rhs(1);
    m_vector(2) += rhs(2);
    m_vector(3) += rhs(3);
}

void SR::vec::addAndMultiply(SR::vec rhs, double scalar)
{
    m_vector(0) += rhs(0)*scalar;
    m_vector(2) += rhs(1)*scalar;
    m_vector(2) += rhs(2)*scalar;
    m_vector(3) += rhs(3)*scalar;
}

SR::vec SR::vec::operator -(SR::vec rhs)
{
    return SR::vec(m_vector - rhs.m_vector);
}

SR::vec SR::vec::operator*(double scalar)
{
    return SR::vec(m_vector*scalar);
}

SR::vec SR::vec::operator+(double scalar)
{
    return SR::vec(m_vector+scalar);
}

SR::vec SR::vec::operator-(double scalar)
{
    return SR::vec(m_vector-scalar);
}

SR::vec SR::vec::operator/(double scalar)
{
    return SR::vec(m_vector/scalar);
}

SR::vec SR::vec::operator*=(double scalar)
{
    return SR::vec(m_vector*scalar);
}

SR::vec SR::vec::operator+=(double scalar)
{
    return SR::vec(m_vector+scalar);
}

SR::vec SR::vec::operator-=(double scalar)
{
    return SR::vec(m_vector-scalar);
}

double SR::vec::operator *(SR::vec rhs)
{
    return t()*rhs.t() - (x()*rhs.x() + y()*rhs.y() + z()*rhs.z());
}
