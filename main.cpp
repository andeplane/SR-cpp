#include <iostream>
#include <armadillo>
#include <cmath>

#include "lorentzgenerator.h"
#include "sr.h"

int main()
{
    std::vector<SR::vec> events;
    events.resize(4);

    SR::vec &eventA = events.at(0);
    SR::vec &eventB = events.at(1);
    SR::vec &eventC = events.at(2);
    SR::vec &eventD = events.at(3);

    double x1 = 1.0;
    double beta = 0.999;

    double dt = 1.0;
    double t0 = 0;
    double t1 = t0 + dt;
    double t2 = t1 + dt;
    double t3 = t2 + dt;

    eventB = SR::vec(t1,t1*beta + x1,0,0) + eventA;
    eventC = SR::vec(t2,t2*beta + x1,0,0) + eventA;
    eventD = SR::vec(t2 + (eventC.x()/beta),0,0,0) + eventA;

    cout << "We have four events described in reference frame S as: " << endl << "A:" << eventA << endl << "B:" << eventB << endl << "C:" << eventC << endl << "D: " << eventD;
    cout << "ds(A,B)^2 = " << (eventB - eventA).lengthSquared() << endl;
    cout << "ds(B,C)^2 = " << (eventC - eventB).lengthSquared() << endl;
    cout << "ds(C,D)^2 = " << (eventD - eventC).lengthSquared() << endl;


    SR::boost(eventA, {1,0,0}, beta);
    SR::boost(eventB, {1,0,0}, beta);
    SR::boost(eventC, {1,0,0}, beta);
    SR::boost(eventD, {1,0,0}, beta);

    cout << "In the reference frame S' moving with a velocity v=0.5 in the x-direction, these events look like: "  << endl << "A:" << eventA << endl << "B:" << eventB << endl << "C:" << eventC << endl << "D: " << eventD;
    cout << "ds(A,B)^2 = " << (eventB - eventA).lengthSquared() << endl;
    cout << "ds(B,C)^2 = " << (eventC - eventB).lengthSquared() << endl;
    cout << "ds(C,D)^2 = " << (eventD - eventC).lengthSquared() << endl;

    return 0;
}

