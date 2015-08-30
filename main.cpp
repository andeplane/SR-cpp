#include <iostream>
#include <armadillo>
#include <cmath>

#include "lorentzgenerator.h"
#include "sr.h"

using namespace arma;
using namespace std;

int main()
{
    SR::vec v(1,2,3,4);
    cout << "This vector: " << endl << v << " has length " << v.length() << endl;
    SR::boost(v, {1,0,0}, 0.2);
    cout << "Boosted vector: " << endl << v << " has length " << v.length() << endl;
    SR::rotate(v, {1,0,0}, M_PI/2);
    cout << "Rotated vector: " << endl << v << " has length " << v.length() << endl;


    return 0;
}

