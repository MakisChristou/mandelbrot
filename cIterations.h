#ifndef CITERATIONS_H
#define CITERATIONS_H

#include "complex.h"

class cIterations
{
    public:
        Complex c;
        int n;

        cIterations()
        {

        }

        cIterations(int n, Complex c)
        {
            this->n = n;
            this->c = c;
        }

};

#endif