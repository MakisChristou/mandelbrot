#ifndef COMPLEX_H
#define COMPLEX_H

class Complex
{
    public:
        long double Re;
        long double Im;

        Complex()
        {

        }

        Complex(long double Re, long double Im)
        {
            this->Re = Re;
            this->Im = Im;
        }
};

#endif
