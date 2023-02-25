#ifndef COLOR_H
#define COLOR_H

#include <cstdint>

class Color
{   
    public:
        uint8_t R;
        uint8_t G;
        uint8_t B;

        Color()
        {
            R = 0;
            G = 0;
            B = 0;
        }

        Color(uint8_t r, uint8_t g, uint8_t b)
        {
            R = r;
            G = g;
            B = b;
        }
};

#endif