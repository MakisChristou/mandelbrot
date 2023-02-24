#ifndef COLOR_H
#define COLOR_H

class Color
{   
    public:
        int R;
        int G;
        int B;

        Color()
        {
            R = 0;
            G = 0;
            B = 0;
        }

        Color(int r, int g, int b)
        {
            R = r;
            G = g;
            B = b;
        }
};

#endif