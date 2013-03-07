#pragma once
#include <cmath>
#include <vector>

//-------------------------------------------------------------------------

template <class T>
inline T saturate(const T &x) { return (x>1) ? 1 : ((x<0) ? 0 : x); };

//-------------------------------------------------------------------------

struct rgbf
{  
    float r, g, b;

    rgbf() { r=0; g=0; b=0; }
    rgbf(float _r, float _g, float _b) { r=_r; g=_g; b=_b; }
    
    inline rgbf operator + (const rgbf& rhs) const
    { return rgbf(r + rhs.r, g + rhs.g, b + rhs.b); } 
    inline rgbf operator - (const rgbf& rhs) const
    { return rgbf(r - rhs.r, g - rhs.g, b - rhs.b); } 
    inline rgbf operator * (const float rhs) const
    { return rgbf(r*rhs, g*rhs, b*rhs); } const
    inline rgbf operator / (const float rhs) const
    { return rgbf(r/rhs, g/rhs, b/rhs); } 
    inline rgbf& operator+=(const rgbf& rhs) 
    { r += rhs.r; g += rhs.g; b += rhs.b; return *this; }
    inline rgbf& operator/=(const float rhs) 
    { r /= rhs; g /= rhs; b /= rhs; return *this; }

    inline float length() const
    { return sqrt(r*r + g*g + b*b); }
    inline rgbf normalize(const float bias=0) const
    {
        float l = this->length();
        return rgbf((r+bias) / l, (g+bias) / l, (b+bias) / l);
    }
    inline float dot(const rgbf& v) const
    { return r*v.r + g*v.g + b*v.b; }
}; 

//-------------------------------------------------------------------------