#ifndef  BASE_INT_VEC3_H
#define  BASE_INT_VEC3_H

/****************************************************************
 * Include files
 ****************************************************************/

//#define WIN32_LEAN_AND_MEAN

#include <stdexcept>
//using std::exception;
using std::domain_error;
//using std::invalid_argument;

/****************************************************************
 * Class definition
 ****************************************************************/

template<typename T> 
struct base_int_vec3
{
    int e[3];

    /* Constructors */
    base_int_vec3<T>(); // Default constructor
    base_int_vec3<T>(int e0, int e1, int e2);

    base_int_vec3<T>& operator+=(const base_int_vec3<T>&);
    base_int_vec3<T>& operator-=(const base_int_vec3<T>&);
    base_int_vec3<T>& operator*=(const int   );
    base_int_vec3<T>& operator/=(const int   );
    int&   operator[](const int   );
    int    operator[](const int   ) const;

    base_int_vec3<T> operator+(const base_int_vec3<T>&) const;
    base_int_vec3<T> operator-(const base_int_vec3<T>&) const;
    int   operator*(const base_int_vec3<T>&) const; // Scalar product
    base_int_vec3<T> operator&(const base_int_vec3<T>&) const; // Crossproduct
    base_int_vec3<T> operator*(const int   ) const;
    base_int_vec3<T> operator/(const int   ) const;

    int sqr_length() const;
};

typedef  base_int_vec3<int>  ivec3;

/****************************************************************
 * Class related functions
 ****************************************************************/

template<typename T>
base_int_vec3<T> operator* (const int, const base_int_vec3<T>&);

/****************************************************************
 * Public functions
 ****************************************************************/

// Default constructor
template<typename T>
base_int_vec3<T>::base_int_vec3()
{
    e[0] = e[1] = e[2] = 0;
}

template<typename T>
base_int_vec3<T>::base_int_vec3(int e0, int e1, int e2)
{
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator+=(const base_int_vec3<T>& rhs)
{
    for (int i = 0; i < 3; i++) e[i] += rhs.e[i];
    return *this;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator-=(const base_int_vec3<T>& rhs)
{
    for (int i = 0; i < 3; i++) e[i] -= rhs.e[i];
    return *this;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator*=(const int k)
{
    for (int i = 0; i < 3; i++) e[i] *= k;
    return *this;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator/=(const int den)
{
    if (!den) throw domain_error("Trying to divide a vector by zero");
    for (int i = 0; i < 3; i++) e[i] /= den;
    return *this;
}

template<typename T>
int& base_int_vec3<T>::operator[](const int i)
{
    return e[i];
}

template<typename T>
int base_int_vec3<T>::operator[](const int i) const
{
    return e[i];
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator+(const base_int_vec3<T>& rhs) const
{
    return base_int_vec3<T>(e[0] + rhs.e[0],
                 e[1] + rhs.e[1],
                 e[2] + rhs.e[2]);
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator-(const base_int_vec3<T>& rhs) const
{
    return base_int_vec3<T>(e[0] - rhs.e[0],
                 e[1] - rhs.e[1],
                 e[2] - rhs.e[2]);
}

template<typename T>
int base_int_vec3<T>::operator*(const base_int_vec3<T>& rhs) const
{
    return e[0]*rhs.e[0] + e[1]*rhs.e[1] + e[2]*rhs.e[2];
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator&(const base_int_vec3<T>& rhs) const
{
    return base_int_vec3<T>(e[1]*rhs.e[2] - e[2]*rhs.e[1],
                 e[2]*rhs.e[0] - e[0]*rhs.e[2],
                 e[0]*rhs.e[1] - e[1]*rhs.e[0]);
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator*(const int k) const
{
    return base_int_vec3<T>(e[0]*k, e[1]*k, e[2]*k);
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator/(const int den) const
{
    if(!den) throw domain_error("Not defined to divide with zero");
    return base_int_vec3<T>(e[0]/den, e[1]/den, e[2]/den);
}

template<typename T>
int base_int_vec3<T>::sqr_length() const
{
    return (e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
}

/****************************************************************
 * Class related functions
 ****************************************************************/

template<typename T>
base_int_vec3<T> operator*(const int lhs, const base_int_vec3<T>& rhs)
{
    return base_int_vec3<T>(lhs*rhs.e[0], lhs*rhs.e[1], lhs*rhs.e[2]);
}

#endif  /* base_int_vec3<T>_H */
