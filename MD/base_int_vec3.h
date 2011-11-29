#ifndef  BASE_INT_VEC3_H
#define  BASE_INT_VEC3_H

/****************************************************************
 * Include files
 ****************************************************************/

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
public:
    T e[3];

    /* Constructors */
    base_int_vec3<T>(); // Default constructor
    base_int_vec3<T>(T e0, T e1, T e2);

    base_int_vec3<T>& operator+=(const base_int_vec3<T>&);
    base_int_vec3<T>& operator-=(const base_int_vec3<T>&);
    //base_int_vec3<T>& operator&=(const base_int_vec3<T>&); // Yet to define
    base_int_vec3<T>& operator*=(const T                );
    base_int_vec3<T>& operator/=(const T                );
    T&                operator[](const T)      ;
    T                 operator[](const T) const;

    base_int_vec3<T> operator-() const;

    base_int_vec3<T> operator+(const base_int_vec3<T>&) const;
    base_int_vec3<T> operator-(const base_int_vec3<T>&) const;
    T                operator*(const base_int_vec3<T>&) const; // Scalar product
    base_int_vec3<T> operator&(const base_int_vec3<T>&) const; // Crossproduct
    base_int_vec3<T> operator*(const T                ) const;
    base_int_vec3<T> operator/(const T                ) const;

    T sqr_length() const;
};

typedef  base_int_vec3<int>  ivec3;

/****************************************************************
 * Class related functions
 ****************************************************************/

template<typename T>
base_int_vec3<T> operator* (const T, const base_int_vec3<T>&);

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
base_int_vec3<T>::base_int_vec3(T e0, T e1, T e2)
{
    e[0] = e0;
    e[1] = e1;
    e[2] = e2;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator+=(const base_int_vec3<T>& rhs)
{
    for (T i = 0; i < 3; i++) e[i] += rhs.e[i];
    return *this;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator-=(const base_int_vec3<T>& rhs)
{
    for (T i = 0; i < 3; i++) e[i] -= rhs.e[i];
    return *this;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator*=(const T k)
{
    for (T i = 0; i < 3; i++) e[i] *= k;
    return *this;
}

template<typename T>
base_int_vec3<T>& base_int_vec3<T>::operator/=(const T den)
{
    if (!den) throw domain_error("Trying to divide a vector by zero");
    for (T i = 0; i < 3; i++) e[i] /= den;
    return *this;
}

template<typename T>
T& base_int_vec3<T>::operator[](const T i)
{
    return e[i];
}

template<typename T>
T base_int_vec3<T>::operator[](const T i) const
{
    return e[i];
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator-() const
{
    return base_int_vec3<T>(-e[0],
                            -e[1],
                            -e[2]);
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
T base_int_vec3<T>::operator*(const base_int_vec3<T>& rhs) const
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
base_int_vec3<T> base_int_vec3<T>::operator*(const T k) const
{
    return base_int_vec3<T>(e[0]*k, e[1]*k, e[2]*k);
}

template<typename T>
base_int_vec3<T> base_int_vec3<T>::operator/(const T den) const
{
    if(!den) throw domain_error("Not defined to divide with zero");
    return base_int_vec3<T>(e[0]/den, e[1]/den, e[2]/den);
}

template<typename T>
T base_int_vec3<T>::sqr_length() const
{
    return (e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
}

/****************************************************************
 * Class related functions
 ****************************************************************/

template<typename T>
base_int_vec3<T> operator*(const T lhs, const base_int_vec3<T>& rhs)
{
    return base_int_vec3<T>(lhs*rhs.e[0], lhs*rhs.e[1], lhs*rhs.e[2]);
}

#endif  /* BASE_INT_VEC3_H */
