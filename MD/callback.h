#ifndef  CALLBACK_H
#define  CALLBACK_H

////////////////////////////////////////////////////////////////
// CLASS
////////////////////////////////////////////////////////////////

template<typename T>
struct callback
{
public:
    // Variables
    T func;
    void *param;

    // Constructor
    callback();
    callback(T function);
    callback(T function, void *parameter);
};

////////////////////////////////////////////////////////////////
// PUBLIC MEMBER FUNCTIONS
////////////////////////////////////////////////////////////////

// Constructor
template<typename T>
callback<T>::callback()
{
    func  = 0;
    param = 0;
}

template<typename T>
callback<T>::callback(T function)
{
    func  = function;
    param = 0       ;
}

template<typename T>
callback<T>::callback(T function, void *parameter)
{
    func  = function ;
    param = parameter;
}

#endif  /* CALLBACK_H */
