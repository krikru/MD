#ifndef PREPROCESSING_H
#define PREPROCESSING_H

////////////////////////////////////////////////////////////////
// STATIC ASSERT
////////////////////////////////////////////////////////////////

#define JOIN2(x, y) x##y
#define JOIN(x, y) JOIN2(x, y)

template<bool>  struct COMPILE_TIME_ERROR;
template<>      struct COMPILE_TIME_ERROR<false>{};
template<int x> struct STATIC_ASSERT_TEST{};
#define STATIC_ASSERT(exp)                                         \
    typedef STATIC_ASSERT_TEST<sizeof(COMPILE_TIME_ERROR<!(exp)>)> \
    JOIN(static_assert_typedef_at_line_, __LINE__)
/*
 * #define STATIC_ASSERT(exp)                    \
 *     do {                                      \
 *         enum { assert_static__ = 1/(exp) };   \
 *     } while (0)
 */
/*
 * #define STATIC_ASSERT(exp) ((void)sizeof(char[1 - 2*!(exp)]))
 */

////////////////////////////////////////////////////////////////
// META FUNCTIONS
////////////////////////////////////////////////////////////////

/*********/
/* LOG_B */
/*********/
/* Calculates floor(log_b(n)) */
template<uint64 b, uint64 n>
struct LOG_B { /* Rounds down */
    /* Check for arguments outside definition set */
    STATIC_ASSERT(b != 0 && b != 1);
    STATIC_ASSERT(n != 0);
    const static uint value = LOG_B<b, n/b + (n < b)>::value + (n >= b);
};
template<uint64 b>
struct LOG_B<b, 1> {STATIC_ASSERT(b > 1); const static uint value = 0;};

/*******/
/* POW */
/*******/
/* Calculates b^n */
template<uint64 b, uint n>
struct POW {
private:
    const static uint64 pre_pow = b > 1 ? POW<b, n/2>::value : b;
    const static uint64 factor  = pre_pow * (n & 1 ? b : 1);
    /* Check for overflow */
    STATIC_ASSERT(pre_pow*factor/(b ? factor : 1) == pre_pow);
public: 
    const static uint64 value = pre_pow * factor;
};
template<uint64 b> struct POW<b, 0> {const static uint64 value = 1;};

#endif /* PREPROCESSING_H */
