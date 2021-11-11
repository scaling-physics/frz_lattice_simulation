#ifndef definitions_H_INCLUDED
#define definitions_H_INCLUDED
# include <vector>		// std::vector<>
# include <array>
# include <memory>

#include <algorithm>
#include <functional>

# define NDEBUG //comment out to turn on assert. costs ~10%
# include <assert.h>	// for assert()



template <typename T>
inline std::vector<T> operator+(const std::vector<T> &A, const std::vector<T> &B)//add vectors, returns new vector
{
    assert(A.size()==B.size());

    std::vector<T> result;
    result.reserve( A.size()); // preallocate memory

    std::transform(A.begin(), A.end(), B.begin(), std::back_inserter(result), std::plus<T>());
    return result;
}

template <typename T>
inline std::vector<T> & operator+=(std::vector<T> &A, const std::vector<T> &B)//+= add vectors, returns by reference
{
    assert(A.size()==B.size());

    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::plus<T>());
    return A;
}


template <typename T>
inline std::vector<T> operator+(const std::vector<T> &A, T &B)//add scalar, returns new vector
{
    std::vector<T> AB;
    AB.reserve( A.size()); // preallocate memory
    AB=A;
    for(T& d : AB)
    {
        d += B;
    }
    return AB;
}

template <typename T>
inline std::vector<T> & operator+=(std::vector<T> &A, T &B)//+= add scalar, returns by reference
{
    for(T& d : A)
    {
        d += B;
    }
    return A;
}

template <typename T, size_t N>
inline T operator*(const std::array<T,N> &A, const std::array<T,N> &B)//dot product
{
    assert(A.size()==B.size());

    T out=0;
    for (int i=0; i<N; ++i)
    {
        out += A[i]*B[i];
    }
    return out;
}

inline static float fast_exp (float x)
{
    //if(x<-1000) return 0.0;
    volatile union
    {
        float f;
        unsigned int i;
    } cvt;

    /* exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1 */
    float t = x * 1.442695041f;
    float fi = floorf (t);
    float f = t - fi;
    int i = (int)fi;
    cvt.f = (0.3371894346f * f + 0.657636276f) * f + 1.00172476f; /* compute 2^f */
    cvt.i += (i << 23);                                          /* scale by 2^i */
    return cvt.f;
}

#endif // definitions_H_INCLUDED

