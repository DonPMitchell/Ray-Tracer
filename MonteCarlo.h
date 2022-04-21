//
//  Pseudo-random number generators
//  D.P. Mitchell  2018/09/11.
//
#pragma once

extern "C"
{
   unsigned __int64 __emulu(unsigned int a, unsigned int b);
}

#pragma intrinsic(__emulu)

extern unsigned     ML_nMarsagliaX, ML_nMarsagliaC;
extern unsigned     RandomLessThanN(unsigned N);

inline unsigned
RandomUnsigned()
{
    unsigned __int64 n;

    n = __emulu(1965537969, ML_nMarsagliaX);
    n = n + ML_nMarsagliaC;
    ML_nMarsagliaX = unsigned(n & 0xFFFFFFFF);
    ML_nMarsagliaC = unsigned(n >> 32);
    return ML_nMarsagliaX;
}

inline float
RandomFloat()
{
    return float(RandomUnsigned()) * (1.0f/4294967296.0f);
}

inline double
RandomDouble()
{
    return (double(RandomUnsigned()) + double(RandomUnsigned())*(1.0/4294967296.0))*(1.0/4294967296.0);
}

inline void
ML_InitializeRandom(unsigned nX = 886459, unsigned nC = 361290869)
{
    ML_nMarsagliaX = nX;
    ML_nMarsagliaC = nC;
}
//
//  Randomly shuffle an array
//
template <class T>
void
RandomShuffle(T rgtItem[], int nItems)
{
    int i, j;
    T tTemp;

    for (j = nItems - 1; j > 0; --j) {
        i = RandomLessThanN(j + 1);
        tTemp = rgtItem[i];
        rgtItem[i] = rgtItem[j];
        rgtItem[j] = tTemp;
    }
}
//
//  True random numbers
//
extern unsigned TrueRandomUnsigned();
//
//  Random point inside disk
//
extern void RandomInsideDisk(double rgf[2]);