#include <cassert>

#include "permutation.hh"

#define MAX_PERM_LEN 20U

static void validate_perm_length(const unsigned permLen)
{
    // Here, we are assuming that unsigned long is at least a 64-bit number
    assert(sizeof(unsigned long) >= 8U);

    // 20! is less than 2^64 and 21! is larger
    assert(permLen <= MAX_PERM_LEN);
}

unsigned long factorial(const unsigned n)
{
    validate_perm_length(n);
    unsigned long f = 1UL;
    for (unsigned i=2U; i<=n; ++i)
        f *= i;
    return f;
}

unsigned long permutationNumber(const unsigned *permutation,
                                const unsigned permLen)
{
    assert(permLen);
    assert(permutation);

    if (permLen == 1U)
    {
        assert(permutation[0] == 0U);
        return 0UL;
    }
    validate_perm_length(permLen);
    unsigned long stride = 1UL;
    unsigned long permNum = 0UL;
    for (unsigned inum=0; inum<permLen-1; ++inum)
    {
        unsigned i = 0, icnt = 0;
        for (; inum != permutation[i] && i < permLen; ++i)
            if (inum < permutation[i])
                ++icnt;
        if (i >= permLen)
            assert(!"In permutationNumber: invalid permuted index");
        permNum += icnt*stride;
        stride *= (permLen - inum);
    }
    return permNum;
}

void orderedPermutation(unsigned long permutationNumber,
                        unsigned *permutation, const unsigned permLen)
{
    if (!permLen)
        assert(!"In orderedPermutation: empty permutation array");
    assert(permutation);
    if (permLen == 1U)
    {
        assert(permutationNumber == 0UL);
        permutation[0] = 0UL;
        return;
    }
    validate_perm_length(permLen);
    unsigned long strides[MAX_PERM_LEN];
    unsigned long stride = 1UL;
    for (unsigned inum=0; inum<permLen; ++inum)
    {
        strides[inum] = stride;
        stride *= (permLen - inum);
    }
    if (permutationNumber >= stride)
        assert(!"In orderedPermutation: permutation number is out of range");

    unsigned posnum[MAX_PERM_LEN];
    const unsigned pmLm1 = permLen-1;
    for (unsigned itmp=pmLm1; itmp>0; --itmp)
    {
        const unsigned inum = itmp - 1;
        posnum[inum] = permutationNumber / strides[inum];
        permutationNumber %= strides[inum];
    }

    for (unsigned i=0; i<permLen; ++i)
        permutation[i] = pmLm1;

    permutation[posnum[0]] = 0U;
    for (unsigned ipos=1; ipos<pmLm1; ++ipos)
    {
        const unsigned newpos = posnum[ipos];
        unsigned cnt = 0U;
        for (unsigned j=0; j<permLen; ++j)
            if (permutation[j] > ipos)
                if (cnt++ == newpos)
                {
                    permutation[j] = ipos;
                    break;
                }
    }
}
