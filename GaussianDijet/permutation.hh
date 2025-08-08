#ifndef PERMUTATION_HH_
#define PERMUTATION_HH_

// Simple factorial. Will generate a run-time error if n!
// is larger than the largest unsigned long.
unsigned long factorial(unsigned n);

// On output, array "permutation" will be filled by the permutation
// of numbers from 0 to permLen-1 which correspond to the given
// permutation number. The input permutation number must be
// less than factorial(permLen).
void orderedPermutation(unsigned long permutationNumber,
                        unsigned *permutation, unsigned permLen);

// A mapping from a permuted set into a linear sequence.
// factorial(permLen) should be less than the largest unsigned long.
unsigned long permutationNumber(const unsigned *permutation,
                                unsigned permLen);

#endif // PERMUTATION_HH_
