RK-LCE:  a in-place LCE data structure
===
Author: Nicola Prezza (nicolapr@gmail.com)

The class rk_lce.hpp contains a LCE data structure taking exactly the same size of the plain text T[1,...,n] (plus a constant number of memory words). The structure supports constant-time access to any text character, O(log n)-time LCE queries with high probability and is constructed in O(n) time. 

The structure is based on Rabin-Karp fingerprinting with Mersenne primes. Since for efficiency reasons in this implementation the prime is fixed (2^127^-1), the structure returns a wrong result with low probability (<2^-120^) only if the text is perturbed by some noise; for real-case texts this can be assumed to be true, but an artificial text can be created that causes the structure to always fail.

TODO: serialization/deserialization.
