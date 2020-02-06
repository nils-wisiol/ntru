# Toy Implementation of NTRU Crypto Systems

[![Build Status](https://travis-ci.org/nils-wisiol/ntru.svg?branch=master)](https://travis-ci.org/nils-wisiol/ntru)
[![Build Status](https://codecov.io/gh/nils-wisiol/ntru/branch/master/graphs/badge.svg)](https://codecov.io/gh/nils-wisiol/ntru/branch/master)

This repository contains a study-toy implementation of the Streamlined NTRU Prime public key crypto system.
Do not use in production.

## API

The implementation mainly consists of two parts.

### polymod

The `polymod` package ships the `Poly` class which models polynomials from ℤ[X]/p for polynomials p and ℤ_N[X]/p for polynomials p and natural numbers N.
It provides basic implementations of polynomial arithmetic, including the Extended Euclidean Algorithm to find inverses of given polynomials.

For polynomial multiplication, a very basic implementation of the sparse multiplication method is used. Note however, while this dramatically increases performance,
polynomials are still stored in a non-sparse way and hence the asymptotic run time of the multiplication is still O(n²).

`Poly` tries to implement a human-readable `__str__` to ease debugging.

### ciphers

The `ciphers` package ships the `StreamlinedNTRUPrime` class which implements mainly three functions,

1. `generate_keys()` returns a tuple `(pk, sk)` with public and secret key, where the public key is a `Poly` instance and the secret key is a tuple of two `Poly` instances;
1. `encrypt(m, pk)` encrypts message `m` (must be from a particular subset of `Poly` instances) with public key `pk`;  
1. `decrypt(c, sk)` attempts to decrypt cipher text `c` using secret key `sk`.

### Known Limitations

Again, this is code for studying, not for security.
Note that decrypt is not secure against chosen-cipher-text-attacks; the behavior of this implementation on invalid cipher text-secret key pairs is undefined.
For more detailed information, check TODOs in the code.  

## Usage

I made extensive use of python doc tests to verify the correct behavior of most sub-routines.

To get started, install python3 and numpy. Then run the example:

    apt-get install python3-numpy  # or change to your package manager
    git clone  # insert repository URL here
    cat example.py  # look at the code
    python3 -m example

