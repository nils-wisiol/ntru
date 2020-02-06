from math import gcd, ceil, sqrt

from numpy import around, int32, zeros
from numpy.random import RandomState

from polymod import Poly


def isprime(n):
    """
    >>> isprime(1109)
    True

    >>> isprime(10**6)
    False
    """
    # this is slow and basic, but saves the sympy dependency
    # alternatively, from sympy import isprime
    for i in range(2, ceil(sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True


class StreamlinedNTRUPrime:
    """
    Training implementation of Streamlined NTRU Prime. Do not use in production.
    
    Following https://ntruprime.cr.yp.to/nist/ntruprime-20190330.pdf
    
    # Used Notation (in Comments and Variable Names)
    R = Z[x]/(x^p - x - 1)
    R/3 = (Z/3)[x]/(x^p - x - 1)
    R/q = (Z/q)[x]/(x^p - x - 1)
    Short = set of small weight-w elements of R

    >>> p, q, w = 23, 113, 5
    >>> cipher = StreamlinedNTRUPrime(p, q, w, seed=1337)
    >>> pk, sk = cipher.generate_keys()
    >>> message = cipher.random_small_poly(w, None, cipher.modulus_r, RandomState(42))
    >>> message == cipher.decrypt(cipher.encrypt(message, pk), sk)
    True
    """
    
    @staticmethod
    def co_prime(a, b):
        return gcd(a, b) == 1

    @staticmethod
    def random_small_poly(weight, N, modulus, random_state=RandomState()):
        """
        Returns a uniformly randomly chosen small polynomial from ℤ_N[x]/modulus (if N not None) or
        ℤ[x]/modulus (if N is None) of weight w that is invertible in ℤ_inv_N[x]/modulus.
        Small polynomials only have coefficients in {-1,0,1}.
        Weight-w polynomials have exactly w non-zero coefficients.

        >>> StreamlinedNTRUPrime.random_small_poly(2, None, Poly([0, 0, 0, 1]), RandomState(42)).deg < 4
        True

        >>> StreamlinedNTRUPrime.random_small_poly(13, None, Poly.xn1(n=44), RandomState(42)).weight
        13
        """
        c_values = 2 * random_state.randint(2, size=weight) - 1  # "small" non-zero coefficients, i.e. -1, 1
        c_pos = random_state.choice(modulus.deg, weight, replace=False)
        cs = zeros((modulus.deg + 1), dtype=int32)
        for i in range(weight):
            cs[c_pos[i]] = c_values[i]
    
        return Poly(cs, N=N, modulus=modulus)  # TODO confirm uniformity on the ring
    
    def _random_small_poly_invertible(self, weight, N, modulus, inv_N):
        """
        Returns a uniformly randomly chosen small polynomial from ℤ_N[x]/modulus (if N not None) or
        ℤ[x]/modulus (if N is None) of weight w that is invertible in ℤ_inv_N[x]/modulus.
        Small polynomials only have coefficients in {-1,0,1}.
        Weight-w polynomials have exactly w non-zero coefficients.

        >>> modulus = Poly.xn1(n=34, N=113)
        >>> p, inv = StreamlinedNTRUPrime(23, 113, 5, seed=27182)._random_small_poly_invertible(12, None, modulus, inv_N=113)
        >>> p.in_ring(N=113, modulus=modulus) * inv
        1

        >>> modulus = Poly.xn1(n=34, N=113)
        >>> p, inv = StreamlinedNTRUPrime(23, 113, 5, seed=31415)._random_small_poly_invertible(12, 113, modulus, inv_N=113)
        >>> p.in_ring(N=113, modulus=modulus) * inv
        1
        """
        while True:
            poly = self.random_small_poly(weight, N, modulus, self.random)
            try:
                inverse = poly.in_ring(N=inv_N, modulus=poly.modulus.in_ring(N=inv_N)).inv()
                return poly, inverse
            except ValueError:
                pass
    
    def __init__(self, p, q, w, seed=None):
        # check parameters (non-exhaustive)
        assert q >= 17
        assert w <= p
        assert isprime(p)
        assert isprime(q)
        assert q > p
        assert self.co_prime(p, q)
        assert 2 * p >= 3 * w
        assert q >= 16 * w + 1
        # TODO assert that x^p - x - 1 is irreducible in (Z/q)[x]

        # initialize object state
        self.p, self.q, self.w = p, q, w
        self.random = RandomState() if seed is None else RandomState(seed)
        self.modulus_r = Poly.monomial(n=p) - Poly.monomial(n=1) - Poly.monomial(n=0)
        self.modulus_r3 = self.modulus_r.in_ring(N=3)
        self.modulus_rq = self.modulus_r.in_ring(N=q)

    def generate_keys(self):
        """
        >>> h, (f, g_r3_inv) = StreamlinedNTRUPrime(23, 113, 5, seed=3141).generate_keys()
        >>> f.coeffs.min(), f.coeffs.max()
        (-1, 1)
        >>> f.weight
        5
        >>> h.N, f.N, g_r3_inv.N
        (113, None, 3)
        """
        # Generate a uniform random small element g in R, repeat until invertible in R/3
        g, g_r3_inv = self._random_small_poly_invertible(self.w, None, self.modulus_r, 3)
        
        # Generate a uniform random f in Short
        f = self.random_small_poly(self.w, None, self.modulus_r, self.random)
        
        # Compute h = g/(3f) in R/q
        g_rq = g.in_ring(N=self.q, modulus=self.modulus_rq)
        three_f = (3*f).in_ring(N=self.q, modulus=self.modulus_rq)
        h = g_rq * three_f.inv()

        return h, (f, g_r3_inv)

    def encrypt(self, plain_text: Poly, public_key: Poly):
        # notation from NIST submission
        r = plain_text
        h = public_key

        # compute h*r in R/q
        hr = h.in_ring(N=self.q, modulus=self.modulus_rq) * r.in_ring(N=self.q, modulus=self.modulus_rq)

        # Round(h*r)
        return Poly(around(hr.coeffs * 3) // 3, N=self.q, modulus=self.modulus_rq)
    
    def decrypt(self, cipher_text: Poly, secret_key):
        # notation from NIST submission
        f, v = secret_key
        c = cipher_text

        # compute 3fc in R/q
        three_f_c = 3 * f.in_ring(N=self.q, modulus=self.modulus_rq) * c.in_ring(N=self.q, modulus=self.modulus_rq)

        # reduce modulo three to obtain e in R/3
        e = Poly(three_f_c.coeffs % 3, N=3, modulus=self.modulus_r3)

        # lift e*v in R/3 to a small polynomial in R
        return (e*v).in_ring(N=None, modulus=self.modulus_r)
