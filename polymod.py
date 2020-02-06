from math import ceil

from numpy import array, zeros, nonzero, trim_zeros, flatnonzero


def Zn_inv(a, n):
    """
    Computes the inverse of a in ℤ_n using the Extended Euclidean algorithm.

    >>> Zn_inv(2, 7)
    4

    >>> Zn_inv(2, 3)
    2

    >>> Zn_inv(2, 256)
    Traceback (most recent call last):
      ...
    ValueError: 2 has no inverse in ℤ_256.
    """
    if a < 0:
        a += n

    t = 0
    r = n
    new_t = 1
    new_r = a

    while new_r != 0:
        quotient = r // new_r
        t, new_t = new_t, t - quotient * new_t
        r, new_r = new_r, r - quotient * new_r

    if r > 1:
        raise ValueError(f'{a} has no inverse in ℤ_{n}.')
    if t < 0:
        t += n

    assert (t * a) % n == 1  # sanity check

    return t


class Poly:
    """
    Represents polynomials from ℤ[X], ℤ[X]/p, ℤ_N[X], and ℤ[X]_N/p for any p in ℤ[X] and N in ℕ.

    >>> Poly([0, 1], N=23)
    X

    >>> p = Poly([0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1], N=2, modulus=Poly([1, 1, 0, 1, 1, 0, 0, 0, 1], N=2))
    >>> p % p.modulus
    1
    """

    @classmethod
    def monomial(cls, n, *, N=None, modulus=None):
        """
        >>> Poly.monomial(16)
        X¹⁶
        """
        return cls([0] * n + [1], N=N, modulus=modulus)
    
    @classmethod
    def xn1(cls, *, n, N=None):
        """
        >>> Poly.xn1(n=25)
        X²⁵ - 1

        >>> 4 * Poly.xn1(n=1000,N=3)
        X¹⁰⁰⁰ - 1
        """
        return cls.monomial(n, N=N) - cls.monomial(n=0, N=N)

    def __init__(self, coeffs=None, *, N=None, modulus=None):
        if coeffs is None:
            coeffs = []
        coeffs = array(coeffs)
        if N:
            coeffs = coeffs[:] % N
            coeffs[coeffs > ceil((N-1)/2)] -= N
        self.coeffs = trim_zeros(coeffs, trim='b')
        self.N = N
        self.modulus = modulus

    @property
    def deg(self):
        """
        Degree of the polynomial.

        >>> Poly.monomial(n=635).deg
        635

        >>> (Poly([1, 1], N=2) * Poly.xn1(n=12, N=2)).deg
        13
        """
        try:
            return nonzero(self.coeffs)[0].max()
        except ValueError:
            return 0  # TODO

    @property
    def weight(self):
        """
        Number of non-zero coefficients.

        >>> Poly([]).weight
        0

        >>> Poly.xn1(n=1337).weight
        2
        """
        return sum(self.coeffs != 0)

    @property
    def lead_c(self):
        """
        >>> Poly().lead_c
        0

        >>> Poly([1, 2, 3, 4]).lead_c
        4

        >>> Poly(list(range(100)), N=2).lead_c
        1
        """
        if len(self.coeffs) == 0:
            return 0
        return self.coeffs[self.deg]

    def __add__(self, other):
        """
        >>> Poly([0, 1, 23], N=23) + Poly([0, 1, 2], N=23)
        2X² + 2X
        >>> Poly([-42, 1, 0, 0, 0, 1, -1, -1]) + Poly([0] * 1234 + [4321])
        4321X¹²³⁴ - X⁷ - X⁶ + X⁵ + X - 42
        """
        self.assert_compatibility_with(other)
        self_c, other_c = self.align_coefficients(self, other)
        return Poly(self_c + other_c, N=self.N, modulus=self.modulus)

    def __sub__(self, other):
        """
        >>> Poly([0, 1]) - Poly([0, 1])
        0
        >>> Poly([0, 1]) - Poly([0, 1, -1, -1])
        X³ + X²
        """
        return self.__add__(-other)

    def __neg__(self):
        return -1 * self

    def __mul__(self, other):
        """
        >>> 3 * Poly([0, 1, 3])
        9X² + 3X

        >>> Poly([0, 1, 3]) * Poly([4, 1])
        3X³ + 13X² + 4X

        >>> Poly([0, 1, 3]) * Poly([4, 1, 2])
        6X⁴ + 5X³ + 13X² + 4X

        >>> Poly([1, 2, 3], N=23) * Poly([5, 7, 11], N=23)
        10X⁴ - 3X³ - 6X² - 6X + 5

        >>> Poly([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1], N=11)
        -X¹⁰ + X⁹ + X⁶ - X⁴ + X² + X - 1

        >>> Poly([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1], N=11) * \
            Poly([0, -1], N=11)
        X¹¹ - X¹⁰ - X⁷ + X⁵ - X³ - X² + X

        >>> Poly([1], N=2) * Poly([2], N=3)
        Traceback (most recent call last):
          ...
        ValueError: Arithmetic for operands undefined. For arithmetic operations, Poly class operands must have equal N and modulus, but given were
        N: 2 and 3
        modulus: None and None.
        """
        # scalar multiplication
        if isinstance(other, (int, float)):
            return Poly(other * self.coeffs, N=self.N, modulus=self.modulus)

        # polynomial multiplication
        self.assert_compatibility_with(other)
        self_c, other_c = self.align_coefficients(self, other)
        length = len(self_c)
        mul_c = zeros((2*length))
        # Sparse Multiplication
        for i in flatnonzero(self_c):
            for j in flatnonzero(other_c):
                mul_c[i + j] += self_c[i] * other_c[j]
        # Operand Scanning
        # for i in range(length):
        #     for j in range(length):
        #         mul_c[i + j] += (self_c[i] * other_c[j])
        if self.N:
            mul_c = mul_c % self.N
        p = Poly(mul_c, N=self.N, modulus=self.modulus)
        return p % self.modulus if self.modulus else p

    __rmul__ = __mul__

    def euclidean_division(self, other):
        if self.N is None:
            raise ValueError('Euclidean division not supported for ℤ[X]')
        # computations in ℤ_N[X], hence modulus=None
        a = Poly(self.coeffs, N=self.N)
        b = Poly(other.coeffs, N=other.N)
        quotient = Poly([], N=self.N)
        remainder = a
        while remainder.deg >= b.deg and remainder != Poly([]):
            c = (remainder.lead_c * Zn_inv(b.lead_c, a.N)) % a.N
            if c == 0: raise ValueError  # sanity check
            divide_by = c * self.monomial(remainder.deg - b.deg, N=a.N)
            quotient += divide_by
            remainder -= b * divide_by
        if not b * quotient + remainder == a: raise ValueError  # sanity check
        return (
            quotient.in_ring(N=self.N, modulus=self.modulus),
            remainder.in_ring(N=self.N, modulus=self.modulus),
        )

    def inv(self):
        """
        >>> Poly([1, 1, 0, 0, 1, 0, 1], N=2, modulus=Poly([1, 1, 0, 1, 1, 0, 0, 0, 1], N=2)).inv()
        X⁷ + X⁶ + X³ + X

        >>> Poly([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1], N=3, modulus=Poly.xn1(n=11, N=3)).inv()
        -X⁹ + X⁸ - X⁷ + X⁵ - X⁴ - X³ - X + 1

        Note that inv() sometimes cannot find an inverse even if one exists, as in this case:
        >>> Poly([-1, 1, 1, 0, -1, 0, 1, 0, 0, 1, -1], N=32, modulus=Poly.xn1(n=11, N=32)).inv()
        Traceback (most recent call last):
          ...
        ValueError: Extended Euclidean algorithm could not find an inverse for -X¹⁰ + X⁹ + X⁶ - X⁴ + X² + X - 1 in ℤ_32[X]/(X¹¹ - 1).

        >>> Poly([1, 1], modulus=Poly.xn1(n=3)).inv()
        Traceback (most recent call last):
          ...
        ValueError: Extended Euclidean algorithm could not find an inverse for X + 1 in ℤ[X]/(X³ - 1).

        >>> Poly([0, 1, 1], N=3, modulus=Poly([-1, -1, 0, 1], N=3)).inv()
        X - 1
        """
        t = p0 = Poly([], N=self.N)
        r = self.modulus
        new_t = Poly([1], N=self.N)
        new_r = self.in_ring(N=self.N, modulus=None)

        try:
            while new_r != p0:
                q = r // new_r
                t, new_t = new_t, t - q * new_t
                r, new_r = new_r, r - q * new_r

        except ValueError:
            raise ValueError(f'Extended Euclidean algorithm could not find an inverse for {self} in '
                             f'{self.ring_name()}.')

        inverse = t / r  # TODO This fails sometimes with remainder != 0. Bug?
        inverse = inverse.in_ring(N=self.N, modulus=self.modulus)

        assert inverse * self % self.modulus == Poly([1], N=self.N, modulus=self.modulus), \
            f'{self * t} was not 1!\nself: {self}\nt: {t}'
        return inverse

    def __mod__(self, other):
        """
        >>> Poly([3, 1, 1], N=6)
        X² + X + 3
        >>> Poly([-1, 1], N=6)
        X - 1
        >>> Poly([3, 1, 1], N=6) % Poly([-1, 1], N=6)
        -1
        """
        _, remainder = self.euclidean_division(other)
        return remainder

    def __truediv__(self, other):
        """
        >>> Poly([3, 1, 1], N=23) / Poly([3, 1, 1], N=23)
        1
        >>> Poly([6, 2, 2], N=23) / Poly([3, 1, 1], N=23)
        2
        >>> Poly([0, 1, 1], N=3) / Poly([2], N=3)
        -X² - X
        """
        result, remainder = self.euclidean_division(other)
        if remainder != Poly([], N=self.N):
            raise ValueError(
                f'Could not divide {self} by {other}, as the remainder was not zero, but {remainder}. '
                f'I.e., we have ({result}) * ({other}) + ({remainder}) = {self}.'
            )
        return result

    def __floordiv__(self, other):
        """
        >>> Poly([3, 1, 1], N=23) // Poly([3, 1, 1], N=23)
        1
        >>> Poly([6, 2, 2], N=23) // Poly([2, 1, 1], N=23)
        2
        >>> Poly([0, 1, 2], N=23) // Poly([0, 1, 1], N=23)
        2
        """
        result, remainder = self.euclidean_division(other)
        return result

    def __eq__(self, other):
        """
        >>> Poly([1, 1], N=23) == Poly([1, 1], N=23)
        True
        >>> Poly([1, 1], N=23) == Poly([1, 1], N=3)
        True
        >>> Poly([1, 1], N=23) == Poly([1, 1, 0], N=23)
        True
        >>> Poly([1, 1], N=23) == Poly([1], N=23)
        False
        >>> Poly() == 0
        False
        """
        if not isinstance(other, self.__class__):
            return False
        self_c, other_c = self.align_coefficients(self, other)
        return all(self_c == other_c)

    def __str__(self):
        """
        Returns self string that represents this polynomial.
        >>> Poly([0, 1])
        X
        >>> Poly([0, 1, 42, 0, -1], N=23)
        -X⁴ - 4X² + X
        >>> Poly([-42, 1.0, 0, 0, 0, 1, -1, -1], N=23)
        -X⁷ - X⁶ + X⁵ + X + 4
        >>> Poly([0] * 1234 + [4321])
        4321X¹²³⁴
        >>> Poly([], N=23)
        0
        >>> Poly([-1, 1])
        X - 1
        >>> Poly([.3, .5])
        0.500000X + 0.300000
        >>> Poly([1, -1])
        -X + 1
        """
        def superscript(digit):
            offset = {1: 184, 2: 176, 3: 176}  # What a mess.
            return chr(offset.get(digit, 0x2070) + digit)

        def monomial_repr(idx, var='X'):
            if idx == 0:
                return ''
            if idx == 1:
                return var
            return var + ''.join([superscript(int(digit)) for digit in str(idx)])

        def sgn(c):
            return '+' if c >= 0 else '-'

        def coeff_repr(c, force=False):
            if c == 1 and not force:
                return '+ '
            if c == -1 and not force:
                return '- '
            if c % 1 == 0:
                return f'{sgn(c)} {int(abs(c))}'
            return f'{sgn(c)} {abs(c):f}'

        r = ' '.join(reversed([
            (coeff_repr(c, force=idx == 0) + monomial_repr(idx)).strip()
            for idx, c in enumerate(self.coeffs[:self.deg+1])
            if c != 0
        ])).strip().lstrip('+ ')

        if r == '':
            return '0'
        if r.startswith('- '):
            return '-' + r[2:]
        return r

    __repr__ = __str__

    def in_ring(self, *, N=None, modulus=None):
        """
        Converts the polynomial into a polynomial with same coefficients in a different ring.
        """
        return Poly(self.coeffs, N=N, modulus=modulus)

    def assert_compatibility_with(self, other):
        """
        If given `other` is arithmetically incompatible with this polynomial, a ValueError is raised.
        """
        if self.N != other.N or self.modulus != other.modulus:
            raise ValueError('Arithmetic for operands undefined. For arithmetic operations, Poly class operands must '
                             'have equal N and modulus, but given were\n'
                             f'N: {self.N} and {other.N}\n'
                             f'modulus: {self.modulus} and {other.modulus}.')

    @staticmethod
    def align_coefficients(p1, p2):
        """
        Brings the coefficient arrays of p1 and p2 into same length, padding with zero where needed.
        """
        length = max(len(p1.coeffs), len(p2.coeffs))
        p1_c = zeros((length, ))
        p1_c[:len(p1.coeffs)] = p1.coeffs
        p2_c = zeros((length, ))
        p2_c[:len(p2.coeffs)] = p2.coeffs
        return p1_c, p2_c

    def ring_name(self):
        """
        Returns a human-readable name for the ring this polynomial is in.
        """
        base = f'ℤ_{self.N}[X]' if self.N is not None else 'ℤ[X]'
        modulus = f'/({self.modulus})' if self.modulus else ''
        return base + modulus
