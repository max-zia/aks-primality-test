"""
Failed attempt at implementing aks primality test.
"""

from math import log, sqrt, floor
from numpy.polynomial import polynomial as P
from numpy import poly1d

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

def aks_test(n):
    # Check if n is a perfect power. If so, return composite.
    if is_perfect_power(n):
        return "composite"
    
    # Find the smallest r such that the multiplicative order of n modulo r
    # is greater than log(n, 2)^2. If r and n are not coprime, skip this r.
    r = 2
    while True:
        if gcd(r, n) != 1:
            r += 1
        elif ord(n, r) > pow(log(n, 2), 2):
            break
        else:
            r += 1

    # If 1 < gcd(a, n) < n for some a <= r, return composite
    for a in range(1, r):
        if gcd(a, n) > 1 and gcd(a, n) < n:
            return "composite"

    # If n <= r, return prime
    if n <= r:
        return "prime"

    # If (x + a)^n mod (x^r - 1, n) != (x^n + a) mod (x^r - 1, n), output congruent
    for a in range(1, floor(sqrt(phi(r)) * log(n, 2))):
        print(f"\n\nPASS {a}")

        # Expand the polynomial (x + a)^n
        polynomial = poly1d([1, a])
        polynomial = pow(polynomial, n)

        # Compute (x + a)^n mod (x^r - 1, n)
        modulus = get_poly_array(r, -1)
        x1 = (polynomial/modulus)[1]
        x1 = reduce_by_modulo(x1, n)

        print(x1)

        # Compute (x^n + a) mod (x^r - 1, n)
        polynomial = get_poly_array(n, a)
        x2 = (polynomial/modulus)[1]
        x2 = reduce_by_modulo(x2, n)

        print(x2)

        if x1 != x2:
            return "composite"
    
    return "prime"


def ord(a, n):
    """
    Computes the multiplicative order of a modulo n, namely the smallest
    number k such that a^k is congruent with 1 (mod n). The multiplicative
    order only exists when a and n are coprime. 
    """
    k = 2

    while True:
        if (pow(a, k) % n) == 1:
            break
        else:
            k += 1
    
    return k


def is_perfect_power(n):
    """
    Returns True if n is a perfect power, that is, if there exist integers
    a and b such that n = a^b. This algorithm was first proposed in
    https://link.springer.com/article/10.1007/BF01228507
    """
    # If n is a perfect bth power, then b <= log(n). Therefore, compute an
    # integer approximation x of n^(1/b) for every such b, starting with b = 2.
    # Only prime values of b need to be used.
    b = 2
    while b <= (log(n, 2) + 1):
        if b in primes:
            # Establish upper limit for search of bases x to prime b
            lim = pow(2, (log(n, 2) / b) + 1)

            # Do a search between 2 and lim
            # BINARY SEARCH would be an optimisation
            x = 2
            while x <= int(lim):
                if (x**b) == n:
                    return True
                else:
                    x += 1
            b += 1
        else:
            b += 1
    
    return False


def gcd(a, b):
    """
    Euclidean algorithm for computing the greatest common divisor of two integers.
    """
    while b > 0:
        r = a % b
        a = b
        b = r

    return a


def phi(n):
    """
    Counts the positive integers up to given integer n that are coprime with n.
    Also known as Euler's totient (or phi) function.
    """
    return len([x for x in range(1, n) if gcd(x, n) == 1])


def get_poly_array(exponent, constant):
    """
    To operate on polynomials, numpy's polynomial module needs array-like
    objects (e.g., [1, 2, 3] would correspond to x^2 + 2x + 3). Therefore,
    this function returns an array [1, 0, ..., 0, constant] corresponding to the
    polynomial (x^exponent + constant). This is needed to compute the modular
    division of polynomials that features in Step 5 of the AKS test. 
    """
    array = [1]
    for i in range(exponent - 1):
        array.append(0)
    array.append(constant)

    return poly1d(array)


def reduce_by_modulo(array, modulus):
    """
    Takes a polynomial array (e.g., [1, 2, 3] would correspond to x^2 + 2x + 3)
    and returns another polynomial array in which each of the coefficients in
    the original have been subjected to modulo division by modulus. 
    """ 
    return poly1d([(x % modulus) for x in array])


n = 31
r = 29
a = 3

# Expand the polynomial (x + a)^n
p = poly1d([1, a])
print(p)
p = p**n

# Get (x^r - 1)
m = get_poly_array(r, -1)

# Compute (x + a)^n mod (x^r - 1)
p = (p/m)[1]

# Reduce by mod n
p = reduce_by_modulo(p, n)


