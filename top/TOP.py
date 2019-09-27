###############################
# Written by Paul Rheinberger #
###############################

##########################################################################################
#
# Imports
#
##########################################################################################
import collections
import hashlib
import random
import binascii
import sys



##########################################################################################
#
# Elliptic Curve
#
##########################################################################################
EllipticCurve = collections.namedtuple('EllipticCurve', 'name p a b g n h')

curve = EllipticCurve(
    'secp256k1',
    # Field characteristic.
    p=115792089210356248762697446949407573530086143415290314195533631308867097853951,
    # Curve coefficients.
    a=-3,
    b=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
    # Base point.
    g=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
       0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5),
    # Subgroup order.
    n=115792089210356248762697446949407573529996955224135760342422259061068512044369,
    # Subgroup cofactor.
    h=1,
)



##########################################################################################
#
# Helper Functions
#
##########################################################################################
def inverse_mod(k, p):
    """Returns the inverse of k modulo p.
    This function returns the only integer x such that (x * k) % p == 1.
    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_mod(-k, p)

    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p



def is_on_curve(point):
    """Returns True if the given point lies on the elliptic curve."""
    if point is None:
        # None represents the point at infinity.
        return True

    x, y = point

    return (y * y - x * x * x - curve.a * x - curve.b) % curve.p == 0



def point_add(point1, point2):
    """Returns the result of point1 + point2 according to the group law."""
    assert is_on_curve(point1)
    assert is_on_curve(point2)

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2 * y1, curve.p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_mod(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p,
              -y3 % curve.p)

    assert is_on_curve(result)

    return result



def scalar_mult(k, point):
    """Returns k * point computed using the double and point_add algorithm."""
    assert is_on_curve(point)

    if k % curve.n == 0 or point is None:
        return None

    if k < 0:
        # k * point = -k * (-point)
        return scalar_mult(-k, point_neg(point))

    result = None
    addend = point

    while k:
        if k & 1:
            # Add.
            result = point_add(result, addend)

        # Double.
        addend = point_add(addend, addend)

        k >>= 1

    assert is_on_curve(result)

    return result



def exp1(e,g,n):
    """ This calculates g^e (mod n) for integers e,g, and n """
    t = 1
    sq = g
    e1 = e
    while(e1!=0):
        if (e1%2)==1:
            t = (sq*t)%n
            e1 = (e1-1)//2
        else:
            e1 = e1//2
        sq = (sq*sq)%n
    return(t)



def mult2(a,b,q,n):
    """ This multiplies two polynomials in GF(p^2)
    that is it calculates c(x) = a(x)*b(x) <mod q(x)>
    where a(x) = a[0]x + a[1] and b(x) = b[0]x + b[1] and
    q(x) = x^2 + q[0]x + q[1]
    """
    t1 = (a[0]*b[1])%n
    t2 = (a[1]*b[0])%n
    t1 = (t1+t2)%n
    t2 = (a[1]*b[1])%n
    t3 = ((n-1)*q[0])%n
    t4 = ((n-1)*q[1])%n
    t5 = (a[0]*b[0])%n
    t3 = (t5*t3)%n
    t4 = (t5*t4)%n
    c = [(t1+t3)%n , (t2+t4)%n]
    return(c)



def exp2(e,g,q,n):
    """ exp2 exponentiates in GF(p^2) by calculating (g(x))^e mod <q(x)> """
    t = [0,1]
    sq = g
    e1 = e
    while(e1!=0):
        if (e1%2)==1:
            t = mult2(sq,t,q,n)
            e1 = (e1-1)//2
        else:
            e1 = e1//2
        sq = mult2(sq,sq,q,n)
    return(t)
            


def CL(c,b,p):
    """ This is the main function used in the Cipolla-Lehmer algorithm for
    calculating square roots mod a prime p """
    t1 = (b*b)%p
    t2 = (4*c)%p
    t2 = (p-t2)%p
    g = (t1+t2)%p
    e = (p-1)//2
    h = exp1(e,g,p)
    s = 1
    if ((h==0) or (h==1)):
        s = 0
    e = e+1
    t1 = ((p-1)*b)%p
    t2 = c%p
    q = [t1,t2]
    a = [1,0]
    b = exp2(e,a,q,p)
    t = s*b[1]
    return(t)



def sqrt1(c,p):
    """ This is the Cipolla-Lehmer algorithm. It calculates a square root
    of c mod p assuming c is a quadratic residue mod p where p is a prime > 2"""
    m1 = 50
    t = 0
    c = c%p
    for i in range(m1):
        y = CL(c,((i+1)%p),p)
        t1 = (y*y)%p
        if (t1==c):
            t = y
            break
    return(t)



def getPoint(x):
    """ Returns a point on the elliptic curve with x as the x coordinate """
    while True:
        y2 = (x * x * x + curve.a * x + curve.b) % curve.p
        y = sqrt1(y2, curve.p)
        if pow(y, 2, curve.p) == y2:
            return (x, y)
        x = x + 1 % curve.n



def setup():
    return curve

##########################################################################################
#
# Test
#
##########################################################################################
def test():
    pwd_input = "Password123!"
    pwd_hash = int(hashlib.sha256(pwd_input.encode('utf-8')).hexdigest(), 16)


    # Sign Up
    pwd_point = getPoint(pwd_hash)

    k = random.randrange(1, curve.n)

    h = scalar_mult(k, pwd_point)

    # Request 
    rho = random.randrange(1, curve.n)
    c = scalar_mult(rho, pwd_point)

    # Respond
    z = scalar_mult(k, c)

    # Finalize
    rho_inverse = inverse_mod(rho, curve.n)
    h_prime = scalar_mult(rho_inverse, z)

    print("h       ", h[0])
    print("c       ", c[0])
    print("z       ", z[0])
    print("h_prime ", h_prime[0])
