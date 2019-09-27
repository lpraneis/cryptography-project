############################
# Written by Logan Praneis #
############################
from math import factorial
from .znUtil import *
from .znUtil import DEBUG as DEBUG

def genshare(n, t, pp):
    """ Generate shares
    :param: n - number of parties
    :param: t - threshold
    :param: pp - public parameters
    :returns: (shares, vki) tuple, where shares, vk are individual lists
    """
    delta, G = pp
    # Generate Polynomaial as [d, d1, d2, ..., dt] where di in Zm

    # test polynomial for standard...
    if (DEBUG):
        poly = [17347086731773102978075690000237900743020218807929968231283538408081081615556425667988651206910109992090780923918455344212810165982683318676177484690666399305361183299557216807105056173518954507482609272584859403513512094325866621284363322891456654771554103467475819278865119664335946033292565600378034305395, 25730005957616554973438949889364004061220096693122195502203481336522068918498861353263727424064461781510724849217619899580306650491885582987382308183389625142906864311862562507207157836636312425246112360773257254354520538255110297459891342595486379890537920548920426896098758441569751548772269281829083087729, 4079927350961859416026048100337255946148027614268264709943066941571424451619927998264086050714984340064738363855499326546880241144214283312952964144788212883498050598172730085166192070120422126541186075610113875534270383412806800814540440493229239095239001319361145004089845387412094965163522250451471642978, 16973282165949417687412212895036635721665896339987157685828800266096559176796149636865549537878758430932757803404054196347938654216398913851877662663073596324536934499846570778417631609094012040206399955363935563371933002737818833275154697485266660638178116631613715787492321832706722219092736605284185868081]
    else:
        poly = [G.d] +  [G.random(G.m) for i in range(1, t)]

    # Evaluate the polynomial at n points (not including 0, as that would be d
    shares = [(i, int(evaluate(poly, i, G.m))) for i in range(1, n+1)]
    # Generate verification keys
    vki = [pow(G.v, s[1], G.N) for s in shares]
    vk = (G.N, G.e)
    # return shares, vki
    return shares, vk


def setup(k, n, t):
    ''' Setup Values
    :param:k - length of primes in bits {512, 1024, 2048}
    :param:n - number of authorities
    :param:t - threshold parameter

    :returns: (shares, vk, pp) tuple
    '''
    G = Zn(n, k)
    delta = factorial(n)
    pp = delta, G
    shares, vk = genshare(n, t, pp)
    return shares, vk, pp


    

def partEval(sk_i, vk, x, pp):
    '''
    Evaluate a single share 
    :param: sk_i - (i, sk) tuple for share
    :param: vk - verification key
    :param: x - bytestring message (password)
    :param: pp - public parameters
    :returns: (i, y_i ) tuple 
    '''
    delta, G = pp
    w = G.hash(x)
    i, ski = sk_i

    # Using 4 * delta * ski now to find in Qn
    # Also, this saves from scaling lambda coefficient later
    signVal = 4 * delta * ski
    return (i, pow(w, signVal,  G.N))
    

def combine(sd, pp, t, x):
    '''
    Combine all shares
    :param: sd - share list of tuples (i, share)
    :param: pp - public parameters
    :param: t - threshold
    :param: x - password 
    :returns: 
        False - token invalid, not enough parties
        tk - token

    '''
    if len(sd) < t:
        # not threshold
        return False
    elif len(set(sd)) != len(sd):
        # Not unique
        return False


    delta, G = pp
    inds, vals = zip(*sd)
    inds, vals = list(inds), list(vals)

    w = 1
    for i in range(len(inds)):
        l = lamb_coeff(inds[i], inds, delta)
        if l < 0:
            w *=  pow(inverse(vals[i], G.N), -1 * l,   G.N)
        else: 
            w *= pow(vals[i], l, G.N)
    w = w % G.N
    
    eprime = 4 * (delta ** 2)
    a, b = egcd(eprime,  G.e) 
    # assert( (a * 4 * delta**2 + b* G.e ) % G.N == 1)

    hashed = G.hash(x)

    # Equality on bottom on p.215  
    lhs = (pow(w, G.e, G.N))
    rhs = (pow(hashed, eprime, G.N))
    if (DEBUG):
         print('Equal?', lhs == rhs)

    if a < 0:
        tk = (pow(inverse(w, G.N), -1 *a, G.N) * pow(hashed, b, G.N))
    else:
        tk = (pow(w, a, G.N) * pow(inverse(hashed, G.N),  -1 * b, G.N))
            
    
    return tk % G.N




def verify(vk, x, tk, pp):
    ''' 
    :param: vk - verification key
    :param: x - password
    :param: tk - token
    :param: pp - public parameters
    :returns:
        True - tk is valid for x
        False - tk is not valid for x
    '''
    delta, G = pp
    lhs = pow(tk, G.e, G.N) 
    rhs = G.hash(x) % G.N
    if (DEBUG):
        print('lhs', lhs)
        print('rhs', rhs)
    return (lhs == rhs)

def test():

    # testing utility, should output true!
    t = 4
    shares, vk, pp = setup(512, 5, t)
    parts = [partEval(shares[i], vk, b'testing', pp) for i in range(len(shares))]
    thresh_parts = [parts[3], parts[2], parts[4], parts[0]] # Share 1,2,4,0
    tk = combine(thresh_parts, pp, t, b'testing')
    v = verify(vk, b'testing', tk, pp)
    print('Verified?', v)

