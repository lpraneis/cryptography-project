############################
# Written by Logan Praneis #
############################
from Crypto.Util.number import getPrime, inverse, isPrime, GCD
from Crypto.Hash import BLAKE2b, SHAKE256
from binascii import hexlify
from Crypto.Random.random import randrange


DEBUG = False

class Zn:
    def __init__(self, n, bits=512):
        """
        :param: n - threshold
        :param: `optional` bits - # bits in prime, default 512
        """
        self.bits = bits
        if self.bits == 512:
            self.static_512_key_parameters()
        elif self.bits == 1024:
            self.static_1024_key_parameters()
        elif self.bits == 2048:
            self.static_2048_key_parameters()
        else:
            print('{} not accepted, using 512 bit keys'.format(bits))
            self.static_512_key_parameters()

        # self.e = getPrime(512)
        self.N = self.p * self.q
        self.totient = (self.p -1) * (self.q -1 ) # not needed ?

        if (DEBUG):
            self.v= 30051438525904603932038634613000506355731333729926882770827749382442527506143269391357699435109939444343116292747632257699116276797824029401181860856620430423376587332544138525001568919623373685708978958507709568035016483077413186364672948764771908340105872500533783665268325774837963538579805171956190567494
            self.e = 53
        else:
            self.v = pow(self.random(self.N), 2, self.N) # choose random v in Qn
            self.e = getPrime(2 * n.bit_length())

        # de = 1 mod m
        self.d = inverse(self.e, self.m)

        # Assertions
        # assert((self.d * self.e) % self.m == 1)


    def static_512_key_parameters(self):
        
        p =14803816792137194961464740500861084314498273884502062848972100269359945855258103338378053276309214964458559920121618145541357486576747004424971736904458523
        pp = 7401908396068597480732370250430542157249136942251031424486050134679972927629051669189026638154607482229279960060809072770678743288373502212485868452229261
        q =13801179485494374276973210326755346111973952066482811157739430718935024900211415657056158753382064093168952785842599167158882530769356959883177037592830167 
        qp = 6900589742747187138486605163377673055986976033241405578869715359467512450105707828528079376691032046584476392921299583579441265384678479941588518796415083

        #Assertions
        # assert(2*pp +1 == p)
        # assert(2*qp +1 == q)
        # assert(isPrime(p))
        # assert(isPrime(pp))
        # assert(isPrime(q))
        # assert(isPrime(qp))


        m = qp * pp
        self.p, self.q, self.m = p, q,  m

    def static_1024_key_parameters(self):
        
        p =246534071970375761184459786036510271602126530601097845003447340565610927982656140909521434213025854892807632901834242304162377103942910006414608086604746070849623005950199278734033772557908558473447995487898954880321914920982943884224020529574357687644637824533207433752051195886250335306153316258231085724679
        pp =123267035985187880592229893018255135801063265300548922501723670282805463991328070454760717106512927446403816450917121152081188551971455003207304043302373035424811502975099639367016886278954279236723997743949477440160957460491471942112010264787178843822318912266603716876025597943125167653076658129115542862339
        q = 237193845588964867838364573394095192553085949412310151962074555573946706610754100537895268206835971266985568809190215487397061161349330302956027041152068154756620520765841644599875281035967259216754074318524466124917071268637733103686000319889308168832600272588994356381064734163808418631163495322774840449899
        qp = 118596922794482433919182286697047596276542974706155075981037277786973353305377050268947634103417985633492784404595107743698530580674665151478013520576034077378310260382920822299937640517983629608377037159262233062458535634318866551843000159944654084416300136294497178190532367081904209315581747661387420224949
        m = qp * pp
        self.p, self.q, self.m = p, q,  m

    def static_2048_key_parameters(self):

        qp = 16519407997707361939917278697292696807145958283342873227055639664808929544449075946731097026777747814338050234368865619420077223796876172267179868234961469872193993413646149696202117888992614387011205348292009713138510643085946000304647985461593515788877790828839468033150495281257467551315858749285485655689570516682601242804566215506106014228740190462622374199399493053814275141339015060712169985692176523505269258927780362254850166113643449875803648913127383151568057512383727881000059264971486794322583022429959296323800845078028893137191759277261981883867946026586310316121476993251483310083056324436500493378121 
        q = 33038815995414723879834557394585393614291916566685746454111279329617859088898151893462194053555495628676100468737731238840154447593752344534359736469922939744387986827292299392404235777985228774022410696584019426277021286171892000609295970923187031577755581657678936066300990562514935102631717498570971311379141033365202485609132431012212028457480380925244748398798986107628550282678030121424339971384353047010538517855560724509700332227286899751607297826254766303136115024767455762000118529942973588645166044859918592647601690156057786274383518554523963767735892053172620632242953986502966620166112648873000986756243 
        p = 47861407752790525997131098361108406122043048720868999328281323903371925202315991865537967768731143464033692447855088627006073694264317365582380595334207133298537924589942253341827745672450433376941500326295068787878081775716133283627670531493685135607319289868027506571347598049080446281591129864228995139148236090883284601639122653328795607705391929418876441742949286444426345072361176594702993097458512980608512112173855855323247693697686504690215496470594493009783878291699542825026419274715235678738119520853777035524370581369370802600982137629966848151596047919239187925248438366319331990827805128321189037230763 
        pp = 23930703876395262998565549180554203061021524360434499664140661951685962601157995932768983884365571732016846223927544313503036847132158682791190297667103566649268962294971126670913872836225216688470750163147534393939040887858066641813835265746842567803659644934013753285673799024540223140795564932114497569574118045441642300819561326664397803852695964709438220871474643222213172536180588297351496548729256490304256056086927927661623846848843252345107748235297246504891939145849771412513209637357617839369059760426888517762185290684685401300491068814983424075798023959619593962624219183159665995413902564160594518615381
        m = qp * pp
        self.p, self.q, self.m = p, q,  m

    @staticmethod
    def random(x):
        """
        :param: x - integer
        :returns: Random (secure) integer in range (0, x-1)
        """
        return randrange(0, x)

    def hash(self, x):
        """
        :param: x (bytestring) : value to hash
        :returns: H(x): x hashed into Zn*
        """

        b = SHAKE256.new()
        b.update(x)

        # reading the bitlength of parameter +1 from hash
        length = int.bit_length(self.N) + 1
        h = int(hexlify(b.read(length)), 16) % self.N
        # calculaitng h in Zn
        return h

    def hashSign(self , m, vi,  l1, delta, si):
        # Not used currently, but could be used to verify token parts are valid
        """
        H'(m)
        :param: m  - message (plaintext)
        :param: vi  - verification key share
        :param: l1 - Secondary security parameter
        :param: delta  - n factorial 
        :param: si - secret share
        :returns: (z,c) - pair
        """
        x = self.hash(m)
        r = randrange(0, pow(2, int.bit_length(self.N) + 2*l1))
        vprime = pow(self.v, r, self.N)
        xtilde = pow(x, 4*delta, self.N)
        xprime = pow(xtilde, r, self.N)
        xi = pow(x, 4 * delta * si, self.N)

        b = SHAKE256.new()
        # Via Shoup p.215, 
        b.update(bytes(self.v % self.N))
        b.update(bytes(xtilde))
        b.update(bytes(vi % self.N))
        b.update(bytes(xi % self.N))
        b.update(bytes(vprime))
        b.update(bytes(xprime))
        c =int(hexlify(b.read(l1)), 16) % self.N
        z = c * si + r

        return (z, c)




def evaluate(poly, x, p=None):
    """ Evaluate a polynomial in Zp*
    :param: poly - list of polynomail coefficients
    :param: x - value to evaluate polynomail at
    :param: p `optional` - prime
    :returns: p(x), integer in Zp*
    """
    evaluated = ((poly[j] * pow(x, j)) for j in range(len(poly)))
    if p:
        return sum(evaluated) % p
    else:
        return sum(evaluated) 



def lamb_coeff(i,  S, delta):
    """ Compute lagrage interpolation, as shown in Shoup
    :param: i - point in S
    :param: S - set of k points in {0, 1, ..., n}
    :param: delta - n!

    :returns: lagrange interpolation at 0, scaled
    """
    v = delta
    # Construct Numerator
    for e in S:
        if e  != i:
            v*= e


    # Construct Denominator
    for e in S:
        if e != i:
            v //= (e - i)

    return v


def egcd(a, b):
    """ Compute Extened GCD
    :param: a - int
    :param: b -int
    :returns:
        (x, y) - such that ax + by = 1
    """
    x0, x1, y0, y1 = 0, 1, 1, 0

    while a != 0:
        q, b, a = b // a, a, b % a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return x0, y0