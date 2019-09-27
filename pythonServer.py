#############################
# Written by Bennett Larson #
#############################
from flask import Flask
from flask import jsonify
from flask import request
from Crypto.Random.random import randrange, getrandbits
from Crypto.Hash import SHA256
import collections
app = Flask(__name__)
from pasta import signUpUser, requestTk, verifyTk, CryptoError

global secrets, shares, vk, pp, x

@app.route('/login', methods=['GET'])
def login():
    SecParam = collections.namedtuple('SecParam', 'pwd rho')
    global shares, vk, pp, x
    username = request.args.get('username', type = str)
    password = request.args.get('password', type = str)
    valid = False
    try:
        C = int(SHA256.new(bytes(username, 'utf-8')).hexdigest(), 16)
        secrets = SecParam(pwd=password, rho=randrange(1, pp.opp.n))
        T = [1,2,3,4]
        tk = requestTk(C, secrets, shares, vk, pp, T, x)
        valid = verifyTk(vk, C, x, tk, pp)
    except CryptoError:
        valid = False
    return jsonify(valid)

@app.route('/register', methods=['GET'])
def register():
    global shares, vk, pp, x
    result = False
    try:
        username = request.args.get('username', type = str)
        password = request.args.get('password', type = str)
        x = getrandbits(128)
        _, _, shares, vk, pp = signUpUser(username, password, 512, 4, 5, x)
        result = True
    except Exception(e):
        result = False
    return jsonify(True)
