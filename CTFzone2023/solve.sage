from pwn import *
from sage.all import *
import random
from hashlib import blake2b
from multiprocessing import Pool
import time

q = 21888242871839275222246405745257275088696311157297823662689037894645226208583
beta = 2203960485148121921418603742825762020974279258880205651966
nprimes = 10000000
#nprimes = 10000
pr = Primes()[:nprimes]
def trivial_fact(n):
    pairs = []
    for elem in pr:
        if n % elem == 0:
            count = 1
            n //= elem
            while n % elem == 0:
                n //= elem
                count += 1
            pairs.append((elem, count))

    return pairs


def recov_lambda(E, p, P):
    Q = E(P.xy()[0]*beta, P.xy()[1])
    lam = discrete_log(Q, P, p, operation='+')
    return lam

def query(P):
    payl = f'ecdh ({P.xy()[0]},{P.xy()[1]})'
    io.send(payl)
    io.recvuntil('session check: ')
    hsh = io.recvline(False).decode()
    return hsh

def hashcoll(hsh, beg, end, P, E):

    id0 = E(0)
    if beg == 0:
        orig = id0
    else:
        orig = ZZ(beg) * P
    for guess in range(beg, end):
        if orig == id0:
            session_key = tuple([0, 0])
        else:
            session_key = tuple(orig.xy())

        hash_check = blake2b((str(session_key) + "0").encode()).hexdigest()
        if hash_check == hsh:
            return (k1 + lam*k2 - guess, p)

        orig += P

    return None

def endo(P, E):
    return E(P.xy()[0]*beta, P.xy()[1])

def nonsquarefree(hsh, p, P, E):
    GG = E.gens()
    for G in GG:
        G_ord = G.order()
        if G_ord % p == 0:
            P = (G_ord//p) * G
            break

    hsh = query(P)
    eP = endo(P, E)
    id0 = E(0)

    for i in range(p):
        if i == 0:
            P1 = id0
        else:
            P1 = ZZ(i) * P
        for j in range(p):
            if j == 0:
                P2 = id0
            else:
                P2 = ZZ(j) * eP
            orig = P1 + P2
            if orig == id0:
                session_key = tuple([0, 0])
            else:
                session_key = tuple(orig.xy())
            hash_check = blake2b((str(session_key) + "0").encode()).hexdigest()
            if hash_check == hsh:
                return [(k1 - i, p), (k2 - j, p)]
    return None


Fq = GF(q)
b = 715007269373311787545362300121277285376352080834761325372306735370567249410
P.<k1, k2> = PolynomialRing(ZZ)
nthreads = 8

dic = dict()
for ntrial in range(500):
    b = random.randint(1, q)
    E = EllipticCurve(Fq, [0, b])
    order = E.order()
    if order not in dic.keys():
        dic[order] = b

# io = remote('localhost',int(30449))
io = remote('wise-sage.ctfz.one', int(1337))
io.recvuntil('help - this info\n>')

remainders = []
moduli = []
    
for order in dic.keys():
    b = dic[order]
    E = EllipticCurve(GF(q), [0, b])
    fact = trivial_fact(order)
    if not fact:
        continue

    print(fact)


    GG = E.gens()
    for p, _ in fact:
        for G in GG:
            G_ord = G.order()
            if G_ord % p == 0:
                P = (G_ord//p) * G
                break
        try:
            lam = recov_lambda(E, p, P)
        except Exception as e:
            if p > 1000:
                continue

            out = nonsquarefree(query(P), p, P, E)
            if out:
                for elem in out:
                    remainders.append(elem[0])
                    moduli.append(elem[1])
                continue

        if p < 200:
            out = hashcoll(query(P), 0, p, P, E)
            if out:
                remainders.append(out[0])
                moduli.append(out[1])
                continue
            
        delta = p//nthreads
        intervals = [[x, x+delta] for x in range(0, p-delta, delta)]
        hsh = query(P)
        with Pool(nthreads) as poo:
            res = poo.starmap(hashcoll, [(hsh, beg, end, P, E) for beg, end in intervals])
            for elem in res:
                if elem:
                    remainders.append(elem[0])
                    moduli.append(elem[1])
                    break 

    print(remainders)
    print(moduli)
    

def handcrt(res, md):
    r0, m0 = res[0], md[0]
    for ri, mi in zip(res[1:], md[1:]):
        if m0 % mi == 0:
            continue
        
        g, x, y = xgcd(m0, mi)
        r0 = r0 * y * mi + ri * x * m0
        m0 = lcm(m0, mi)

    return r0, m0


def LLLCRT(res, md):
    A, _ = Sequence(res).coefficient_matrix()
    A = A.dense_matrix().T
    Q = diagonal_matrix(md)
    M = block_matrix([[A, 1],
                      [Q, 0]])

    weights = [1 for _ in range(A.ncols())] + [2^128]*2 + [1] 
    weights = [2^150 // w for w in weights]
    W = diagonal_matrix(weights)

    M = M*W
    M = M.change_ring(ZZ)
    M = M.BKZ(block_size = 10)
    M /= W

    for row in M:
        print(row)
        if set(row[:A.nrows()-1]) == set([0]) and row[-1] == 1:
            return (int(row[-3]), int(row[-2]))
        if set(row[:A.nrows()-1]) == set([0]) and row[-1] == -1:
            return (-int(row[-3]), -int(row[-2]))


# res, mod = handcrt(remainders, moduli)
# print(int(mod).bit_length())

# A, _= Sequence([res]).coefficient_matrix()
# A = A.dense_matrix().T.augment(identity_matrix(3))

# M = A.stack(vector([mod, 0, 0, 0]))

# weights = [1, 2^128, 2^128, 1]
# weights = [2^150 // w for w in weights]
# W = diagonal_matrix(weights)


# M = M*W
# M = M.change_ring(ZZ)
# M = M.BKZ(block_size=4)
# M /= W

# for row in M:
#     if row[0] == 0 and row[-1] == 1:
#         k1 = int(row[1])
#         k2 = int(row[2])
#         break

k1, k2 = LLLCRT(remainders, moduli)

print(k1, k2)


curve_order = 21888242871839275222246405745257275088548364400416034343698204186575808495617
a1 = 9931322734385697763 
b1 = -147946756881789319000765030803803410728 
a2 = 147946756881789319010696353538189108491 
b2 = 9931322734385697763

P.<eps1, eps2, s> = PolynomialRing(ZZ)
f1 = curve_order * k1 - curve_order * s + (b2 * s + eps1) * a1 + (-b1 * s + eps2) * a2
f2 = curve_order * k2 + (b2 * s + eps1) * b1 + (-b1 * s + eps2) * b2

A ,_ = Sequence([f1, f2]).coefficient_matrix()

A = A.dense_matrix().T

M = A.augment(identity_matrix(3)).change_ring(ZZ)

wt = [1, 1, 2^254, 2^254, 1]
wt = [2^256//w for w in wt]
W = diagonal_matrix(wt)

M = M * W
M = M.change_ring(ZZ)
M = M.BKZ(block_size=4)
M /= W

row = M[0]
e1 = ZZ(row[2])
e2 = row[3]
out = int((-e1 * pow(b2, -1, curve_order)) % curve_order)

payl = f'ecdh ({0},{0})'
io.send(payl)
io.recvuntil('session check: ')
print(io.recvline())
session_key = tuple([0, 0])
hash_check = blake2b((str(session_key) + "0").encode()).hexdigest()
print(hash_check)

from Crypto.Cipher import AES

key = blake2b((str(session_key) + "1").encode()).digest()[
    : AES.key_size[-1]
]

import os

iv = os.urandom(16)
cipher = AES.new(key, AES.MODE_CBC, iv=iv)
out = '0' * (16 - len(str(out)) % 16) + str(out)
print(out)
msg = str(out).encode()
ct = cipher.encrypt(msg)
payl = f'answer {(iv+ct).hex()}'
io.send(payl)
io.interactive()
