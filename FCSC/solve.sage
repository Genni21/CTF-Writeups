import json
from Crypto.Hash import SHA256
from Crypto.Cipher import AES
from Crypto.Util.Padding import unpad

Q=8383489
n=512
P.<y>=PolynomialRing(Integers(Q))
S.<x>=QuotientRing(P, y^n+1)

def sw(arr):
    p=0
    for elem in range(len(arr)):
        p=p+ arr[elem]*x^elem

    return p




f=open("output.txt", "r")
tup=json.loads(f.readline())
collective=json.loads(f.readline())
actual=collective[0]
z1=sw(actual['z1'])
y1=sw(actual['y1'])
c=sw(actual['c'])
s1=(z1 - y1) /c
t=sw(tup['t'])
a=sw(tup['a'])
s2=t - (s1*a)

key = SHA256.new(str(s1.list()).encode() + str(s2.list()).encode()).digest()

last=json.loads(f.readline())
from binascii import unhexlify

iv= unhexlify(last['iv'])
enc=unhexlify(last['enc'])

E = AES.new(key, AES.MODE_CBC, iv)

print(unpad(E.decrypt(enc), 16)[:-1].decode())
