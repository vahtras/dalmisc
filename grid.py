def readfirst(quadfile):
    with open(quadfile, 'rb') as qf:
        for f in quaditer(qf):
            return f

def quaditer(qf):
    import struct
    isize = 4
    dsize = 8
    while True:
        ibytes = qf.read(2*isize)
        if not ibytes: raise StopIteration
        n, nb = struct.unpack(2*'i', ibytes)
        ibytes = qf.read(2*nb*isize)
        fbytes = qf.read(4*n*dsize)
        xyzw = struct.unpack(4*n*'d', fbytes)
        for i in range(n):
            yield xyzw[4*i: 4*(i+1)]
        

if __name__ == "__main__":
    import sys
    with open(sys.argv[1], 'rb') as qf:
        for x,y,z,w in quaditer(qf):
            print x,y,z,w
    
