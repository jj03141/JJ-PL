from math import *

class transformacje:
    def __init__(self, model: str = 'wgs84'):
        if model == 'wgs84':
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == 'grs80':
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == 'krasowski':
            self.a = 6378245.000
            self.b = 6356863.018773
        else:
            raise NotImplementedError(f'{model} model not implemented')
        self.flat = (self.a - self.b) / self.a
        self.e = sqrt(2 * self.flat - self.flat ** 2)
        self.e2 = (2 * self.flat - self.flat ** 2)
    
    def dms(self, x):
        sig = ' '
        if x < 0:
            sig = ' - '
            x = abs(x)
        x = x * 180/pi
        d = int(x)
        m = int(60 * (x - d))
        s = (x - d - m/60) * 3600
        return sig, "%3d %2d %7.5f" % (d,m,s)
    
    def hirvonen(self, X, Y, Z, output='dec_degree'):
        p = sqrt(X**2 + Y**2)
        f = atan(Z/(p*(1 - self.e2)))
        while True:
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)
            h = (p / cos(f)) - N
            fp = f
            f = atan(Z/(p*(1 - self.e2 * (N / (N + h)))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = atan2(Y, X)
        if output == 'dec_degree':
            return degrees(f), degrees(l), h
        elif output == 'dms':
            f = self.dms(degrees(f))
            l = self.dms(degrees(l))
            return f, l, h
        else:
            raise NotImplementedError(f'{output} - output format not defined')

if __name__ == '__main__':
    geo = transformacje(model='wgs84')
    X = 3664940.500
    Y = 1409153.590
    Z = 5009571.170
    fi, lam, h = geo.hirvonen(X, Y, Z, output='dms')
    print(fi, lam, h)
