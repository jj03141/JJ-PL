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
        f_dg = degrees(f)
        l_dg = degrees(l)
        while True:
            f_dg = f_dg - 360
            if f_dg - 360 < 0:
                break
        while True:
            l_dg = l_dg - 360
            if l_dg - 360 < 0:
                break
        if output == 'dec_degree':
            return f_dg, l_dg, h
        elif output == 'dms':
            f = self.dms(f_dg)
            l = self.dms(l_dg)
            return f, l, h
        else:
            raise NotImplementedError(f'{output} - output format not defined')
            
    def flh2XYZ(self, f, l, h):
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        X = (N + h) * cos(f) * cos(l)
        Y = (N + h) * cos(f) * sin(l)
        Z = (N * (1 - self.e2) + h) * sin(f)
        return X, Y, Z
    
    def BL292(self, f, l):
        l0 = radians(19)
        a = self.a
        e2 = self.e2
        a2 = a**2
        b2 = a2 * (1 - e2)
        e_2 = (a2 - b2)/b2
        dl = l - l0
        dl2 = dl**2
        dl4 = dl**4
        t = tan(f)
        t2 = t**2
        t4 = t**4
        n2 = e_2 * (cos(f)**2)
        n4 = n2 ** 2
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        e4 = e2**2
        e6 = e2**3
        A0 = 1 - (e2/4) - ((3*e4)/64) - ((5*e6)/256)
        A2 = (3/8) * (e2 + e4/4 + (15*e6)/128)
        A4 = (15/256) * (e4 + (3*e6)/4)
        A6 = (35*e6)/3072
        sigma = a * ((A0 * f) - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        xgk = sigma + ((dl**2)/2) * N * sin(f) * cos(f) * (1 + ((dl**2)/12)*(cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
        ygk = dl * N * cos(f) * (1 + (dl2/6) * (cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
        x_92 = xgk * 0.9993 - 5300000
        y_92 = ygk * 0.9993 + 500000
        return x_92, y_92
    
    def BL200(self, f, l):
        a = self.a
        e2 = self.e2
        if abs(l - 15) <= 1.5:
            l0_deg = 15
        elif abs(l - 18) < 1.5:
            l0_deg = 18
        elif abs(l - 21) <= 1.5:
            l0_deg = 21
        else:
            l0_deg = 24
        l0 = radians(l0_deg)
        a2 = a**2
        b2 = a2 * (1 - e2)
        e_2 = (a2 - b2)/b2
        dl = l - l0
        dl2 = dl**2
        dl4 = dl**4
        t = tan(f)
        t2 = t**2
        t4 = t**4
        n2 = e_2 * (cos(f)**2)
        n4 = n2 ** 2
        N = Np(f,a,e2)
        e4 = e2**2
        e6 = e2**3
        A0 = 1 - (e2/4) - ((3*e4)/64) - ((5*e6)/256)
        A2 = (3/8) * (e2 + e4/4 + (15*e6)/128)
        A4 = (15/256) * (e4 + (3*e6)/4)
        A6 = (35*e6)/3072
        sigma = a * ((A0 * f) - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6*f))
        xgk = sigma + ((dl**2)/2) * N * sin(f) * cos(f) * (1 + ((dl**2)/12)*(cos(f)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (cos(f)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
        ygk = dl * N * cos(f) * (1 + (dl2/6) * (cos(f)**2) * (1 - t2 + n2) + (dl4/120) * (cos(f)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
        strefa = int(l0 * 180/pi)/3
        x_00 = xgk * 0.999923
        y_00 = ygk * 0.999923 + strefa * 1000000 + 500000
        return x_00, y_00
    
    def xyz2neu(self, x, y, z):
        a = self.a
        e2 = self.e2
        f = self.hirvonen(x, y, z)[0]
        l = self.hirvonen(x, y, z)[1]
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        
        T = np.array([[-sin(l)         , cos(l),                0],
                      [-sin(f) * cos(l), -sin(f) * sin(l), cos(f)],
                      [cos(f) * cos(l) , cos(f) * sin(l), sin(f)]])
        '''
        N = [-T[0][0] * X - T[1][0] * Y - T[2][0] * Z,
             -T[0][1] * X - T[1][1] * Y - T[2][1] * Z,
             -T[0][2] * X - T[1][2] * Y - T[2][2] * Z]

        E = [-T[0][0] * X - T[1][0] * Y - T[2][0] * Z,
             -T[0][1] * X - T[1][1] * Y - T[2][1] * Z,
             -T[0][2] * X - T[1][2] * Y - T[2][2] * Z]
        
        U = X * sin(f) - Y * sin(l) * cos(f) - Z * cos(l) * cos(f)
        '''
        N = -sin(f) * cos(l) * x - sin(f) * sin(l) * y + cos(f) * z
        E = -sin(l) * x + cos(l) * y
        U = cos(f) * cos(l) * x + cos(f) * sin(l) * y  + sin(f) * z

        return N, E, U

        



if __name__ == '__main__':
    geo = transformacje(model='wgs84')
    X = 3664940.500
    Y = 1409153.590
    Z = 5009571.170
    fi, lam, h = geo.xyz2neu(X, Y, Z)
    print(fi, lam, h)


