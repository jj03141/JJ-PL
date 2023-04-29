from math import *
import argparse as arg
import numpy as np

class transformacje:
    def __init__(self, model: str = 'wgs84'):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            e2 - mimośród^2
        """
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
    
    
    def hirvonen(self, X, Y, Z):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (X, Y, Z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (fi, lambda, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 mm.     
        Parametery
        ----------
        X, Y, Z : FLOAT
            współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        f
            [stopnie dziesiętne] - szerokość geodezyjna
        l
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        """
        a = self.a
        e2 = self.e2
        p = sqrt(X**2 + Y**2)
        f = atan(Z/(p*(1 - e2)))
        while True:
            N = self.a / sqrt(1 - self.e2 * sin(f)**2)
            h = (p / cos(f)) - N
            fp = f
            f = atan(Z/(p*(1 - e2 * (N / (N + h)))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = atan2(Y,X)
        return(f,l,h)
            
    def flh2XYZ(self, f, l, h):
        """
        Odwrotny algorytm Hirvonena - algorytm transformacji współrzędnych geodezyjnych 
        długość szerokość i wysokośc elipsoidalna(fi, lambda, h) na współrzędne ortokartezjańskie  (X, Y, Z).

        Parametery
        ----------
        f, l, h : FLOAT
            [dec_degree] współrzędne geodezyjne, 

        Returns
        -------
        X, Y, Z : FLOAT
            [metry] współrzędne ortokartezjańskie
        """
        f = radians(f)
        l = radians(l)
        N = self.a / sqrt(1 - self.e2 * sin(f)**2)
        X = (N + h) * cos(f) * cos(l)
        Y = (N + h) * cos(f) * sin(l)
        Z = (N * (1 - self.e2) + h) * sin(f)
        return X, Y, Z
    
    def BL292(self, x, y, z):
        """
        Transformacja f,l -> PL-1992 - algorytm transformacji współrzędnych elipsoidalnych (f, l)
        na współrzędne płaskie w układzie odniesienia PL-1992 (x_92, y_92). 
   
        Parametery
        ----------
        X, Y, Z : FLOAT
            [metry] współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        x_92, y_92 : FLOAT
            [metry] współrzędne płaskie w układzie PL-1992
        """
        f = self.hirvonen(x, y, z)[0]
        l = self.hirvonen(x, y, z)[1]
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
    
    def BL200(self, x, y, z):
        """
        Transformacja f,l -> PL-2000 - algorytm transformacji współrzędnych elipsoidalnych (f, l)
        na współrzędne płaskie w układzie odniesienia PL-2000 (x_00, y_00). 

        Parametery
        ----------
        X, Y, Z : FLOAT
            [metry] współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        x_00, y_00 : FLOAT
            [metry] współrzędne płaskie w układzie PL-2000
        """
        a = self.a
        e2 = self.e2
        f = self.hirvonen(x, y, z)[0]
        l = self.hirvonen(x, y, z)[1]
        if abs(degrees(l) - 15) <= 1.5:
            l0_deg = 15
        elif abs(degrees(l) - 18) < 1.5:
            l0_deg = 18
        elif abs(degrees(l) - 21) <= 1.5:
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
        strefa = int(l0 * 180/pi)/3
        x_00 = xgk * 0.999923
        y_00 = ygk * 0.999923 + strefa * 1000000 + 500000
        return x_00, y_00
    
    def Rneu(self, f, l):
        """
        Macierz R w transformacji współrzędnych XYZ na NEU jest macierzą rotacji, która pozwala przeliczyć 
        współrzędne z układu kartezjańskiego na współrzędne związanego z Ziemią układu współrzędnych geodezyjnych NEU.
        Wykorzystujemy bibliotekę numpy.
        
        Parametery
        ----------
        fi, lam: FLOAT
            [dec_degree] współrzędne fi, lambda w układzie geodezyjnym, 
       
        Returns
        -------
        R : array
            [niemianowane] macierz rotacji
        """
        R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                      [-np.sin(f) * np.sin(l), np.cos(l), np.cos(f) * np.sin(l)],
                      [np.cos(f), 0, np.sin(f)]])
        return(R)
    
    def xyz2neu(self, xa, ya, za, xb, yb, zb):
        """
        Transformacja XYZ -> NEU - algorytm transformacji współrzędnych wektora pomiędzy dwoma punktami w układzie współrzędnych 
        ortokartezjańskich (X, Y, Z) na współrzędne wektora pomiędzy dwoma punktami w układzie NEU: North, East, Up (N, E, U). 
        Wykorzystujemy bibliotekę numpy.

        Parametery
        ----------
        Xa, Ya, Za : FLOAT
            [metry] współrzędne w układzie orto-kartezjańskim, 
        
        Xb, Yb, Zb : FLOAT
            [metry] współrzędne punktu referencyjnego w układzie orto-kartezjańskim, 

        Returns
        -------
        N, E, U : FLOAT
            [metry] współrzędne w układzie NEU
        """
        
        a = self.a
        e2 = self.e2
        dxyz = np.array([xb, yb, zb]) - np.array([xa, ya, za])
        f = self.hirvonen(xa, ya, za)[0]
        l = self.hirvonen(xa, ya, za)[1]
        R = self.Rneu(f, l)
        dneu = -np.linalg.solve(R, dxyz)
        N = dneu[0]
        E = dneu[1]
        U = dneu[2]

        return N, E, U


X = []
Y = []
Z = []
F = []
L = []
H = []
X92 = []
Y92 = []
X00 = []
Y00 = []
N = []
E = []
U = []

with open('wsp_inp.txt', 'r') as plik:
    '''
    Wczytanie pliku i wyodrębnienie podanych w nim współrzędnych 
    za pomocą pętli for. Odczytane dane dodajemy do list, a następnie 
    transformujemy je do innych układów współrzędnych.
    '''
    lines = plik.readlines()
    t = 0
    for i in lines:
        t = t+1
        if t > 4:
            x = i.split(',')
            X.append(float(x[0]))
            Y.append(float(x[1]))
            Z.append(float(x[2]))
            
if __name__ =='__main__':
    '''
    Wykonujemy transformacje wczesniej wyodrębnionych współrzędnych i 
    przypisujemy je do specjalnie utworzonych dla każdego układu list.
    '''
    geo = transformacje(model = 'wgs84')
    for a, b, c in zip(X, Y, Z):
        f, l, h = geo.hirvonen(a, b, c)
        F.append(degrees(f))
        L.append(degrees(l))
        H.append(h)
        x92, y92 = geo.BL292(a, b, c)
        X92.append(x92)
        Y92.append(y92)
        x00, y00 = geo.BL200(a, b, c)
        X00.append(x00)
        Y00.append(y00)
        n, e, u = geo.xyz2neu(a, b, c, 100, 100, 100)
        N.append(n)
        E.append(e)
        U.append(u)

'''
Tworzymy plik, w którym zapisujemy uzyskane dla każdego układu wyniki,
wykorzystujemy do tego metodę f-string.
'''        
plik=open("wsp_out.txt","w")
plik.write(f'Współrzędne BLH, PL-1992, PL-2000, NEU stacji permanentnej GNSS \n')
plik.write(f'Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu \n')
plik.write(f'# ----------------------------------------------------- \n')
plik.write(f'# BLH-------------------------------------------------- \n')
plik.write(f'  B[d]         L[d]         H[m] \n')
plik.write(f'# ----------------------------------------------------- \n')
for a,b,c in zip(F,L,H):
    a = f'{a:7.4f}'
    b = f'{b:7.4f}'
    c = f'{c:7.4f}'
    plik.write(f'{a},      {b},      {c} \n')

plik.write(f'# ----------------------------------------------------- \n')
plik.write(f'# PL-1992---------------------------------------------- \n')
plik.write(f'  X[m]         Y[m] \n')
plik.write(f'# ----------------------------------------------------- \n')
for a,b in zip(X92,Y92):
    a = f'{a:7.3f}'
    b = f'{b:7.3f}'
    plik.write(f'{a},   {b} \n')
    
plik.write(f'# ----------------------------------------------------- \n')
plik.write(f'# PL-2000---------------------------------------------- \n')
plik.write(f'  X[m]         Y[m] \n')
plik.write(f'# ----------------------------------------------------- \n')
for a,b in zip(X00,Y00):
    a = f'{a:7.3f}'
    b = f'{b:7.3f}'
    plik.write(f'{a},   {b} \n')

plik.write(f'# ----------------------------------------------------- \n')
plik.write(f'# NEU---------------------------------------------- \n')
plik.write(f'  N[m]         E[m]         U[m] \n')
plik.write(f'# ----------------------------------------------------- \n')

for a,b,c in zip(N,E,U):
    a = f'{a:7.3f}'
    b = f'{b:7.3f}'
    c = f'{c:7.3f}'
    plik.write(f'{a},   {b},      {c} \n')
plik.close()

if __name__ == '__main__':
    '''
    Utworzenie programu do transformacji współrzędnych geodezyjnych między 
    różnymi układami odniesienia i elipsoidami, który odczytuje argumenty z wiersza poleceń 
    i zapisuje wyniki do pliku tekstowego. W zależności od wybranego układu i modelu elipsoidy,
    program wywołuje odpowiednie funkcje transformacji i zapisuje wyniki do konsoli i pliku.
    Wykorzystujemy bibliotekę argparse.
    '''
    parser = arg.ArgumentParser(description = 'XYZ -> PL-1992')
    parser.add_argument('--dane', type = str, choices = ['XYZ', 'BLH'], default = 'XYZ', help = 'Typ wprowadzanych współrzędnych (BLH lub XYZ), domyslnie: XYZ' )
    parser.add_argument('x', type = float, help = 'Współrzędna X')
    parser.add_argument('y', type = float, help = 'Współrzędna Y')
    parser.add_argument('z', type = float, help = 'Współrzędna Z')
    parser.add_argument('-x_ref', type = float, help = 'Współrzędna X punktu referencyjnego', default = 100.00)
    parser.add_argument('-y_ref', type = float, help = 'Współrzędna Y punktu referencyjnego', default = 100.00)
    parser.add_argument('-z_ref', type = float, help = 'Współrzędna Z punktu referencyjnego', default = 100.00)
    parser.add_argument('--model', type = str, choices = ['wgs84', 'grs80', 'krasowski'], default = 'wgs84', help = 'Model elipsoidy (wgs84, grs80 lub krasowski), domyslnie: wgs84')
    parser.add_argument('--uklad', type = str, choices = ['PL-1992', 'PL-2000', 'BLH', 'NEU'], default = 'BLH', help= 'System współrzędnych (PL-1992, PL-2000, BLH, NEU), domyslnie: BLH')
    parser.add_argument('--output', type = str, default = 'output.txt', help = 'Nazwa pliku z wynikami, domyslnie: output.txt')
    args = parser.parse_args()
    
    geo = transformacje(model=args.model)
    if args.dane == 'XYZ':
        if args.uklad == 'PL-1992':
            x92, y92 = geo.BL292(args.x, args.y, args.z)
            wynik = f'X = {x92:.3f} [m]; Y = {y92:.3f} [m] | {args.uklad} | {args.model}'
            print(wynik)
        elif args.uklad == 'PL-2000':
            x00, y00 = geo.BL200(args.x, args.y, args.z)
            wynik = f'X = {x00:.3f} [m]; Y = {y00:.3f} [m] | {args.uklad} | {args.model}'
            print(wynik)
        elif args.uklad == 'BLH':
            f, l, h = geo.hirvonen(args.x, args.y, args.z)
            wynik = f'fi = {degrees(f):.4f} [deg]; lam = {degrees(l):.4f} [deg]; h = {h:.3f} [m] | {args.uklad} | {args.model}'
            print(wynik)
        elif args.uklad == 'NEU':
            n, e, u = geo.xyz2neu(args.x, args.y, args.z, args.x_ref, args.y_ref, args.z_ref)
            wynik = f'N = {n:.3f} [m]; E = {e:.3f} [m]; U = {u:.3f} [m] | {args.uklad} | {args.model}'
            print(wynik)
    elif args.dane == 'BLH':
        X, Y, Z = geo.flh2XYZ(args.x, args.y, args.z)
        wynik = f'X = {X:.3f} [m]; Y = {Y:.3f} [m]; Z = {Z:.3f} [m] | wsp. ortokartezjańskie |'
        print(wynik)
        
    with open(args.output, 'a') as file:
        file.write(wynik + '\n')

        



