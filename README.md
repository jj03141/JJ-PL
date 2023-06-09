# Kalkulator Geodezyjny
## Opis
Program umożliwia przeliczanie współrzędnych geodezyjnych między układami XYZ (współrzędne ortokartezjańskie), BLH (współrzędne geodezyjne) oraz PL-1992 i PL-2000. Program obsługuje elipsoidy: WGS84, GRS80 oraz elipsoidę Krasowskiego.
## Wymagania
- Python 3.9
- Biblioteki: math, argparse, numpy
## System operacyjny
Program został napisany dla systemu operacyjnego Windows.
## Instrukcja obsługi
Program umożliwia transformację XYZ -> BLH, BLH -> XYZ, XYZ -> PL-1992, XYZ -> PL-2000 oraz XYZ -> NEU dla wcześniej podanych elipsoid.
Dane wejściowe i wyjściowe programu są obsługiwane w formacie float.
Istnieją dwie możliwości wprowadzenia współrzędnych: wprowadzanie ręczne (--input cmd), czyli opcja domyślna lub wczytanie pliku z rozszerzeniem .txt (--input txt).
Aby uzyskać współrzędne przeliczone na wybrany przez nas układ musimy otworzyć w oknie **cmd** ścieżkę do folderu z naszym plikiem (przykład ścieżki: *C:\Users\user\OneDrive\Pulpit\informatyka\projekt1*), \
a następnie wpisać "python" oraz nazwę naszego pliku. Dalej w tej samej linijce możemy wpisywać współrzędne XYZ (wartośći współrzędnych wyrażone w metrach) lub BLH (wartości B oraz L wprowadzamy w stopniach dziesiętnych, natomiast H w metrach) dla transformacji BLH -> XYZ. 

**Przykład poniżej:** \
*python "infa_projekt.py" 3664940.500 1409153.590 5009571.170*

Mamy także możliwość wyboru układu wprowadzanych przez nas współrzędnych,  docelowego układu do którego chcemy przeliczyć nasze współrzędne oraz wyboru elipsoidy (spośród trzech wyżej wymienionych), a także wybrania pliku docelowego, do którego chcemy zapisać otrzymane wyniki.
W przypadku braku ich wybrania, zostaną automatycznie wybrane domyślne ustawienia tj. dla układu współrzędnych wejściowych "XYZ" dla modelu "WGS84", dla układu "BLH" oraz dla pliku "output.txt".

**Przykład:** \
*python infa_projekt.py --dane XYZ 3664940.500 1409153.590 5009571.170 --model grs80 --uklad BLH --output wyniki.txt*

**Przykład otrzymanych wyników:** \
*fi = 52.0973 [deg]; lam = 21.0315 [deg]; h = 141.399 [m] | BLH | grs80*

W przypadku transformacji XYZ -> NEU musimy wprowadzić także współrzędne punktu referencyjnego, współrzędne punktu referencyjnego przyjmują domyślne wartości (100,100,100) w sytuacji gdy nie zostaną one podane. Współrzędne te są współrzędnymi w układzie XYZ.

**Przykład:** \
*python infa_projekt.py 3664940.500 1409153.590 5009571.170 -x_ref 3664941.500 -y_ref 1409152.590 -z_ref 5009570.170 --uklad NEU*

Jeśli chcemy przetransformować współrzędne zawarte w pliku tekstowym, wówczas  musimy podać jego nazwę (-txt wsp_wejściowe.txt), nazwę pliku wyjściowego (-txt_out wsp_wyjściowe.txt) oraz wybrać opcję --input txt. Ważne jest aby dane w pliku były rozdzielone przecinkiem, a także aby nie zawierał on spacji a współrzędne każdego punktu zaczynały się od nowego wiersza. Separatorem rozwinięcia dziesiętnego liczby powinna być kropka. Należy pamiętać, że plik powinien znajdować się w tym samym folderze roboczym co nasz program.

**Fragment pliku .txt gotowego do wczytania:** \
*3664940.500,1409153.590,5009571.170* \
*3664940.510,1409153.580,5009571.167* \
*3664940.520,1409153.570,5009571.167* \
*3664940.530,1409153.560,5009571.168* \
*3664940.520,1409153.590,5009571.170*

**Przykład:** \
*python infa_projekt.py --input txt --model grs80 -txt wsp_input.txt -txt_out wsp_output.txt*

Podobnie jak w przypadku współrzędnych wprowadzanych ręcznie mamy możliwość wybrania układu wyjściowego dla transformowanych współzędnych, a także typu współrzędnych wejściowych, czy modelu elipsoidy. 

## Znane błędy i nietypowe zachowania
- Program zwraca błąd w przypadku podania niepoprawnego modelu elipsoidy lub systemu współrzędnych.
- Program zwraca błąd dla transformacji XYZ -> BLH  w przypadku podania współrzędnych X=0 Y=0, dla których nie jest możliwe jednoznaczne określenie współrzędnych w układzie BLH. 
- Nieobsługiwane są niektóre formaty danych wejściowych (np. tekstowe).
- Transformacja XYZ -> PL-1992 oraz XYZ ->PL-2000 dla elipsoidy Krasowskiego zwraca błędne wyniki, nie powinna więc być używana.
