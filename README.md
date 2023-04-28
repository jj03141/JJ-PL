# Kalkulator Geodezyjny
## Opis
Program umożliwia przeliczanie współrzędnych geodezyjnych między układami XYZ (współrzędne ortokartezjańskie), BLH (współrzędne geodezyjne) oraz PL-1992 i PL-2000. Program obsługuje elipsoidy: WGS84, GRS80 oraz elipsoidę Krasowskiego.
## Wymagania
- Python 3.9
- Biblioteki: math, argparse
## System operacyjny
Program został napisany dla systemu operacyjnego Windows.
## Instrukcja obsługi
Program umożliwia transformację XYZ -> BLH, BLH -> XYZ, XYZ -> PL-1992, XYZ -> PL-2000 oraz XYZ -> NEU dla wcześniej podanych elipsoid.
Dane wejściowe i wyjściowe programu są obsługiwane w formacie float.
Aby uzyskać współrzędne przeliczone na wybrany przez nas układ musimy otworzyć w oknie **cmd** ścieżkę do folderu z naszym plikiem (przykład ścieżki: *C:\Users\user\OneDrive\Pulpit\informatyka\projekt1*), \
a następnie wpisać "python" oraz nazwę naszego pliku. Dalej w tej samej linijce możemy wpisywać współrzędne XYZ. 

**Przykład poniżej:** \
*python "infa_projekt.py" 3664940.500 1409153.590 5009571.170*

Mamy także możliwość wyboru docelowego układu do którego chcemy przeliczyć nasze współrzędne oraz wyboru elipsoidy (spośród trzech wyżej wymienionych), a także wybrania pliku docelowego, do którego chcemy zapisać otrzymane wyniki.
W przypadku braku ich wybrania, zostaną automatycznie wybrane domyślne ustawienia tj. dla modelu "WGS84", dla układu "BLH" oraz dla pliku "output.txt".

**Przykład:** \
*3664940.500 1409153.590 5009571.170 --model grs80 --uklad BLH --output wyniki.txt*

**Przykład otrzymanych wyników:** \
*fi = 52.0973 [deg]; lam = 21.0315 [deg]; h = 141.399 [m] | BLH | grs80*
## Znane błędy i nietypowe zachowania
- Program zwraca błąd w przypadku podania niepoprawnego modelu elipsoidy lub systemu współrzędnych.
- Nieobsługiwane są niektóre formaty danych wejściowych (np. tekstowe).
