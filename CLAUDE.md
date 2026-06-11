# Spec: geometría alt/az de los cerros de Bogotá y alineaciones solares

Genera el código Python descrito abajo. Usa `numpy`, `pandas`, `matplotlib` y
`astropy`. Organiza el resultado en scripts reproducibles que guarden los CSV y
las figuras como archivos.

## Datos de entrada (WGS84)

| Punto | Latitud | Longitud | Elevación |
|---|---|---|---|
| Plaza de Bolívar (origen) | 4.59806 | −74.07611 | 2600 m (nivel del suelo) |
| Monserrate (cima) | 4.60583 | −74.05639 | 3152 m |
| Guadalupe (cima) | 4.59194 | −74.05417 | 3360 m |

Constante: `R = 6_371_000` m (radio terrestre medio).

## Marco de referencia

Sistema cartesiano local ENU centrado en la Plaza de Bolívar:
`x = Este`, `y = Norte` (metros), observador siempre a nivel del suelo (2600 m).
Convierte cada cima a ENU relativo a la plaza con la aproximación esférica
(`E = Δlon·deg2rad·R·cos(lat0)`, `N = Δlat·deg2rad·R`, `U = h − 2600`).

## 1. Función azimut/altitud

Desde un observador en `(x, y, 0)` hacia una cima en ENU `(E, N, U)`:

- `dE = E − x`, `dN = N − y`, `dU = U`
- `azimut = atan2(dE, dN)` en grados, normalizado a [0, 360); convención desde
  el norte verdadero, sentido horario (N=0, E=90).
- `horiz = hypot(dE, dN)`
- corrección por curvatura: `curv = horiz² / (2R)`
- `altitud = atan2(dU − curv, horiz)` en grados.

## 2. Grilla y CSV

Grilla cada 100 m dentro de un radio de 1 km de la plaza
(`arange(-1000, 1001, 100)` en x e y, máscara `x²+y² ≤ 1000²`). Para cada punto
calcula az y alt a Monserrate y Guadalupe. Exporta un CSV con columnas:
`x_east_m, y_north_m, lon, lat, az_monserrate, alt_monserrate, az_guadalupe, alt_guadalupe`.

## 3. Figuras de campo y mapa

- **Panel 2×2 de contornos** sobre la grilla: azimut y altitud de cada cima,
  enmascarando fuera del círculo de 1 km, con flecha roja desde el origen hacia
  cada cima y el círculo de 1 km punteado.
- **Mapa geográfico** en lon/lat: grilla, plaza, ambas cimas, líneas de visión,
  círculo de 1 km, flecha norte y barra de escala de 500 m. Usa
  `aspect = 1/cos(lat0)` para distancias correctas.
- Nota: los tiles de OpenStreetMap pueden estar bloqueados. Incluye además, en
  comentarios, el código con `contextily` + `geopandas`
  (`to_crs(3857)`, `cx.add_basemap(..., source=cx.providers.OpenStreetMap.Mapnik)`)
  para añadir el basemap real al ejecutarlo con internet.

## 4. Intersección de azimuts fijos y paralelogramo de incertidumbre

Cada condición «azimut fijo a una cima» es una recta de posiciones del
observador (recta que pasa por la cima):
`E − tan(β)·N = E_cima − tan(β)·N_cima`.

- Resuelve la intersección exacta (sistema 2×2) para
  `az(Monserrate) = 66°` y `az(Guadalupe) = 115°`. Reporta E/N desde la plaza,
  lon/lat, distancia y rumbo.
- Con incertidumbre **±1°** en ambos azimuts, calcula los 4 vértices del
  paralelogramo de incertidumbre (intersecciones de las rectas a 65/67° con
  114/116°). Reporta extensión E-O y N-S, diagonal máxima, área, semiejes
  principales (vía SVD) y la dilución geométrica `1/sin(ángulo de cruce)`.
- Superpón el paralelogramo (con recuadro de zoom y vértices etiquetados por su
  par de azimuts) sobre el panel 2×2 de campos.

## 5. Sol 2026 y «Sol detrás de la cima» (astropy)

Configura `EarthLocation` en la plaza, frame `AltAz`, `get_sun`, sin refracción,
zona horaria UTC−5. Desactiva descarga IERS
(`iers.conf.auto_download = False`, `iers.conf.auto_max_age = None`).

- **Efemérides** para solsticios y equinoccios de 2026 (≈ 20-mar, 21-jun,
  23-sep, 21-dic): hora de salida local, azimut de salida (alt cruza 0°) y
  altitud al mediodía.
- **Lugar donde el Sol queda detrás de la cima**: para cada instante del
  recorrido solar matutino con `(az_sol, alt_sol)`, mapea a un punto del
  observador:
  - `horiz = U_cima / tan(alt_sol)`
  - `x = E_cima − horiz·sin(az_sol)`
  - `y = N_cima − horiz·cos(az_sol)`
  Cada solución es una **recta/locus** de puntos (un punto por cada altura del
  Sol). Hazlo para:
  - **Solsticio de junio → detrás de Monserrate**
  - **Solsticio de diciembre → detrás de Guadalupe**
  Reporta puntos representativos (altura del Sol, hora local, distancia y rumbo
  desde la plaza, lon/lat) y el punto más cercano a la plaza. Dibuja ambos loci
  en metros E/N con marcas de altura solar; marca el punto más cercano a la
  plaza. Recorta las anotaciones al marco visible.
- **Intersección de los dos loci**: calcula el punto (E/N desde la plaza) donde
  se cruzan las dos rectas (intersección de polilíneas), su lon/lat, distancia y
  rumbo, y la hora local + altura del Sol de cada evento desde ese punto.

## Notas / supuestos

- Geometría sin refracción atmosférica (la refracción eleva el Sol bajo
  ~0.3–0.5°; afecta tiempos por 1–2 min cerca del horizonte).
- Suelo plano a 2600 m alrededor de la plaza (sabana de Bogotá).
- La elevación de Guadalupe es el dato menos cierto (3260–3360 m); se usa
  3360 m.
- Un desplazamiento del observador de ~1 km cambia la posición del Sol en
  <0.02°; el Sol se calcula en la plaza.
