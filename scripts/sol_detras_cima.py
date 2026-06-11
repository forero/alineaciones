"""«Sol detrás de la cima» en Bogotá, solsticios 2026 (astropy).

¿Qué hace este script?
Imagina que estás parado cerca de la Plaza de Bolívar una mañana y miras
hacia los cerros orientales. En algún momento el Sol, que va subiendo,
queda EXACTAMENTE detrás de la cima de Monserrate (o de Guadalupe).
¿Desde qué lugares de la ciudad se ve esa alineación? Ese conjunto de
lugares es una línea sobre el mapa (un "locus"), y este script la calcula
para dos fechas especiales:

  - el solsticio de junio (21-jun-2026)    -> Sol detrás de Monserrate
  - el solsticio de diciembre (21-dic-2026) -> Sol detrás de Guadalupe

Además calcula dónde se cruzan las dos líneas (un punto desde el cual
ocurren AMBAS alineaciones, una en junio y otra en diciembre) y dibuja
todo sobre un mapa con las calles reales de Bogotá.

Supuestos: geometría sin refracción atmosférica, suelo plano a 2600 m,
hora local fija UTC-5.

Salidas (archivos que produce):
  data/efemerides_2026.csv             (salidas del Sol y mediodías)
  data/locus_junio_monserrate.csv      (línea de junio, punto por punto)
  data/locus_diciembre_guadalupe.csv   (línea de diciembre)
  data/interseccion_loci.csv           (el punto donde se cruzan)
  data/calles_osm.json                 (caché de calles de OpenStreetMap)
  figures/loci_sol_detras_cima.png     (la figura final)
"""

# ============================ IMPORTACIONES ============================
# "import" trae herramientas (librerías) ya escritas por otras personas.

import json                  # leer/escribir datos en formato JSON (texto)
import urllib.parse          # codificar texto para enviarlo por internet
import urllib.request        # descargar datos de una página de internet
from pathlib import Path     # manejar rutas de archivos y carpetas

import matplotlib            # librería para hacer gráficas

# "Agg" le dice a matplotlib que dibuje en memoria y guarde a archivo,
# sin intentar abrir una ventana (útil en un supercomputador sin pantalla).
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # interfaz de matplotlib para graficar
import numpy as np               # cálculos numéricos con arreglos (vectores)
import pandas as pd              # tablas de datos (como hojas de cálculo)

from astropy.utils import iers   # astropy: librería de astronomía

# astropy normalmente descarga de internet unas tablas de la rotación
# terrestre (IERS). Aquí lo desactivamos para que el script funcione
# igual con o sin internet (la precisión que perdemos es despreciable).
iers.conf.auto_download = False
iers.conf.auto_max_age = None

import astropy.units as u                                  # unidades físicas (grados, metros, horas)
from astropy.coordinates import AltAz, EarthLocation, get_sun  # posiciones astronómicas
from astropy.time import Time                              # manejo de fechas/horas astronómicas

# ======================== CARPETAS DE SALIDA ===========================
# __file__ es la ruta de ESTE archivo; .parents[1] sube dos niveles:
# de scripts/sol_detras_cima.py a la carpeta raíz del proyecto.
RAIZ = Path(__file__).resolve().parents[1]
DATA = RAIZ / "data"      # aquí guardaremos los CSV
FIGS = RAIZ / "figures"   # aquí guardaremos la figura
DATA.mkdir(exist_ok=True)  # crear la carpeta si no existe todavía
FIGS.mkdir(exist_ok=True)

# ============================ GEOMETRÍA ================================
# Trabajamos en un sistema de coordenadas "local": un plano con origen en
# la Plaza de Bolívar donde x apunta al Este y y apunta al Norte, ambos
# en metros. Se llama sistema ENU (East, North, Up).

R_TIERRA = 6_371_000.0  # radio medio de la Tierra en metros

LAT0 = 4.59806    # latitud de la Plaza de Bolívar (grados)
LON0 = -74.07611  # longitud de la plaza (negativa = oeste de Greenwich)
H0 = 2600.0       # altura del suelo en Bogotá (metros sobre el mar)

# Diccionario con la posición de cada cerro: latitud, longitud y altura.
CIMAS_GEO = {
    # Monserrate: centroide del edificio de la basílica según OSM
    "monserrate": {"lat": 4.60563, "lon": -74.05542, "h": 3152.0},
    "guadalupe":  {"lat": 4.59194, "lon": -74.05417, "h": 3360.0},
}

# Coseno de la latitud de la plaza. Se necesita porque los "grados de
# longitud" miden menos metros a medida que uno se aleja del ecuador.
_COS_LAT0 = np.cos(np.radians(LAT0))


def geo_a_enu(lat, lon, h):
    """Convierte (latitud, longitud, altura) a metros (E, N, U) desde la plaza.

    La idea: una diferencia de latitud o longitud (en radianes),
    multiplicada por el radio de la Tierra, da una distancia en metros
    (longitud de arco). Es la aproximación esférica, válida para
    distancias cortas como las de este problema.
    """
    # Diferencia de longitud -> metros hacia el Este (corregida por cos(lat))
    E = np.radians(lon - LON0) * R_TIERRA * _COS_LAT0
    # Diferencia de latitud -> metros hacia el Norte
    N = np.radians(lat - LAT0) * R_TIERRA
    # Altura sobre el nivel del suelo del observador (2600 m)
    U = h - H0
    return E, N, U


def enu_a_geo(x, y):
    """Operación inversa: de metros (x=Este, y=Norte) a (longitud, latitud)."""
    # Metros al este -> grados de longitud (deshaciendo la fórmula anterior)
    lon = LON0 + np.degrees(x / (R_TIERRA * _COS_LAT0))
    # Metros al norte -> grados de latitud
    lat = LAT0 + np.degrees(y / R_TIERRA)
    return lon, lat


def distancia_rumbo(x, y):
    """Distancia en línea recta y rumbo desde la plaza hasta el punto (x, y).

    El rumbo es el ángulo medido desde el Norte girando hacia el Este
    (Norte = 0°, Este = 90°, Sur = 180°, Oeste = 270°).
    """
    # np.hypot calcula sqrt(x² + y²), o sea el teorema de Pitágoras
    d = np.hypot(x, y)
    # arctan2(x, y) da el ángulo del punto; el % 360 lo deja entre 0 y 360
    rumbo = np.degrees(np.arctan2(x, y)) % 360.0
    return d, rumbo


# Convertimos de una vez las dos cimas a coordenadas ENU (metros).
# Por ejemplo Monserrate queda en E ≈ 2293 m, N ≈ 842 m, U ≈ 552 m.
CIMAS_ENU = {n: geo_a_enu(c["lat"], c["lon"], c["h"]) for n, c in CIMAS_GEO.items()}

# ============================== EL SOL =================================

# Bogotá usa todo el año la hora UTC-5 (5 horas detrás del meridiano de
# Greenwich). Lo guardamos como una cantidad con unidades de "horas".
UTC_OFFSET = -5 * u.hour

# Le decimos a astropy DÓNDE está el observador sobre la Tierra.
PLAZA = EarthLocation(lat=LAT0 * u.deg, lon=LON0 * u.deg, height=H0 * u.m)


def sol_altaz(tiempos_utc):
    """Para una lista de instantes, devuelve dónde está el Sol en el cielo.

    El resultado son dos ángulos por instante:
      az  (azimut):  dirección horizontal (Norte = 0°, Este = 90°, ...)
      alt (altitud): qué tan alto sobre el horizonte (0° = horizonte,
                     90° = justo encima de la cabeza).
    Sin refracción: como pressure no se especifica, astropy usa 0 y NO
    simula cómo la atmósfera "levanta" la imagen del Sol.
    """
    # El "marco" AltAz describe el cielo visto desde la plaza en esos instantes
    marco = AltAz(obstime=tiempos_utc, location=PLAZA)
    # get_sun da la posición del Sol; transform_to la traduce a az/alt
    s = get_sun(tiempos_utc).transform_to(marco)
    # .deg extrae los valores numéricos en grados
    return s.az.deg, s.alt.deg


def hora_local(t_utc):
    """Convierte un instante UTC a texto con la hora local de Bogotá (HH:MM:SS)."""
    # Sumar el desfase (-5 h) da la hora local; .iso es texto tipo
    # "2026-06-21 06:48:00.000" y [11:19] recorta solo "06:48:00".
    return (t_utc + UTC_OFFSET).iso[11:19]


def muestrear_dia(fecha, h_ini=5.0, h_fin=14.0, paso_s=30.0):
    """Genera instantes cada 30 s entre las 5:00 y las 14:00 hora local.

    Así "filmamos" la mañana entera y podemos seguir el movimiento del Sol
    paso a paso.
    """
    # Lista de horas: 5.0, 5.00833..., ... hasta 14.0 (en horas decimales)
    horas = np.arange(h_ini, h_fin, paso_s / 3600.0)
    # Medianoche local convertida a UTC (restar -5 h equivale a sumar 5 h)
    t0 = Time(f"{fecha} 00:00:00") - UTC_OFFSET
    # A la medianoche le sumamos cada hora de la lista -> instantes del día
    return t0 + horas * u.hour


# ========================== EFEMÉRIDES 2026 ============================
# "Efemérides" = tabla con datos astronómicos de fechas concretas.
# Aquí: ¿a qué hora sale el Sol, hacia dónde sale, y qué tan alto llega?

FECHAS = {"equinoccio marzo": "2026-03-20", "solsticio junio": "2026-06-21",
          "equinoccio septiembre": "2026-09-23",
          "solsticio diciembre": "2026-12-21"}

print("=== Efemérides solares 2026 en la Plaza de Bolívar (UTC-5, sin "
      "refracción) ===")
filas = []  # aquí iremos acumulando una fila de resultados por fecha
for evento, fecha in FECHAS.items():
    # Posición del Sol cada 30 segundos durante la mañana de esa fecha
    t = muestrear_dia(fecha)
    az, alt = sol_altaz(t)

    # SALIDA DEL SOL: buscamos el momento en que la altitud pasa de
    # negativa (bajo el horizonte) a positiva (sobre el horizonte).
    # alt[:-1] son todos los valores menos el último; alt[1:] todos menos
    # el primero. Comparándolos encontramos el "cruce" entre dos muestras.
    i = np.where((alt[:-1] < 0) & (alt[1:] >= 0))[0][0]
    # Interpolación lineal: ¿en qué fracción del intervalo se cruzó el 0?
    frac = -alt[i] / (alt[i + 1] - alt[i])
    # Con esa fracción estimamos la hora exacta y el azimut de salida
    t_salida = t[i] + frac * (t[i + 1] - t[i])
    az_salida = az[i] + frac * (az[i + 1] - az[i])

    # MEDIODÍA SOLAR: el instante en que el Sol alcanza su altura máxima.
    j = np.argmax(alt)        # posición del valor más grande de alt
    alt_mediodia = alt[j]
    t_mediodia = t[j]

    # Guardamos los resultados de esta fecha como un diccionario (una fila)
    filas.append({"evento": evento, "fecha": fecha,
                  "hora_salida_local": hora_local(t_salida),
                  "az_salida_deg": az_salida,
                  "hora_mediodia_local": hora_local(t_mediodia),
                  "alt_mediodia_deg": alt_mediodia})
    print(f"  {evento:22s} {fecha}: salida {hora_local(t_salida)} local, "
          f"az salida = {az_salida:6.2f}°, alt mediodía = "
          f"{alt_mediodia:5.2f}° ({hora_local(t_mediodia)})")

# Convertimos la lista de filas en una tabla y la guardamos como CSV
df_efe = pd.DataFrame(filas)
df_efe.to_csv(DATA / "efemerides_2026.csv", index=False, float_format="%.4f")
print(f"CSV escrito en {DATA / 'efemerides_2026.csv'}\n")


# ================= LOCI «SOL DETRÁS DE LA CIMA» ========================
# La idea geométrica clave:
# Si el Sol está a una altitud alt_sol y yo lo veo justo detrás de una
# cima cuya cumbre está U metros por encima de mí, entonces mi distancia
# horizontal a la cima tiene que ser exactamente
#     horiz = U / tan(alt_sol)
# (un triángulo rectángulo: cateto vertical U, ángulo alt_sol).
# Y además tengo que estar en la dirección OPUESTA al Sol vista desde la
# cima. Conociendo el azimut del Sol, eso fija mi posición (x, y).
# Como el Sol se mueve durante la mañana, cada instante da un punto
# distinto: todos juntos forman una línea sobre el mapa, el "locus".

def locus_detras_de(fecha, cima, alt_min=2.0, alt_max=60.0):
    """Calcula el locus de observadores que ven el Sol detrás de `cima`."""
    # Posición de la cima en metros desde la plaza (E, N) y su altura U
    E, N, U = CIMAS_ENU[cima]
    # Trayectoria del Sol durante esa mañana
    t = muestrear_dia(fecha)
    az, alt = sol_altaz(t)

    # Nos quedamos solo con el tramo MATUTINO: desde el comienzo del
    # muestreo hasta el momento de altura máxima (mediodía solar).
    j = np.argmax(alt)
    sel = slice(0, j + 1)          # "rebanada" de índices 0, 1, ..., j
    t, az, alt = t[sel], az[sel], alt[sel]
    # Y descartamos alturas demasiado bajas (puntos lejísimos) o muy
    # altas (puntos pegados a la cima): nos quedamos entre 2° y 60°.
    ok = (alt > alt_min) & (alt < alt_max)
    t, az, alt = t[ok], az[ok], alt[ok]

    # Las funciones trigonométricas de numpy trabajan en radianes
    az_r, alt_r = np.radians(az), np.radians(alt)
    # Distancia horizontal a la que debo estar de la cima (triángulo)
    horiz = U / np.tan(alt_r)
    # Retroceder desde la cima esa distancia, en dirección contraria al
    # Sol: sin(az) da la componente Este y cos(az) la componente Norte.
    x = E - horiz * np.sin(az_r)
    y = N - horiz * np.cos(az_r)

    # Para cada punto calculamos también su lon/lat y su distancia/rumbo
    # desde la plaza, y armamos una tabla con todo.
    lon, lat = enu_a_geo(x, y)
    d, rumbo = distancia_rumbo(x, y)
    return pd.DataFrame({"hora_local": [hora_local(ti) for ti in t],
                         "alt_sol_deg": alt, "az_sol_deg": az,
                         "x_east_m": x, "y_north_m": y, "lon": lon,
                         "lat": lat, "dist_plaza_m": d, "rumbo_deg": rumbo})


def reportar_locus(df, titulo, alts_repr=(5, 10, 15, 20, 30, 45)):
    """Imprime puntos representativos del locus y el más cercano a la plaza."""
    print(f"=== {titulo} ===")
    print("  puntos representativos (interpolados en altitud solar):")
    a = df["alt_sol_deg"].values   # columna de altitudes como arreglo
    for alt0 in alts_repr:         # para cada altitud "redonda" (5°, 10°, ...)
        # Si esa altitud no ocurre en este locus, la saltamos
        if not (a.min() <= alt0 <= a.max()):
            continue
        # searchsorted encuentra entre cuáles dos muestras cae alt0
        k = np.searchsorted(a, alt0)
        # Fracción del camino entre la muestra k-1 y la muestra k
        f = (alt0 - a[k - 1]) / (a[k] - a[k - 1])
        # Interpolamos linealmente las columnas numéricas en esa fracción
        num = df[["dist_plaza_m", "rumbo_deg", "lon", "lat"]]
        fila = num.iloc[k - 1] + f * (num.iloc[k] - num.iloc[k - 1])
        # la hora no se interpola numéricamente; tomar la más cercana
        hora = df.iloc[k if f > 0.5 else k - 1]["hora_local"]
        print(f"    alt_sol = {alt0:4.1f}° ({hora} local): "
              f"d = {fila['dist_plaza_m']:7.1f} m, rumbo = "
              f"{fila['rumbo_deg']:6.2f}°, lon = {fila['lon']:.6f}°, "
              f"lat = {fila['lat']:.6f}°")
    # idxmin da la fila cuya distancia a la plaza es la más pequeña
    imin = df["dist_plaza_m"].idxmin()
    c = df.loc[imin]
    print(f"  punto más cercano a la plaza: d = {c['dist_plaza_m']:.1f} m, "
          f"rumbo = {c['rumbo_deg']:.2f}°, alt_sol = {c['alt_sol_deg']:.2f}°,"
          f" {c['hora_local']} local ({c['lon']:.6f}°, {c['lat']:.6f}°)\n")
    return imin


# Calculamos los dos loci: junio/Monserrate y diciembre/Guadalupe
locus_jun = locus_detras_de("2026-06-21", "monserrate")
locus_dic = locus_detras_de("2026-12-21", "guadalupe")
# Guardamos cada locus completo (punto por punto) como CSV
locus_jun.to_csv(DATA / "locus_junio_monserrate.csv", index=False,
                 float_format="%.6f")
locus_dic.to_csv(DATA / "locus_diciembre_guadalupe.csv", index=False,
                 float_format="%.6f")

# Imprimimos el resumen de cada locus y recordamos cuál fila era la más
# cercana a la plaza (la usaremos para marcarla en la figura).
imin_jun = reportar_locus(locus_jun,
                          "Solsticio de junio: Sol detrás de Monserrate")
imin_dic = reportar_locus(locus_dic,
                          "Solsticio de diciembre: Sol detrás de Guadalupe")


# ================== INTERSECCIÓN DE LOS DOS LOCI =======================
# Cada locus es una "polilínea": muchos puntos unidos por segmentos
# rectos. Para hallar dónde se cruzan las dos, revisamos segmento contra
# segmento usando la fórmula paramétrica de la recta:
#   punto del segmento 1:  p + t·r   (con t entre 0 y 1)
#   punto del segmento 2:  q + u·s   (con u entre 0 y 1)
# Igualando se despejan t y u; si ambos quedan entre 0 y 1, los segmentos
# de verdad se tocan.

def interseccion_polilineas(P, Q):
    """Primer cruce entre las polilíneas P y Q (arreglos de puntos Nx2).

    Devuelve (punto, iP, fP, iQ, fQ): el punto de cruce, en qué segmento
    de cada polilínea ocurre (iP, iQ) y en qué fracción de ese segmento
    (fP, fQ), para luego poder interpolar la hora y la altitud del Sol.
    """
    # Recorremos cada segmento de la primera polilínea...
    for i in range(len(P) - 1):
        p, r = P[i], P[i + 1] - P[i]   # inicio y vector del segmento i
        # ...y lo comparamos contra TODOS los segmentos de la segunda a la
        # vez (numpy opera sobre arreglos completos, sin otro bucle).
        q = Q[:-1]                     # inicios de los segmentos de Q
        s = Q[1:] - Q[:-1]             # vectores de los segmentos de Q
        # "den" es el producto cruzado r×s; si es 0, son paralelos
        den = r[0] * s[:, 1] - r[1] * s[:, 0]
        con = np.abs(den) > 1e-12      # True donde NO son paralelos
        dq = q - p                     # vector de p al inicio de cada segmento
        # Parámetros t y u de la intersección (np.nan donde son paralelos)
        tt = np.where(con, (dq[:, 0] * s[:, 1] - dq[:, 1] * s[:, 0]) / den,
                      np.nan)
        uu = np.where(con, (dq[:, 0] * r[1] - dq[:, 1] * r[0]) / den, np.nan)
        # Hay cruce real solo si ambos parámetros caen entre 0 y 1
        hit = con & (tt >= 0) & (tt <= 1) & (uu >= 0) & (uu <= 1)
        if hit.any():
            k = np.argmax(hit)         # índice del primer True
            # El punto de cruce es p + t·r
            return p + tt[k] * r, i, tt[k], k, uu[k]
    return None  # las polilíneas no se cruzan


def evento_en(df, i, f):
    """Hora local y altitud del Sol en el segmento i, fracción f, del locus."""
    # Altitud: interpolación lineal entre las dos muestras del segmento
    alt = (1 - f) * df.iloc[i]["alt_sol_deg"] + f * df.iloc[i + 1]["alt_sol_deg"]
    # Hora: tomamos la de la muestra más cercana (el paso es de solo 30 s)
    hora = df.iloc[i if f < 0.5 else i + 1]["hora_local"]
    return hora, alt


# Extraemos de cada tabla solo las columnas (x, y) como arreglos Nx2
P = locus_jun[["x_east_m", "y_north_m"]].values
Q = locus_dic[["x_east_m", "y_north_m"]].values
res = interseccion_polilineas(P, Q)

print("=== Intersección de los dos loci ===")
if res is None:
    print("  (no se cruzan dentro del rango calculado)")
else:
    # Desempaquetamos el resultado: punto de cruce e índices/fracciones
    (xi, yi), iP, fP, iQ, fQ = res
    # Lo describimos de todas las formas útiles: lon/lat, distancia, rumbo
    lon_i, lat_i = enu_a_geo(xi, yi)
    d_i, rumbo_i = distancia_rumbo(xi, yi)
    # ¿A qué hora y con qué altura del Sol ocurre cada alineación allí?
    hora_jun, alt_jun = evento_en(locus_jun, iP, fP)
    hora_dic, alt_dic = evento_en(locus_dic, iQ, fQ)
    print(f"  E = {xi:+.1f} m, N = {yi:+.1f} m (desde la plaza)")
    print(f"  lon = {lon_i:.6f}°, lat = {lat_i:.6f}°")
    print(f"  distancia = {d_i:.1f} m, rumbo = {rumbo_i:.2f}°")
    print(f"  evento junio (detrás de Monserrate):   {hora_jun} local, "
          f"alt_sol = {alt_jun:.2f}°")
    print(f"  evento diciembre (detrás de Guadalupe): {hora_dic} local, "
          f"alt_sol = {alt_dic:.2f}°")
    # Guardamos el punto de intersección con todos sus datos en un CSV
    pd.DataFrame([{"x_east_m": xi, "y_north_m": yi, "lon": lon_i,
                   "lat": lat_i, "dist_plaza_m": d_i, "rumbo_deg": rumbo_i,
                   "hora_local_junio": hora_jun, "alt_sol_junio_deg": alt_jun,
                   "hora_local_diciembre": hora_dic,
                   "alt_sol_diciembre_deg": alt_dic}]).to_csv(
        DATA / "interseccion_loci.csv", index=False, float_format="%.6f")
    print(f"  CSV escrito en {DATA / 'interseccion_loci.csv'}")


# =================== CALLES DE BOGOTÁ (OpenStreetMap) ==================
# Para que la figura tenga contexto urbano, descargamos las calles del
# centro de Bogotá desde OpenStreetMap usando su API "Overpass".

def cargar_calles(x_rango=(-1100.0, 1100.0), y_rango=(-1100.0, 1100.0)):
    """Descarga (o lee de caché) las vías de OSM alrededor de la plaza.

    `x_rango`/`y_rango` son los límites en metros del recuadro que
    necesitamos cubrir. La primera vez descarga de internet y guarda el
    resultado en data/calles_osm.json; las siguientes veces lo lee de ahí
    (y solo vuelve a descargar si el recuadro pedido es más grande que el
    guardado). Si no hay internet ni caché, devuelve None.
    """
    # Convertimos las esquinas del recuadro de metros a lon/lat,
    # porque la API de OSM trabaja con coordenadas geográficas.
    lon_min, lat_min = enu_a_geo(x_rango[0], y_rango[0])
    lon_max, lat_max = enu_a_geo(x_rango[1], y_rango[1])
    # bbox = "bounding box", el rectángulo (sur, oeste, norte, este)
    bbox = [round(v, 5) for v in (lat_min, lon_min, lat_max, lon_max)]

    cache = DATA / "calles_osm.json"
    if cache.exists():
        # Ya hay datos guardados: los leemos
        guardado = json.loads(cache.read_text())
        b = guardado.get("bbox")
        # ¿El recuadro guardado cubre completamente el que pedimos ahora?
        if b and b[0] <= bbox[0] and b[1] <= bbox[1] \
                and b[2] >= bbox[2] and b[3] >= bbox[3]:
            return guardado["osm"]   # sí: usamos la caché, sin internet

    # Consulta en el lenguaje de Overpass: "dame todas las vías
    # (way[highway]) dentro de este rectángulo, con su geometría".
    consulta = (f"[out:json][timeout:60];way[highway]"
                f"({bbox[0]:.5f},{bbox[1]:.5f},{bbox[2]:.5f},{bbox[3]:.5f});"
                f"out geom;")
    # Preparamos la petición HTTP (el User-Agent identifica quién pregunta)
    req = urllib.request.Request(
        "https://overpass-api.de/api/interpreter",
        data=urllib.parse.urlencode({"data": consulta}).encode(),
        headers={"User-Agent": "alineaciones-bogota/0.1 (analisis academico)",
                 "Content-Type": "application/x-www-form-urlencoded"})
    try:
        # Enviamos la petición y leemos la respuesta (JSON con las vías)
        with urllib.request.urlopen(req, timeout=90) as r:
            datos = json.loads(r.read())
    except OSError as exc:
        # Sin internet (o el servidor falló): avisamos y seguimos sin calles
        print(f"AVISO: no se pudieron descargar las calles OSM ({exc}); "
              "la figura saldrá sin overlay de calles")
        return None
    # Guardamos la respuesta junto con su bbox para las próximas corridas
    cache.write_text(json.dumps({"bbox": bbox, "osm": datos}))
    print(f"Calles OSM: {len(datos['elements'])} vías, caché en {cache}")
    return datos


def dibujar_calles(ax, calles):
    """Dibuja las calles sobre los ejes `ax` (las principales más oscuras)."""
    principales = {"trunk", "primary", "secondary", "tertiary"}
    # Cada "elemento" de la respuesta es una vía (calle, carrera, sendero...)
    for via in calles["elements"]:
        geom = via.get("geometry")   # lista de puntos lat/lon de la vía
        if not geom:
            continue                 # algunas vías vienen sin geometría
        # Pasamos los puntos a arreglos de numpy...
        lat = np.array([p["lat"] for p in geom])
        lon = np.array([p["lon"] for p in geom])
        # ...y los convertimos a metros ENU para dibujarlos en nuestro mapa
        x, y, _ = geo_a_enu(lat, lon, H0)
        # Las vías importantes se dibujan más gruesas y oscuras.
        # zorder=1 las pone DEBAJO de los loci (que se dibujan después).
        if via.get("tags", {}).get("highway") in principales:
            ax.plot(x, y, color="dimgray", lw=1.0, alpha=0.8, zorder=1)
        else:
            ax.plot(x, y, color="gray", lw=0.4, alpha=0.5, zorder=1)


# ============================== FIGURA =================================

# Límites del mapa: corte en -500 m al oeste; al este lo justo para que
# se vean ambos cerros (la cima más oriental más un margen de 170 m para
# que quepa su etiqueta); ±1000 m en norte/sur.
XLIM = (-500, max(E for E, _, _ in CIMAS_ENU.values()) + 170)
YLIM = (-1000, 1000)

# Creamos la figura (el "lienzo") y unos ejes (el área de dibujo)
fig, ax = plt.subplots(figsize=(13, 9))
# Pedimos calles para un recuadro un poco más grande que el mapa (margen
# de 100 m) para que no queden bordes vacíos.
calles = cargar_calles(x_rango=(XLIM[0] - 100, XLIM[1] + 100),
                       y_rango=(YLIM[0] - 100, YLIM[1] + 100))
if calles is not None:
    dibujar_calles(ax, calles)

# Los loci se DIBUJAN solo hasta una altitud solar de 25° (más allá los
# puntos se pegan a los cerros y no aportan); los CSV sí tienen todo.
ALT_MAX_TRAZO = 25.0
trazo_jun = locus_jun[locus_jun["alt_sol_deg"] <= ALT_MAX_TRAZO]
trazo_dic = locus_dic[locus_dic["alt_sol_deg"] <= ALT_MAX_TRAZO]

# Dibujamos los dos loci con el mismo código, recorriendo una lista de
# tuplas (tabla, color, texto de la leyenda, fila más cercana a la plaza).
for df, color, etiqueta, imin in (
        (trazo_jun, "tab:orange",
         "21-jun: Sol detrás de Monserrate", imin_jun),
        (trazo_dic, "tab:blue",
         "21-dic: Sol detrás de Guadalupe", imin_dic)):
    # La línea del locus
    ax.plot(df["x_east_m"], df["y_north_m"], "-", color=color, lw=2,
            label=etiqueta)
    # Marcas sobre la línea cada 5° de altura del Sol, con su hora local
    a = df["alt_sol_deg"].values
    ultima_etiqueta = None   # recordamos dónde pusimos el último texto
    for alt0 in range(5, 61, 5):
        if not (a.min() <= alt0 <= a.max()):
            continue                       # esa altitud no está en el trazo
        k = np.searchsorted(a, alt0)       # muestra más cercana a alt0
        x0, y0 = df.iloc[k][["x_east_m", "y_north_m"]]
        # Solo dibujamos si el punto cae dentro del mapa
        if XLIM[0] < x0 < XLIM[1] and YLIM[0] < y0 < YLIM[1]:
            ax.plot(x0, y0, "o", color=color, ms=4)   # el puntito
            # El texto solo si no está a menos de 250 m del anterior
            # (para que no se amontonen las etiquetas)
            if (ultima_etiqueta is None
                    or np.hypot(x0 - ultima_etiqueta[0],
                                y0 - ultima_etiqueta[1]) > 250):
                ax.annotate(f"{alt0}° {df.iloc[k]['hora_local'][:5]}",
                            (x0, y0), xytext=(6, -4),
                            textcoords="offset points", fontsize=7,
                            color=color, annotation_clip=True)
                ultima_etiqueta = (x0, y0)
    # Cuadrito hueco en el punto del locus más cercano a la plaza
    c = df.loc[imin]
    ax.plot(c["x_east_m"], c["y_north_m"], "s", color=color, ms=9, mfc="none",
            mew=2)
    ax.annotate(f"más cercano: {c['dist_plaza_m']:.0f} m",
                (c["x_east_m"], c["y_north_m"]), xytext=(8, 8),
                textcoords="offset points", fontsize=8, color=color,
                annotation_clip=True)

if res is not None:
    # Estrella negra en la intersección de los dos loci, con sus
    # coordenadas tanto en metros E/N como en latitud/longitud
    ax.plot(xi, yi, "k*", ms=16, zorder=5)
    ax.annotate(f"intersección\nE/N: ({xi:+.0f}, {yi:+.0f}) m\n"
                f"lat: {lat_i:.5f}°\nlon: {lon_i:.5f}°",
                (xi, yi), xytext=(12, -40), textcoords="offset points",
                fontsize=9, annotation_clip=True,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray",
                          alpha=0.85))

    # RECUADRO DE ZOOM de 250 m de lado centrado en la intersección.
    # inset_axes crea un mini-mapa dentro del mapa grande; los cuatro
    # números son su posición y tamaño como fracción de los ejes grandes.
    axz = ax.inset_axes([0.50, 0.06, 0.20, 0.31],
                        xlim=(xi - 125.0, xi + 125.0),
                        ylim=(yi - 125.0, yi + 125.0))
    # Dentro del mini-mapa repetimos: calles, los dos loci y la estrella
    if calles is not None:
        dibujar_calles(axz, calles)
    axz.plot(trazo_jun["x_east_m"], trazo_jun["y_north_m"], "-",
             color="tab:orange", lw=2)
    axz.plot(trazo_dic["x_east_m"], trazo_dic["y_north_m"], "-",
             color="tab:blue", lw=2)
    axz.plot(xi, yi, "k*", ms=14, zorder=5)
    axz.set_aspect("equal")     # misma escala en x y en y
    axz.set_xticks([])          # sin números en los ejes del mini-mapa
    axz.set_yticks([])
    axz.set_title("zoom 250 m × 250 m", fontsize=8)
    # Dibuja el rectángulo en el mapa grande y las líneas que lo conectan
    # con el mini-mapa
    ax.indicate_inset_zoom(axz, edgecolor="black")

# Triángulo negro en la Plaza de Bolívar (el origen de coordenadas)
ax.plot(0, 0, "k^", ms=10)
ax.annotate("Plaza de Bolívar", (0, 0), xytext=(8, 6),
            textcoords="offset points", fontsize=9)
# Triángulos cafés en las dos cimas, con su nombre
for nombre in ("monserrate", "guadalupe"):
    E, N, _ = CIMAS_ENU[nombre]
    ax.plot(E, N, "^", color="saddlebrown", ms=11)
    ax.annotate(nombre.capitalize(), (E, N), xytext=(8, 6),
                textcoords="offset points", fontsize=9, color="saddlebrown",
                annotation_clip=True)

# Toques finales: límites, escala igual en ambos ejes, títulos, leyenda
ax.set_xlim(*XLIM)
ax.set_ylim(*YLIM)
ax.set_aspect("equal")   # 1 m en x mide lo mismo que 1 m en y
ax.set_xlabel("x = Este [m]")
ax.set_ylabel("y = Norte [m]")
ax.set_title("Loci de observadores con el Sol detrás de la cima "
             "(solsticios 2026, sin refracción)")
ax.legend(loc="lower left", fontsize=9)
ax.grid(alpha=0.3)       # cuadrícula suave de fondo

# Guardamos la figura como imagen PNG y cerramos para liberar memoria
salida = FIGS / "loci_sol_detras_cima.png"
fig.savefig(salida, dpi=150)
plt.close(fig)
print(f"\nFigura escrita en {salida}")
