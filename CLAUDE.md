# Spec: «Sol detrás de la cima» — cerros de Bogotá, solsticios 2026

El proyecto consiste en **un único script reproducible**,
`scripts/sol_detras_cima.py`, que calcula desde dónde se ve el Sol salir
exactamente detrás de Monserrate (solsticio de junio) y de Guadalupe
(solsticio de diciembre), y genera una figura y los CSV asociados. Usa
`numpy`, `pandas`, `matplotlib` y `astropy`.

El código debe estar **comentado línea a línea**, en español, de manera que
un estudiante de primer semestre de universidad pueda entender cada paso.

## Datos de entrada (WGS84)

| Punto | Latitud | Longitud | Elevación |
|---|---|---|---|
| Plaza de Bolívar (origen) | 4.59806 | −74.07611 | 2600 m (nivel del suelo) |
| Monserrate (basílica) | 4.60563 | −74.05542 | 3152 m |
| Guadalupe (cima) | 4.59194 | −74.05417 | 3360 m |

Constante: `R = 6_371_000` m (radio terrestre medio).

## Marco de referencia

Sistema cartesiano local ENU centrado en la Plaza de Bolívar:
`x = Este`, `y = Norte` (metros), observador siempre a nivel del suelo
(2600 m). Cada cima se convierte a ENU relativo a la plaza con la
aproximación esférica (`E = Δlon·deg2rad·R·cos(lat0)`,
`N = Δlat·deg2rad·R`, `U = h − 2600`). Azimut desde el norte verdadero,
sentido horario (N=0, E=90).

## Sol con astropy

`EarthLocation` en la plaza, frame `AltAz`, `get_sun`, sin refracción
(`pressure = 0`, el default), zona horaria fija UTC−5. Descarga IERS
desactivada (`iers.conf.auto_download = False`,
`iers.conf.auto_max_age = None`).

### Efemérides 2026

Para solsticios y equinoccios de 2026 (20-mar, 21-jun, 23-sep, 21-dic):
hora de salida local (alt cruza 0° interpolando linealmente), azimut de
salida y altitud máxima (mediodía solar). → `data/efemerides_2026.csv`.

### Loci «Sol detrás de la cima»

Para cada instante del recorrido solar matutino con `(az_sol, alt_sol)`,
el observador que ve el Sol exactamente detrás de la cima está en:

- `horiz = U_cima / tan(alt_sol)`
- `x = E_cima − horiz·sin(az_sol)`
- `y = N_cima − horiz·cos(az_sol)`

Cada solución es una **recta/locus** de puntos (uno por cada altura del
Sol). Se calcula para:

- **Solsticio de junio (21-jun) → detrás de Monserrate**
- **Solsticio de diciembre (21-dic) → detrás de Guadalupe**

Se reportan puntos representativos (altura del Sol, hora local, distancia
y rumbo desde la plaza, lon/lat) y el punto más cercano a la plaza.
→ `data/locus_junio_monserrate.csv`, `data/locus_diciembre_guadalupe.csv`.

### Intersección de los dos loci

Punto (E/N desde la plaza) donde se cruzan las dos polilíneas, su lon/lat,
distancia y rumbo, y la hora local + altura del Sol de cada evento visto
desde ese punto. → `data/interseccion_loci.csv`.

## Figura (`figures/loci_sol_detras_cima.png`)

Ambos loci en metros E/N con:

- **Overlay de las cuadras y calles de Bogotá** descargadas de
  OpenStreetMap vía la API de Overpass (`way[highway]` en el bbox del
  marco), con caché local en `data/calles_osm.json` que guarda su bbox y
  se renueva sola si no cubre el recuadro pedido. Vías principales
  (trunk/primary/secondary/tertiary) en gris oscuro, el resto en gris
  claro. Si no hay red ni caché, la figura sale sin calles (con aviso).
- **Marco**: corte en −500 m al oeste; al este lo justo para que se vean
  ambos cerros (E de la cima más oriental + margen para su etiqueta);
  ±1000 m en norte/sur. `aspect` igual.
- Las líneas de los loci se dibujan **solo hasta una altitud solar de
  25°** (los CSV sí guardan el locus completo, hasta 60°).
- Marcas de altura solar cada 5° con hora local (sin etiquetar puntos
  amontonados), punto más cercano a la plaza de cada locus, plaza y cimas
  marcadas, intersección con estrella.
- **Recuadro de zoom de 250 m × 250 m** centrado en la intersección, con
  calles y ambos loci, conectado al mapa con `indicate_inset_zoom`.
- Recorta las anotaciones al marco visible (`annotation_clip`).

## Cómo correr (NERSC Perlmutter)

```bash
module load python
python3 scripts/sol_detras_cima.py
```

## Notas / supuestos

- Geometría sin refracción atmosférica (la refracción eleva el Sol bajo
  ~0.3–0.5°; afecta tiempos por 1–2 min cerca del horizonte).
- Suelo plano a 2600 m alrededor de la plaza (sabana de Bogotá).
- La elevación de Guadalupe es el dato menos cierto (3260–3360 m); se usa
  3360 m.
- La posición de Monserrate es el centroide del edificio de la basílica
  según OSM (4.60563, −74.05542); la cota usada (3152 m) es la citada para
  la explanada del santuario.
- Un desplazamiento del observador de ~1 km cambia la posición del Sol en
  <0.02°; el Sol se calcula en la plaza.
