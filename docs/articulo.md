> **Nota:** borrador v1 en español, conservado como referencia. La versión vigente del artículo (v2, en inglés, con la matriz del estado del arte y los tres casos de estudio) es [`paper.md`](paper.md).

# NURBS-Maxima: un paquete de álgebra computacional para el modelado exacto de curvas y superficies NURBS y su geometría diferencial

**Borrador v0.1**. Autores: *[completar]*. Afiliación: *[completar]*. Correspondencia: *[completar]*.

> Nota para los autores: este borrador sigue la estructura habitual de una revista de software
> científico (por ejemplo, *SoftwareX* o *Journal of Open Source Software*) y de revistas de
> modelado geométrico. Las cifras de la Sección 5 proceden de `tests/test_nurbs.mac` ejecutado
> con Maxima 5.46.0 (SBCL) y geomdl 5.5.1 en un Intel Xeon a 2,1 GHz. Antes de enviar el artículo,
> repita las mediciones en el equipo definitivo y compruebe las referencias.

## Resumen

Las B-splines racionales no uniformes (NURBS) son el estándar de facto para la representación
de geometría en diseño asistido por computador (CAD), gráficos y análisis isogeométrico.
Las bibliotecas existentes están escritas en lenguajes de cálculo numérico (C++, Python) y trabajan
exclusivamente en coma flotante. Presentamos **NURBS-Maxima**, un paquete de código abierto para
el sistema de álgebra computacional Maxima. Implementa la evaluación, la diferenciación y la
manipulación de curvas y superficies NURBS y el cálculo de sus invariantes de geometría diferencial,
en aritmética exacta o en coma flotante, con el mismo código. El paquete cubre funciones base
B-spline y sus derivadas, curvas y superficies racionales con derivadas de orden arbitrario,
transformaciones afines y proyectivas, inserción de nudos, interpolación global, primitivas cuádricas
exactas, curvaturas gaussiana, media y principales, detección de singularidades, integrales globales
(área, volumen y curvatura total) y exportación a formatos de visualización y de intercambio.
La validación combina pruebas simbólicas (por ejemplo, la identidad `x² + y² = 1` sobre la
circunferencia NURBS), contrastes con fórmulas analíticas y con el teorema de Gauss-Bonnet, y una
comparación con la biblioteca NURBS-Python (geomdl). Los resultados coinciden con geomdl hasta
3·10⁻¹⁴ en puntos y derivadas de hasta segundo orden. Las 93 pruebas son reproducibles con un único
comando.

**Palabras clave:** NURBS, B-splines, álgebra computacional, Maxima, geometría diferencial,
curvatura gaussiana, software reproducible.

## 1. Introducción

Las NURBS combinan la flexibilidad local de las B-splines con la capacidad de representar
exactamente cónicas y cuádricas mediante pesos racionales [1, 2, 5]. Las bibliotecas abiertas más
difundidas, como openNURBS, SISL, verb y NURBS-Python/geomdl [3], están orientadas al cálculo
numérico. Para la docencia y la investigación en geometría computacional resulta valioso disponer
de una implementación dentro de un sistema de álgebra computacional, por tres razones:

1. **Verificación simbólica**: comprobar identidades de forma exacta y no solo hasta la tolerancia
   de máquina, por ejemplo que la circunferencia racional cuadrática tiene curvatura constante o
   que un parche reproduce una cuádrica.
2. **Derivación de fórmulas**: obtener expresiones cerradas de curvaturas en un parche como funciones
   racionales de los parámetros, útiles en análisis isogeométrico [6] y en la enseñanza.
3. **Transparencia didáctica**: los algoritmos de Piegl y Tiller [1] se pueden leer y modificar
   línea a línea en un lenguaje de alto nivel, sin compilación.

Maxima [7] es un sistema de álgebra computacional libre, descendiente de Macsyma, con una interfaz
de cuaderno (wxMaxima) muy extendida en la enseñanza universitaria. Hasta donde sabemos *[verificar
con una búsqueda bibliográfica]*, no existe un paquete NURBS completo para Maxima.

Las contribuciones de este trabajo son:

* una implementación en Maxima puro de los algoritmos fundamentales de NURBS, válida tanto
  en aritmética exacta como en coma flotante;
* un módulo de geometría diferencial de curvas y superficies NURBS que incluye clasificación de
  puntos, detección de singularidades e integrales globales;
* una metodología de validación reproducible con verificación simbólica, invariantes geométricos
  y comparación con una biblioteca de referencia independiente.

## 2. Fundamentos

### 2.1 Funciones base B-spline

Dado un vector de nudos no decreciente `U = {u_0, …, u_m}`, las funciones base de grado `p` se
definen mediante la recurrencia de Cox-de Boor [8, 9]:

```
N_{i,0}(u) = 1 si u_i ≤ u < u_{i+1}, y 0 en otro caso
N_{i,p}(u) = (u − u_i)/(u_{i+p} − u_i) · N_{i,p−1}(u) + (u_{i+p+1} − u)/(u_{i+p+1} − u_{i+1}) · N_{i+1,p−1}(u)
```

El paquete implementa la búsqueda del tramo (A2.1), la evaluación de las `p+1` funciones no nulas
(A2.2) y sus derivadas (A2.3) [1]. A2.2 solo realiza operaciones de cuerpo sobre `u` y los nudos,
así que funciona igual con `u` simbólico. Así se obtienen los polinomios exactos de cada tramo sin
código adicional.

### 2.2 Curvas y superficies NURBS

Una curva NURBS de grado `p` es

```
C(u) = Σ_i N_{i,p}(u) w_i P_i / Σ_i N_{i,p}(u) w_i ,
```

y una superficie es el producto tensorial análogo con pesos `w_{ij}`. La implementación trabaja en
coordenadas homogéneas `P_i^w = (w_i P_i, w_i)`: evalúa la B-spline polinómica en R⁴ y proyecta el
resultado. Las derivadas racionales se obtienen con la regla de Leibniz (A4.2 y A4.4 [1]):

```
C^(k) = ( A^(k) − Σ_{i=1..k} C(k,i) w^(i) C^(k−i) ) / w .
```

La invariancia proyectiva de las NURBS permite aplicar cualquier matriz homogénea 4×4 a los puntos
de control ponderados. Así, rotaciones, escalados y transformaciones perspectivas son exactos.

### 2.3 Geometría diferencial

Para una superficie regular `S(u,v)` con normal `n = S_u × S_v / |S_u × S_v|`, las formas
fundamentales son `E = S_u·S_u`, `F = S_u·S_v`, `G = S_v·S_v`, `L = S_uu·n`, `M = S_uv·n` y
`N = S_vv·n`, y

```
K = (LN − M²)/(EG − F²),     H = (EN − 2FM + GL)/(2(EG − F²)),     κ_{1,2} = H ± √(H² − K) [4].
```

Un punto es singular cuando `|S_u × S_v| = 0`, como ocurre en los polos de la esfera obtenida por
revolución. El paquete detecta estos puntos con una tolerancia relativa a `|S_u||S_v|`. Clasifica
los puntos regulares como elípticos, hiperbólicos, parabólicos, planos o umbilicales a partir de
`K` y de `H² − K`. En este último caso evita la raíz cuadrada, que amplificaría el error de redondeo.

Las integrales globales `∬ dA`, `∬ K dA` y `⅓∬ S·n dA` se calculan con cuadratura de
Gauss-Legendre de 5×5 puntos sobre una subdivisión de cada parche no vacío.

## 3. Diseño e implementación

### 3.1 Representación

Curvas y superficies se representan con estructuras de Maxima (`defstruct`):
`nurbscrv(degree, knots, ctrlpts, weights)` y
`nurbssrf(degree_u, degree_v, knots_u, knots_v, ctrlpts, weights)`. Los índices matemáticos
empiezan en 0, como en la literatura, y las listas de Maxima en 1. Los constructores validan el
vector de nudos (longitud, monotonía y multiplicidad) y la positividad de los pesos, e informan del
error con mensajes en español.

### 3.2 Aritmética exacta y en coma flotante

Ninguna rutina fuerza la conversión a coma flotante. Si el parámetro de evaluación es flotante,
el resultado se convierte al final con `float`. Si es exacto, el resultado conserva radicales y
racionales. La función `nurbs_float` convierte un objeto completo para los muestreos densos.
Las funciones `nurbs_curve_span_expr` y `nurbs_surface_span_expr` devuelven la expresión racional
cerrada de un tramo, a la que se pueden aplicar `diff`, `ratsimp` o `integrate`.

### 3.3 Organización por etapas

| Etapa | Módulo | Funciones principales |
|---|---|---|
| 1 | Núcleo | `nurbs_find_span`, `nurbs_basis_funs`, `nurbs_ders_basis_funs`, `bspline_basis`, `knots_*` |
| 2 | Curvas | `nurbs_curve`, `nurbs_curve_point`, `nurbs_curve_derivs`, `nurbs_transform`, `nurbs_curve_insert_knot`, `nurbs_interpolate_curve` |
| 3 | Superficies | `nurbs_surface`, `nurbs_surface_point`, `nurbs_surface_derivs`, `nurbs_surface_normal`, `nurbs_revolve`, `nurbs_surface_mesh` |
| 4 | Geometría diferencial | `nurbs_gaussian_curvature`, `nurbs_mean_curvature`, `nurbs_principal_*`, `nurbs_classify_point`, `nurbs_singular_points`, `nurbs_surface_area` |
| 5 | Visualización | `nurbs_plot`, `nurbs_draw_*`, `nurbs_export_{plot,obj,stl,vtk,json,curve_csv}` |

La visualización usa el paquete `draw` de Maxima (gnuplot) y dentro de wxMaxima dibuja en línea.
La exportación VTK incluye los campos escalares `K` y `H` para inspeccionar la curvatura en
ParaView. La exportación JSON sigue el esquema de geomdl, lo que permite intercambiar geometría
exacta entre ambos sistemas.

## 4. Metodología de validación

Las pruebas se organizan en cuatro niveles:

1. **Exactitud simbólica.** Se comprueban identidades en aritmética exacta: la expresión de los cuatro
   tramos de la circunferencia cumple `x² + y² − 1 = 0`; `C(1/8) = (1/√2, 1/√2)`; un parche de la
   esfera evaluado en parámetros racionales está sobre la esfera; las bases de Bézier coinciden con
   los polinomios de Bernstein; la curvatura gaussiana de `z = xy` es `−1/(1 + x² + y²)²`.
2. **Valores de libro.** Se reproduce el ejemplo 2.3 de [1] (tramo, bases `1/8, 3/4, 1/8` y sus
   derivadas).
3. **Invariantes geométricos.** Se verifican la partición de la unidad, la invariancia de la curva
   frente a la inserción de nudos y al escalado de pesos, la invariancia afín y proyectiva, el teorema
   de Gauss-Bonnet (`∬K dA = 4π` en la esfera y `0` en el toro), el área y el volumen del toro, la
   curvatura y la torsión de curvas planas y la ortogonalidad del triedro de Frenet y de las direcciones
   principales.
4. **Comparación con una implementación independiente.** El guion `tests/compare_geomdl.py`
   construye los mismos objetos en geomdl [3] y exporta puntos, derivadas, normales, interpolación
   e inserción de nudos a `tests/reference_geomdl.mac`. La batería de Maxima los compara con
   tolerancias de 10⁻⁸ a 10⁻¹².

## 5. Resultados

### 5.1 Precisión

**Tabla 1.** Error máximo absoluto en las pruebas representativas (de 93 pruebas, 93 correctas).

| Prueba | Error máximo |
|---|---|
| Circunferencia: `x² + y² = 1` (simbólico, 4 tramos) | 0 (exacto) |
| Curvatura de la circunferencia `r = 3` (14 puntos) | 1,7·10⁻¹⁶ |
| Longitud de la circunferencia frente a `2π` | 3,6·10⁻¹⁵ |
| Inserción de nudos (12 puntos) | 2,7·10⁻¹⁵ |
| Invariancia proyectiva | 4,4·10⁻¹⁶ |
| Toro: ecuación implícita (169 puntos) | < 10⁻¹³ |
| Toro: `K` frente a la expresión analítica (8 puntos) | 2,2·10⁻¹⁵ |
| Toro: `H` frente a la expresión analítica (8 puntos) | 2,4·10⁻¹⁵ |
| Esfera: `K = 1` | 8,9·10⁻¹⁶ |
| Gauss-Bonnet en la esfera (`4π`) | 3,7·10⁻¹² |
| Gauss-Bonnet en el toro (`0`) | 6,6·10⁻¹⁶ |
| Área del toro (`4π²Rr`) | 1,4·10⁻¹² |
| Volumen del toro (`2π²Rr²`) | 4,5·10⁻¹³ |

**Tabla 2.** Comparación con geomdl 5.5.1: error máximo en puntos y derivadas hasta segundo orden.

| Objeto | Magnitudes comparadas | Error máximo |
|---|---|---|
| Bézier cúbica plana | `C, C', C''` (9 parámetros) | 0 |
| Circunferencia racional | `C, C', C''` | 1,1·10⁻¹⁴ |
| B-spline cúbica 3D | `C, C', C''` | 7,1·10⁻¹⁵ |
| NURBS cúbica con pesos | `C, C', C''` | 1,4·10⁻¹⁴ |
| Bézier bicúbica | `S, S_u, S_v, S_uu, S_uv, S_vv` (6 puntos) | 4,4·10⁻¹⁶ |
| Toro racional | ídem | 3,2·10⁻¹⁴ |
| Esfera racional | ídem | 1,4·10⁻¹⁴ |
| Normales (toro, Bézier) | `n` | 2,5·10⁻¹⁶ |
| Interpolación global | nudos / puntos de control | 1,1·10⁻¹⁶ / 5,3·10⁻¹⁵ |
| Inserción de nudos | puntos homogéneos | 2,2·10⁻¹⁶ |

La única discrepancia apreciable aparece en la longitud de arco (1,75·10⁻³). geomdl la aproxima
con una poligonal de puntos evaluados. NURBS-Maxima usa cuadratura de Gauss-Legendre, y su valor
coincide hasta 1,1·10⁻¹² con `quad_qags` aplicado a la expresión exacta por tramos. La diferencia
se debe, por tanto, al método aproximado de la referencia.

### 5.2 Rendimiento

**Tabla 3.** Tiempo medio por llamada (Maxima 5.46/SBCL frente a geomdl 5.5.1/CPython 3.13,
mismo equipo, objetos en coma flotante).

| Operación | NURBS-Maxima | geomdl |
|---|---|---|
| Punto de la circunferencia | 0,25 ms | 0,007 ms |
| Punto del toro | 0,60 ms | 0,015 ms |
| Derivadas del toro (orden 2) | 2,8 ms | 0,10 ms |
| Curvatura gaussiana del toro | 2,9 ms | n/d |
| Malla 31×31 del toro con normales | 2,5 s | n/d |

NURBS-Maxima es entre 30 y 40 veces más lento que geomdl. La penalización es aceptable en
aplicaciones docentes y de verificación, y es el precio de la generalidad simbólica.

### 5.3 Ejemplos de uso

Los cuatro ejemplos de `examples/` ilustran:

* (01) que el peso central de una Bézier cuadrática racional determina la cónica (elipse, parábola
  o hipérbola), y la subdivisión de De Casteljau;
* (02) la circunferencia exacta, la no uniformidad de su parametrización (`|C'|` varía entre 5,66
  y 6,63) y sus transformaciones proyectivas;
* (03) las curvaturas de un parche bicúbico y las expresiones exactas para el paraboloide
  hiperbólico;
* (04) la distribución de curvatura del toro, con zonas elípticas (`K > 0`) en el exterior,
  hiperbólicas (`K < 0`) en el interior y parabólicas en las circunferencias superior e inferior,
  junto con la verificación de Gauss-Bonnet y la exportación VTK.

*[Insertar aquí las figuras generadas por los ejemplos: 01_bezier_curve.png, 03_bezier_surface.png
y 04_torus.png, más una figura de ParaView con el campo de K sobre 04_torus.vtk.]*

## 6. Discusión y limitaciones

La principal fortaleza del paquete es la coexistencia de evaluación numérica y simbólica en un
mismo código. Las limitaciones actuales son tres:

* **Rendimiento**: al ser interpretado, no es adecuado para mallas muy densas ni para optimización.
* **Crecimiento de expresiones**: las curvaturas simbólicas de superficies con pesos irracionales,
  como el toro o la esfera con `√2/2`, producen expresiones muy grandes. Se puede mitigar
  reparametrizando con pesos racionales o evaluando en puntos exactos.
* **Cobertura**: faltan la elevación y reducción de grado, la inserción de nudos en superficies,
  las superficies recortadas, las intersecciones y la aproximación por mínimos cuadrados.

## 7. Conclusiones

NURBS-Maxima pone a disposición de la comunidad de Maxima y wxMaxima un conjunto completo y validado
de herramientas para curvas y superficies NURBS y su geometría diferencial. Coincide con una
biblioteca de referencia hasta la precisión de máquina y permite, además, verificaciones exactas.
El trabajo futuro incluye las operaciones de la Sección 6, la compilación selectiva de rutinas
críticas con `translate`/`compile` y la integración con análisis isogeométrico.

## Disponibilidad y reproducibilidad

* Código: *[URL del repositorio]*, licencia MIT.
* Requisitos: Maxima ≥ 5.40, gnuplot (opcional) y Python 3 con geomdl (solo para regenerar las
  referencias).
* Reproducción: `python3 tests/compare_geomdl.py && maxima -b tests/test_nurbs.mac`.

## Referencias

*[Verificar todas las referencias antes del envío.]*

1. L. Piegl, W. Tiller. *The NURBS Book*, 2.ª ed. Springer, Berlín, 1997.
2. G. Farin. *Curves and Surfaces for CAGD: A Practical Guide*, 5.ª ed. Morgan Kaufmann, 2002.
3. O. R. Bingol, A. Krishnamurthy. NURBS-Python: An open-source object-oriented NURBS modeling
   framework in Python. *SoftwareX* 9 (2019) 85–94.
4. M. P. do Carmo. *Differential Geometry of Curves and Surfaces*. Prentice-Hall, 1976.
5. W. Tiller. Rational B-splines for curve and surface representation. *IEEE Computer Graphics and
   Applications* 3(6) (1983) 61–69.
6. T. J. R. Hughes, J. A. Cottrell, Y. Bazilevs. Isogeometric analysis: CAD, finite elements, NURBS,
   exact geometry and mesh refinement. *Computer Methods in Applied Mechanics and Engineering* 194
   (2005) 4135–4195.
7. Maxima, a Computer Algebra System. Versión 5.46.0, https://maxima.sourceforge.io/
8. M. G. Cox. The numerical evaluation of B-splines. *IMA Journal of Applied Mathematics* 10 (1972)
   134–149.
9. C. de Boor. On calculating with B-splines. *Journal of Approximation Theory* 6 (1972) 50–62.
