# APP_Lab3

## Introducción

En esta memoria se describe el proceso de paralelización del código `stencil.c`, tanto usando comunicaciones punto a punto como collectivas de vecindad para la comunicación de las fronteras entre procesos. 

En primer lugar se discutirá la estrategia general de paralelización, común a las versiones del código, profundizando en los detalles comunes a las dos versiones.

A continuación nos centraremos en los patrones de comunicación, y se explicarán, tanto detalles comunes relacionados con la comunicación, como peculiaridades e implementaciones concretas de las comunicaciones punto a punto y las colectivas de vecindad. 

Finalmente, se comparará el rendimiento de las dos versiones para comparar el rendimiento de los patrones de comunicación.

## Implementación paralela, detalles generales

### Distribución de procesos
-----

La estrategía seguida para paralelizar el programa `stencil.c` consiste en distribuir los procesos de forma análoga a los datos sobre los que vamos a trabajar. 

Para ello se ha usado un topologia virtual en forma de malla en dos dimensiones, usando la función `MPI_Cart_create()`.

```c++
#define DIMENSION 2
int periods[DIMENSIONS] = {0,0};
MPI_Comm MPI_COMM_MESH;
MPI_Cart_create( MPI_COMM_WORLD , DIMENSIONS, dims , periods , 1 , &MPI_COMM_MESH);
```

Nótese que esta función necesita como parámetro las dimensiones de la malla. Para que esto sea escalable y el programa funcione con cualquier número de procesos, dejamos que MPI decida esa distribución por nosotros usando la función `MPI_Dims_create()`, con todas las dimensiones iniciadas a 0, de forma que MPI puede elegir como distribuir los procesos sin restricciones.

```c++
int dims[DIMENSIONS] = {0,0};
MPI_Dims_create( numprocs , DIMENSIONS , dims);
```

### Distribución de datos
-----

Igual que los procesos, los datos se han distribuido entre los procesos en bloques bidimensionales. Es decir, si antes todos los datos estaban en la memoria del mismo proceso, ahora cada proceso solo reserva memoria para un bloque, y no para todos los datos. El tamaño de esos bloques depende de las dimensiones obtenidas con `MPI_Dims_create()`.

![Distribución bidimensional de datos sobre procesos][dataDist]

[dataDist]: https://i.stack.imgur.com/2JO0M.jpg "Data Distribution over processes"

### Reindexación de las fuentes de calor
-----

Otro punto que tratar es la gestión de las fuentes de calor "`sources`", que son los puntos en los que en cada iteracion se genera calor. Estos puntos inicialmente tienen unos índices respecto  toda la matriz, pero ahora que cada proceso solo tiene un bloque de datos, puede que esas fuentes no se encuentren en los datos que el proceso controla, o en caso de que si se correspondan con estos datos, puede ser necesario modificar esos indices para que sean correctos dentro del bloque local del proceso. 

Para esto se ha implementado una función `reindex_source()` que cumple estas dos funciones, comprueba si el `source` pertenece a este proceso, y en caso afirmativo, recalibra sus indices para adaptarlo al bloque de datos local.

### Comunicación de los datos frontera
-----

El bucle principal permanece práctimente sin cambios, la única modificación es que se ha añadido la comunicación de los datos frontera con otros procesos. Es decir, ahora al comienzo de cada iteración los procesos vecinos se comunican los datos necesarios para hacer todos los cálculos.

### Salida de los resultados
-----

El programa devuelve dos salidas.

Por un lado el calor total del sistema al final de las iteraciones y el tiempo de cálculo en segundos.

Por otro lado, un archivo bmp con la salida visual del mapa de calor.

Respecto al calor total y el tiempo de cálculo, esos datos se recolectan en un proceso usando `MPI_Reduce()`.

Respecto al mapa de calor, al acabar las iteraciones, los bloques finales se recolectan en el proceso principal, que es el encargado de generar ese archivo de salida. Esto se ha implementado de forma muy básica y poco eficiente, ya que no es el objetivo de la práctica.

## Detalles de la comunicación 

### Comunicación punto a punto
-----

### Comunicación con colectivas de vecindad
-----