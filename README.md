# APP_Lab3

## Introducción

En esta memoria se describe el proceso de paralelización del código `stencil.c`, tanto usando comunicaciones punto a punto como colectivas de vecindad para la comunicación de las fronteras entre procesos. 

En primer lugar se discutirá la estrategia general de paralelización, común a las versiones del código, profundizando en los detalles comunes a las dos versiones.

A continuación nos centraremos en los patrones de comunicación, y se explicarán, tanto detalles comunes relacionados con la comunicación, como peculiaridades e implementaciones concretas de las comunicaciones punto a punto y las colectivas de vecindad. 

Finalmente, se comparará el rendimiento de las dos versiones para comparar el rendimiento de los patrones de comunicación.

## Implementación paralela, detalles generales

### Distribución de procesos
-----

La estrategía seguida para paralelizar el programa `stencil.c` consiste en distribuir los procesos de forma análoga a los datos sobre los que vamos a trabajar. 

Para ello se ha usado un topología virtual en forma de malla en dos dimensiones, usando la función `MPI_Cart_create()`.

```c
#define DIMENSION 2
int periods[DIMENSIONS] = {0,0};
MPI_Comm MPI_COMM_MESH;
MPI_Cart_create( MPI_COMM_WORLD , DIMENSIONS, dims , periods , 1 , &MPI_COMM_MESH);
```

Nótese que esta función necesita como parámetro las dimensiones de la malla. Para que esto sea escalable y el programa funcione con cualquier número de procesos, dejamos que MPI decida esa distribución por nosotros usando la función `MPI_Dims_create()`, con todas las dimensiones iniciadas a 0, de forma que MPI puede elegir cómo distribuir los procesos sin restricciones.

```c
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

Otro punto que tratar es la gestión de las fuentes de calor "`sources`", que son los puntos en los que en cada iteración se genera calor. Estos puntos inicialmente tienen unos índices respecto  toda la matriz, pero ahora que cada proceso solo tiene un bloque de datos, puede que esas fuentes no se encuentren en los datos que el proceso controla, o en caso de que si se correspondan con estos datos, puede ser necesario modificar esos índices para que sean correctos dentro del bloque local del proceso. 

Para esto se ha implementado una función `reindex_source()` que cumple estas dos funciones, comprueba si el `source` pertenece a este proceso, y en caso afirmativo, recalibra sus índices para adaptarlo al bloque de datos local.

### Comunicación de los datos frontera
-----

El bucle principal permanece prácticamente sin cambios, la única modificación es que se ha añadido la comunicación de los datos frontera con otros procesos. Es decir, ahora al comienzo de cada iteración los procesos vecinos se comunican los datos necesarios para hacer todos los cálculos.

### Salida de los resultados
-----

El programa devuelve dos salidas.

Por un lado el calor total del sistema al final de las iteraciones y el tiempo de cálculo en segundos.

Por otro lado, un archivo bmp con la salida visual del mapa de calor.

Respecto al calor total y el tiempo de cálculo, esos datos se recolectan en un proceso usando `MPI_Reduce()`.

Respecto al mapa de calor, al acabar las iteraciones, los bloques finales se recolectan en el proceso principal, que es el encargado de generar ese archivo de salida. Esto se ha implementado de forma muy básica y poco eficiente, ya que no es el objetivo de la práctica.

## Detalles de la comunicación 

### Tipos de datos de las comunicaciones
-----

En este código stencil, un proceso necesita recibir su halo de los procesos vecinos. Es decir, recibirá de los vecinos de su misma fila las columnas laterales de su halo, y de los vecinos de su misma columna las filas superiores e inferiores de su halo.

Por este motivo, para representar esta información en las comunicaciones se han creado los tipos derivados `MPI_ROW` y `MPI_COLUMN`.

```c
// Create Datatypes
MPI_Type_contiguous( width, MPI_DOUBLE , &MPI_ROW);
MPI_Type_commit( &MPI_ROW);
MPI_Type_vector( height , 1 , width+2 , MPI_DOUBLE , &MPI_COLUMN);
MPI_Type_commit( &MPI_COLUMN);
```

### Comunicación punto a punto
-----

Para las comunicaciones punto a punto lo primero que necesitamos es identificar a nuestros vecinos inmediatos. Para ello, utilizamos `MPI_Cart_Shift()`.

```c
// Get neighbours
MPI_Cart_shift( MPI_COMM_MESH, 0, 1, &top , &bot);
MPI_Cart_shift( MPI_COMM_MESH, 1, 1, &left , &right);
```

Con esta información, procedemos a mover los datos a nuestros cuatro vecinos. En nuestra implementación hemos decidido seguir el siguiente orden: 

1. Izquierda
2. Derecha
3. Arriba
4. Abajo

Este movimiento de datos se ha implementado utilizando la función `MPI_Sendrecv`, de forma que, por ejemplo, en el caso en que movemos la información a la izquierda, un proceso simultáneamente envía a su vecino de la izquierda y recibe de su vecino de la derecha.

```c
// Send to the left
MPI_Sendrecv( buff+1+width+2 , 1 , MPI_COLUMN, left, 10, 
                buff+width+2+width+1 , 1 , MPI_COLUMN , right , 10 , 
                MPI_COMM_MESH , MPI_STATUS_IGNORE);
  
// Send to the right
MPI_Sendrecv( buff+width+2+width , 1 , MPI_COLUMN , right , 20 , 
                buff+width+2 , 1 , MPI_COLUMN , left , 20 , 
                MPI_COMM_MESH , MPI_STATUS_IGNORE);

// Send to the top
MPI_Sendrecv( buff+1+width+2, 1 , MPI_ROW , top , 30 , 
                buff+1+((width+2) * (height+1)) , 1 , MPI_ROW , bot , 30 , 
                MPI_COMM_MESH , MPI_STATUS_IGNORE);
  
//Send to the bottom
MPI_Sendrecv( buff+1+((width+2) * height) , 1 , MPI_ROW , bot , 40 , 
                buff+1 , 1 , MPI_ROW , top , 40 , 
                MPI_COMM_MESH , MPI_STATUS_IGNORE);
```


### Comunicación con colectivas de vecindad
-----

Como queremos mandar datos diferentes a cada vecino se ha decidido implementar esta comunicación usando `MPI_Neighbor_alltoall()`. En concreto, como queremos enviar y recibir filas o columnas, dependiendo del vecino con el que nos comunicamos, se ha decidido usar la versión `w` de esta colectiva.

```c
// Order: top -> bot -> left -> right 
MPI_Neighbor_alltoallw( buff , sendcounts , sdispls , sendtypes , buff , recvcounts , rdispls , recvtypes , MPI_COMM_MESH);
```

## Comprobación de resultados

Para asegurar que el programa paralelo funciona correctamente se han ejecutado ambas versiones con distintas configuraciones de procesos y se ha comparado que las salidas coinciden con las de la versión de referencia. En concreto se ha comprobado que:

 - El calor total coincide con el de la versión de referencia.
 - El fichero de salida `heat.svg` es exactamente igual que el generado por la versión de referencia. Para esto se ha utilizado el comando `diff`.

### Evaluación del rendimiento
-----

```sh
srun stencil 2048 1 1024
```

| .....        | 1 Process | 2 Processes | 4 Processes | 8 Processes | 16 Processes |
| :-----------:|----------:| -----------:| -----------:| -----------:| ------------:|
|Punto a punto | 54.846 s  | 28.145 s    | 14.885 s    | 8.367 s     | 4.512 s      |
|  Colectivas  | 54.315 s  | 27.899 s    | 14.406 s    | 8.190 s     | 4.441 s      |

El código escala de forma similar en ambas implementaciones, aunque el uso de las colectivas de vecindad si ayuda a obtener un rendimiento ligeramente superior para cualquier cantidad de procesos. 