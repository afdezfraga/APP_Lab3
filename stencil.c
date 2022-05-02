// Compile with mpicc -o stencil stencil.c -lm
// Example for run: mpirun -np 1 ./stencil 200 1 500
// Output: heat.svg

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// row-major order
#define ind(i,j) (j)*(width+2)+i

#define DIMENSIONS 2

void
communicate_borders(int top, int bot, int left, int right, 
                  int width, int heigth, double* buff, 
                  MPI_Datatype MPI_ROW, MPI_Datatype MPI_COLUMN,
                  MPI_Comm MPI_COMM_MESH, int* coords);

void
reindex_source(int* sources, int* my_sources, 
              int n, int* coords, int* dims);

// Since this is not important it has been done the easy way
void 
gather_final_result(double* main_buff, double* local_buff, 
                    int rank, int numprocs, int width, int height, int n, 
                    int procs_per_line, MPI_Comm MPI_COMM_MESH);

void printarr(double *a, int n) {
  // does nothing right now, should record each "frame" as image
  FILE *fp = fopen("heat.svg", "w");
  const int size = 5;

  int width = n;

  fprintf(fp, "<html>\n<body>\n<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">");

  fprintf(fp, "\n<rect x=\"0\" y=\"0\" width=\"%i\" height=\"%i\" style=\"stroke-width:1;fill:rgb(0,0,0);stroke:rgb(0,0,0)\"/>", size*n, size*n);
  for(int i=1; i<n+1; ++i)
    for(int j=1; j<n+1; ++j) {
      int rgb = (a[ind(i,j)] > 0) ? rgb = (int)round(255.0*a[ind(i,j)]) : 0.0;
      if(rgb>255) rgb=255;
      if(rgb) fprintf(fp, "\n<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"stroke-width:1;fill:rgb(%i,0,0);stroke:rgb(%i,0,0)\"/>", size*(i-1), size*(j-1), size, size, rgb, rgb);
    }
  fprintf(fp, "</svg>\n</body>\n</html>");


  fclose(fp);
}

int main(int argc, char **argv) {

  int n = atoi(argv[1]); // nxn grid
  int energy = atoi(argv[2]); // energy to be injected per iteration
  int niters = atoi(argv[3]); // number of iterations
  double *tmp;

  int rank, cart_rank, numprocs;
  int dims[DIMENSIONS] = {0,0}, periods[DIMENSIONS] = {0,0}, coords[DIMENSIONS];
  int top, bot, right, left;
  int height, width;

  MPI_Comm MPI_COMM_MESH;
  MPI_Datatype MPI_ROW, MPI_COLUMN;

  //MPI_Init(NULL, NULL);
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);	
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Init the mesh 
  MPI_Dims_create( numprocs , DIMENSIONS , dims);
  MPI_Cart_create( MPI_COMM_WORLD , DIMENSIONS, dims , periods , 1 , &MPI_COMM_MESH);
  MPI_Comm_rank(MPI_COMM_MESH, &cart_rank);
  MPI_Cart_coords( MPI_COMM_MESH , cart_rank , DIMENSIONS , coords);

  // Get my size
  width = (coords[1] < n % dims[1]) ? (n / dims[1]) + 1 : n / dims[1];
  height = (coords[0] < n % dims[0]) ? (n / dims[0]) + 1 : n / dims[0];

  // Init memory
  double *aold = (double*)calloc(1,(width+2)*(height+2)*sizeof(double)); // 1-wide halo zones!
  double *anew = (double*)calloc(1,(width+2)*(height+2)*sizeof(double)); // 1-wide halo-zones!
  double *afull;
  if(cart_rank == 0)
    afull = (double*)calloc(1,(n+2)*(n+2)*sizeof(double)); // 1-wide halo-zones!

  // Get neighbours
  MPI_Cart_shift( MPI_COMM_MESH, 0, 1, &top , &bot);
  MPI_Cart_shift( MPI_COMM_MESH, 1, 1, &left , &right);

  // Create Datatypes
  MPI_Type_contiguous( width, MPI_DOUBLE , &MPI_ROW);
  MPI_Type_commit( &MPI_ROW);
  MPI_Type_vector( height , 1 , width+2 , MPI_DOUBLE , &MPI_COLUMN);
  MPI_Type_commit( &MPI_COLUMN);

  #define nsources 3
  int sources[nsources][2] = {{n/2,n/2}, {n/3,n/3}, {n*4/5,n*8/9}};

  int my_sources[nsources][2] = {{-1,-1}, {-1,-1}, {-1,-1}};

  for (int i = 0; i < nsources; i++)
  {
    reindex_source(sources[i], my_sources[i], n, coords, dims);
  }
  
  double heat=0.0, last_heat, max_t; // total heat in system
  double t=-MPI_Wtime();
  for(int iter=0; iter<niters; ++iter) {
    heat = 0.0;
    communicate_borders(top, bot, left, right, width, height, aold, MPI_ROW, MPI_COLUMN, MPI_COMM_MESH, coords);
    for(int j=1; j<height+1; ++j) {
      for(int i=1; i<width+1; ++i) {
        anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
        heat += anew[ind(i,j)];
      }
    }
    for(int i=0; i<nsources; ++i) {
      if(my_sources[i][0] != -1)
        anew[ind(my_sources[i][0], my_sources[i][1])] += energy; // heat source
    }
    tmp=anew; anew=aold; aold=tmp; // swap arrays
  }
  t+=MPI_Wtime();
  MPI_Reduce( &heat , &last_heat , 1 , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_MESH);
  MPI_Reduce( &t, &max_t , 1 , MPI_DOUBLE , MPI_MAX , 0 , MPI_COMM_MESH);
  gather_final_result(afull, anew, 
                    cart_rank, numprocs, width, height, n, 
                    dims[1], MPI_COMM_MESH);
  if (cart_rank == 0){
    printarr(afull, n);
    printf("last heat: %f time: %f\n", last_heat, max_t);
  }
  MPI_Finalize();
}


// Since this is not important it has been done the easy way
void 
gather_final_result(double* main_buff, double* local_buff, 
                    int rank, int numprocs, int width, int height, int n, 
                    int procs_per_line, MPI_Comm MPI_COMM_MESH) {

  int widths[numprocs];
  int heights[numprocs];

  MPI_Gather( &width , 1, MPI_INT , widths , 1 , MPI_INT , 0 , MPI_COMM_MESH);
  MPI_Gather( &height , 1, MPI_INT , heights , 1 , MPI_INT , 0 , MPI_COMM_MESH);

  int j_start = 1;
  int i_start = 1;
  
  if (rank == 0){
    // Write tile [0,0] to main buffer
    for (int j = j_start; j < height+1; j++)
    {
      for (int i = j_start; i < width+1; i++)
      {
        main_buff[j*(n+2)+i] = local_buff[j*(width+2)+i];
      }
      
    }

    //Update where next local_buffer will be written in main_buffer
    if(1 % procs_per_line == 0){
      i_start=1;
      j_start=j_start+height;
    } else {
      i_start += width;
    }
    
    // Receive tiles from other procs and write them to main buffer
    for (int p = 1; p < numprocs; p++)
    {
      MPI_Recv( local_buff , heights[p]*(widths[p]+2) , MPI_DOUBLE , p , MPI_ANY_TAG , MPI_COMM_MESH , MPI_STATUS_IGNORE);
      for (int j = 0; j < heights[p]; j++)
      {
        for (int i = 0; i < widths[p]; i++)
        {
          main_buff[(j_start+j)*(n+2)+(i_start+i)] = local_buff[j*(widths[p]+2)+i+1];
        }
        
      }
      //Update where next local_buffer will be written in main_buffer
      if((p+1) % procs_per_line == 0){
        i_start=1;
        j_start=j_start+heights[p];
      } else {
        i_start += widths[p];
      }
      
    }
  } else {
    // Send to P0
    MPI_Send( local_buff+width+2 , height*(width+2) , MPI_DOUBLE , 0 , rank , MPI_COMM_MESH);
  }
    

}


void
communicate_borders(int top, int bot, int left, int right, 
                  int width, int heigth, double* buff, 
                  MPI_Datatype MPI_ROW, MPI_Datatype MPI_COLUMN,
                  MPI_Comm MPI_COMM_MESH, int* coords){
             
  // Send to the left
  MPI_Sendrecv( buff+1+width+2 , 1 , MPI_COLUMN, left, 10, buff+width+2+width+1 , 1 , MPI_COLUMN , right , 10 , MPI_COMM_MESH , MPI_STATUS_IGNORE);
  
  // Send to the right
  MPI_Sendrecv( buff+width+2+width , 1 , MPI_COLUMN , right , 20 , buff+width+2 , 1 , MPI_COLUMN , left , 20 , MPI_COMM_MESH , MPI_STATUS_IGNORE);

  // Send to the top
  MPI_Sendrecv( buff+1+width+2, 1 , MPI_ROW , top , 30 , buff+1+((width+2) * (heigth+1)) , 1 , MPI_ROW , bot , 30 , MPI_COMM_MESH , MPI_STATUS_IGNORE);
  
  //Send to the bottom
  MPI_Sendrecv( buff+1+((width+2) * heigth) , 1 , MPI_ROW , bot , 40 , buff+1 , 1 , MPI_ROW , top , 40 , MPI_COMM_MESH , MPI_STATUS_IGNORE);
}

void
reindex_source(int* sources, int* my_sources, 
              int n, int* coords, int* dims){

  int row_pos = sources[1] - 1;
  int col_pos = sources[0] - 1;

  int row_start, row_end, col_start, col_end;

  int row_local = n / dims[0];
  int row_resto = n % dims[0];

  int col_local = n / dims[1];
  int col_resto = n % dims[1];

  if (coords[0] < row_resto){
    row_start = coords[0] * (row_local + 1); 
    row_end = (coords[0] + 1) * (row_local + 1);
  } else {
    row_start = (coords[0] * row_local) + row_resto;
    row_end = ((coords[0] + 1) * row_local) + row_resto;
  }

    if (coords[1] < col_resto){
    col_start = coords[1] * (col_local + 1); 
    col_end = (coords[1] + 1) * (col_local + 1);
  } else {
    col_start = (coords[1] * col_local) + col_resto;
    col_end = ((coords[1] + 1) * col_local) + col_resto;
  }

  if(row_pos >= row_start && row_pos < row_end && col_pos >= col_start && col_pos < col_end){
    my_sources[0] = sources[0] - col_start; 
    my_sources[1] = sources[1] - row_start;
  }
}