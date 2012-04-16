#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <cutil.h>
#include "plot.h"

#ifdef FPS
#include <time.h>
#endif

#define pi      3.141592653589
#define BLACK   0
#define RED     1
#define EPSILON DBL_EPSILON

int keys[256];
int plot_sizex;  
int plot_sizey;
int win;

typedef unsigned long long int ullong;
#define FPS
#define RIC   0
#define ff    8 

#define N        (int)(512*ff*ff)
#define L        20.0
#define radius   (0.2*2.5/ff)
#define Npercell  36

#define epsilon  46.0
#define speed    0.2

#define vhappy_black  0.0
#define vhappy_red    0.2
#define damp_coeff    1.0

#define dt  1e-1
#define R   2*radius 
#define R2  R*R
#define FR  2*R
#define FR2 FR*FR

#define TILE    16 
#define TILEX   TILE
#define TILEY   TILE
#define THREADS ((int)(N/(TILE*TILE)))

#define ERROR_CHECK { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
    printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__);}}

//-----------------------------------------------------------
// some defines and what not 
//------------------------------------------------------------
void   ran_seed(long j);
float  ran_ran2();
ullong vseed;
ullong vran;

void simulate(float a, float e, int s);
void init_circle(float *x, float *v, int *t);

float eps = EPSILON; 

__device__ float mymod(float a, float b){
  return a - b*(int)(a/b) + b*(a<0);
}

__device__ void coords_to_index(float *x, int *size, int *index){   
    index[0] = (int)(x[0]/L  * size[0]);
    index[1] = (int)(x[1]/L  * size[1]);
}

__device__ int mod_rvec(int a, int b, int *image){
    *image = 1;
    if (b==0) {if (a==0) *image=0; return 0;}
    if (a>b)  return a-b-1;
    if (a<0)  return a+b+1;
    *image = 0;
    return a;
}

void coords_to_index2(float *x, int *size, int *index){   
    index[0] = (int)(x[0]/L  * size[0]);
    index[1] = (int)(x[1]/L  * size[1]);
}

//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    float alpha_in = 0.1; 
    float eta_in   = 0.1;
    int seed_in     = 0;
    
    CUT_DEVICE_INIT(argc, argv);

    //if (argc == 1) 
        simulate(alpha_in, eta_in, seed_in);
    /*else if (argc == 4){
        alpha_in = atof(argv[1]);
        eta_in   = atof(argv[2]);
        seed_in  = atoi(argv[3]);
        simulate(alpha_in, eta_in, seed_in);
    }
    else {
        printf("usage:\n");
        printf("\t./entbody [alpha] [eta] [seed]\n");
    }*/
    return 0;
}


__global__ void step(float *x, float *v, int *type, 
                     unsigned int *cells, unsigned int *count, 
                     float *col, int *size, int size_total){
    //int bx = blockIdx.x;     
    //int by = blockIdx.y;     
    int tx = threadIdx.x;    

    //int dx = blockDim.x;
    //int dy = blockDim.y;
 
    int n = blockDim.x*blockIdx.x + tx;

    int index[2];
    if (n < size_total)
        count[n] = 0;
    if (n < size_total*Npercell)
        cells[n] = 0;

    __syncthreads();

    coords_to_index(&x[2*n], size, index);
    int t = index[0] + index[1]*size[0];
    unsigned int pos = atomicInc(&count[t], 0xffffffff);
    cells[Npercell*t + pos] = n;

    __syncthreads();

    int j,k;
    int tt[2];
    int tix[2];
    int image[2];
    float dx[2];

    float fx = 0.0;
    float fy = 0.0;
    float wx = 0.0;
    float wy = 0.0;
    
    for (tt[0]=-1; tt[0]<=1; tt[0]++){
    for (tt[1]=-1; tt[1]<=1; tt[1]++){
        tix[0] = mod_rvec(index[0]+tt[0],size[0]-1,&image[0]);
        tix[1] = mod_rvec(index[1]+tt[1],size[1]-1,&image[1]);

        int ind = tix[0] + tix[1]*size[0]; 

        for (j=0; j<count[ind]; j++){
            int tn = cells[Npercell*ind+j];

            float dist = 0.0;
            for (k=0; k<2; k++){
                dx[k] = x[2*tn+k] - x[2*n+k];
        
                if (image[k])
                    dx[k] += L*tt[k];
                dist += dx[k]*dx[k];
            }

            //===============================================
            // force calculation - hertz
            if (dist > 1e-6 && dist < R2){
                float r0 = R; 
                float l  = sqrt(dist);
                float co = epsilon * (1-l/r0)*(1-l/r0) * (l<r0);
                fx += - dx[k] * co;
                fy += - dx[k] * co;
                col[n] += co*co*dx[0]*dx[0]; 
                col[n] += co*co*dx[1]*dx[1]; 
            }
            //===============================================
            // add up the neighbor veocities
            if (dist > 1e-6 && dist < FR2 && type[n] == RED && type[tn] == RED){
                wx += v[2*n+0];
                wy += v[2*n+1];
            }                           
        }
    } } 

    //=====================================
    // flocking force 
    float wlen = wx*wx + wy*wy;
    if (type[n] == RED && wlen > 1e-6){
        fx += speed * wx / wlen; 
        fy += speed * wy / wlen;
    }

    //====================================
    // self-propulsion
    float vlen = v[2*n+0]*v[2*n+0] + v[2*n+1]*v[2*n+1];
    float vhappy = type[n]==RED?vhappy_red:vhappy_black;
    if (vlen > 1e-6){
        fx += damp_coeff*(vhappy - vlen)*v[2*n+0]/vlen;
        fy += damp_coeff*(vhappy - vlen)*v[2*n+1]/vlen;
    }
    
    // Newton-Stomer-Verlet
    v[2*n+0] += fx * dt;
    v[2*n+1] += fy * dt;

    x[2*n+0] += v[2*n+0] * dt;
    x[2*n+1] += v[2*n+1] * dt;
    
    // boundary conditions 
    if (x[2*n+0] >= L-EPSILON || x[2*n+0] < 0)
        x[2*n+0] = mymod(x[2*n+0], L);
    if (x[2*n+1] >= L-EPSILON || x[2*n+1] < 0)
        x[2*n+1] = mymod(x[2*n+1], L);

    col[n] = col[n]/4; 
}

//==================================================
// simulation
//==================================================
void simulate(float alpha, float eta, int seed){
    printf("Simulating %i particles\n", N);

    ran_seed(seed);
    int i;
    //int exit = 0;

    printf("Freeing local memory...\n");
    int *type   = (int*)malloc(sizeof(int)*N);
    float *rad = (float*)malloc(sizeof(float)*N); 
    float *col = (float*)malloc(sizeof(float)*N); 
    for (i=0; i<N; i++){ type[i] = 0; rad[i] = 0.0;}

    float *x = (float*)malloc(sizeof(float)*2*N);
    float *v = (float*)malloc(sizeof(float)*2*N);
    float *f = (float*)malloc(sizeof(float)*2*N);
    float *w = (float*)malloc(sizeof(float)*2*N);
    for (i=0; i<2*N; i++){x[i] = v[i] = f[i] = w[i] = 0.0;}

    float time_end = 1e2;

    #ifdef PLOT
    plot_init(); 
    plot_clear_screen();
    #endif

    //-------------------------------------------------
    // initialize
    printf("Initializing...\n");
    if (RIC){
        for (i=0; i<N; i++){
            float t = 2*pi*ran_ran2();
    
            rad[i] = radius;
            x[2*i+0] = L*ran_ran2();
            x[2*i+1] = L*ran_ran2();
     
            if (ran_ran2() > 0.3){
                v[2*i+0] = 0.0;
                v[2*i+1] = 0.0;
                type[i] = BLACK;
            }
            else {
                v[2*i+0] = vhappy_red * sin(t);
                v[2*i+1] = vhappy_red * cos(t);
                type[i] = RED;
            } 
        }
    }
    else {
        for (i=0; i<N; i++)
            rad[i] = radius;
        init_circle(x, v, type);
    }

    //-------------------------------------------------------
    // make boxes for the neighborlist
    int size[2];
    int size_total = 1;
    for (i=0; i<2; i++){
        size[i] = (int)(L / (FR)); 
        size_total *= size[i];
    }

    unsigned int *count  = (unsigned int*)malloc(sizeof(unsigned int)*size_total);
    unsigned int *cells  = (unsigned int*)malloc(sizeof(unsigned int*)*Npercell*size_total);
    unsigned int *count2  = (unsigned int*)malloc(sizeof(unsigned int)*size_total);
    unsigned int *cells2  = (unsigned int*)malloc(sizeof(unsigned int*)*Npercell*size_total);
    for (i=0; i<size_total; i++)
        count[i] = 0;
    for (i=0; i<size_total*Npercell; i++)
        cells[i] = 0;

    int index[2];
    for (i=0; i<N; i++){
        coords_to_index2(&x[2*i], size, index);
        int t = index[0] + index[1]*size[0];
        cells[Npercell*t + count[t]] = i;
        count[t]++; 
    }


    //==========================================================
    // where the magic happens
    //==========================================================
    printf("Freeing device memory...\n");
    int mem_size2 = sizeof(int)*2;
    int imem_size = sizeof(int)*N;
    int fmem_size = sizeof(float)*N;
    int fmem_siz2 = sizeof(float)*N*2;
    int mem_cell  = sizeof(unsigned int)*size_total;
    int mem_cell2 = sizeof(unsigned int)*size_total*Npercell;

    unsigned int *cu_count  = NULL;
    unsigned int *cu_cells  = NULL;    
    int *cu_size   = NULL;
    int *cu_type   = NULL; 

    float *cu_rad  = NULL;
    float *cu_col  = NULL;
    float *cu_x    = NULL;
    float *cu_v    = NULL; 
 
    cudaMalloc((void**) &cu_count, mem_cell);
    cudaMalloc((void**) &cu_cells, mem_cell2);
    cudaMalloc((void**) &cu_size,  mem_size2); 
 
    cudaMalloc((void**) &cu_type,  imem_size);
    cudaMalloc((void**) &cu_rad,   fmem_size);
    cudaMalloc((void**) &cu_col,   fmem_size);
    cudaMalloc((void**) &cu_x,     fmem_siz2);
    cudaMalloc((void**) &cu_v,     fmem_siz2);
    
    printf("Copying problem...\n");
    cudaMemcpy(cu_size,  size,  mem_size2, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_type,  type,  imem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_rad,   rad,   fmem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_col,   col,   fmem_size, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_x,     x,     fmem_siz2, cudaMemcpyHostToDevice);
    cudaMemcpy(cu_v,     v,     fmem_siz2, cudaMemcpyHostToDevice);
    
    cudaMemset(cu_count, 0, mem_cell);
    cudaMemset(cu_cells, 0, mem_cell2);
    ERROR_CHECK

    // ================================================
    // Initialize the block and grid dimensions here
    // ================================================
    #ifdef FPS
    struct timespec start;
    clock_gettime(CLOCK_REALTIME, &start);
    #endif

    int frames = 0;
    float t=0.0;

    dim3 grid(TILEX, TILEY, 1);
    dim3 thrd(THREADS, 1, 1);
    
    struct timespec startt;
    clock_gettime(CLOCK_REALTIME, &startt);

    int blocks = 256;
    printf("Simulating %ix%i = %i...\n", blocks, N/blocks);
    for (t=0.0; t<time_end; t+=dt){
 
        //step<<<grid, thrd>>>(cu_x, cu_v, cu_type, cu_cells, cu_count, cu_col, cu_size, size_total);
        step<<<blocks,N/blocks>>>(cu_x, cu_v, cu_type, cu_cells, cu_count, cu_col, cu_size, size_total);
        cudaThreadSynchronize();

        frames++;
        if (frames % 100 == 0){
            struct timespec end;
            clock_gettime(CLOCK_REALTIME, &end);
            printf("%i : %f\n", frames, (100)/((end.tv_sec - startt.tv_sec) + (end.tv_nsec - startt.tv_nsec)/1e9));
            ERROR_CHECK
            clock_gettime(CLOCK_REALTIME, &startt);
        }

        #ifdef PLOT
        plot_clear_screen();
        plot_render_particles(x,col,etc  );
        #endif
    }

    #ifdef FPS
    struct timespec end;
    clock_gettime(CLOCK_REALTIME, &end);
    printf("frames = %i\n", frames);
    printf("time   = %f\n",(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9 );
    printf("fps = %f\n", frames/((end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9));
    #endif

    cudaMemcpy(cells2, cu_cells, mem_cell2, cudaMemcpyDeviceToHost);
    cudaMemcpy(count2, cu_count, mem_cell, cudaMemcpyDeviceToHost);
    cudaMemcpy(type, cu_type, imem_size, cudaMemcpyDeviceToHost);

    int tcc = 0;
    int tcc2 = 0;
    for (i=0; i<size_total; i++){
        tcc += count[i];
        tcc2 += count2[i];
    }
    printf("-- %i %i --\n", tcc, tcc2);

    cudaFree(cu_count);
    cudaFree(cu_cells);
    cudaFree(cu_type);
    cudaFree(cu_rad);
    cudaFree(cu_col);
    cudaFree(cu_x);
    cudaFree(cu_v);
    ERROR_CHECK
  
    free(cells);
    free(count);
 
    free(x);
    free(v);
    free(f);
    free(w);
    free(rad);
    free(type);

    #ifdef PLOT
    plot_clean(); 
    #endif
}




//=================================================
// extra stuff
//=================================================
void ran_seed(long j){
  vseed = j;  vran = 4101842887655102017LL;
  vran ^= vseed; 
  vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
  vran = vran * 2685821657736338717LL;
}

float ran_ran2(){
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    ullong t = vran * 2685821657736338717LL;
    return 5.42101086242752217e-20*t;
}

void init_circle(float *x, float *v, int *type){
    int i;
    for (i=0; i<N; i++){
        float tx = L*ran_ran2();
        float ty = L*ran_ran2();
        float tt = 2*pi*ran_ran2();

        x[2*i+0] = tx;
        x[2*i+1] = ty;
        
        float dd = sqrt((tx-L/2)*(tx-L/2) + (ty-L/2)*(ty-L/2));

        // the radius for which 30% of the particles are red on avg
        float rad = sqrt(0.15*L*L / pi);
        if (dd < rad)
            type[i] = RED;

        if (type[i] == RED){
            v[2*i+0] = vhappy_red*cos(tt);
            v[2*i+1] = vhappy_red*sin(tt);
        }
        else {
            v[2*i+0] = 0.0;
            v[2*i+1] = 0.0;
        }
    }   
} 


#ifdef PLOT
//========================================================
// all of the plotting functionality
//========================================================
void key_down(unsigned char key, int x, int y){
    keys[key] = 1;
}
void key_up(unsigned char key, int x, int y){
    keys[key] = 0;
}

void plot_init(){
  plot_sizex = 640;
  plot_sizey = 640;
  win = 0;
  plot_init_opengl();
    int i;
    for (i=0; i<256; i++)
        keys[i] = 0;
}

void plot_clean(){
  plot_end_opengl();
}

//=============================================================
// OpenGL functionality
// http://www.andyofniall.net/2d-graphics-with-opengl/
//=============================================================
void plot_init_opengl(){
  int argc = 1;
  char *argv = (char*)malloc(sizeof(char)*42);
  sprintf(argv, "./entbody");

  glutInit(&argc, &argv);	         
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(plot_sizex, plot_sizey);
  glutInitWindowPosition(100,100);
  win = glutCreateWindow("EntBody Simulation");	

  glDisable(GL_DEPTH_TEST);
  glClearColor(1.0, 1.0, 1.0, 0.0);	/* set background to white */
  glutKeyboardFunc(key_down);
  glutKeyboardUpFunc(key_up);
  glViewport(0,0,plot_sizex, plot_sizey);

  glutMainLoopEvent();
  free(argv);
}

void plot_end_opengl(){
  glutDestroyWindow(win);
}

int plot_clear_screen(){
  glClear(GL_COLOR_BUFFER_BIT);
  return 1;
}


int plot_render_particles(double *x, double *rad, int *type, double *shade){
    // focus on the part of scene where we draw nice
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, L, L, 0, 0, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // lets draw our viewport just in case its not square
    glBegin(GL_LINE_LOOP);
      glVertex2f(0, 0);
      glVertex2f(0, L);
      glVertex2f(L, L);
      glVertex2f(L, 0);
    glEnd();
    
    glDisable(GL_POINT_SMOOTH);
    glPointSize(3);

    #ifdef POINTS 
    glBegin(GL_POINTS);
    #else
    double t=0;
    #endif

    int i; 
    for (i=0; i<N; i++){
        float tx = (float)x[2*i+0];
        float ty = (float)x[2*i+1];

        double c = fabs(shade[i]);
        if (c < 0) c = 0.0;
        if (c > 1.0) c = 1.0;

        float cr = c;
        float cg = c;
        float cb = c;
        float ca = 1.0;

        if (type[i] == 1) {
            cr = 0.9;//if (cr < 0.2) cr = 0.2;
            cg = 0.05;
            cb = 0.05;
        }
        
        #ifdef POINTS
        plot_set_draw_color(cr,cg,cb,ca);
        glVertex2f(tx, ty);
        #else
        double rx = rad[i];
        uint secs = 15;
        plot_set_draw_color(cr,cg,cb,ca);
        glBegin(GL_POLYGON);
        for (t=0; t<2*pi; t+=2*pi/secs)
          glVertex2f(tx + rx*cos(t), ty + rx*sin(t));
        glEnd();
        plot_set_draw_color(0.0,0.0,0.0,1.0);
        glBegin(GL_LINE_LOOP);
        for (t=0; t<2*pi; t+=2*pi/secs)
          glVertex2f(tx + rx*cos(t), ty + rx*sin(t));
        glEnd();
        #endif
    }
    #ifdef POINTS
    glEnd();
    #endif

    glutSwapBuffers();
    glutMainLoopEvent();

   return 0;
}

void plot_set_draw_color(float cr, float cg, float cb, float ca){
  glColor4f(cr, cg, cb, ca);
}
#endif
