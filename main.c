#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef PLOT
#include "plot.h"
#endif

#ifdef FPS
#include <time.h>
#endif

//===========================================
// do we want to record various measurements?
//#define VORTICITY_TIMESERIES
//#define VELOCITY_DISTRIBUTION
//#define TEMPERATURE_BINS
//===========================================

//-----------------------------------------------------------
// some defines and helper functions for NBL
//------------------------------------------------------------
#define pi      3.141592653589
#define EPSILON DBL_EPSILON
#define BLACK   0
#define RED     1
#define RADS    10
#define BINS    50

void simulate(double a, double e, int s);

void   init_circle(double *x, double *v, int *t, double s, long N, double L);
void   temperature(double *x, double *v, int *t, int N, double L, int bins[RADS][BINS]);
double angularmom(double *x, double *v, int *t, int N);

void   coords_to_index(double *x, int *size, int *index, double L);
int    mod_rvec(int a, int b, int p, int *image);
double mymod(double a, double b);

void   ran_seed(long j);
double ran_ran2();
unsigned long long int vseed;
unsigned long long int vran;


//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    double alpha_in = 0.1; 
    double eta_in   = 0.1;
    int seed_in     = 0;

    if (argc == 1) 
        simulate(alpha_in, eta_in, seed_in);
    else if (argc == 4){
        alpha_in = atof(argv[1]);
        eta_in   = atof(argv[2]);
        seed_in  = atoi(argv[3]);
        simulate(alpha_in, eta_in, seed_in);
    }
    else {
        printf("usage:\n");
        printf("\t./entbody [alpha] [eta] [seed]\n");
    }
    return 0;
}



//==================================================
// simulation
//==================================================
void simulate(double alpha, double eta, int seed){
    ran_seed(seed);
    int  RIC  = 0;

    int    NMAX    = 50;
    int    N       = 1000;
    double radius  = 1.0;
    double L       = 1.05*sqrt(pi*radius*radius*N);

    int pbc[] = {1,1};

    double epsilon = 100.0;
    double T       = eta;
    double speed   = alpha;

    double vhappy_black = 0.0;
    double vhappy_red   = 1.0;
    double damp_coeff   = 1.0;
    double Tglobal      = 0.0;

    double dt  = 1e-1;
    double t   = 0.0;
    double R   = 2*radius; 
    double R2  = R*R;
    double FR  = 2*R;
    double FR2 = FR*FR;

    int i, j, k;

    int *type   = (int*)malloc(sizeof(int)*N);
    int *neigh  = (int*)malloc(sizeof(int)*N);
    double *rad = (double*)malloc(sizeof(double)*N); 
    double *col = (double*)malloc(sizeof(double)*N); 
    for (i=0; i<N; i++){ type[i] = neigh[i] = rad[i] = 0;}

    double *x = (double*)malloc(sizeof(double)*2*N);
    double *v = (double*)malloc(sizeof(double)*2*N);
    double *f = (double*)malloc(sizeof(double)*2*N);
    double *w = (double*)malloc(sizeof(double)*2*N);
    double *o = (double*)malloc(sizeof(double)*2*N);
    for (i=0; i<2*N; i++){o[i] = x[i] = v[i] = f[i] = w[i] = 0.0;}

    #ifdef PLOT 
    double time_end = 1e20;
    #else
    double time_end = 2e3;
    #endif

    #ifdef PLOT 
        int *key;
        double kickforce = 2.0;
        plot_init(); 
        plot_clear_screen();
        key = plot_render_particles(x, rad, type, N, L,col);
    #endif

    //-------------------------------------------------
    // initialize
    if (RIC){
        for (i=0; i<N; i++){
            double t = 2*pi*ran_ran2();
    
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
        init_circle(x, v, type, vhappy_red, N, L);
    }

    //-------------------------------------------------------
    // make boxes for the neighborlist
    int size[2];
    int size_total = 1;
    for (i=0; i<2; i++){
        size[i] = (int)(L / (FR)); 
        size_total *= size[i];
    }

    int *count = (int*)malloc(sizeof(int)*size_total);
    int *cells = (int*)malloc(sizeof(int)*size_total*NMAX);
    for (i=0; i<size_total; i++)
        count[i] = 0;
    for (i=0; i<size_total*NMAX; i++)
        cells[i] = 0;

    //==========================================================
    // where the magic happens
    //==========================================================
    int frames = 0;

    #ifdef FPS
    struct timespec start;
    clock_gettime(CLOCK_REALTIME, &start);
    #endif

    double angularmom_avg    = 0.0;
    double angularmom_std    = 0.0;
    double angularmom_sq_avg = 0.0;
    double angularmom_sq_std = 0.0;
    int angularmom_count = 0;

    double momentumx_avg = 0.0;
    double momentumx_std = 0.0;
    double momentumy_avg = 0.0;
    double momentumy_std = 0.0;
    int momentum_count = 0;
    

    #ifdef VORTICITY_TIMESERIES
    FILE *file1 = fopen("angularmom.txt", "wb");
    #endif

    #ifdef VELOCITY_DISTRIBUTION
    FILE *file2 = fopen("velocities.txt", "wb");
    #endif

    #ifdef TEMPERATURE_BINS 
    int bins[RADS][BINS];
    char name[80];
    sprintf(name, "temp_%0.2f.txt", damp_coeff);
    FILE *file3 = fopen(name, "w");
    for (i=0; i<RADS; i++){
        for (j=0; j<BINS; j++){
            bins[i][j] = 0;
        }
    }
    #endif

    int showplot = 1;
    for (t=0.0; t<time_end; t+=dt){

        int index[2];
        for (i=0; i<size_total; i++)
            count[i] = 0;

        for (i=0; i<N; i++){
            coords_to_index(&x[2*i], size, index, L);
            int t = index[0] + index[1]*size[0];
            cells[NMAX*t + count[t]] = i;
            count[t]++; 
        }

        int tt[2];
        int tix[2];
        int image[2];
        double dx[2];
        int goodcell, ind, n;
        double r0, l, co, dist;
        double wlen, vlen, vhappy;

        #ifdef OPENMP
        #pragma omp parallel for private(i,dx,index,tt,goodcell,tix,ind,j,n,image,k,dist,r0,l,co,wlen,vlen,vhappy)
        #endif 
        for (i=0; i<N; i++){
            f[2*i+0] = 0.0;
            f[2*i+1] = 0.0;
            w[2*i+0] = 0.0;
            w[2*i+1] = 0.0;
            neigh[i] = 0;
            
            coords_to_index(&x[2*i], size, index, L);

            for (tt[0]=-1; tt[0]<=1; tt[0]++){
            for (tt[1]=-1; tt[1]<=1; tt[1]++){
                goodcell = 1;    
                for (j=0; j<2; j++){
                    tix[j] = mod_rvec(index[j]+tt[j],size[j]-1,pbc[j],&image[j]);
                    if (pbc[j] < image[j])
                        goodcell=0;
                }

                if (goodcell){
                    ind = tix[0] + tix[1]*size[0]; 

                    for (j=0; j<count[ind]; j++){
                        n = cells[NMAX*ind+j];

                        dist = 0.0;
                        for (k=0; k<2; k++){
                            dx[k] = x[2*n+k] - x[2*i+k];
                    
                            if (image[k])
                                dx[k] += L*tt[k];
                            dist += dx[k]*dx[k];
                        }

                        //===============================================
                        // force calculation - hertz
                        if (dist > 1e-10 && dist < R2){
                            r0 = R; 
                            l  = sqrt(dist);
                            co = epsilon * (1-l/r0)*(1-l/r0) * (l<r0);
                            for (k=0; k<2; k++){
                                f[2*i+k] += - dx[k] * co;
                                col[i] += 0.0;//co*co*dx[k]*dx[k]; 
                            }
                        }
                        //===============================================
                        // add up the neighbor veocities
                        if (dist > 1e-10 && dist < FR2 && type[n] == RED && type[i] == RED){
                             for (k=0; k<2; k++)
                                w[2*i+k] += v[2*n+k];
                            neigh[i]++;
                        }                           
                    }
                }
            } } 

            //=====================================
            // flocking force 
            wlen = w[2*i+0]*w[2*i+0] + w[2*i+1]*w[2*i+1];
            if (type[i] == RED && neigh[i] > 0 && wlen > 1e-6){
                f[2*i+0] += speed * w[2*i+0] / wlen; 
                f[2*i+1] += speed * w[2*i+1] / wlen;
            }

            //====================================
            // self-propulsion
            vlen = v[2*i+0]*v[2*i+0] + v[2*i+1]*v[2*i+1];
            vhappy = type[i]==RED?vhappy_red:vhappy_black;
            if (vlen > 1e-6){
                f[2*i+0] += damp_coeff*(vhappy - vlen)*v[2*i+0]/vlen;
                f[2*i+1] += damp_coeff*(vhappy - vlen)*v[2*i+1]/vlen;
            }
        
            //=======================================
            // noise term
            if (type[i] == RED){
                double rt = 2*pi*ran_ran2();
                f[2*i+0] += T * cos(rt);
                f[2*i+1] += T * sin(rt);
            }

            f[2*i+0] += Tglobal*(ran_ran2()-0.5);
            f[2*i+1] += Tglobal*(ran_ran2()-0.5);

            //=====================================
            // kick force
            f[2*i+0] += o[2*i+0]; o[2*i+0] = 0.0;
            f[2*i+1] += o[2*i+1]; o[2*i+1] = 0.0;
        }
        #ifdef OPENMP
        #pragma omp barrier
        #endif

        // now integrate the forces since we have found them
        #ifdef OPENMP
        #pragma omp parallel for private(j)
        #endif 
        for (i=0; i<N;i++){
            // Newton-Stomer-Verlet
            #ifdef PLOT
            if (key['h'] != 1){
            #endif
            v[2*i+0] += f[2*i+0] * dt;
            v[2*i+1] += f[2*i+1] * dt;

            x[2*i+0] += v[2*i+0] * dt;
            x[2*i+1] += v[2*i+1] * dt;
            #ifdef PLOT
            }   
            #endif

            #ifdef VELOCITY_DISTRIBUTION
            double ttt = sqrt(v[2*i+0]*v[2*i+0] + v[2*i+1]*v[2*i+1]) ;
            fwrite(&ttt, sizeof(double), 1, file2);
            #endif

            // boundary conditions 
            for (j=0; j<2; j++){
                if (pbc[j] == 1){
                    if (x[2*i+j] >= L-EPSILON || x[2*i+j] < 0)
                        x[2*i+j] = mymod(x[2*i+j], L);
                }
                else {
                    const double restoration = 1.0;
                    if (x[2*i+j] >= L){x[2*i+j] = 2*L-x[2*i+j]; v[2*i+j] *= -restoration;}
                    if (x[2*i+j] < 0) {x[2*i+j] = -x[2*i+j];    v[2*i+j] *= -restoration;}
                    if (x[2*i+j] >= L-EPSILON || x[2*i+j] < 0){x[2*i+j] = mymod(x[2*i+j], L);}
                }
            }

            // just check for errors
            if (x[2*i+0] >= L || x[2*i+0] < 0.0 ||
                x[2*i+1] >= L || x[2*i+1] < 0.0)
                printf("out of bounds\n");
            
            col[i] = col[i]/4; 
        }
        #ifdef OPENMP
        #pragma omp barrier
        #endif

        #ifdef PLOT 
        if (frames % 100 == 0){
            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L,col);
        }
        #endif
        frames++;

        #ifdef TEMPERATURE_BINS
        temperature(x, v, type, N, L, bins);
        #endif

        angularmom_count++;
        
        double vtemp     = angularmom(x,v,type,N);
        double delta     = vtemp    - angularmom_avg;
        angularmom_avg    = angularmom_avg    + delta    / angularmom_count; 
        angularmom_std    = angularmom_std    + delta    * (vtemp    - angularmom_avg);
        
        double vtemp_sq  = vtemp*vtemp;
        double delta_sq  = vtemp_sq - angularmom_sq_avg;
        angularmom_sq_avg = angularmom_sq_avg + delta_sq / angularmom_count;
        angularmom_sq_std = angularmom_sq_std + delta_sq * (vtemp_sq - angularmom_sq_avg);


        double linearmomx = 0.0;
        double linearmomy = 0.0;
        int linearmomc = 0;
        for (i=0; i<N; i++){
            if (type[i] == RED){
                linearmomx += v[2*i+0];
                linearmomy += v[2*i+1];
                linearmomc++;
            }
        }
        linearmomx /= linearmomc;
        linearmomy /= linearmomc;
        momentum_count++;

        double deltax = linearmomx - momentumx_avg;
        double deltay = linearmomy - momentumy_avg;
        momentumx_avg = momentumx_avg + deltax / momentum_count;
        momentumy_avg = momentumy_avg + deltay / momentum_count;
        momentumx_std = momentumx_std + deltax * (linearmomx - momentumx_avg);
        momentumy_std = momentumy_std + deltay * (linearmomy - momentumy_avg);

        #ifdef VORTICITY_TIMESERIES
        fwrite(&vtemp, sizeof(double), 1, file1);
        #endif

        #ifdef PLOT
        if (key['f'] == 1)
            showplot = !showplot;
        if (key['k'] == 1)
            vhappy_red = 0.0;
        if (key['q'] == 1)
            break;
        if (key['w'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+1] = -kickforce;
            }
        }
        if (key['s'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+1] = kickforce;
            }
        }
        if (key['a'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+0] = -kickforce;
            }
        }
        if (key['d'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+0] = kickforce;
            }
        }
        if (key['9'] == 1)
            Tglobal -= 0.01;
        if (key['0'] == 1)
            Tglobal += 0.01;
        if (key['8'] == 1)
            Tglobal = 0.0;
        #endif
    }
    // end of the magic, cleanup
    //----------------------------------------------
    #ifdef FPS
    struct timespec end;
    clock_gettime(CLOCK_REALTIME, &end);
    printf("fps = %f\n", frames/((end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9));
    #endif

    #ifdef VORTICITY_TIMESERIES
    fclose(file1);
    #endif

    #ifdef VELOCITY_DISTRIBUTION
    fclose(file2);
    #endif

    #ifdef TEMPERATURE_BINS
    for (i=0; i<RADS; i++){
        for (j=0; j<BINS; j++){
            fprintf(file3, "%i ", bins[i][j]);
        }
        fprintf(file3, "\n");
    }
    fclose(file3);
    #endif

    printf("tend = %f\n", t);
    angularmom_std    = angularmom_std    / (angularmom_count - 1);
    angularmom_sq_std = angularmom_sq_std / (angularmom_count - 1);
 
    momentumx_std = momentumx_std / (momentum_count - 1);
    momentumy_std = momentumy_std / (momentum_count - 1);
 
    printf("%f %f %f %f %f %f %f %f\n", angularmom_avg, sqrt(angularmom_std), 
                                        angularmom_sq_avg, sqrt(angularmom_sq_std), 
                                        momentumx_avg, sqrt(momentumx_std), 
                                        momentumy_avg, sqrt(momentumy_std));

    free(cells);
    free(count);
 
    free(x);
    free(v);
    free(f);
    free(w);
    free(o);
    free(neigh);
    free(rad);
    free(type);
    free(col);

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

double ran_ran2(){
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    unsigned long long int t = vran * 2685821657736338717LL;
    return 5.42101086242752217e-20*t;
}

void init_circle(double *x, double *v, 
                 int *type, double speed, long N, double L){
    int i;
    for (i=0; i<N; i++){
        double tx = L*ran_ran2();
        double ty = L*ran_ran2();
        double tt = 2*pi*ran_ran2();

        x[2*i+0] = tx;
        x[2*i+1] = ty;
        
        // the radius for which 30% of the particles are red on avg
        double dd = sqrt((tx-L/2)*(tx-L/2) + (ty-L/2)*(ty-L/2));
        double rad = sqrt(0.15*L*L / pi);

        //if (i<0.15*N)
        if (dd < rad)
            type[i] = RED;
        else
            type[i] = BLACK;

        if (type[i] == RED){
            v[2*i+0] = speed*cos(tt);
            v[2*i+1] = speed*sin(tt);
        }
        else {
            v[2*i+0] = 0.0;
            v[2*i+1] = 0.0;
        }
    }   
} 

//=======================================
// NBL - neighborlist helper functions
//=======================================
inline double mymod(double a, double b){
  return a - b*(int)(a/b) + b*(a<0);
}

inline void coords_to_index(double *x, int *size, int *index, double L){   
    index[0] = (int)(x[0]/L  * size[0]);
    index[1] = (int)(x[1]/L  * size[1]);
}

inline int mod_rvec(int a, int b, int p, int *image){
    *image = 1;
    if (b==0) {if (a==0) *image=0; return 0;}
    if (p != 0){
        if (a>b)  return a-b-1;
        if (a<0)  return a+b+1;
    } else {
        if (a>b)  return b;
        if (a<0)  return 0;
    }
    *image = 0;
    return a;
}



//==========================================
// measurement functions
//=========================================
double angularmom(double *x, double *v, int *t, int N){
    int i=0;
    double ang = 0.0;
    double cmx = 0.0;
    double cmy = 0.0;
    int count = 0;

    for (i=0; i<N; i++){
        if (t[i] == RED){
            cmx += x[2*i+0];
            cmy += x[2*i+1];    
            count++;
        }
    }
    cmx /= count;
    cmy /= count;

    for (i=0; i<N; i++){
        if (t[i] == RED){
            double tx = x[2*i+0] - cmx;
            double ty = x[2*i+1] - cmy;
            double vx = v[2*i+0];
            double vy = v[2*i+1];
            double tv = vx*ty - vy*tx;
            ang += tv;
        }
    }

    return ang/count;
}


void temperature(double *x, double *v, int *t, int N, double L, int bins[RADS][BINS]){
    int i=0;
    double cmx = 0.0;
    double cmy = 0.0;
    int count = 0;

    for (i=0; i<N; i++){
        if (t[i] == RED){
            cmx += x[2*i+0];
            cmy += x[2*i+1];    
            count++;
        }
    }
    cmx /= count;
    cmy /= count;

    for (i=0; i<N; i++){
        if (t[i] == RED){
            double dx = x[2*i+0] - cmx;
            double dy = x[2*i+1] - cmy;
            double rr = sqrt(dx*dx + dy*dy);
            double vv = sqrt(v[2*i+0]*v[2*i+0] + v[2*i+1]*v[2*i+1]);

            int rad = RADS * 2*rr/L;
            int bin = BINS * vv/3;
            if (rad < RADS && bin < BINS){
                bins[rad][bin]++;
            }
        }
    }
}

