#include <stdio.h>
#include <math.h>
#include "plot.h"

#define pi 3.14159265358

int keys[256];
int plot_sizex;  
int plot_sizey;
int win;

void key_down(unsigned char key, int x, int y){
  keys[key] = 1;
}

void key_up(unsigned char key, int x, int y){
  keys[key] = 0;
}

void plot_init(){
  plot_sizex = 680;
  plot_sizey = 680;
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


int *plot_render_particles(double *x, double *rad, int *type, long N, double L, double *shade){
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
    glPointSize(1);

    #ifdef POINTS 
    glBegin(GL_POINTS);
    #else
    double t=0;
    #endif

    int i;
    float tx, ty, cr, cg, cb, ca;
    double c, rx;
    uint secs;

    #ifdef OPENMP
    //#pragma omp parallel for private(tx,ty, c, cr, cg, cb, ca, rx, t, secs)
    #endif
    for (i=0; i<N; i++){
        tx = (float)x[2*i+0];
        ty = (float)x[2*i+1];

        c = fabs(shade[i]);
        if (c < 0) c = 0.0;
        if (c > 1.0) c = 1.0;

        cr = c;
        cg = c;
        cb = c;
        ca = 1.0;

        if (type[i] == 1) {
            cr = 0.9;//if (cr < 0.2) cr = 0.2;
            cg = 0.05;
            cb = 0.05;
        }
        
        #ifdef POINTS
        plot_set_draw_color(cr,cg,cb,ca);
        glVertex2f(tx, ty);
        #else
        rx = rad[i];
        secs = 15;
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
    #ifdef OPENMP 
    //#pragma omp barrier
    #endif

    #ifdef POINTS
    glEnd();
    #endif

    glutSwapBuffers();
    glutMainLoopEvent();

    return keys;
}

void plot_set_draw_color(float cr, float cg, float cb, float ca){
  glColor4f(cr, cg, cb, ca);
}

