#include <stdio.h>
#include <math.h>
#include "plot.h"

#define pi 3.14159265358

int keys[256];
int plot_sizex;  
int plot_sizey;
int win;

void exit_key(unsigned char key, int x, int y){
  if (key == 'q')
    commence_exit = 1;
}


void plot_init(){
  plot_sizex = 640;
  plot_sizey = 640;
  win = 0;
  plot_init_opengl();
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
  glutKeyboardFunc(exit_key);
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


int plot_render_particles(double *x, double *rad, int *type, long N, double L, double *shade){
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

    if (commence_exit == 1){
      commence_exit = 0;
      return 1;
    }
    return 0;
}

int plot_exit_func(){
    glutMainLoopEvent();
    if (commence_exit == 1){
        commence_exit = 0;
        return 1;
    }
    return 0;
}

void plot_set_draw_color(float cr, float cg, float cb, float ca){
  glColor4f(cr, cg, cb, ca);
}

