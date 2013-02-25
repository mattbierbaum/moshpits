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
  plot_sizex = 720;
  plot_sizey = 720;
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

#ifdef OPENIL
ILuint img;
void plot_initialize_canvas(){
  ilInit();
  ilGenImages(1, &img);
  ilBindImage(img);
  ilTexImage(plot_sizex, plot_sizey, 0, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);
}

void plot_saveimage(const char* name){
  ILubyte *data = ilGetData();

  int viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  int x      = viewport[0];
  int y      = viewport[1];
  int width  = viewport[2];
  int height = viewport[3];

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);
  
  ilEnable(IL_FILE_OVERWRITE);
  ilSaveImage(name);
}
#endif

void draw_arrow(double px, double py, double vx, double vy, int inverse){
    GLfloat savedLineWidth = 1.0f;
    glGetFloatv(GL_LINE_WIDTH, &savedLineWidth);
    //glLineWidth(1.8*f*f*s);//sqrt(vx*vx+vy*vy)/20.);

    float t;
    float f = 1.6; //2.2
    float h = 0.7; //1.1
    float s = sqrt(vx*vx+vy*vy);
    double pt1x, pt1y, pt2x, pt2y, theta;
    pt1x = px-f*vx/2;
    pt1y = py-f*vy/2;
    pt2x = px+f*vx/2;
    pt2y = py+f*vy/2;
    theta = atan2(pt2y-pt1y, pt2x-pt1x);

    double perpx, perpy;
    perpx = -s/10.* f*vy/2;
    perpy = s/10.* f*vx/2;

    plot_set_draw_color(1.0,0.0,0.0,0.0);
    if (inverse)
        plot_set_draw_color(0.0,0.0,0.0,0.0);
    glBegin(GL_POLYGON);
    for (t=theta; t<theta+2*pi; t+=2*pi/3)
      glVertex2f(pt2x + h*cos(t)*s/f, pt2y + h*sin(t)*s/f);
    glEnd();

    plot_set_draw_color(1.0,0.0,0.0,0.0);
    if (inverse)
        plot_set_draw_color(0.0,0.0,0.0,0.0);
    glBegin(GL_POLYGON);
    glVertex2f(pt1x+perpx, pt1y+perpy);
    glVertex2f(pt1x-perpx, pt1y-perpy);
    glVertex2f(pt2x-perpx, pt2y-perpy);
    glVertex2f(pt2x+perpx, pt2y+perpy);
    glEnd();

    glLineWidth(savedLineWidth);
}

int *plot_render_particles(double *x, double *rad, int *type, long N, double L, double *shade, int forces,
                           double cmx, double cmy, int docom, int *pbc, double *v, int doarrows){
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

        if (forces){
            c = fabs(shade[i]);
            if (c < 0) c = 0.0; 
            if (c > 1.0) c = 1.0;
            cr = cg = cb = c;
            ca = 1.0;
            if (type[i] == 1) {
                cr = 0.9;
                cg = 0.05;
                cb = 0.05;
            }
        } else {
            cr = 1.0; 
            cg = 1.0;
            cb = 1.0;
            ca = 1.0;
            if (type[i] == 1) {
                cr = 0.00; 
                cg = 0.00;
                cb = 0.00;
            }
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

    if (doarrows)
        for (i=0; i<N; i++)
            if (type[i] == 1) draw_arrow(x[2*i+0], x[2*i+1], v[2*i+0], v[2*i+1], forces);
    #ifdef OPENMP 
    //#pragma omp barrier
    #endif

    if (docom == 1){
        double rx = 2;
        int secs = 15;
        double t;
        glColor4f(1.0,1.0,1.0,1.0);
        glBegin(GL_POLYGON);
        for (t=0; t<2*pi; t+=2*pi/secs)
          glVertex2f(cmx + rx*cos(t), cmy + rx*sin(t));
        glEnd();
        glColor4f(0.0,0.0,0.0,1.0);
        glBegin(GL_LINE_LOOP);
        for (t=0; t<2*pi; t+=2*pi/secs)
          glVertex2f(cmx + rx*cos(t), cmy + rx*sin(t));
        glEnd();

        for (i=0; i<N; i++){
            if (type[i] == 1){
                double tx = x[2*i+0] - cmx;
                double ty = x[2*i+1] - cmy;
                
                if (pbc[0] && tx > L/2)  tx -= L;
                if (pbc[1] && ty > L/2)  ty -= L;
                if (pbc[0] && tx < -L/2) tx += L;
                if (pbc[1] && ty < -L/2) ty += L;
                glBegin(GL_LINE_LOOP);
                glVertex2f(x[2*i+0], x[2*i+1]);
                glVertex2f(x[2*i+0]-tx, x[2*i+1]-ty);
                glEnd();
            }
        }
    }
 
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

