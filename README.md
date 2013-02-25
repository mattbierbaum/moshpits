Collective motion at heavy metal concerts
=========================================

A simple n-body code to simulate a 2D set of particles
using a cell neighbor locator.  
Uses OpenGL for display capabilities and is reasonably fast.

There are options in the Makefile so that it works easily
with other systems (edit as necessary):
 - DOPLOT - display with opengl
 - FPS    - calculate frames per second (not portable, POSIX only)
 - POINTS - make the opengl display points, not fancy circles (a speed issue)
 - IMAGES - use OpenIL to save screen shots to disk

To compile, simply `make`.

There are several dependencies required to use all features:
 - freeglut - used for simple OpenGL bindings.  This is different than regular glut and not compatible.
 - OpenIL - open image library used to save screenshots to various image formats.

To install these on a debian system, use the command:

    sudo apt-get install freeglut3-dev libdevil-dev    

There are several helper scripts which generate several of the plots in the paper. 
They are written in Python and located in the scripts directory.  Each one
requires different files to be present in order to run.  The generation of these
are controlled by the top `#define` statements in `main.c`.

Check out the related project at <a href="http://github.com/mattbierbaum/moshpits.js">Moshpits.js</a>
