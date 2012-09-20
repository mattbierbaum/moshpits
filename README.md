Simple nbody code to simulate a 2D set of particles
using a cell neighbor locator.  Has OpenGL display capabilities and is reasonably fast.

Options in the Makefile so that it works easily
with other systems (edit as necessary):
 - DOPLOT - display with opengl
 - FPS    - calculate frames per second (not portable)
 - POINTS - make the opengl display points, not fancy circles

To compile, simply do:
 - make

NOTE: helper scripts should not necessarily be run in one directory.  Have a look first, they are not complicated.

Check out the related project at <a href="http://github.com/mattbierbaum/moshpits.js">Moshpits.js</a>
