"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

x_data = np.loadtxt('datos.dat')
y_data = np.loadtxt('datosevol.dat')
x = x_data[1:-1,0]
n_y = len(y_data[:,0])

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0,100),ylim=(-1,1))
line, = ax.plot([],[],'-r',lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([],[])
    return line,

# animation function.  This is called sequentially
def animate(j):
    y = y_data[j,:]
    line.set_data(x,y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames = n_y ,interval=1,init_func=init,blit=True)
plt.xlabel("posicion en $x$")
plt.ylabel("posicion en $y$")
plt.show()

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
