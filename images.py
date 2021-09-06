import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import HTML
from IPython.display import display


# for step in range(100):
#     df = pd.read_csv(f'data/data_{step}.csv', header=None)
#
#     fig, ax = plt.subplots(figsize=(3, 3))
#     ax.set_xlim([0, 25])
#     ax.set_ylim([0, 25])
#
#     for i in range(len(df)):
#         ax.plot(df.iloc[i, 1]+df.iloc[i, 0]*0.2, df.iloc[i, 2]+df.iloc[i, 0]*0.1, '.', color='r', alpha=1-i/len(df))
#
#     plt.savefig(f'images/{step}.jpg', dpi=100)
#     plt.close()
#     print(f"{step} step done")


n_particles = 100
dim = 3
n_steps = 1000
skip = 100
# T = 1
# plot_scale = 150*np.sqrt(n_steps*(T/n_steps))

fig = plt.figure(figsize=(4, 4))
ax = fig.gca(projection='3d')

color = 'blue'
line = ax.plot([], [], [], '-', c=color)[0]
points = [None]*n_particles*2
for i in range(n_particles-1):
    points[i] = ax.plot([], [], [], '.', c=color)[0]
    points[i+n_particles] = ax.plot([], [], [], '.', c='grey', alpha=0.4)[0]
points[n_particles-1] = ax.plot([], [], [], '.', c='tab:red')[0]
points[-1] = ax.plot([], [], [], '.', c='darkred', alpha=0.4)[0]

ax.set_xlim3d(0, 25)
ax.set_ylim3d(0, 25)
ax.set_zlim3d(0, 25)
plt.suptitle(fr'{n_particles} particles sim', fontsize=20)


def animate(i):
    print(f"Processing {i+1}/{n_steps}")
    df = pd.read_csv(f'data/data_{i*skip}.csv', header=None)
    # line.set_data(result[0][0][:i], result[0][1][:i])
    # line.set_3d_properties(result[0][2][:i])

    for k in range(n_particles):
        points[k].set_data(df.iloc[k, 0], df.iloc[k, 1])
        points[k+n_particles].set_data(df.iloc[k, 0], df.iloc[k, 1])
        points[k].set_3d_properties(df.iloc[k, 2])
        points[k+n_particles].set_3d_properties(0)

    ax.view_init(30, 20+i*0.5)

    return (line, *points)


# call the animator. blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                               frames=n_steps, interval=30, blit=True)

# use if you want to save it into mp4
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=600, extra_args=['-vcodec', 'libx264'])
anim.save('3d-animated.mp4', writer=writer)

# display in IPython
# display(HTML(anim.to_jshtml()))
plt.close()
