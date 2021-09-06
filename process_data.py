import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

steps = 10**5
particles = 100
stop = 10**5

data = np.zeros(shape=(particles, 6, steps))

actual = 0
for step in range(0, steps, 100):
    df = pd.read_csv(f'data/data_{step}.csv', header=None)
    data[:, :, actual] = df.values

    actual += 1
    if step > stop:
        break

data = data[:, :, 0:actual]

# Potential = np.zeros(shape=actual)
# sigma = 0.86
# epsilon = 0.0001
# cutoff = 5*sigma
# for i in range(particles):
#     dr = data[:, :3, :]-data[i, :3, :]
#     dr2 = (dr**2).sum(axis=1)
#     dr2_mask = ~((dr2 > cutoff**2)|(dr2 == 0))
#     dr2_inv = sigma**2*dr2_mask*(1/np.clip(dr2, a_min=0.00001, a_max=None))
#     u = (4*epsilon*(dr2_inv**6-dr2_inv**3)).sum(axis=0)
#     Potential += u
# Potential /= 2

v2 = data[:, 3:6, :]**2
# KE = v2.sum(axis=1).sum(axis=0)
v = v2.sum(axis=1)**(1/2)
# Momentum = (data.sum(axis=0)**2).sum(axis=0)**(1/2)


# fig, axs = plt.subplots(2, 3, figsize=(15, 10))
# # plt.plot(KE/KE.mean()-1, label='Kinetic Energy')
# # plt.plot(np.clip(Potential, None, 0.01), label='Potential Energy')
# # plt.hist(KE+Potential, label='Full Energy')
# axs[0, 0].plot(KE)
# axs[0, 0].set_title('Kinetic Energy')
# axs[1, 0].plot(KE, '.')
# axs[1, 0].set_ylim([9515.1, 9515.3])
# axs[0, 1].plot(Potential)
# axs[0, 1].set_title('Potential Energy')
# axs[1, 1].plot(Potential, '.')
# axs[1, 1].set_ylim([-0.05, 0.1])
# axs[0, 2].plot(Potential*2+KE)
# axs[0, 2].set_title(r'Kinetic + Potential $\times$ 2')
# axs[1, 2].plot(Potential*2+KE, '.')
# axs[1, 2].set_ylim([9515.1, 9515.5])
# plt.suptitle(f'Корреляция {np.corrcoef(Potential, KE)[0][1]}', fontsize=20)
# plt.savefig('Full_Energy_100000.pdf')

# fig, ax = plt.subplots(figsize=(5,5))
# plt.plot(Momentum)
# plt.savefig('Momentum_100000.pdf')

# fig, ax = plt.subplots(figsize=(5, 5))
# plt.plot(data[0, 3, :])
# plt.plot(data[1, 3, :])
# plt.plot(data[2, 3, :])
# plt.plot(data[3, 3, :])
# plt.plot(data[4, 3, :])
# plt.savefig('x_velocity_100000.pdf')


fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
axs[0].hist(v[:, 0], bins=25)
axs[0].set_xlim([0, 20])
axs[0].set_ylim([0, 12])
axs[1].hist(v[:, 50], bins=25)
axs[1].set_xlim([0, 20])
axs[1].set_ylim([0, 12])
axs[2].hist(v[:, 99], bins=25)
axs[2].set_xlim([0, 20])
axs[2].set_ylim([0, 12])
plt.savefig('hist_velocities_10000.pdf')
#
# for i in range(100):
#     fig, ax = plt.subplots(figsize=(4, 4))
#     plt.hist(v[:, i], bins=25)
#     plt.xlim([0, 20])
#     plt.ylim([0, 12])
#     plt.savefig(f'images/{i}.jpg')
#     plt.close()

# MSD = ((data[:, :3, :]-data[:, :3, 0:1])**2).sum(axis=1).mean(axis=0)
#
# fig, ax = plt.subplots(figsize=(4, 4))
# plt.plot(MSD)
# plt.savefig('Mean_squared_displacement_100000.pdf')
#

