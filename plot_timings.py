#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

fname = 'build/timings.dat'
df = pd.read_csv(fname, sep='\s+', header=None)

iters = df[0]
alloc_times = df[1]
dealloc_times = df[2]
kernel_times = df[3]
#total_times = df[4]

plt.xlabel('# iterations')
plt.ylabel('Time (s)')
plt.title('CLOUDSC memory pinning overhead')

plt.plot(iters, alloc_times, label='Init')
plt.plot(iters, dealloc_times, label='Final')
plt.plot(iters, kernel_times, label='Kernel')
#plt.plot(iters, total_times, label='Total')

fname = 'build/timings_packed.dat'
df = pd.read_csv(fname, sep='\s+', header=None)

alloc_times = df[1]
dealloc_times = df[2]
plt.plot(iters, alloc_times, color='tab:blue', linestyle='dashed', label='Init packed')
plt.plot(iters, dealloc_times, color='tab:orange', linestyle='dashed', label='Final packed')

fname = 'build/timings_static.dat'
df = pd.read_csv(fname, sep='\s+', header=None)

alloc_times = df[1]
dealloc_times = df[2]
plt.plot(iters, alloc_times[:158], color='tab:blue', linestyle='dashdot', label='Init baseline')
plt.plot(iters, dealloc_times[:158], color='tab:orange', linestyle='dashdot', label='Final baseline')

plt.xrange = [0, 160]

plt.legend(loc=(0.7, 0.125))
plt.savefig('build/pinning_overhead.png', bbox_inches='tight', dpi=300)
