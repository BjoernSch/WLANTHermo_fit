#!/usr/bin/env python3
import os
from scipy.optimize import curve_fit
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot
import math

file = open('datei.csv')
temperaturen = []
widerstaende = []

rn = 1

for line in file:
    (t, r) = line.split(';')
    temperaturen.append(float(t.replace(',','.')))
    widerstaende.append(float(r.replace(',','.')))

def r2t_wlt(rt, a, b, c):
    v = math.log(rt/rn)
    t = (1/(a + b*v + c*v*v)) - 273
    return t

def r2t_wlt_plot(x, a, b, c):
    t = []
    for rt in x:
        v = math.log(rt/rn)
        t.append((1/(a + b*v + c*v**2)) - 273)
    return t

def diff_list(x, y):
    diff = []
    for a, b in zip(x, y):
        diff.append(a-b)
    return diff

(popt, pcov) = curve_fit(r2t_wlt_plot, np.array(widerstaende), np.array(temperaturen), maxfev=2000)

(a, b, c) = popt
perr = np.sqrt(np.diag(pcov))
(err_a, err_b, err_c) = perr

temperaturen_curve = r2t_wlt_plot(widerstaende, *popt)
fig=pyplot.figure(1, figsize=(10.91,7.48))
pyplot.subplot(2, 1, 1)
pyplot.plot(temperaturen, widerstaende, 's', label='Messpunkte')
pyplot.plot(temperaturen_curve, widerstaende, 'k:', label='Kurve')
pyplot.ylabel('Widerstand (k\u03a9)')
pyplot.subplot(2, 1, 2)
pyplot.plot(temperaturen, diff_list(temperaturen, temperaturen_curve), 'r-', label='Differenz')
pyplot.xlabel('Temperatur (°C)')
pyplot.ylabel('Abweichung (°C)')
pyplot.grid(True)
fig.subplots_adjust(top=0.8)
pyplot.figtext(0.3, 0.85, 'a: {0:1.7e} \u00b1 {3:1.7e}\nb: {1:1.7e} \u00b1 {4:1.7e}\nc: {2:1.7e} \u00b1 {5:1.7e}\nRn = {6} k\u03a9'.format(a,b,c,err_a,err_b,err_c,rn), bbox=dict(facecolor='white'))
pyplot.savefig("test.png")

print('Werte für sensor.conf')
print('a = {0:1.7e}'.format(a))
print('b = {0:1.7e}'.format(b))
print('c = {0:1.7e}'.format(c))
print('Rn = {0}'.format(rn))
