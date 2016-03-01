#!/usr/bin/env python3
import os
from scipy.optimize import curve_fit
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot
import math

file = open('datei.csv')
temperaturen = []
widerstaende = []

for line in file:
    print(line)
    (t, r) = line.split(';')
    temperaturen.append(float(t.replace(',','.')))
    widerstaende.append(float(r.replace(',','.')))

def r2t_wlt(rt, a, b, c):
    v = math.log(rt/1)
    t = (1/(a + b*v + c*v*v)) - 273
    return t

def r2t_wlt_plot(x, a, b, c):
    t = []
    for rt in x:
        v = math.log(rt/256)
        t.append((1/(a + b*v + c*v**2)) - 273)
    return t

def diff_list(x, y):
    diff = []
    for a, b in zip(x, y):
        diff.append(a-b)
    return diff

popt, pcov = curve_fit(r2t_wlt_plot, widerstaende, temperaturen)

(a, b, c) = popt
print(pcov)
temperaturen_curve = r2t_wlt_plot(widerstaende, *popt)
fig=pyplot.figure(1, figsize=(10.91,7.48))
pyplot.subplot(2, 1, 1)
pyplot.plot(temperaturen, widerstaende, 's', label='Messpunkte')
pyplot.plot(temperaturen_curve, widerstaende, 'k:', label='Kurve')
pyplot.ylabel('Widerstand (kOhm)')
pyplot.subplot(2, 1, 2)
pyplot.plot(temperaturen, diff_list(temperaturen, temperaturen_curve), 'r-', label='Differenz')
pyplot.xlabel('Temperatur (°C)')
pyplot.ylabel('Abweichung (°C)')
pyplot.grid(True)
fig.subplots_adjust(top=0.8)
pyplot.figtext(0.93, 0.5, 'a: ' + str(a) + '\nb: ' + str(b) + '\nc: '+ str(c), bbox=dict(facecolor='white'))
pyplot.savefig("test.png")

