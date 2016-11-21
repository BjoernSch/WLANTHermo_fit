#!/usr/bin/env python3
# coding=utf-8

import os
from scipy.optimize import curve_fit
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot
import math
import sys
import csv
import json
import matplotlib.transforms as mtransforms

class wlt_fit:
    def __init__(self, filename):
        self.filename = filename
        self.name = os.path.basename(filename[:-4])
        self.directory = os.path.dirname(filename)
        self.id = 0
        self.rn = None
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.start_a = 3.35e-3
        self.start_b = 2.5e-4
        self.start_c = 3e-6
        self.err_a = 0.0
        self.err_b = 0.0
        self.err_c = 0.0
        self.beta25_85 = None
        self.rt_table = []
        self.r0 = 0
        self.r25 = 0
        self.r85 = 0
        self.temp_adc20 = 0
        self.temp_adc4075 = 0

        self.rmess = 47
        self.uref = 3.3
        self.adc_max = 2**12

        self.report = []
        self.mintemp = -40
        self.maxtemp = 315
        self.steptemp = 5
        self.csvmode = 'de'
        self.highres = 0.1

    def do_fit(self):
        self.read_csv()
        if self.rn is None:
            self.rn = self.search_rn()
        else:
            real_rn = self.search_rn()
            factor =  real_rn / self.rn
            self.err_a = self.err_a * factor
            self.err_b = self.err_b * factor
            self.err_c = self.err_c * factor
        self.fit()
        self.write_config()

    def do_fitreport(self):
        self.calc_error()
        self.write_fitcsv()
        self.plot_fit()

    def do_report(self):
        self.calc_report()
        self.write_reportcsv()
        self.write_reportdata()
        self.plot_report()

    def read_csv(self):
        self.rt_table = []
        if self.csvmode == 'de':
            delimiter = ';'
            dec_delimiter = ','
        else:
            delimiter = ','
            dec_delimiter = '.'

        with open(self.filename, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=delimiter)
            for line in reader:
                if not line:
                    continue
                (t, r) = line
                if t.startswith('temp'):
                    continue
                if not (t == '' or r == ''):
                    self.rt_table.append({'temp_measured': float(t.replace(',','.')),
                                  'r_ntc': float(r.replace(',','.'))})

    def r2t_wlt(self, rt, a, b, c):
        v = math.log(rt/self.rn)
        t = (1/(a + b*v + c*v*v)) - 273
        return t

    def r2t_wlt_plot(self,x, a, b, c):
        t = []
        for rt in x:
            v = math.log(rt/self.rn)
            t.append((1/(a + b*v + c*v**2)) - 273)
        return t

    def diff_list(self, x, y):
        diff = []
        for a, b in zip(x, y):
            diff.append(a-b)
        return diff

    def fit(self):
        (popt, pcov) = curve_fit(self.r2t_wlt_plot, np.array([x['r_ntc'] for x in self.rt_table]), np.array([x['temp_measured'] for x in self.rt_table]), (self.start_a, self.start_b, self.start_c), maxfev=2000)
        (self.a, self.b, self.c) = popt
        perr = np.sqrt(np.diag(pcov))
        (self.err_a, self.err_b, self.err_c) = perr

    def search_rn(self):
        rn = 1
        lastdiff = None
        for line in self.rt_table:
            diff = abs(line['temp_measured'] - 25)
            if lastdiff is None or diff < lastdiff:
                lastdiff = diff
                rn = line['r_ntc']
        return rn

    def calc_error(self):
        for line in self.rt_table:
            line['temp_calc'] = self.r2t_wlt(line['r_ntc'], self.a, self.b, self.c)
            line['temp_diff'] = line['temp_measured'] - line['temp_calc']

    def write_fitcsv(self):
        if self.csvmode == 'de':
            delimiter = ';'
            dec_delimiter = ','
        else:
            delimiter = ','
            dec_delimiter = '.'
        with open('{}_curvefit.csv'.format(os.path.join(self.directory, self.name)),'w') as csvfile:
            fieldnames = ['temp_measured', 'r_ntc', 'temp_calc', 'temp_diff']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            for line in self.rt_table:
                conv_line = {}
                for key, value in line.items():
                    conv_line[key] = str(value).replace('.',dec_delimiter)
                writer.writerow(conv_line)

    def plot_fit(self):
        print('Creating fitting plot...')
        fig = pyplot.figure(1, figsize=(10.91,7.48))
        pyplot.subplot(2, 1, 1)
        pyplot.plot([x['temp_measured'] for x in self.rt_table], [x['r_ntc'] for x in self.rt_table], 's', label='Measurement points')
        pyplot.plot([x['temp_calc'] for x in self.rt_table], [x['r_ntc'] for x in self.rt_table], 'k:', label='Curve')
        pyplot.ylabel('Resistance (k\u03a9)')
        pyplot.subplot(2, 1, 2)
        pyplot.plot([x['temp_measured'] for x in self.rt_table], [x['temp_diff'] for x in self.rt_table], 'r-', label='Deviation')
        pyplot.xlabel('Temperature (°C)')
        pyplot.ylabel('Deviation (°C)')
        pyplot.grid(True)
        fig.subplots_adjust(top=0.8)
        pyplot.figtext(0.3, 0.85, 'a: {0:1.7e} \u00b1 {3:1.7e}\nb: {1:1.7e} \u00b1 {4:1.7e}\nc: {2:1.7e} \u00b1 {5:1.7e}\nRn = {6} k\u03a9'.format(self.a,self.b,self.c,self.err_a,self.err_b,self.err_c,self.rn), bbox=dict(facecolor='white'))
        pyplot.savefig('{}_curve.png'.format(os.path.join(self.directory, self.name)),dpi=200)
        pyplot.close()

    def write_config(self):
        with open('{}.conf'.format(os.path.join(self.directory, self.name)),'w') as conffile:
            conffile.write('[{}]\n'.format(self.id))
            conffile.write('number = {}\n'.format(self.id))
            conffile.write('name = {}\n'.format(self.name))
            conffile.write('a = {0:1.7e}\n'.format(self.a))
            conffile.write('b = {0:1.7e}\n'.format(self.b))
            conffile.write('c = {0:1.7e}\n'.format(self.c))
            conffile.write('Rn = {0}\n'.format(self.rn))

    def t2r(self, temp):
        r = self.rn * np.exp((np.sqrt((self.b*273+self.b*temp)**2-4*(self.a*temp+273*self.a-1)*(self.c*temp+273*self.c))+self.b*(-1*temp)-273*self.b)/(2*(self.c*temp+273*self.c)))
        return r

    def get_adc(self, u):
        adc = int(round(u/(self.uref/(self.adc_max+1))))
        return adc

    def get_u_from_adc(self, adc):
        u = adc * (self.uref/(self.adc_max+1))
        return u

    def get_u(self, r_ntc):
        u = (self.uref/(r_ntc+self.rmess))*r_ntc
        return u

    def get_beta(self, r1, t1, r2, t2):
        beta = math.log(r1 / r2) / (1/(t1 + 273.15) - 1/(t2 + 273.15))
        return beta

    def calc_report(self):
        self.report = []
        lastval = None
        peak_res = None
        peak_res_temp = None
        highres_min = None
        highres_max = None
        self.report = []
        for adc_value in range(1, self.adc_max - 1):
            adc_value = self.adc_max - 1 - adc_value
            u_adc = self.get_u_from_adc(adc_value)
            r_ntc = self.rmess * ((4096/adc_value) - 1)
            temp = self.r2t_wlt(r_ntc, self.a, self.b, self.c)
            if lastval is None:
                resolution = False
                lastval = temp
            else:
                resolution = lastval - temp
                if peak_res is None or peak_res > resolution:
                    peak_res = resolution
                    peak_res_temp = temp

                if resolution <= 0.1:
                    if highres_min is None:
                        highres_min = temp
                    highres_max = temp
                lastval = temp

                self.report.append({'temp': temp,
                               'r_ntc': r_ntc,
                               'u_adc':u_adc,
                               'adc_value': adc_value,
                               'resolution': resolution})

        self.temp_adc20 = self.r2t_wlt(self.rmess * ((4096/20) - 1), self.a, self.b, self.c)
        self.temp_adc4075 = self.r2t_wlt(self.rmess * ((4096/4076) - 1), self.a, self.b, self.c)
        self.r0 = self.t2r(0)
        self.r25 = self.t2r(25)
        self.r85 = self.t2r(85)
        self.beta25_85 = self.get_beta(self.r25, 25, self.r85, 85)

        self.highres_min = highres_min
        self.highres_max = highres_max
        self.highres_area = highres_max - highres_min

        self.peak_res = peak_res
        self.peak_res_temp = peak_res_temp

    def write_reportcsv(self):
        if self.csvmode == 'de':
            delimiter = ';'
            dec_delimiter = ','
        else:
            delimiter = ','
            dec_delimiter = '.'

        with open('{}_report.csv'.format(os.path.join(self.directory, self.name)),'w') as csvfile:
            fieldnames = ['temp', 'r_ntc', 'u_adc', 'adc_value', 'resolution']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            for line in self.report:
                conv_line = {}
                for key, value in line.items():
                    conv_line[key] = str(value).replace('.',dec_delimiter)
                writer.writerow(conv_line)

    def write_reportdata(self):
        print('Writing report data...')
        fields = ['name', 'rn', 'a', 'b', 'c', 'err_a', 'err_b', 'err_c', 'rmess', 'r0', 'r25', 'r85', 'beta25_85', 'highres_min', 'highres_max', 'highres_area', 'peak_res', 'peak_res_temp', 'temp_adc4075', 'temp_adc20']
        reportdata = {}
        for field in fields:
            reportdata[field] = getattr(self, field)

        report_template = '''
## {name}
### Probe performance data

Values based on {rmess:.0f}k measurement resistor.
![Sensor performance chart]({name}_resolution.png)

Property | Symbol | Value
-------- | -------- | --------
Resistance at 0°C | R<sub>25</sub> | {r0:.2f}k
Resistance at 25°C | R<sub>25</sub> | {r25:.2f}k
Resistance at 85°C | R<sub>25</sub> | {r85:.2f}k
Beta 25°C to 85°C | B<sub>25/85</sub>| {beta25_85:.0f}K
Minimum measurable temperature | | {temp_adc4075:.1f}°C
Minimum high-res temperature | | {highres_min:.1f}°C
Highest resolution || {peak_res:1.2e}°C/step at {peak_res_temp:.1f}°C
Maximum high-res temperature | | {highres_max:.1f}°C
Maximum measurable temperature | | {temp_adc20:.1f}°C

### Probe curve data
![Probe fit chart]({name}_curve.png)

Property | Symbol | Value
-------- | -------- | --------
Resistance near 25°C | R<sub>25</sub><sup>1</sup> | {rn:.2f}k
Steinhart-Hart coefficient | a | {a:1.7e} ± {err_a:1.7e}
Steinhart-Hart coefficient | b | {b:1.7e} ± {err_b:1.7e}
Steinhart-Hart coefficient | c | {c:1.7e} ± {err_c:1.7e}

<sup>1</sup>: The deviation between this R<sub>25</sub> and the R<sub>25</sub> shown above is not relevant, this R<sub>25</sub> is taken from the original data point which is closest to 25°C. The value taken as a factor into the calculation of the final value and serves only a scaling purpose to the Steinhart-Hart coefficients.
'''

        with open('{}.md'.format(os.path.join(self.directory, self.name)),'w') as mdfile:
            mdfile.write(report_template.format(**reportdata))

    def plot_report(self):
        print('Creating resolution plot...')
        fig = pyplot.figure(1, figsize=(10.91,7.48), dpi=400)
        x_data = np.asanyarray([x['temp'] for x in self.report if x['resolution'] < 10])
        y_data = np.asanyarray([x['resolution'] for x in self.report if x['resolution'] < 10])
        fig, ax = pyplot.subplots()
        ax.plot(x_data, y_data, '-', label='Resolution')
        ax.set_yscale('log')
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        ax.axhline(0.1, color='green', lw=2, alpha=0.5)
        ax.axhline(1, color='yellow', lw=2, alpha=0.8)
        ax.fill_between(x_data, 0, 1, where=y_data <= self.highres, facecolor='green', alpha=0.3, transform=trans, linewidth=0.0)
        ax.fill_between(x_data, 0, 1, where=(y_data <= 1) & (~ (y_data <= self.highres)), facecolor='yellow', alpha=0.3, transform=trans, linewidth=0.0)
        pyplot.axvline(self.temp_adc20, color='red')
        pyplot.axvline(self.temp_adc4075, color='red')
        pyplot.xlabel('Temperature (°C)')
        pyplot.ylabel('Resolution (°C/step)')
        pyplot.grid(True)
        pyplot.savefig('{}_resolution.png'.format(os.path.join(self.directory, self.name)),dpi=200)
        pyplot.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.exit('usage: python3 %s file.csv' % sys.argv[0])

    for file_name in sys.argv[1:]:

        if not os.path.isfile(file_name):
            sys.exit('File "{}" not found'.format(file_name))

        if not file_name.endswith('.csv'):
            sys.exit('File "{}" is not a .csv file'.format(file_name))

        print('Starting with %s' % file_name)
        fitter = wlt_fit(file_name)

        print('Curve fitting...')
        fitter.do_fit()
        print('Creating fitting report...')
        fitter.do_fitreport()
        print('Creating probe report...')
        fitter.do_report()
