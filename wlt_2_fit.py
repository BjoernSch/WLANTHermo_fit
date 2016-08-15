#!/usr/bin/env python3
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

class wlt_fit:
    def __init__(self, name):
        self.name = name
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
        self.rt_table = []
        
        self.rmess = 47
        self.uref = 3.3
        self.adc_max = 2**12

        self.report = []
        self.mintemp = -40
        self.maxtemp = 315
        self.steptemp = 5
        self.csvmode = 'de'
        
        self.beta_t1 = 25
        self.beta_t2 = 85

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
        self.calc_beta()
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
        
        with open('{}.csv'.format(self.name), 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=delimiter)
            for line in reader:
                (t, r) = line
                if t.startswith('temp'):
                    continue
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
        with open('{}_curvefit.csv'.format(self.name),'w') as csvfile:
            fieldnames = ['temp_measured', 'r_ntc', 'temp_calc', 'temp_diff']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            for line in self.rt_table:
                conv_line = {}
                for key, value in line.items():
                    conv_line[key] = str(value).replace('.',dec_delimiter)
                writer.writerow(conv_line)
    
    def plot_fit(self):
        fig = pyplot.figure(1, figsize=(10.91,7.48))
        pyplot.subplot(2, 1, 1)
        pyplot.plot([x['temp_measured'] for x in self.rt_table], [x['r_ntc'] for x in self.rt_table], 's', label='Messpunkte')
        pyplot.plot([x['temp_calc'] for x in self.rt_table], [x['r_ntc'] for x in self.rt_table], 'k:', label='Kurve')
        pyplot.ylabel('Widerstand (k\u03a9)')
        pyplot.subplot(2, 1, 2)
        pyplot.plot([x['temp_measured'] for x in self.rt_table], [x['temp_diff'] for x in self.rt_table], 'r-', label='Differenz')
        pyplot.xlabel('Temperatur (°C)')
        pyplot.ylabel('Abweichung (°C)')
        pyplot.grid(True)
        fig.subplots_adjust(top=0.8)
        pyplot.figtext(0.3, 0.85, 'a: {0:1.7e} \u00b1 {3:1.7e}\nb: {1:1.7e} \u00b1 {4:1.7e}\nc: {2:1.7e} \u00b1 {5:1.7e}\nRn = {6} k\u03a9'.format(self.a,self.b,self.c,self.err_a,self.err_b,self.err_c,self.rn), bbox=dict(facecolor='white'))
        pyplot.savefig('{}_curve.png'.format(self.name),dpi=200)
        pyplot.close()
        
    def write_config(self):
        with open('{}.conf'.format(self.name),'w') as conffile:
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
    
    def get_u(self, r_ntc):
        u = (self.uref/(r_ntc+self.rmess))*r_ntc
        return u
    
    
    def calc_report(self):
        self.report = []
        lastval = None
        peak_res = 0
        peak_res_temp = None
        highres_min = None
        highres_max = None
        for temp in range(self.mintemp, self.maxtemp + self.steptemp, self.steptemp):
            r_ntc = self.t2r(temp)
            u_adc = self.get_u(r_ntc)
            adc = self.get_adc(u_adc)
            if lastval is None:
                resolution = 0
            else:
                resolution = (lastval - adc)/self.steptemp
                if resolution >= 10:
                    if highres_min is None:
                        highres_min = temp
                    highres_max = temp
                if resolution > peak_res:
                    peak_res = resolution
                    peak_res_temp = temp
            lastval = adc
            self.report.append({'temp': temp,
                               'r_ntc': r_ntc,
                               'u_adc':u_adc,
                               'adc': adc,
                               'resolution': resolution})
        self.highres_min = highres_min
        self.highres_max = highres_max
        self.highres_area = highres_max - highres_min
        self.peak_res = peak_res
        self.peak_res_temp = peak_res_temp
        
        
    def calc_beta(self):
        rt1 = self.t2r(self.beta_t1)
        rt2 = self.t2r(self.beta_t2)
        self.beta = round(math.log(rt1 / rt2) / (1/(self.beta_t1 + 273.15) - 1/(elf.beta_t2 + 273.15)))

    def write_reportcsv(self):
       if self.csvmode == 'de':
            delimiter = ';'
            dec_delimiter = ','
       else:
            delimiter = ','
            dec_delimiter = '.'
    
       with open('{}_report.csv'.format(self.name),'w') as csvfile:
            fieldnames = ['temp', 'r_ntc', 'u_adc', 'adc', 'resolution']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            for line in self.report:
                conv_line = {}
                for key, value in line.items():
                    conv_line[key] = str(value).replace('.',dec_delimiter)
                writer.writerow(conv_line)
    
    def write_reportdata(self):
        fields = ['name', 'rn', 'a', 'b', 'c', 'err_a', 'err_b', 'err_c', 'rmess', 'highres_min', 'highres_max', 'highres_area', 'peak_res', 'peak_res_temp', 'beta_t1', 'beta_t2', 'beta']
        reportdata = {}
        for field in fields:
            reportdata[field] = getattr(self, field)
        with open('{}_report.json'.format(self.name),'w') as jsonfile:
            jsonfile.write(json.dumps(reportdata, sort_keys=True, indent=4, separators=(',', ': ')))
    
    def plot_report(self):
        fig = pyplot.figure(1, figsize=(10.91,7.48), dpi=400)
        pyplot.plot([x['temp'] for x in self.report], [x['resolution'] for x in self.report], '-', label='Auflösung')
        pyplot.xlabel('Temperatur (°C)')
        pyplot.ylabel('Auflösung (Stufen/°C)')
        pyplot.grid(True)
        pyplot.savefig('{}_resolution.png'.format(self.name),dpi=200)
        pyplot.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit('usage: python3 %s file.csv' % sys.argv[0])

    file_name = sys.argv[1]

    if not os.path.isfile(file_name):
        sys.exit('File "{}" not found'.format(file_name))

    if not file_name.endswith('.csv'):
        sys.exit('File "{}" is not a .csv file'.format(file_name))

    name = os.path.splitext(file_name)[0]

    fitter = wlt_fit(name)

    fitter.do_fit()
    fitter.do_fitreport()
    fitter.do_report()
    
