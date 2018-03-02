# this module will be imported in the into your flowgraph
import sys
import os
import time
import ephem
import operator
import math
import numpy



lastt = time.time()
lasttpt = time.time()
fft_buffer = [-900.0]*2048
fft2_buffer = [0.0]*2048
baseline_buffer=[0.0]*2048
freq_mask=[1.0]*2048
freq_mask_processed=False
tpwra=-99
tpwrb=-99
first_time = 0
pacet = time.time()
corr_sin = 1.0e-12
corr_cos = 1.0e-12
dpwr=0.0

import copy
def fft_log(p,p2,corr,frq,bw,longitude,normalize,prefix,decln,flist,again,ffa,mode,zt):
    global fft_buffer
    global first_time
    global lastt
    global lasttpt
    global baseline_buffer
    global tpwra, tpwrb
    global dpwr
    global freq_mask
    global freq_mask_processed
    global pacet
    global corr_cos, corr_sin

    preflist=prefix.split(",")
    
    if ((time.time() - pacet) < 0.05):
        return False
    
    pacet = time.time()
    
    if (len(p) != len(freq_mask)):
        fft_buffer = [-900.0]*len(p)
        fft2_buffer = [0.0]*len(p)
        baseline_buffer=[0.0]*len(p)
        freq_mask=[1.0]*len(p)
    
    fplen = float(len(p))
    if flist != "" and freq_mask_processed == False:
        nmsks = 0
        freq_mask_processed = True
        stf = frq-(bw/2.0)
        incr = bw/fplen
        flist = flist.split(",")
        for i in range(0,len(freq_mask)):
            for f in flist:
                if math.fabs(stf-float(f)) <= (incr*2.5):
                    freq_mask[i] = 1.0e-15
                    print "Setting mask at F %f position %d (%f)" % (float(f), i, incr*float(i))
                    nmsks += 1
            stf += incr
        print "Total mask percentage %f" % ((float(nmsks)/fplen)*100.0)
                
    if fft_buffer[10] < -800:
        for i in range(0,len(fft_buffer)):
            fft_buffer[i] = p[i]

    pwra = 0.0
    pwrb = 0.0
    diff = 0.0
    a = 0.250 * ffa
    
    #
    # Use numpy to process the log10-FFT buffers
    #
    ta = numpy.divide(p,[10.0]*len(p))
    ta = numpy.power([10.0]*len(ta),ta)
    ta = numpy.multiply(ta,freq_mask)
    pwra = numpy.sum(ta)
    
    tb = numpy.divide(p2,[10.0]*len(p))
    tb = numpy.power([10.0]*len(tb),tb)
    tb = numpy.multiply(tb,freq_mask)
    pwrb = numpy.sum(tb)
    
    #
    # Differential FFT
    #
    td = numpy.subtract(tb,ta)
    tf = ta
    if (mode == "diff" or mode == "differential"):
        tf = td
    
    #
    # To guard against log10 blowing up
    #
    tf = numpy.add([1.0e-20]*len(tf), tf)
    
    #
    # Turn linears back into logs
    #
    tf = numpy.log10(tf)
    
    #
    # Integrate as a vector
    #
    tf = numpy.multiply([a]*len(tf),tf)
    tf2 = numpy.multiply([1.0-a]*len(tf), fft_buffer)
    tf = numpy.add(tf2,tf)
    
    fft_buffer = tf
    
    atp=0.02
    #
    # Gain settings
    #
    corr_a = corr.real * again
    corr_b = corr.imag * again
    pwra *= again
    pwrb *= again
    
    #
    # Prime the integrator
    #
    if (tpwra < -10):
        tpwra = pwra
        tpwrb = pwrb
        
    if (corr_cos < -1.0e-10):
        corr_cos = corr_a
        corr_sin = corr_b
    
    #
    # Integration
    #
    corr_cos = (atp * corr_a) + ((1.0 - atp)*corr_cos)
    corr_sin = (atp * corr_b) + ((1.0 - atp)*corr_sin)

    tpwra = (atp * pwra) + ((1.0 - atp)*tpwra)
    tpwrb = (atp * pwrb) + ((1.0 - atp)*tpwrb)
    
    diff = tpwra - tpwrb

    if (first_time == 0):
        first_time = int(time.time())
    
    dpwr = tpwra    
    if (mode == "differential" or mode == "diff"):
        dpwr = diff
    
    if (mode == "interferometer" or mode == "corr" or mode == "correlator"):
        dpwr = corr_cos

    #
    # Allow integrators to settle, etc, so don't write "ramp up" data
    #
    if ((time.time() - first_time) > 90):
        if (time.time() - lasttpt) >= 2:
            lasttpt = time.time()
            
            ltp = time.gmtime()
            for prefix in preflist:
                fn="%s%04d%02d%02d-tp.csv" % (prefix, ltp.tm_year, ltp.tm_mon, ltp.tm_mday)
                f = open (fn, "a")
                f.write ("%02d,%02d,%02d," % (ltp.tm_hour, ltp.tm_min, ltp.tm_sec))
                f.write ("%s," % cur_sidereal(longitude))
                f.write("%5.5f," % (frq/1.0e6))
                f.write("%5.5f," % bw)
                f.write("%5.1f," % decln)
                f.write("%10.7f," % tpwra)
                f.write("%10.7f," % tpwrb)
                f.write("%10.7f," % diff)
                f.write("%10.7f,%10.7f\n" %  (corr_cos, corr_sin))
                f.close()
        
        if (time.time() - lastt) >= 20:
            lastt = time.time()
            ltp = time.gmtime()
            tm = smooth(fft_buffer)
            bb = smooth(baseline_buffer)
            if normalize:
                tm = norm(smooth(fft_buffer))
                bb = norm(baseline_buffer)
            tm = map(operator.sub, tm, bb)
            
            for prefix in preflist:
                fn = "%s%04d%02d%02d-spec.csv" % (prefix, ltp.tm_year, ltp.tm_mon, ltp.tm_mday)
                f = open (fn, "a")
                f.write ("%02d,%02d,%02d," % (ltp.tm_hour, ltp.tm_min, ltp.tm_sec))
                f.write ("%s," % cur_sidereal(longitude))
                f.write("%5.5f," % (frq/1.0e6))
                f.write("%5.5f," % bw)
                f.write("%5.1f," % decln)
                
                st = cur_sidereal(longitude)
                st = st.split(",")
                st_h = float(st[0])
                st_h += float(st[1])/60.0
                st_h += float(st[2])/3600.0
                    
                for i in range(0,len(fft_buffer)):
                    f.write("%6.2f," % tm[i])
                f.write ("\n")
                f.close()
                
                if (math.fabs(st_h - zt) < (60.0/3600.0)):
                    baseline_setter(1)
            return 1
        
    return 0
 
def smooth(vect,a=0.4):
    ovect=[0]*len(vect)
    val=vect[0]
    for i in range(0,len(vect)):
        val = (a*vect[i]) + ((1.0 - a)*val)
        ovect[i] = val
    
    return ovect   
    
def avgvect(v1, v2):
    q = map(operator.add, v1, v2)
    return map(operator.div, q, [2]*len(v1))

def cur_sidereal(longitude):
    longstr = "%02d" % int(longitude)
    longstr = longstr + ":"
    longitude = abs(longitude)
    frac = longitude - int(longitude)
    frac *= 60
    mins = int(frac)
    longstr += "%02d" % mins
    longstr += ":00"
    x = ephem.Observer()
    x.date = ephem.now()
    x.long = longstr
    jdate = ephem.julian_date(x)
    tokens=str(x.sidereal_time()).split(":")
    hours=int(tokens[0])
    minutes=int(tokens[1])
    seconds=int(float(tokens[2]))
    sidt = "%02d,%02d,%02d" % (hours, minutes, seconds)
    return (sidt)
    
def curr_findx(pace):
    global  sched_index
    global schedule
    global baseline_count
    global MAX_BASELINE
    
    if baseline_count > MAX_BASELINE:
        return 0

    cur=schedule[sched_index]
    sched_index = sched_index + 1
    if sched_index >= len(schedule):
        sched_index = 0
    return cur



def curr_diff(pace,normalize):
    global fft_buffer
    global baseline_buffer
    if normalize:
        x = map(operator.sub, norm(smooth(fft_buffer)), norm(smooth(baseline_buffer)))
    else:
        x = map(operator.sub, smooth(fft_buffer), smooth(baseline_buffer))
    return x
    
def norm(vect):
    m=min(vect)
    s=[m]*len(vect)
    return map(operator.sub, vect, s)

def baseline_setter(thing):
    global baseline_buffer
    global fft_buffer
    
    if thing != 0:
        baseline_buffer = smooth(fft_buffer,a=0.75)

def baseline_clearer(thing):
    global baseline_buffer
    
    if thing != 0:
        baseline_buffer = [0.0]*len(fft_buffer)

def lmst_string(pacer,longitude):
    return cur_sidereal(longitude).replace(",", ":")


TPLEN=3600
tp_vect=[0.0]*TPLEN
def get_tp_vect(pacer):
    global tp_vect
    global dpwr
    global TPLEN
    
    
    tp_vect=tp_vect[0:(TPLEN-1)]
    tp_vect=[dpwr]+tp_vect
    
    if len(tp_vect) != TPLEN:
        print "Blarf!!! short TP_VECT"
    
    return (tp_vect)
            

tpwr = 0 
def lmst_hours(pacer,longitude):
    x = lmst_string(1,longitude)
    
    x = x.split(":")
    hours = float(x[0])
    hours += float(x[1])/60.0
    hours += float(x[2])/3600.0
    
    return hours

