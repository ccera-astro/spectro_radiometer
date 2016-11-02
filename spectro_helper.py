# this module will be imported in the into your flowgraph
import sys
import os
import time
import ephem
import operator
import math



lastt = time.time()
lasttpt = time.time()
fft_buffer = [-900.0]*2048
fft2_buffer = [0.0]*2048
baseline_buffer=[0.0]*2048
freq_mask=[1.0]*2048
freq_mask_processed=False
tpwr=-99
first_time = 0
pacet = time.time()
import copy
def fft_log(p,p2,corr,frq,bw,longitude,normalize,prefix,decln,flist,again,ffa,mode,zt):
    global fft_buffer
    global first_time
    global lastt
    global lasttpt
    global baseline_buffer
    global tpwr
    global freq_mask
    global freq_mask_processed
    global pacet
    
    preflist=prefix.split(",")
    
    if ((time.time() - pacet) < 0.05):
        return False
    
    pacet = time.time()
    
    
    if flist != "" and freq_mask_processed == False:
        freq_mask_processed = True
        stf = frq-(bw/2)
        incr = bw/2048.0
        flist = flist.split(",")
        for i in range(0,len(freq_mask)):
            for f in flist:
                if math.abs(stf-float(f)) < incr:
                    freq_mask[i] = 1.0e-12
            stf += incr
                
    
    if fft_buffer[10] < -800:
        for i in range(0,len(fft_buffer)):
            fft_buffer[i] = p[i]

    pwr = 0.0
    a = 0.250 * ffa
    for i in range(0,len(fft_buffer)):
        pwr += (math.pow(10.0,p[i]/10.0)*freq_mask[i])
        q = p[i]
        if (mode == "diff" or mode == "differential"):
            pwr -= math.pow(10.0,p2[i]/10.0)
            q = math.pow(10.0,p[i]/10.0)
            q -= math.pow(10.0,p2[i]/10.0)
            if (q <= 0.0):
                q = 1.0e-12
            q = math.log10(q)
   
        fft_buffer[i] = (a * q) + ((1.0 - a) * fft_buffer[i])
    
    if (mode != "corr" and mode != "correlator" and mode != "interferometer"):
        pwr = pwr * again
    else:
        pwr = corr.real * again
    
    if (tpwr < -10):
        tpwr = pwr
    
    atp=0.1
    tpwr = (atp * pwr) + ((1.0 - atp)*tpwr)

    if (first_time == 0):
        first_time = int(time.time())
    
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
                f.write("%10.7f," % tpwr)
                f.write("%10.7f,%10.7f\n" %  (corr.real, corr.imag))
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
                
                if (math.abs(st_h - zt) < (60.0/3600.0)):
                    baseline_setter(1)
                


                for i in range(0,len(fft_buffer)):
                    f.write("%6.2f," % tm[i])
                f.write ("\n")
                f.close()
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
        baseline_buffer = [0.0]*2048

def lmst_string(pacer,longitude):
    return cur_sidereal(longitude).replace(",", ":")


TPLEN=3600
tp_vect=[0.0]*TPLEN
def get_tp_vect(pacer):
    global tp_vect
    global tpwr
    global TPLEN
    
    
    tp_vect=tp_vect[0:(TPLEN-1)]
    tp_vect=[tpwr]+tp_vect
    
    if len(tp_vect) != TPLEN:
        print "Blarf!!! short TP_VECT"
    
    return (tp_vect)
            

def lmst_hours(pacer,longitude):
    x = lmst_string(1,longitude)
    
    x = x.split(":")
    hours = float(x[0])
    hours += float(x[1])/60.0
    hours += float(x[2])/3600.0
    
    return hours

