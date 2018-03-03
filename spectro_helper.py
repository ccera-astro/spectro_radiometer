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

def do_annotation(ra,dec,baseline,annotation,bw,abw,freq,srate,prefix):
    ltp = time.gmtime()
    if len(annotation) <= 0:
		return True
    fn = "%s%04d%02d%02d-annotation.txt" % (prefix, ltp.tm_year, ltp.tm_mon, ltp.tm_mday)
    fp = open(fn, "a")
    fp.write ("Annotation update: %04d%02d%02d-%02d:%02d:%02d\n" % 
        (ltp.tm_year, ltp.tm_mon, ltp.tm_mday, ltp.tm_hour, ltp.tm_min, ltp.tm_sec))
    fp.write ("RA: %6.2f\n" % ra)
    fp.write ("DEC: %6.2f\n" % dec)
    fp.write ("Baseline length: %6.2fm\n" % baseline)
    fp.write ("Annotation: %s\n" % annotation)
    fp.write ("Input Filter Bandwidth: %6.2fMHz\n" % (bw/1.0e6))
    fp.write ("Analog Bandwidth: %6.2fMHz\n" % (abw/1.0e6))
    fp.write ("Frequency: %8.4fMHz\n" % (freq/1.0e6))
    fp.write ("Sample rate: %6.2fMsps\n\n" % (srate/1.0e6))
    fp.close()
    return True
    
    
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

#
#
# Calculations for automated fringe-stopping
#
# We are called on a regular basis to produce a complex rotation value that is
#  applied to one leg of the interferometer, to reduce the fringe frequency to
#  close-to zero, thus allowing longer integration times.
#
def ha (ra,longit):
    #
    # First get current sidereal time as as HH:MM:SS string
    #
    lmst = cur_sidereal (longit)
    parts = lmst.split(",")
    
    #
    # Re-express lmst (string) as a decimal hours
    #
    lmst = float(parts[0])
    lmst += float(parts[1]) / 60.0
    lmst += float(parts[2]) / 3600.0
    
    #
    # Now compute relative hour-angle between us, and the object in question (at some RA)
    #   When we're done, we have the relative hour-angle in radians
    #
    h = ra - lmst
    return (h)
#
# Phase accumulator
#
phase_accum = 0.0

#
# For time-delta calcs
#
last_time_phase = -99.0

#
# Earth rotation in radians/second
#
earth = math.radians(360.0 / 86400.0)

#
# Since part of the calculation only needs to be done when baseline or
#  frequency changes, we keep 'em global
gbaseline = 0.0
gfreq = 0.0
basefreqrot = 0.0

#
# Another part only changes when dec/latit changes
#
gdec = -99.0
glatit = -99.0
declatrot = 0.0

#
# The combination of the two
#
baserot = 0.0

#
# Speed of light, in case it wasn't obvious :)
#
C = 299792000.0

#
#
# Compute phase-rotation value based either on manual input, or automatic
#   fringe rotation
#
# pacer   - ignored, simply used as a "dummy" variable in the flow-graph to
#           make sure that we're called regularly.
# ra       - The RA in fractional hours
# dec      - The DEC in fractional degrees
# longit   - Our geographic longitude
# latit    - Our geographic latitude
# baseline - Baseline length, in meters
# ena      - Whether automagic or manual is in use
# manval   - The manual value
# freq     - The sky frequency
#
#
# These calculations assume a meridian-transit type setup, with a strictly east-west baseline
#
def fringe_stop (pacer, ra, dec, longit, latit, baseline, ena, manval, freq):
    global phase_accum
    global last_time_phase
    global earth
    global gbaseline
    global gfreq
    global basefreqrot
    global declatrot
    global gdec
    global glatit
    global baserot
 
    #
    # Just compute the fixed-value from the degrees input by the UI
    #
    if ena == False:
        radians = math.radians(manval)
        return (complex(math.cos(radians),math.sin(radians)))

   
    #
    # Handle initialization of our time-delta calculator
    #
    if last_time_phase < 0:
        last_time_phase = time.time()

    #
    # Compute time delta compared to last-time
    #
    now = time.time()
    tdelt = now - last_time_phase
    
    #
    # Remember "now" for next time
    #
    last_time_phase = now
    
    changed = False
    
    #
    # If baseline or freq parameters have changed
    #
    if baseline != gbaseline or freq != gfreq:
        gbaseline = baseline
        gfreq = freq
        basefreqrot = (baseline * earth)/(C / freq)
        changed = True
    
    #
    # If dec or latit changes
    #
    if dec != gdec or latit != glatit:
        gdec = dec
        glatit = latit
        declatrot = math.cos(math.radians(dec))*math.cos(math.radians(latit))
        changed = True
    
    #
    # Update the combined "baserot"
    #
    if changed == True:
        baserot = basefreqrot * declatrot
        
    #
    # First get current sidereal time as as HH:MM:SS string
    #
    lmst = cur_sidereal (longit)
    parts = lmst.split(",")
    
    #
    # Re-express lmst (string) as a decimal hours
    #
    lmst = float(parts[0])
    lmst += float(parts[1]) / 60.0
    lmst += float(parts[2]) / 3600.0
    
    #
    # Now compute relative hour-angle between us, and the object in question (at some RA)
    #   When we're done, we have the relative hour-angle in radians
    #
    h = ra - lmst
    h *= math.radians((360.0 / 24.0))
        
    #
    # Compute fringe rate in radians/hour
    #
    # h       - relative hour-angle, in radians between us and source
    #
    # Other parameters to the equation don't change with time, so they're
    #   efficiently pre-computed above, and re-computed when the input
    #   parameters change
    #
    F = baserot*math.sin(h)
    
    #
    # Because I get lost in the units, the above may be in radians/hour
    #   or radians/second so we can just adjust this scaling for the calculation below
    #
    scale = 1.0/3600.0
    
    #
    # Update the phase accumulator (stored in radians)
    # Handle phase-wrap gracefully
    # (although, the Python transcendentals handle this just fine themselves)
    #
    phase_accum += tdelt*(F * scale)
    if phase_accum > 2.0*math.pi:
        phase_accum -= (2.0*math.pi)
    
    #
    # We return a value that can be used by the phase rotator in the flow-graph
    #
    rval = complex(math.cos(phase_accum),math.sin(phase_accum))
    return (rval)

