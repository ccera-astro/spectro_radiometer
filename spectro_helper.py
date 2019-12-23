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
fft_buffer = [1.0e-25]*2048
fft2_buffer = [0.0]*2048
baseline_buffer=[0.0]*2048
baseline_buffer2=[0.0]*2048
freq_mask=[1.0]*2048
freq_mask_processed=False
tpwra=-99
tpwrb=-99
first_time = 0
pacet = time.time()
corr_sin = 1.0e-12
corr_cos = 1.0e-12
dpwr=0.0
apwr=0.0
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
    
    
def fft_log(p,p2,corr,frq,bw,longitude,normalize,prefix,decln,flist,again,ffa,mode,zt,decfile,tpi,spi):
    global fft_buffer
    global first_time
    global lastt
    global lasttpt
    global baseline_buffer
    global baseline_buffer2
    global fft2_buffer
    global tpwra, tpwrb
    global dpwr
    global apwr
    global freq_mask
    global freq_mask_processed
    global pacet
    global corr_cos, corr_sin

    preflist=prefix.split(",")
    
    if ((time.time() - pacet) < 0.05):
        return False
    
    pacet = time.time()
    
    if (len(p) != len(freq_mask)):
        fft_buffer = [1.0e-25]*len(p)
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
                
    if fft_buffer[10] < 1.0e-24:
        fft_buffer = numpy.copy(p)

    pwra = 0.0
    pwrb = 0.0
    diff = 0.0
    a = 0.05 * ffa

    
    #
    # Use numpy to process the log10-FFT buffers and turn them
    #  into linear values
    #
    ta = numpy.divide(p,[10.0]*len(p))
    ta = numpy.power([10.0]*len(ta),ta)
    ta = numpy.multiply(ta,freq_mask)
    pwra = numpy.sum(ta)
    
    tb = numpy.divide(p2,[10.0]*len(p))
    tb = numpy.power([10.0]*len(tb),tb)
    tb = numpy.multiply(tb,freq_mask)
    pwrb = numpy.sum(tb)
    
    ffts = [ta, tb]
    
    w = 0
    for  tf in ffts:
        #
        # To guard against log10 blowing up
        #
        tf = numpy.add([1.0e-22]*len(tf), tf)

        #
        # Integrate as a vector
        #
        tf = numpy.multiply([a]*len(tf),tf)
        tf2 = numpy.multiply([1.0-a]*len(tf), fft_buffer if w == 0 else fft2_buffer)
        tf = numpy.add(tf2,tf)

        #
        # Update it BEFORE we convert to log10, dummy
        #
        if (w == 0):
            fft_buffer = tf
        else:
            fft2_buffer = tf
        
        w += 1
    
    
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
    added = tpwra + tpwrb

    if (first_time == 0):
        first_time = int(time.time())
    
    dpwr = tpwra   
    apwr = 0.0
    if (mode == "differential" or mode == "diff"):
        dpwr = diff
        apwr = added
    
    if (mode == "interferometer" or mode == "corr" or mode == "correlator"):
        dpwr = corr_cos
        apwr = corr_sin

    #
    # Allow integrators to settle, etc, so don't write "ramp up" data
    #
    if ((time.time() - first_time) > 20):
            
        if (time.time() - lasttpt) >= tpi:
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
                f.write("%10.7f," % added)
                f.write("%10.7f,%10.7f\n" %  (corr_cos, corr_sin))
                f.close()
        
        if (time.time() - lastt) >= spi:
            lastt = time.time()
            ltp = time.gmtime()
            
            ffts = [fft_buffer, fft2_buffer]
            for w in [0,1]:
                for prefix in preflist:
                    fn = "%s%04d%02d%02d-%d-spec.csv" % (prefix, ltp.tm_year, ltp.tm_mon, ltp.tm_mday, w)
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
                        
                    for i in range(0,len(ffts[w])):
                        try:
                            f.write("%6.2f," % (math.log10(ffts[w][i])*10.0))
                        except:
                            f.write("???,")
                    f.write ("\n")
                    f.close()
                    
                    if (math.fabs(st_h - zt) < (35.0/3600.0)):
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


#
# Get (baselined) FFT buffer, whichever is selected for display
#
def curr_diff(pace,normalize,expected,which):
    global fft_buffer
    global baseline_buffer
    global fft2_buffer
    global baseline_buffer2
    
    inbuf = [fft_buffer, fft2_buffer]
    basebuf = [baseline_buffer, baseline_buffer2]
    
    if (len(fft_buffer) != expected):
        return ([-120.0]*expected)
    
    if (numpy.sum(basebuf[which]) == 0.0):
        x = smooth(inbuf[which])
    else:
        y = map(operator.add, basebuf[which], [1.0e-24]*len(basebuf[which]))
        x = map(operator.truediv, smooth(inbuf[which]), smooth(y))
    x = numpy.add(x, [1.0e-24]*len(x))
    x = numpy.abs(x)
    return numpy.multiply(numpy.log10(x),[10.0]*len(x))
    
def not_the_norm(vect):
    m=min(vect)
    s=[m]*len(vect)
    return map(operator.sub, vect, s)

def baseline_setter(thing):
    global baseline_buffer
    global baseline_buffer2
    global fft_buffer
    global fft_buffer2
    
    if thing != 0:
        baseline_buffer = smooth(fft_buffer,a=0.5)
        baseline_buffer = numpy.multiply(baseline_buffer,[0.95]*len(baseline_buffer))
        baseline_buffer2 = smooth(fft2_buffer,a=0.5)
        baseline_buffer2 = numpy.multiply(baseline_buffer2, [0.95]*len(baseline_buffer2))

def baseline_clearer(thing):
    global baseline_buffer
    global fft_buffer
    global baseline_buffer2
    global fft2_buffer
    
    if thing != 0:
        baseline_buffer = [0.0]*len(fft_buffer)
        baseline_buffer2 = [0.0]*len(fft2_buffer)

def lmst_string(pacer,longitude):
    return cur_sidereal(longitude).replace(",", ":")


TPLEN=3600
tp_vect=[0.0]*TPLEN
tp_vect2=[0.0]*TPLEN
def get_tp_vect(pacer):
    global tp_vect
    global tp_vect2
    global dpwr
    global TPLEN
    
    
    tp_vect=tp_vect[0:(TPLEN-1)]
    tp_vect=[dpwr]+tp_vect
    
    tp_vect2=tp_vect2[0:(TPLEN-1)]
    tp_vect2=[apwr]+tp_vect2
    
    if len(tp_vect) != TPLEN:
        print "Blarf!!! short TP_VECT"
    
    return ([tp_vect,tp_vect2])
            

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

def doppler_start(ifreq, dfreq, bw):
    if (dfreq == 0.0):
        stf = ifreq - bw/2.0
        dop = (ifreq - stf) / ifreq
        dop *= 299792.0
    else:
        stf = ifreq - bw/2.0
        dop = (dfreq - stf) / dfreq
        dop *= 299792.0
    return dop

ui_decln = [time.time(), None]
f_decln = [time.time(), None]
def get_decln(decln,decfile,dummy):
    global ui_decln
    global f_decln
    try:
        newdec = open(decfile,"r").readline().strip("\n")
        newdec = float(newdec)
        if (newdec != f_decln):
            f_decln = [time.time(), none]
    except:
        pass
    
    if (decln != ui_decln):
        ui_decln = [time.time(), decln]

    if (f_decln[0] > ui_decln[0]):
        return(f_decln[1])
    else:
        return(ui_decln[1])
 
def plotlabel(mode, which):
    w1 = {"total" : "TP", "tp" : "TP", "diff" : "Diff", "differential" : "Diff", "interf" : "Cos", "interferometer" : "Cos",
        "correlator" : "Cos"}
    w2 = {"total" : "Zero", "tp" : "Zero", "diff" : "Sum", "differential" : "Sum", "interf" : "Sin", "interferometer" : "Sin",
        "correlator" : "Sin"}
    
    ary = [w1,w2]
    
    di = ary[which]
    return (di[mode])

def get_spec_labels(mode):
    w1 = {"total" : "SKY", "tp" : "SKY", "diff" : "SKY", "differential" : "SKY", "interf" : "East", "interferometer" : "East",
        "correlator" : "East"}
    w2 = {"total" : "Zero", "tp" : "Zero", "diff" : "REF", "differential" : "REF", "interf" : "West", "interferometer" : "West",
        "correlator" : "West"}
    
    return ([w1[mode],w2[mode]])
    
