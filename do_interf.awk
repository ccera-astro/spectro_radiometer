BEGIN {GRAN=20
       FS=","
        }
/./ {h=$4
     m=$5
     s=$6
     
     seconds = (h * 3600) + (m * 60) + s
     bin = int(seconds/GRAN)
     
     bindata[bin] += $13
     bin2data[bin] += $14
     bin3data[bin] += ($10+$11)
     bincount[bin] += 1
    }
    
END {
       ab4 = -9999
       a = 0.1
       for (i = 1; i < (86400/GRAN); i++)
       {
            if (bincount[i] > 0)
            {
				hour = (i * GRAN)
				hour /= 3600.0
				b = bindata[i]/bincount[i]
				b2 = bin2data[i]/bincount[i]
				b3 = bin3data[i]/bincount[i]
				b4 = (b*b)+(b2*b2)
				if (ab4 < -1000)
				{
				    ab4 = b4
				}
				ab4 = (a * b4) + ((1.0-a) * ab4)
				printf ("%5.3f %9.4f %9.4f %9.4f %9.4f\n", hour, b, b2, b3, sqrt(ab4))
			}
		}
	}
