spectro_radiometer.py: spectro_radiometer.grc
	grcc -d . spectro_radiometer.grc

install: spectro_radiometer.py spectro_helper.py
	cp spectro_radiometer.py /usr/local/bin
	cp spectro_helper.py /usr/local/bin
	chmod 755 /usr/local/bin/spectro_*.py
	ln -s /usr/local/bin/spectro_radiometer.py /usr/local/bin/spectro_radiometer

clean:
	rm -f spectro_radiometer.py
