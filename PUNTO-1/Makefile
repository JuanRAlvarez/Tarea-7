animacion.py: datosevol.dat
	python animacion.py

datosevol.dat: waveq.x
	./waveq.x 5000

waveq.x: waveq.c
	cc waveq.c -o waveq.x -lm

clean:
	rm *.dat *.x