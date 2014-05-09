animacion.py: datosevol.dat
	python animacion.py

datosevol.dat: schrod.x
	./schrod.x 5000

schrod.x: schrod.c
	cc schrod.c -o schrod.x -lm

clean:
	rm *.dat *.x