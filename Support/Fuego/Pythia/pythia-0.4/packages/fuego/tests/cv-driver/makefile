FF = g77
CC = gcc
COPT = -O2
FOPT = -O

MECHANISM = HydrogenOxygen.ck2

cv: driver.o ddebdf.o burn_cv.o fuego.o
	$(CC) -o cv driver.o ddebdf.o burn_cv.o fuego.o -lg2c -lm

driver.o: driver.c
	$(CC) $(COPT) -c driver.c -o driver.o

fuego.o: fuego.py
	python fuego.py --file=fuego.c --mechanism=$(MECHANISM) --output=c
	$(CC) $(COPT) -c fuego.c -o fuego.o
 
ddebdf.o: ddebdf.f
	$(FF) $(FOPT) -c ddebdf.f

burn_cv.o: burn_cv.f
	$(FF) $(FOPT) -c burn_cv.f

clean:
	rm -f *.o
