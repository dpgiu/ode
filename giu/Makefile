FC = gfortran
FLAGS = -g -Wall
TARGET = main

# $@ sta per il nome del file da produrre (a sinistra dei due punti)
# $< sta per la prima dipendenza (a destra dei due punti)
# $^ sta per tutte le dipendenze (tutte le cose a destra dei due punti)

$(TARGET): eulero.o solvers.o
	$(FC) $^ -o $(TARGET)

solvers.o: solvers.f90
	$(FC) $(FLAGS) -c $< -o $@

eulero.o: eulero.f90 solvers.o
	$(FC) $(FLAGS) -c $< -o $@

clean:
	rm *.o *.mod $(TARGET)


