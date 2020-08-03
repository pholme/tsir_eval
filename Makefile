SRC = .
CFLAGS = -Ofast -march=native
LDFLAGS = 
CC = gcc

OBJ1 = o/tsir.o o/misc.o o/heap.o o/pcg_rnd.o
OBJ2 = o/tsir_ref.o o/pcg_rnd.o

all : tsir tsir_ref

tsir: $(OBJ1)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

tsir_ref: $(OBJ2)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

o/tsir.o : $(SRC)/tsir.c $(SRC)/tsir.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/tsir.c -o $@

o/tsir_ref.o : $(SRC)/tsir_ref.c $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/tsir_ref.c -o $@

o/misc.o : $(SRC)/misc.c $(SRC)/tsir.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/misc.c -o $@

o/heap.o : $(SRC)/heap.c $(SRC)/tsir.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/heap.c -o $@

o/pcg_rnd.o : $(SRC)/pcg_rnd.c $(SRC)/tsir.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/pcg_rnd.c -o $@
