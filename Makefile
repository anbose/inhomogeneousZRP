CC = mpic++

RANDOM = -Lpcg-c-0.94/src 
#-lpcg_random

#Libraries and opt
FLAGS = -Wall -O3

all: PDE

PDE: main.cc

	${CC} -o mainRKMC.exe $^ ${RANDOM} ${FLAGS}

test: test.cc

	${CC} -o test.exe $^ ${RANDOM} ${FLAGS}

clean:
	rm -f *.exe *.o
