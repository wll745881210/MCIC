CC = g++
Option = -O3 -Wall -std=c++11

Objs = main.o rand_base.o rand_gamma.o rand_knscat.o \
	rand_planck.o electron.o

Out_file = mcic

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

rand_base.o : rand_base.cpp rand_base.h
	$(CC) $(Option) -c -o rand_base.o rand_base.cpp

rand_gamma.o : rand_gamma.cpp rand_gamma.h
	$(CC) $(Option) -c -o rand_gamma.o rand_gamma.cpp

rand_knscat.o : rand_knscat.cpp rand_knscat.h
	$(CC) $(Option) -c -o rand_knscat.o rand_knscat.cpp

rand_planck.o : rand_planck.cpp rand_planck.h
	$(CC) $(Option) -c -o rand_planck.o rand_planck.cpp

electron.o : electron.cpp electron.h
	$(CC) $(Option) -c -o electron.o electron.cpp

clean:
	rm *.o *~ clpt.exe
