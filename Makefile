CC = g++-4.8.4
Option = -O3 -Wall -std=c++11 -fopenmp

Objs = main.o rand_base.o rand_gamma.o rand_knscat.o \
	rand_planck.o electron.o photon.o profile.o \
	input.o driver.o seed.o

Out_file = mcic

main : $(Objs)
	$(CC) $(Option) $(Objs) -o $(Out_file)

main.o : main.cpp
	$(CC) $(Option) -c -o main.o main.cpp

input.o : input.cpp input.h
	$(CC) $(Option) -c -o input.o input.cpp

driver.o : driver.cpp driver.h
	$(CC) $(Option) -c -o driver.o driver.cpp

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

photon.o : photon.cpp photon.h
	$(CC) $(Option) -c -o photon.o photon.cpp

profile.o : profile.cpp profile.h
	$(CC) $(Option) -c -o profile.o profile.cpp

seed.o : seed.cpp seed.h
	$(CC) $(Option) -c -o seed.o seed.cpp

clean:
	rm *.o *~ mcic
