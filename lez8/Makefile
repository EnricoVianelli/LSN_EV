CC = g++
CFLAGS = -Wall -O3 -g --std=c++17

main.exe: main.o VMC_simu.o ../RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.exe: %.o VMC_simu.o ../RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.o : %.cpp
	$(CC) -o $@ -c $<  $(CFLAGS)
clean:
	rm *.o *.exe
