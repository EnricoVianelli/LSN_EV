CC = g++
CFLAGS = -Wall -O3 -g --std=c++17

main.exe: main.o GA.o ../RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.exe: %.o GA.o ../RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.o : %.cpp
	$(CC) -o $@ -c $<  $(CFLAGS)
clean:
	find OUTPUT -type f -exec rm {} \;
	rm *.o *.exe
