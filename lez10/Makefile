CC = mpicxx
CFLAGS = -Wall -O3 -g --std=c++17 -fPIC

main.exe: main.o GA.o ../RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.exe: %.o ../RndGen/random.o
	$(CC) -o $@ $^ $(CFLAGS)
%.o : %.cpp
	$(CC) -o $@ -c $<  $(CFLAGS)
clean:
	rm *.exe *.o
remove:
	rm OUTPUT/Indipendent/* OUTPUT/Migration/*

