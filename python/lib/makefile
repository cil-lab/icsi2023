CC = gcc

default: libicsi2023.a

libicsi2023.a: icsi2023.o
	ar rcs $@ $^

icsi2023.o: icsi2023.c icsi2023.h
	$(CC) -c $< -fPIC

clean:
	rm *.o *.a
