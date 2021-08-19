app: qs.o lgz.o
	gcc -o app qs.o lgz.o -lm

qs.o: qs.c
	gcc -c qs.c

lgz.o: lgz.c
	gcc -c lgz.c

clean:
	rm app *.o
