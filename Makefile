app: qs.o lgz.o heap.o
	gcc -o app qs.o lgz.o heap.o -lm

qs.o: qs.c
	gcc -c qs.c

lgz.o: lgz.c
	gcc -c lgz.c

heap.o: heap.c
	gcc -c heap.c

clean:
	rm app *.o
