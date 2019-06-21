default:euler

euler.o: euler.c $(HEADERS)
	gcc -c euler.c -o euler.o -std=c99
euler: euler.o
	gcc euler.o -o euler -lblas -llapack  -lm
clean:
	-rm -f euler.o
	-rm -f euler
