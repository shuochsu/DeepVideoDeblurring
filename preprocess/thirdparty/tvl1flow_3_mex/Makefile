CC=gcc-4.8
CFLAGS=-std=c99 -g 
#OMPFLAGS=-fopenmp

tvl1flow: main.c tvl1flow_lib.c bicubic_interpolation.c mask.c zoom.c iio.o
	$(CC) $(CFLAGS) $(OMPFLAGS) -o tvl1flow main.c iio.o -lpng -ljpeg -ltiff -lm

iio.o: iio.c
	$(CC) $(CFLAGS) -DNDEBUG -D_GNU_SOURCE -c iio.c

clean:
	rm -f iio.o main.o tvl1flow
