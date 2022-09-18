CC = gcc
CFLAGS = -std=c99 -Wall -Ofast
TARGET = Tests/Test_uniform64 Tests/Test_binary Tests/Test_double Tests/Test_noise Tests/Test_gaussian Tests/Test_ginv Tests/Test_lwe Tests/Test_tlwe Tests/Test_tgsw Tests/Test_bootstrapping Tests/Test_fft Tests/Test_sanitize Tests/Test_bootstrapping_fft
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all : $(TARGET)

Tests/Test_uniform64: Tests/Test_uniform64.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_binary: Tests/Test_binary.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_double: Tests/Test_double.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_noise: Tests/Test_noise.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_gaussian: Tests/Test_gaussian.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I fftw3 -L fftw3 -lfftw3l -lm

Tests/Test_ginv: Tests/Test_ginv.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_lwe: Tests/Test_lwe.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_tlwe: Tests/Test_tlwe.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_tgsw: Tests/Test_tgsw.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_bootstrapping: Tests/Test_bootstrapping.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_fft: Tests/Test_fft.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_sanitize: Tests/Test_sanitize.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

Tests/Test_bootstrapping_fft: Tests/Test_bootstrapping_fft.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS) -I/usr/local/include -L/usr/local/lib -lfftw3l -lm

clean:
	rm -f *.o core