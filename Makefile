
CC = gcc
CFLAGS = -std=c99 -Wall -Ofast
LIBS = -I /usr/local/include -L /usr/local/lib -lfftw3 -lm
TARGET = Tests/Test_uniform64 Tests/Test_binary Tests/Test_double Tests/Test_noise Tests/Test_gaussian Tests/Test_ginv Tests/Test_lwe Tests/Test_tlwe Tests/Test_tlwe_fft Tests/Test_tgsw Tests/Test_bootstrapping Tests/Test_fft Tests/Test_sanitize Tests/Test_sanitize_pkc Tests/Test_bootstrapping_fft
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all : $(TARGET)

Tests/Test_uniform64: Tests/Test_uniform64.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_binary: Tests/Test_binary.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_double: Tests/Test_double.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_noise: Tests/Test_noise.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_gaussian: Tests/Test_gaussian.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_ginv: Tests/Test_ginv.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_lwe: Tests/Test_lwe.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_tlwe: Tests/Test_tlwe.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_tlwe_fft: Tests/Test_tlwe_fft.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_tgsw: Tests/Test_tgsw.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_bootstrapping: Tests/Test_bootstrapping.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_fft: Tests/Test_fft.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_sanitize: Tests/Test_sanitize.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_sanitize_pkc: Tests/Test_sanitize_pkc.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

Tests/Test_bootstrapping_fft: Tests/Test_bootstrapping_fft.c $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS) $(LIBS)

clean:
	rm -f *.o core

