CC = gcc
ARCH = arch -arm64

OPENSSL_PREFIX = /opt/homebrew/opt/openssl@3
FFTW_PREFIX = /opt/homebrew/opt/fftw

CFLAGS =-std=c99 -Wall -Ofast -I$(OPENSSL_PREFIX)/include -I$(FFTW_PREFIX)/include
#CFLAGS =-std=c99 -Wall -O3 -pg -g -Wall -Wextra -fsanitize=address,undefined -I$(OPENSSL_PREFIX)/include -I$(FFTW_PREFIX)/include

LDFLAGS = -L$(OPENSSL_PREFIX)/lib -L$(FFTW_PREFIX)/lib -lfftw3 -lm -lssl -lcrypto

SRC = $(filter-out $(wildcard Tests/*.c), $(wildcard *.c))
OBJ = $(SRC:.c=.o)


# Liste des tests à compiler
TESTS = \
	Test_uniform64 \
	Test_binary \
	Test_double \
	Test_noise \
	Test_gaussian \
	Test_small_gaussian \
	Test_ginv \
	Test_lwe \
	Test_tlwe \
 	Test_product \
	Test_tlwe_fft \
	Test_tgsw \
	Test_tgsw_fft \
	Test_fft \
	Test_bootstrapping_fft


TARGETS = $(addprefix Tests/, $(TESTS))

all: $(TARGETS)


Tests/%: Tests/%.c $(OBJ)
	$(ARCH) $(CC) $(CFLAGS) $< $(OBJ) -o $@ $(LDFLAGS)

# Compilation des .o
%.o: %.c
	$(ARCH) $(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o core $(TARGETS)
