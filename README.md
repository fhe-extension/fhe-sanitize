tfhe-sanitize/tfhe-sanitize is a C implementation of the TFHE bootstrapping with an additional sanitization algorithm.

In order to compile the code, you will need:
 - a recent version of a C compiler,
 - FFTW3 installed on your machine.

Depending on the directory, you might need to change the path and the FLAGS in the Makefile.

We provide unitary tests for each function in the Test directory. 

The sanitize algorithm is implemented in the fft.c file. You can test it using ./Tests/Test_bootstrapping_fft
