# FilterLargePolys
Command line tools to filter out {0, 1} polynomials with large minima on the unit circle.

## Contents

[**filter101seq**](https://github.com/jonas-jankauskas/FilterLargePolys/blob/main/filter101seq.c) - a C program to filter binary sequences with large minimal absolute values of their [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform).

**Command line**

`filter101seq N l t`

**Parameters**

- *N* - FFT size (no. of sample points). Must be a positive integer *>l*. Small powers of 2 are preferred (2,4,8,16,32,64,128,256 etc.) for speed.
- *l* - length of a binary sequence. Must be an integer between 4 and 64.
- *t* - minimal absolute FFT value allowed (positive real, double type, i.e. 1.234)

**Description**

The program generates all *{0, 1}* binary sequences of length *l*, avoiding reversals and [palyndromes](https://en.wikipedia.org/wiki/Palindrome/). This is done using and asymetric version of [binary Gray code](https://en.wikipedia.org/wiki/Gray_code/) and evaluates their FFT minimal absolute value, say, *m* (that is, the minimal absolute value of the corresponding {0, 1} polynomial evaluated at every *N*-th root of unity). If that value *m* is less than the allowed threshold value *t*, the sequence is simply discarded. If *m>=t*, they are printed to *stdin* together with the minimal FFT absolute value *m*, which serves as an approximation of the true minima of the corresponding {0, 1} polynomial on the unit circle |*z*|=1.

**Example**

```
user@user_machine:~/FFT$ ./filter101seq 64 10 0.9
#Samples N=64, length=10, threshold t=0.900
#Sequences and their FFT minima:
1000101101   0.934
1001000111   1.000
1011110011   0.916
#Total: 120
#Left:  3
#Time:  0.011927 sec.
```
---

[**doFFT**](https://github.com/jonas-jankauskas/FilterLargePolys/blob/main/doFFT.c) - a simple C program to evaluate minimal absolute FFT value of a given {0,1} binary sequence.

**Command line**

`doFFT N seq101`

**Parameters**
- *N* - FFT size (no. of sample points). Must be a positive integer. Small powers of 2 are preffered (2,4,8,16,32,64,128,256 etc.)
- *seq101* - a string consisting of 0 or 1, without separating symbols like spaces etc.

**Description**

This program simply takes a {0,1}-string and returns the minimal absolute value *m* of its Fourier transform of size *N*. This corresponds to the evaluation of the corresponding {0,1}-polynomial whose coefficients match the given string *seq101* on the unit circle at each of the *N*-th roots of unity. As *N* increases, *m* approximates the true minima of that polynomial on the unit circle |*z*|=1.

**Example**

```
user@user_machine:~/FFT$ ./doFFT 512 1000101101

#FFT Minimal amplitude:
0.9102364331
```
## Installation

You need [FFTW3 library](https://www.fftw.org/) for the Fast-Fourier Transform and its header files installed on your machine in directories located on default include and link paths of the *gcc* compiler (typically in `/usr/local/`). In many distributions, the library is already installed, but you might still need to obtain header files. This can be done using standard package manager.


For instance, on my *Ubuntu 20.24* machine, *FFTW3* library (called *libfftw3-bin*) comes installed by default. Lets check it:
```
user@user_machine:~/FFT$ sudo apt search libfftw3-bin
[sudo] user password: ********** 
Sorting... Done
Full Text Search... Done
libfftw3-bin/focal,now 3.3.8-2ubuntu1 amd64 [installed,automatic]
  Library for computing Fast Fourier Transforms - Tools
```

The header files are contained in 'libfftw3-dev. It can be installed like this:

`sudo apt install libfftw3-dev`

After downloading the .c files from *GitHub*, compile them with:
```
user@user_machine:~/FFT$ gcc -o filter101seq filter101seq.c -lfftw3 -lm
user@user_machine:~/FFT$ gcc -o doFFT doFFT.c -lfftw3 -lm
```

# Important considerations

## Factors affecting FFT accuracy and speed

- The higher the number of sample points *N* in the FFT, the better the value *m* approximates the true minima of a {0,1} polynomial.
- All FFT computations are internaly done and the value *m* is returned using *C* real *double* type, which is typically limited to about 5-6 correct digits after decimal point.
- As *l* and *N* becomes large, rounding errors tend to accumulate in FFT.
- The larger *l* and *N* get, the slower computations become.
- The FFT libraries typically are optimized for the best performance and reasonable accuracy when *N* is a power of 2.

## Example

The successive FFT minimas *m* of the binary sequence 1000101101 of length *l*=10 approximate the true minima of the polynomial as *N* increases from 2^4 to 2^21=2097152. The improvement in accuracy essentially stops for *N*>2^19=524288 due to the limitations of the *C* real *double* type.

```
user@user_machine:~/FFT$ ./doFFT 16 1000101101
DFT Minimal amplitude: 0.9862355235
user@user_machine:~/FFT$ ./doFFT 32 1000101101
DFT Minimal amplitude: 0.9862355235
user@user_machine:~/FFT$ ./doFFT 64 1000101101
DFT Minimal amplitude: 0.9335162464
user@user_machine:~/FFT$ ./doFFT 128 1000101101
DFT Minimal amplitude: 0.9151236997
user@user_machine:~/FFT$ ./doFFT 256 1000101101
DFT Minimal amplitude: 0.9115784013
user@user_machine:~/FFT$ ./doFFT 512 1000101101
DFT Minimal amplitude: 0.9102364331
user@user_machine:~/FFT$ ./doFFT 1024 1000101101
DFT Minimal amplitude: 0.9101108641
user@user_machine:~/FFT$ ./doFFT 2048 1000101101
DFT Minimal amplitude: 0.9099762003
user@user_machine:~/FFT$ ./doFFT 4096 1000101101
DFT Minimal amplitude: 0.9099762003
user@user_machine:~/FFT$ ./doFFT 128 1000101101
DFT Minimal amplitude: 0.9151236997
user@user_machine:~/FFT$ ./doFFT 256 1000101101
DFT Minimal amplitude: 0.9115784013
user@user_machine:~/FFT$ ./doFFT 512 1000101101
DFT Minimal amplitude: 0.9102364331
user@user_machine:~/FFT$ ./doFFT 1024 1000101101
DFT Minimal amplitude: 0.9101108641
user@user_machine:~/FFT$ ./doFFT 2048 1000101101
DFT Minimal amplitude: 0.9099762003
user@user_machine:~/FFT$ ./doFFT 4096 1000101101
DFT Minimal amplitude: 0.9099762003
user@user_machine:~/FFT$ ./doFFT 8192 1000101101
DFT Minimal amplitude: 0.9099726884
user@user_machine:~/FFT$ ./doFFT 16384 1000101101
DFT Minimal amplitude: 0.9099713548
user@user_machine:~/FFT$ ./doFFT 32768 1000101101
DFT Minimal amplitude: 0.9099712487
user@user_machine:~/FFT$ ./doFFT 64536 1000101101
DFT Minimal amplitude: 0.9099711235
user@user_machine:~/FFT$ ./doFFT 131072 1000101101
DFT Minimal amplitude: 0.9099711086
user@user_machine:~/FFT$ ./doFFT 262144 1000101101
DFT Minimal amplitude: 0.9099711074
user@user_machine:~/FFT$ ./doFFT 524288 1000101101
DFT Minimal amplitude: 0.9099711050
user@user_machine:~/FFT$ ./doFFT 1048576 1000101101
DFT Minimal amplitude: 0.9099711050
user@user_machine:~/FFT$ ./doFFT 2097152 1000101101
DFT Minimal amplitude: 0.9099711050
```
Due to the listed factors, one typically gets a lot of "false positive candidates", with true minima actually smaller than the threshold value *t*.

## Workflow with large *l*

When working with the long binary sequences (*l* >= 30), it is important to balance the speed of computation and the accuracy of the reported minima.

- For the initial filtration, set *N* to the smallest power of 2 that is equal to or larger *4l*. Typically,this will be something like *N*=128,256 or 512. Larger *N* means fewer 'false positive' candidates, but at the cost of excessive computation time.

- Keep the threshold value *t* a bit smaller than you actually need (say, 1.95 instead of 2.0) to account for the possible FFT rounding errors.

- Pipe the results into a large *.txt* file. eg.:
`$ ./filter101seq 128 30 1.95 > results.txt`

- After initial coarse filtering is complete, run a *doFFT* on each example in *results.txt* file, each time doubling the *N* value and removing more and more "false positives", until only a manageable small number of examples are left.

- Find the true minimas on the remaining examples with high accuracy (within correct 4-5 decimal digits) by running *doFFT* with a large *N*, say, *N*=4096 etc.

- Use the newly found largest of the *m* values (with a slight negative corection for possible rounding errors) as a threshold value to double-verify your results for the current length *l* and also as an initial threshold for larger *l*.
