# FilterLargePolys
Command line tools to filter out {0, 1} polynomials with large minima on the unit circle.

## Contents

*filter101seq.c* - a C program to filter binary sequences with large minimal absolute values of their [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform).

**Command line**

`filter101seq N l t`

**Parameters**

- *N* - FFT size (no. of sample points). Must be a positive integer. Small powers of 2 are preferred (2,4,8,16,32,64,128,256 etc.)
- *l* - length of a binary sequence. Must be an integer between 4 and 64.
- *t* - minimal absolute FFT value allowed (positive real, double type, i.e. 1.234)

**Description** the program generates all *{0, 1}* binary sequences of length *l*, avoiding reversals and [palyndromes](https://en.wikipedia.org/wiki/Palindrome/). This is done using and asymetric version of [binary Gray code](https://en.wikipedia.org/wiki/Gray_code/) and evaluates their FFT minimal absolute value, say, *m* (that is, the minimal absolute value of the corresponding {0, 1} polynomial evaluated at every *N*-th root of unity). If that value *m* is less than the allowed threshold value *t*, the sequence is simply discarded. If *m>=t*, they are printed to *stdin* together with the minimal FFT absolute value *m*, which serves as an approximation of the true minima of the corresponding *{0, 1}* polynomial on the unit circle *|z|=1*.

**Example**

```
user$./filter101seq 64 10 0.85
#Samples N=10, length=64, threshold t=0.900
#Sequences and their FFT minima:
1000101101   0.934
1001000111   1.000
1011110011   0.916
#Total: 120
#Left:  3
#Time:  0.011508 sec.
```
---

*doFFT.c* - a simple C program to evaluate minimal absolute FFT value of a given {0,1} binary sequence.

**Command line**

`doFFT N seq101`

**Parameters**
- *N* - FFT size (no. of sample points). Must be a positive integer. Small powers of 2 are preffered (2,4,8,16,32,64,128,256 etc.)
- *seq101* - a string consisting of 0 or 1, without separating symbols like spaces etc.

**Description**

This program simply takes a *{0,1}*-string and returns the minimal absolute value *m* of its Fourier transform of size *N*. This corresponds to the evaluation of the corresponding *{0,1}*-polynomial whose coefficients match the given string *seq101* on the unit circle at each of the *N*-th roots of unity. As *N* increases, *m* approximates the true minima of that polynomial on the unit circle *|z|=1*.

**Example**

```
user$ ./doFFT 512 1000101101

#FFT Minimal amplitude:
0.9102364331
```
