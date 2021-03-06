# FCC
[![Actions Status](https://github.com/FrancoisGrondin/fcc/workflows/tests/badge.svg)](https://github.com/FrancoisGrondin/fcc/actions)

Fast Cross-Correlation

## Compile the library

Create a directory for binary files and call cmake:

```
mkdir bin
cd bin
cmake ../
```

Compile the library:

```
make
```

## Run a demo

Suppose you have a wave file called `audio.wav` with multiple channels, a sample rate of 16000 samples/sec, and a sample format of signed 16-bit, you can run the demo as follows:

To run GCC, and store results to a csv file called `tdoas_gcc.csv`:

```
./tdoa -i audio.wav -m gcc -o tdoas_gcc.csv
```

To run FCC, and store results to a csv file called `tdoas_fcc.csv`:

```
./tdoa -i audio.wav -m fcc -o tdoas_fcc.csv
```

## Measure performances

The TDoAs are computed with a sequence of 4 operations:

1) stft: Compute the Short-Time Fourier Transform (STFT) for each microphone
2) scmphat: Compute the Spatial Correlation Matrix and perform Phase Transform
3) gcc: Perform Generalized Cross-Correlation
-or-
3) fcc: Perform Fast Cross-Correlation
4) qi: Compute quadratic interpolation

To measure the execution time of all operations with two microphones for 1,000,000 iterations, we can run the following command:

```
./speed -c 2 -r 1000000
```
