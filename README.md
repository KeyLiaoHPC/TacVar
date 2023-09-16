# FilTAndVKern

## 1 Introduction

VKern generates run time variations which can be configured by given statistical distributions. By using different timing methods to measure the run time of the specific distribution, one can investigate the instability caused by timing fluctuations in parallel short-interval timing. Furthermore, based on the timing fluctuation distribution sample, FilT can help estimate the real distribution of timing results by decoupling the timing fluctuation from measurement results.

This repository is a preliminary code base for validating the methodology. For broaden the usage of the method, it will soon be integrated into our open-source instrumentation-based profiling tool PerfHound.

## 2 A simple use case

### 2.1 Preparing the test environment

### 2.2 Compiling source codes

### 2.3 Running VKern and visulizing timing fluctuations

### 2.4 Sampling timing fluctuations

### 2.5 Filtering the timing fluctuation with FilT

``` bash

$ ./filt.x --help
Usage: filt.x [OPTION...]

  -g, --hz-ns=FREQ           Tick per ns of the clock
  -l, --plow=PROB            Lowest threshold of possibility of a data bin
  -m, --met-file=FILE        Input measurement file
  -n, --nsamp=NSAMP          Number of samples in each optimization step
  -s, --sample-file=FILE     Input timing fluctuation file
  -w, --width=TIME           The least interval of a time bin (ns).
  -x, --cut-x=PROB           Cut the highest probability of met array
  -y, --cut-y=PROB           Cut the highest probability of timing fluctuation
                             array
  -?, --help                 Give this help list
      --usage                Give a short usage message

```