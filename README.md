# TacVar

TacVar is a toolset for tackling variability in parallel running time measurement. It includes several tools to measure, quantify and filter unstable timing readings in multi-core processors and accelerators.

Here is core tools in TacVar:
- ParTES: Measuring the minimal 'measurable' running time of your server under different run-time environment.
- Vkern: Obeserving the deviation and error of your performance measurement results caused by timing fluctuation.
- Filter: Filtering a noisy running time distribution for better measurement accuracy.

Related publications:

Q. Liao and J. Lin, "TacVar: Tackling Variability in Short-Interval Timing Measurements on X86 Processors," 2024 IEEE 24th International Symposium on Cluster, Cloud and Internet Computing (CCGrid), Philadelphia, PA, USA, 2024, pp. 496-506, doi: 10.1109/CCGrid59990.2024.00062. (Best paper award)

Q. Liao, S. Zuo, Y. Wang and J. Lin. Timing Method and Evaluation Metrics for CPU Performance Variation Detections [J]. Chinese Journal of Computers，2024，47：456-472, doi: 10.11897/SP.J.1016.2024.00456. ([Download Paper](http://cjc.ict.ac.cn/online/onlinepaper/lqc-2024229165011.pdf))

LIAO Qiucheng, ZHOU Yang, LIN Xinhua. Metrics and Tools for Evaluating the Deviation in Parallel Timing[J]. Computer Science, 2025, 52(5): 41-49. doi: 10.11896/jsjkx.241200053 ([Download Paper](https://www.jsjkx.com/CN/article/openArticlePDF.jsp?id=23149))







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