# GNSS L1C/A SDR

OS(Open Source) GNSS(Global Navigation Satellite System) L1C/A(L1-Coarse Acquisition) SDR(Software-Defined Radio) CMEX Edit

**The program code will be committed later.**

## Repository Info

OS GNSS L1C/A SDR CMEX Edit is SDR based on [OS GNSS L1C/A SDR](https://github.com/gnsscusdr/CU-SDR-Collection).

OS GNSS L1C/A SDR is configured with 5 main block, initialization, **acquisition**, **tracking**, navigation, civil processing. The project aims to accelerate acquisition and tracking block through CMEX. Table below comparing the tack time performance of the following configurations.

||Reference (s)|CMEX (s)|Proportion (%)|
|-|-|-|-|
|Acquisition|16.296|9.916|60.849|
|Tracking|395.276|35.788|9.042|

## Development Environment
### Software
- MATLAB
- C(MATLAB CMEX)

### Hardware

- MATLAB Benchmark Result

![](./Image/MATLAB%20Bench%201.png)
![](./Image/MATLAB%20Bench%202.png)

## CMEX build

- [cmex_acquisition.c](./Non-Assist/Acquisition/cmex_acquisition.c)
  - acquisition based on CMEX
  - build
    ```console
     mex -v -lfftw3 -lmwblas -I'./Include/fftw' './Non-Assist/Acquisition/cmex_acquisition.c' -L'./Include/fftw' libfftw3.lib -outdir './Non-Assist/Acquisition/mex' COMPFLAGS="/openmp $COMPFLAGS"
    ```
- [tracking_cmex.c](./Non-Assist/Tracking/tracking_cmex.c)
  - tracking based on CMEX
  - build
    ```console
    mex -v './Non-Assist/Tracking/tracking_cmex.c' -outdir './Non-Assist/Tracking/mex'
    ```
