#include "mex.h"
#include "matrix.h"
#include <fftw3.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <windows.h>
#include <process.h>
#include <complex.h>
#include <omp.h>

#define SettingField    43
#define SatelliteCount  40

#define CHANNEL_MAX     12

#define bitValue(reg, N) (int)((reg & (0x01 << (N - 1))) == 0 ?  -1: 1) 

static char legendre[10223]={0};

enum SettingIndex{
    msToProcess                 = 0,
    numberOfChannels,
    skipNumberOfBytes,
    fileName,
    dataType,
    fileType,
    IF,
    samplingFreq,
    codeFreqBasis,
    codeLength,
    skipAcquisition,
    acqSatelliteList,
    acqSearchBand,
    acqNonCohTime,
    acqThreshold,
    acqStep,
    resamplingThreshold,
    resamplingflag,
    dllDampingRatio,
    dllNoiseBandwidth,
    dllCorrelatorSpacing,
    pllDampingRatio,
    pllNoiseBandwidth,
    intTime,
    navSolPeriod,
    elevationMask,
    useTropCorr,
    truePosition,
    plotTracking,
    plotAcquisition,
    plotNavigation,
    c,
    startOffset,
    CNo
};


enum ResultIndex{
    RstIdx_carrFreq         = 0,
    RstIdx_codePhase           ,
    RstIdx_peakMetric          
};

enum ProcessStatus{
    None = 0                    ,
    Success = 1                 ,
    MemoryOut_Data = 2          ,
};

typedef struct _AcqResult{
    mxArray *a_carrFreq                     ;
    mxArray *a_codePhase                    ;
    mxArray *a_peakMetric                   ;
    double *carrFreq                        ;
    double *codePhase                       ;
    double *peakMetric                      ;

    int ThreadResult                        ;
    int ErrorLoopIdx                        ;
}AcqResult;

typedef struct _Setting{
    double msToProcess;
    double numberOfChannel;
    double skipNumberOfBytes;
    char* fileName;
    char* dataType;
    double fileType;
    double IF;
    double samplingFreq;
    double codeFreqBasis;
    double codeLength;
    double skipAcquisition;
    double *acqSatelliteList;
    double  acqSatelliteLength;
    double acqSearchBand;
    double acqNonCohTime;
    double acqThreshold;
    double acqStep;
    double resamplingThreshold;
    double resamplingFlag;
    double dllDampingRatio;
    double dllNoiseBandwidth;
    double dllCorrelatorSpacing;
    double pllDampingRatio;
    double pllNoiseBandwidth;
    double intTime;
    double navSolPeriod;
    double elevationMask;
    double useTropCorr;
    //double truePosition;
    double plotTracking;
    double plotAcquisition;
    double plotNavigation;
    double c;
    double startOffset;
    double maxPRN                       ;
    double PRNLength                    ;
    //double CNo;
}Setting;

typedef struct _ThreadParameter{
    Setting setting                     ;
    AcqResult *result                   ;
    int stx                             ;
    int edx                             ;
    double* recvCode                    ;
}ThreadParameter;

typedef struct _FileInfo{
    const signed char * dataPtr        ;
    unsigned long Length               ;
}DataFileInfo;

int initSetting(Setting *setting, const mxArray *settingArray);
void calcPRNAcquisition(Setting setting, AcqResult* result, int Index, double *recvCode);
double *zeros(int length);
double *inf(int length);
AcqResult initResult(Setting setting);
double *generateCAcode(int PRN);
void gen_legendre_sequence();
void dispose(Setting setting);
void makeCATable(Setting setting, int PRN, fftw_complex *output, int N);
AcqResult acquisition_Norm(Setting setting, DataFileInfo dataInfo);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    Setting setting;
    DataFileInfo dataInfo;
    AcqResult result;
    int maxPRN = 0;
    clock_t start, end;
    start = clock();

    const char *fieldNames[31];

    clock_t start_block, end_block;

    if(nrhs!=5)
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs",
                "Three input required.");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs",
                "Too many output arguments.");
    else if(!mxIsStruct(prhs[0]))
        mexErrMsgIdAndTxt( "MATLAB:inputNotStruct",
                "Input[2] must be a structure.");
    if(initSetting(&setting, prhs[0])){
    }
    dataInfo.dataPtr = mxGetPr(prhs[1]);
    dataInfo.Length = mxGetScalar(prhs[2]);
    setting.PRNLength = mxGetScalar(prhs[3]);
    setting.maxPRN = mxGetScalar(prhs[4]);

    for(int i = 0; i < 3; i++){
        fieldNames[i] = (char*)mxMalloc(20);
    }

    memcpy(fieldNames[RstIdx_carrFreq], "carrFreq", sizeof("carrFreq"));
    memcpy(fieldNames[RstIdx_codePhase], "codePhase", sizeof("codePhase"));
    memcpy(fieldNames[RstIdx_peakMetric], "peakMetric", sizeof("peakMetric"));

    plhs[0] = mxCreateStructMatrix(1, 1, 3, fieldNames);

    mxFree(fieldNames[0]);
    mxFree(fieldNames[1]);
    mxFree(fieldNames[2]);


#ifdef DEBUG_MODE
    // mexPrintf("\n\n---------------------Display Channel Data---------------------\n\n");
    // for(int i = 0; i < 12; i++){
    //     mexPrintf("channel[%d].PRN: %f\n", i, chnl[i].PRN);
    //     mexPrintf("channel[%d].acquiredFreq: %f\n", i, chnl[i].acquiredFreq);
    //     mexPrintf("channel[%d].codePhase: %f\n", i, chnl[i].codePhase);
    //     mexPrintf("channel[%d].status: %s\n", i, chnl[i].status);
    // }
    mexPrintf("\n\n---------------------Display Setting Data---------------------\n\n");
    mexPrintf("fileName: %s\n", setting.fileName);
    mexPrintf("fileType: %f\n", setting.fileType);
    mexPrintf("dataType: %s\n", setting.dataType);
    mexPrintf("IF: %f\n", setting.IF);
    mexPrintf("samplingFreq: %f\n", setting.samplingFreq);
    mexPrintf("acqThreshold: %f\n", setting.acqThreshold);
    mexPrintf("skipNumberOfBytes: %f\n", setting.skipNumberOfBytes);
    mexPrintf("msToProcess: %f\n", setting.msToProcess);    
    for(int i = 0; i < (int)setting.acqSatelliteLength; i++){
        mexPrintf("acqSatelliteList[%d]: %f\n", i + 1, setting.acqSatelliteList[i]);
    }
    mexPrintf("numberOfChannel: %f\n", setting.numberOfChannel);
    mexPrintf("codeFreqBasis: %f\n", setting.codeFreqBasis);

    mexPrintf("codeLength: %f\n", setting.codeLength);
    mexPrintf("skipAcquisition: %f\n", setting.skipAcquisition);
    mexPrintf("acqSearchBand: %f\n", setting.acqSearchBand);
    mexPrintf("acqNonCohTime: %f\n", setting.acqNonCohTime);
    mexPrintf("acqStep: %f\n", setting.acqStep);
    mexPrintf("resamplingThreshold: %f\n", setting.resamplingThreshold);
    mexPrintf("resamplingFlag: %f\n", setting.resamplingFlag);
    mexPrintf("dllDampingRatio: %f\n", setting.dllDampingRatio);
    mexPrintf("dllNoiseBandwidth: %f\n", setting.dllNoiseBandwidth);
    mexPrintf("dllCorrelatorSpacing: %f\n", setting.dllCorrelatorSpacing);
    mexPrintf("pllDampingRatio: %f\n", setting.pllDampingRatio);
    mexPrintf("pllNoiseBandwidth: %f\n", setting.pllNoiseBandwidth);
    mexPrintf("intTime: %f\n", setting.intTime);
    mexPrintf("navSolPeriod: %f\n", setting.navSolPeriod);
    mexPrintf("elevationMask: %f\n", setting.elevationMask);
    mexPrintf("useTropCorr: %f\n", setting.useTropCorr);
    mexPrintf("plotTracking: %f\n", setting.plotTracking);
    mexPrintf("plotAcquisition: %f\n", setting.plotAcquisition);
    mexPrintf("plotNavigation: %f\n", setting.plotNavigation);
    mexPrintf("C: %f\n", setting.c);
    mexPrintf("startOffset: %f\n", setting.c);
#endif
    result = acquisition_Norm(setting, dataInfo);
    dispose(setting);
    
    mxSetFieldByNumber(plhs[0], 0, RstIdx_carrFreq, result.a_carrFreq);
    mxSetFieldByNumber(plhs[0], 0, RstIdx_codePhase, result.a_codePhase);
    mxSetFieldByNumber(plhs[0], 0, RstIdx_peakMetric, result.a_peakMetric);

    end = clock();
    mexPrintf("\nSignal Acquisition: %f ms\n", (double)(end - start) / CLOCKS_PER_SEC); 
    mexPrintf("-----------------------------------------------\n");
}

int initSetting(Setting *setting, const mxArray *settingArray){
    int        ifield, nfields;
    mxArray * tmp;
    size_t buflen;

    nfields = mxGetNumberOfFields(settingArray);

    tmp = mxGetFieldByNumber(settingArray,0, msToProcess);
    setting->msToProcess = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, numberOfChannels);
    setting->numberOfChannel = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, skipNumberOfBytes);
    setting->skipNumberOfBytes = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, fileName);
    buflen = mxGetN(tmp) + 1;
    setting->fileName = mxMalloc(buflen);
    mxGetString(tmp, setting->fileName, (mwSize) buflen);

    tmp = mxGetFieldByNumber(settingArray,0, fileType);
    setting->fileType = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, dataType);
    buflen = mxGetN(tmp) + 1;
    setting->dataType = mxMalloc(buflen);
    mxGetString(tmp, setting->dataType, (mwSize) buflen);

    tmp = mxGetFieldByNumber(settingArray,0, IF);
    setting->IF = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, samplingFreq);
    setting->samplingFreq = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, codeFreqBasis);
    setting->codeFreqBasis = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, codeLength);
    setting->codeLength = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, skipAcquisition);
    setting->skipAcquisition = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, acqSatelliteList);
    setting->acqSatelliteLength = mxGetN(tmp);
    setting->acqSatelliteList = mxGetPr(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, acqSearchBand);
    setting->acqSearchBand = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, acqNonCohTime);
    setting->acqNonCohTime = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, acqThreshold);
    setting->acqThreshold = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, acqStep);
    setting->acqStep = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, resamplingThreshold);
    setting->resamplingThreshold = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, resamplingflag);
    setting->resamplingFlag = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, dllDampingRatio);
    setting->dllDampingRatio = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, dllNoiseBandwidth);
    setting->dllNoiseBandwidth = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, dllCorrelatorSpacing);
    setting->dllCorrelatorSpacing = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, pllDampingRatio);
    setting->pllDampingRatio = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, pllNoiseBandwidth);
    setting->pllNoiseBandwidth = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, intTime);
    setting->intTime = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, numberOfChannels);
    setting->numberOfChannel = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, navSolPeriod);
    setting->navSolPeriod = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, elevationMask);
    setting->elevationMask = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, useTropCorr);
    setting->useTropCorr = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, plotTracking);
    setting->plotTracking = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, plotAcquisition);
    setting->plotAcquisition = mxGetScalar(tmp);
    
    tmp = mxGetFieldByNumber(settingArray,0, plotNavigation);
    setting->plotNavigation = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, c);
    setting->c = mxGetScalar(tmp);

    tmp = mxGetFieldByNumber(settingArray,0, startOffset);
    setting->startOffset = mxGetScalar(tmp);

    return 0;
}

double *zeros(int length){
    double *ptr = (double*)mxMalloc(sizeof(double) * length);
    double *start = ptr;
    for(int i = 0; i < length; i++){
        *(start++) = 0;
    }
    return ptr;
}

double *inf(int length){
    double *ptr = (double*)mxMalloc(sizeof(double) * length);
    double *start = ptr;
    for(int i = 0; i < length; i++){
        *(start++) = INT_MAX;
    }
    return ptr;
}

AcqResult initResult(Setting setting){
    AcqResult result;
    int maxPRN = 0;

    result.a_carrFreq = mxCreateDoubleMatrix(1, setting.maxPRN, 0);
    result.carrFreq   = mxGetPr(result.a_carrFreq);
    result.a_codePhase = mxCreateDoubleMatrix(1, setting.maxPRN, 0);
    result.codePhase   = mxGetPr(result.a_codePhase);
    result.a_peakMetric = mxCreateDoubleMatrix(1, setting.maxPRN, 0);
    result.peakMetric   = mxGetPr(result.a_peakMetric);

    return result;
}

AcqResult acquisition_Norm(Setting setting, DataFileInfo dataInfo){
    clock_t start, end;
    clock_t start_block, end_block;
    struct timespec start_2;
    AcqResult result = initResult(setting);
    const signed char* dataPtr = dataInfo.dataPtr;
    
    int samplesPerCode = (int)round(setting.samplingFreq / (setting.codeFreqBasis / setting.codeLength));

    fftw_complex *temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 2 * samplesPerCode);
    fftw_complex *recvCode_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 2 * samplesPerCode);
    fftw_complex *localOscil = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 2  * samplesPerCode);
    fftw_complex *localOscil_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 2  * samplesPerCode);
    fftw_complex *corr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * 2  * samplesPerCode);
    fftw_plan p1,p2,p3;

    double ts = 1 / setting.samplingFreq;

    double sigPower = 0;

    int initFreq = setting.IF - setting.acqSearchBand;

    const signed char* dataIdx = dataPtr;


    int frequencyBinIndex = -1;
    double MaxPow = -1;
    int MaxPowIdx = -1;
    double maxVal = 0;
    int maxidx = -1;
    
    int numberOfFreqBins = (int)(setting.acqSearchBand * 2 / setting.acqStep) + 1;  
    double* coarseFreqBin = (double *)mxMalloc(sizeof(double) * numberOfFreqBins);
    
    int fineStep = 25;
    int numOfFineBins = (int)((setting.acqStep / fineStep) * 2 + 1);

    p1 = fftw_plan_dft_1d(2*samplesPerCode, localOscil, localOscil_fft, FFTW_FORWARD, FFTW_MEASURE);
    p2 = fftw_plan_dft_1d(2*samplesPerCode, temp, recvCode_fft, FFTW_FORWARD, FFTW_MEASURE);
    p3 = fftw_plan_dft_1d(2*samplesPerCode, temp, corr, FFTW_BACKWARD, FFTW_MEASURE);

    static const int k = 0; // Just the diagonal; 0 super-diagonal bands
    static const double alpha = 1.0;
    static const int lda = 1;
    static const int incx = 1;
    static const double beta = 0.0;
    static const int incy = 1;
    static const int n = 36000;

    double sigSampling[72000];

    for(int sigSample = 0; sigSample < 2*samplesPerCode; sigSample++){
        sigSampling[sigSample] = 2 * M_PI * ts * sigSample;
    }

    omp_lock_t writelock;

    omp_init_lock(&writelock);
    int PRNIdx = 0;
    int MaxIdx = setting.PRNLength;
    
    #pragma omp parallel for
    for(PRNIdx = 0; PRNIdx < MaxIdx; PRNIdx++){
        // start = clock();
        int PRN = setting.acqSatelliteList[PRNIdx];
        double *CACode= (double*)malloc(sizeof(double) * 72000);
        double buffer[72000];
        
        #pragma omp critical
        {
            makeCATable(setting, PRN, localOscil, 2 * samplesPerCode);    
            fftw_execute(p1);
            for(int i = 0; i < 2 * samplesPerCode; i++){
                CACode[2* i] = localOscil_fft[i][0];
                CACode[2* i + 1] = localOscil_fft[i][1];
            }
        }

        double maxMag = 0;
        int maxidx = -1;
        int freq = -1;
        int freqIdx = 0;
        for(freqIdx = 0; freqIdx < numberOfFreqBins; freqIdx++){
            double results[36000] = {0, };
            coarseFreqBin[freqIdx] = setting.IF + setting.acqSearchBand - setting.acqStep * freqIdx;
            for(int nonCohIdx = 0; nonCohIdx < setting.acqNonCohTime; nonCohIdx++){
                const signed char *sigIdx = dataIdx + nonCohIdx * 2 * samplesPerCode;
                int sigSample = 0;
                for(sigSample =  0; sigSample < 2*samplesPerCode; sigSample++){
                    double sigCarrReal = cos(coarseFreqBin[freqIdx] * sigSampling[sigSample]);
                    double sigCarrImag = sin(coarseFreqBin[freqIdx] * sigSampling[sigSample]); 
                    double localCA_real = CACode[2* sigSample];
                    double localCA_img = CACode[2* sigSample + 1];
                    signed char a = sigIdx[2*sigSample];
                    signed char b = sigIdx[2*sigSample + 1];
                    buffer[2*sigSample] = (double)a * sigCarrReal + (double)b * sigCarrImag;
                    buffer[2*sigSample + 1] = (double)b * sigCarrReal - (double)a * sigCarrImag;
                }
                #pragma omp critical
                {
                    for(sigSample = 0; sigSample < 2*samplesPerCode; sigSample++){
                        temp[sigSample][0] = buffer[2*sigSample];
                        temp[sigSample][1] = buffer[2*sigSample + 1];
                    }
                    fftw_execute(p2);
                                                     
                    for(sigSample =  0; sigSample < 2*samplesPerCode; sigSample++){
                        temp[sigSample][0] = CACode[2*sigSample] * recvCode_fft[sigSample][0] + CACode[2*sigSample + 1]*recvCode_fft[sigSample][1];
                        temp[sigSample][1] = CACode[2*sigSample] * recvCode_fft[sigSample][1] - CACode[2*sigSample + 1] * recvCode_fft[sigSample][0];
                    } 
                    fftw_execute(p3);
                    for(int i = 0; i < 2 * samplesPerCode; i++){
                        buffer[2* i] = corr[i][0];
                        buffer[2* i + 1] = corr[i][1];
                    }
                }

                for(sigSample = 0; sigSample <  2*samplesPerCode; sigSample++){
                    results[sigSample] += sqrt(buffer[2* sigSample]/100000000 * buffer[2* sigSample]/100000000 + buffer[2* sigSample + 1]/100000000 * buffer[2* sigSample + 1]/100000000);
                    if(maxMag <= results[sigSample]){
                        maxMag = results[sigSample];
                        maxidx = sigSample;
                        freq = freqIdx;
                    }
                }
            }    
        } 

        double threshold = result.peakMetric[PRN - 1] = maxMag;
        double fineInitFreq = setting.IF + setting.acqSearchBand - setting.acqStep * freq;
        
        if(maxidx + samplesPerCode > 2 * samplesPerCode)
            maxidx = maxidx - samplesPerCode;

        end = clock();

        if(threshold > 80){ //setting.acqThreshold){
            start = clock();
            double maxMag = 0;
            double maxFineFreq = 0;
            
            double *caCode = generateCAcode(PRN);

            double maxFinePow = 0;
            double minFineBin = 0;
            int fineFreqIdx = 0;
            for(fineFreqIdx = 0; fineFreqIdx < numOfFineBins; fineFreqIdx++){

                const signed char *sigIdx = dataIdx + 2 * maxidx;
                double fineFreq = fineInitFreq + setting.acqStep / 2 - fineStep * fineFreqIdx;

                double basebandSigReal_1msSum[40] = {0, };
                double basebandSigImage_1msSum[40] = {0, };

                for(int i = 0; i < 40; i++){
                    for(int j = 0; j < samplesPerCode; j++){
                        int idx = i * samplesPerCode + j;
                        double tt = (ts * idx) / (1/setting.codeFreqBasis);
                        int caIdx = (int)(floor(tt)) % ((int)setting.codeLength);

                        double sigCarrReal = cos(fineFreq * 2 * M_PI * ts * idx);
                        double sigCarrImag = sin(fineFreq * 2 * M_PI * ts * idx); 

                        signed char a = sigIdx[2*idx];
                        signed char b = sigIdx[2*idx + 1];

                        double basebandSigReal = caCode[caIdx] * ((double)a * sigCarrReal + (double)b * sigCarrImag);
                        double basebandSigImage = caCode[caIdx] * ((double)b * sigCarrReal - (double)a * sigCarrImag);

                        basebandSigReal_1msSum[i] += basebandSigReal;
                        basebandSigImage_1msSum[i] += basebandSigImage;
                    }
                }
                for(int i = 0; i < 20; i++){
                    double comReal = 0;
                    double comImag = 0;
                    for(int  j = 0; j < 20; j++){
                        comReal += basebandSigReal_1msSum[i + j];
                        comImag += basebandSigImage_1msSum[i + j];
                    }
                    double pow = sqrt(comReal * comReal + comImag * comImag);
                    if(pow > maxFinePow){
                        maxFinePow = pow;
                        maxFineFreq = fineFreq;
                    }
                }
            }
            end = clock();
            
            result.carrFreq[PRN - 1] = maxFineFreq;
            result.codePhase[PRN - 1] = maxidx + 1;
    
            // Downsampling Recovery - not used
        }        
    }
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_free(temp);
    fftw_free(recvCode_fft);
    fftw_free(localOscil_fft);
    fftw_free(corr);
    
    return result;
}

void makeCATable(Setting setting, int PRN, fftw_complex *output, int N){
    int samplesPerCode = (int)round(setting.samplingFreq / (setting.codeFreqBasis / setting.codeLength));

    double ts = 1 / setting.samplingFreq;
    double tc = 1 / setting.codeFreqBasis;
    int i = 0;
    double *caCode = generateCAcode(PRN);
    
    for(i = 1; i < N + 1; i++){
        if(i < samplesPerCode){
            output[i - 1][0] = caCode[((int)ceil(ts * i / tc) - 1)];
            output[i - 1][1] = 0;
        }
        else{
            output[i - 1][0] = 0;
            output[i - 1][1] = 0;
        }
    }

    output[0][0] = caCode[0];
    output[samplesPerCode - 1][0] = caCode[(int)setting.codeLength - 1];

}

void dispose(Setting setting){
    // mxFree(setting.fileName);
    // mxFree(setting.path);
    // mxFree(setting.dataType);
    // mxFree(setting.SYS_type);
}

double *generateCAcode(int PRN){
    const static short g2s[]={
        5,   6,   7,   8,  17,  18, 139, 140, 141, 251,
        252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
        473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
        861, 862,
                  145, 175,  52,  21, 237, 235, 886, 657,
        634, 762, 355, 1012, 176, 603, 130, 359, 595, 68,
        386
    };
    
    double *CACode = (double *)mxMalloc(sizeof(double) * 1023);

    short g2shift = g2s[PRN - 1];
    double g1[1023];
    double g2[1023];
    int g1_reg = 0x00;
    int g2_reg = 0x00;
    int g1_saveBit = -1;
    int g2_saveBit = -1;
    int i = 0;
    for(i = 0; i < 1023; i++)
    {
        g1[i] = bitValue(g1_reg, 10);
        g2[i] = bitValue(g2_reg, 10);
        g1_saveBit = (bitValue(g1_reg, 3) * bitValue(g1_reg, 10)) == 1 ? 1 : 0;
        g2_saveBit = (bitValue(g2_reg, 2) * bitValue(g2_reg, 3) * bitValue(g2_reg, 6) * bitValue(g2_reg, 8) * bitValue(g2_reg, 9) * bitValue(g2_reg, 10)) == 1 ? 1 : 0;;
        g1_reg = (g1_reg << 1) | g1_saveBit;
        g2_reg = (g2_reg << 1) | g2_saveBit;
    }
    
    for(i = 0; i < 1023; i++)
    {
        int g2_idx = (i + 1023 - g2shift) % 1023;
        CACode[i] = -g1[i] * g2[g2_idx];
    }

    return CACode;
}