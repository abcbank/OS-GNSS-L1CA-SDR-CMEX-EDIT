#include "mex.h"
#include"matrix.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <windows.h>
#include <process.h>

#define CHANNEL_MAX     12

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
    acqSearchStep,
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

enum ChannelIndex{
    PRN = 0                     ,
    acquiredFreq                ,
    codePhase                   ,
    status
};

enum ResultIndex{
    rst_status           = 0,
    rst_PRN                 ,
    rst_absoluteSample      ,
    rst_codeFreq            ,
    rst_carrFreq            ,
    rst_I_P                 ,
    rst_I_E                 ,
    rst_I_L                 ,
    rst_Q_P                 ,
    rst_Q_E                 ,
    rst_Q_L                 , 
    rst_dllDiscr            ,
    rst_dllDiscrFilt        ,
    rst_pllDiscr            ,
    rst_pllDiscrFilt        ,
    rst_remCodePhase        ,
    rst_remCarrPhase        ,
    rst_CNo
};

enum ProcessStatus{
    None = 0,
    Success = 1,
    MemoryOut_Data = 2,
};

typedef struct _Result{
    
    mxArray *status                     ;
    mxArray *PRN                        ;
    mxArray *absoluteSample             ;
    mxArray *codeFreq                   ;
    mxArray *carrFreq                   ;
    mxArray *I_P                        ;
    mxArray *I_E                        ;
    mxArray *I_L                        ;
    mxArray *Q_P                        ;
    mxArray *Q_E                        ;
    mxArray *Q_L                        ;
    mxArray *dllDiscr                   ;
    mxArray *dllDiscrFilt               ;
    mxArray *pllDiscr                   ;
    mxArray *pllDiscrFilt               ;
    mxArray *remCodePhase               ;
    mxArray *remCarrPhase               ;
    mxArray *VSMValue                   ;
    mxArray *VSMIndex                   ;
    char *p_status                      ;
    double *p_PRN                       ;
    double *p_absoluteSample            ;
    double *p_codeFreq                  ;
    double *p_carrFreq                  ;
    double *p_I_P                       ;
    double *p_I_E                       ;
    double *p_I_L                       ;
    double *p_Q_P                       ;
    double *p_Q_E                       ;
    double *p_Q_L                       ; 
    double *p_dllDiscr                  ;
    double *p_dllDiscrFilt              ;
    double *p_pllDiscr                  ;
    double *p_pllDiscrFilt              ;
    double *p_remCodePhase              ;
    double *p_remCarrPhase              ;
    double *p_SMValue                   ;
    double *p_VSMIndex                  ;
    double *p_VSMValue                  ;

    int ThreadResult                    ;
    int ErrorLoopIdx                    ;
}Result;

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
    double acqSearchStep;
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
    //double CNo;
}Setting;

typedef struct _Channel{
    double PRN                          ;
    double acquiredFreq                 ;
    double codePhase                    ;
    char *status                        ;
}Channel;

typedef struct _FileInfo{
    const signed char * dataPtr        ;
    unsigned long Length               ;
}DataFileInfo;

typedef struct _ThreadParameter{
    Result **result                    ;
    int channelIndex                    ;
    Channel channel                       ;
    Setting setting                     ;
    DataFileInfo fileInfo          ;
}ThreadParameter;


int initSetting(Setting *setting, const mxArray *settingArray);
int initChannel(Channel channel[], const mxArray *channelArray);
Result* chn_calc_DP(Channel channel, DataFileInfo FileInfo, Setting setting, Result *result);
double *zeros(int length);
double *inf(int length);
Result initResult(Setting setting);
void calcLoopCoef(double *tau1, double *tau2, double LBW, double zeta, double k);
void calcLoopCoefCarr(double *pf1, double *pf2, double *pf3, Setting setting);
void generateCAcode(double PRN, double *output);
void gen_legendre_sequence();
unsigned int calcChannelData(void *params);
void dispose(Channel channel[], Setting setting);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    Channel chnl[CHANNEL_MAX];
    const signed char *dataPtr;
    Setting setting;
    Result result[CHANNEL_MAX];
    DataFileInfo dataInfo;

    time_t startTimer = time(NULL);
    struct tm* startTime = localtime(&startTimer);

    mexPrintf("\nCMEX Function Start: \t%2d:%2d:%2d\n", startTime->tm_hour, startTime->tm_min, startTime->tm_sec);

    const char *fieldNames[18];

    if(nrhs!=4)
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs",
                "Three input required.");
    else if(nlhs > 1)
        mexErrMsgIdAndTxt( "MATLAB:maxlhs",
                "Too many output arguments.");
    else if(!mxIsStruct(prhs[3]))
        mexErrMsgIdAndTxt( "MATLAB:inputNotStruct",
                "Input[2] must be a structure.");
    initSetting(&setting, prhs[3]);
    initChannel(chnl, prhs[0]);
    
    dataInfo.dataPtr = mxGetPr(prhs[1]);
    dataInfo.Length = mxGetScalar(prhs[2]);

    for(int i = 0; i < 18; i++){
        fieldNames[i] = (char*)mxMalloc(20);
    }

    memcpy(fieldNames[rst_status], "status", sizeof("status"));
    memcpy(fieldNames[rst_PRN], "PRN", sizeof("PRN"));
    memcpy(fieldNames[rst_absoluteSample], "absoluteSample", sizeof("absoluteSample"));
    memcpy(fieldNames[rst_codeFreq], "codeFreq", sizeof("codeFreq"));
    memcpy(fieldNames[rst_carrFreq], "carrFreq", sizeof("carrFreq"));
    memcpy(fieldNames[rst_I_P], "I_P", sizeof("I_P"));
    memcpy(fieldNames[rst_I_E], "I_E", sizeof("I_E"));
    memcpy(fieldNames[rst_I_L], "I_L", sizeof("I_L"));
    memcpy(fieldNames[rst_Q_P], "Q_P", sizeof("Q_P"));
    memcpy(fieldNames[rst_Q_E], "Q_E", sizeof("Q_E"));
    memcpy(fieldNames[rst_Q_L], "Q_L", sizeof("Q_L"));
    memcpy(fieldNames[rst_dllDiscr], "dllDiscr", sizeof("dllDiscr"));
    memcpy(fieldNames[rst_dllDiscrFilt], "dllDiscrFilt", sizeof("dllDiscrFilt"));
    memcpy(fieldNames[rst_pllDiscr], "pllDiscr", sizeof("pllDiscr"));
    memcpy(fieldNames[rst_pllDiscrFilt], "pllDiscrFilt", sizeof("pllDiscrFilt"));
    memcpy(fieldNames[rst_remCodePhase], "remCodePhase", sizeof("remCodePhase"));
    memcpy(fieldNames[rst_remCarrPhase], "remCarrPhase", sizeof("remCarrPhase"));
    memcpy(fieldNames[rst_CNo], "CNo", sizeof("CNo"));

    plhs[0] = mxCreateStructMatrix(1, CHANNEL_MAX, 18, fieldNames);

    mxFree(fieldNames[0]);
    mxFree(fieldNames[1]);
    mxFree(fieldNames[2]);
    mxFree(fieldNames[3]);
    mxFree(fieldNames[4]);
    mxFree(fieldNames[5]);
    mxFree(fieldNames[6]);
    mxFree(fieldNames[7]);
    mxFree(fieldNames[8]);
    mxFree(fieldNames[9]);
    mxFree(fieldNames[10]);
    mxFree(fieldNames[11]);
    mxFree(fieldNames[12]);
    mxFree(fieldNames[13]);
    mxFree(fieldNames[14]);
    mxFree(fieldNames[15]);
    mxFree(fieldNames[16]);
    mxFree(fieldNames[17]);

#ifdef DEBUG_MODE
    mexPrintf("\n\n---------------------Display Channel Data---------------------\n\n");
    for(int i = 0; i < 12; i++){
        mexPrintf("channel[%d].PRN: %f\n", i, chnl[i].PRN);
        mexPrintf("channel[%d].acquiredFreq: %f\n", i, chnl[i].acquiredFreq);
        mexPrintf("channel[%d].codePhase: %f\n", i, chnl[i].codePhase);
        mexPrintf("channel[%d].status: %s\n", i, chnl[i].status);
    }
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
    mexPrintf("acqSearchStep: %f\n", setting.acqSearchStep);
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
    HANDLE thread[CHANNEL_MAX];
    unsigned int threadID[CHANNEL_MAX];
    if(legendre[0] == 0) gen_legendre_sequence();
    
    for(int i = 0; i < CHANNEL_MAX; i++){
        int idx = i;
        result[idx] = initResult(setting);
        ThreadParameter param = {
            &result   ,
            idx       ,
            *(chnl + idx),
            setting ,
            dataInfo
        };
        thread[idx] = (HANDLE)_beginthreadex(NULL, 0, calcChannelData, (void*)&param, 0, threadID + idx);

        Sleep(10);
    }

    WaitForMultipleObjects(CHANNEL_MAX, thread, TRUE, INFINITE);

    fieldNames[0] = (char*)mxMalloc(20);
    fieldNames[1] = (char*)mxMalloc(20);
    for(int i = CHANNEL_MAX - 1; i >= 0; i--){
        if(result[i].ThreadResult == MemoryOut_Data){
            mexPrintf("Channel[%d]/PRN: %d Exception!\n", i, (int)(*result[i].p_PRN));
            mexPrintf("Not able to read the specified number of samples  for tracking, exiting!\n");
            result[i].status = mxCreateString("-");
            mxSetFieldByNumber(plhs[0], i, rst_status, result[i].status);
            mxSetFieldByNumber(plhs[0], i, rst_PRN, result[i].PRN);
            mexPrintf("LoopCnt: %d\n\n", result[i].ErrorLoopIdx);
            continue;
        }

        result[i].status = mxCreateString("T");

        mxSetFieldByNumber(plhs[0], i, rst_status, result[i].status);
        mxSetFieldByNumber(plhs[0], i, rst_PRN, result[i].PRN);
        mxSetFieldByNumber(plhs[0], i, rst_absoluteSample, result[i].absoluteSample);
        mxSetFieldByNumber(plhs[0], i, rst_codeFreq, result[i].codeFreq);
        mxSetFieldByNumber(plhs[0], i, rst_carrFreq, result[i].carrFreq);
        mxSetFieldByNumber(plhs[0], i, rst_I_P, result[i].I_P);
        mxSetFieldByNumber(plhs[0], i, rst_I_E, result[i].I_E);
        mxSetFieldByNumber(plhs[0], i, rst_I_L, result[i].I_L);
        mxSetFieldByNumber(plhs[0], i, rst_Q_E, result[i].Q_E);
        mxSetFieldByNumber(plhs[0], i, rst_Q_P, result[i].Q_P);
        mxSetFieldByNumber(plhs[0], i, rst_Q_L, result[i].Q_L);
        mxSetFieldByNumber(plhs[0], i, rst_dllDiscr, result[i].dllDiscr);
        mxSetFieldByNumber(plhs[0], i, rst_dllDiscrFilt, result[i].dllDiscrFilt);
        mxSetFieldByNumber(plhs[0], i, rst_pllDiscr, result[i].pllDiscr);
        mxSetFieldByNumber(plhs[0], i, rst_pllDiscrFilt, result[i].pllDiscrFilt);
        mxSetFieldByNumber(plhs[0], i, rst_remCodePhase, result[i].remCodePhase);
        mxSetFieldByNumber(plhs[0], i, rst_remCarrPhase, result[i].remCarrPhase);
        
    }
    mxFree(fieldNames[0]);
    mxFree(fieldNames[1]);
    
    dispose(chnl, setting);
    
    time_t endTimer = time(NULL);
    struct tm* endTime = localtime(&endTimer);

    mexPrintf("CMEX Function Finished: %2d:%2d:%2d\n", endTime->tm_hour, endTime->tm_min, endTime->tm_sec);
    mexPrintf("-----------------------------------------------\n");
}

unsigned int calcChannelData(void *params){
    ThreadParameter* param = (ThreadParameter*)params;
    Result *result = param->result;
    int i = param->channelIndex;
    Channel chnl = param->channel;
    DataFileInfo FileInfo = param->fileInfo;
    Setting setting = param->setting;
    chn_calc_DP(chnl, FileInfo, setting, &result[i]);

    return 0;
}

int initChannel(Channel channel[], const mxArray *channelArray){
    int ifield, nfields;
    mxArray *tmp;
    size_t buflen;
    
    for(int i = 0; i < CHANNEL_MAX; i++){
        tmp = mxGetFieldByNumber(channelArray,i, PRN);
        channel[i].PRN = mxGetScalar(tmp);
        
        tmp = mxGetFieldByNumber(channelArray,i, acquiredFreq);
        channel[i].acquiredFreq = mxGetScalar(tmp);
        
        tmp = mxGetFieldByNumber(channelArray,i, codePhase);
        channel[i].codePhase = mxGetScalar(tmp);
        
        tmp = mxGetFieldByNumber(channelArray,i, status);
        buflen = mxGetN(tmp) + 1;
        channel[i].status = mxMalloc(buflen);
        mxGetString(tmp, channel[i].status, (mwSize)buflen);
    }
    return 0;
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

    tmp = mxGetFieldByNumber(settingArray,0, acqSearchStep);
    setting->acqSearchBand = mxGetScalar(tmp);

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

Result initResult(Setting setting){
    Result result;
    int codePeriods = (int)round(setting.msToProcess / 1000 / setting.intTime);

    result.ThreadResult = None;
    result.ErrorLoopIdx = -1;

    result.status = mxCreateString("-");
    result.p_status = mxGetPr(result.status);

    result.absoluteSample = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_absoluteSample   = mxGetPr(result.absoluteSample);

    result.codeFreq = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_codeFreq   = mxGetPr(result.codeFreq);
    result.carrFreq = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_carrFreq   = mxGetPr(result.carrFreq);

    result.I_P = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.I_E = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.I_L = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.Q_P = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.Q_E = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.Q_L = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_I_P   = mxGetPr(result.I_P);
    result.p_I_E   = mxGetPr(result.I_E);
    result.p_I_L   = mxGetPr(result.I_L);
    result.p_Q_P   = mxGetPr(result.Q_P);
    result.p_Q_E   = mxGetPr(result.Q_E);
    result.p_Q_L   = mxGetPr(result.Q_L);

    result.dllDiscr = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.dllDiscrFilt = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.pllDiscr = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.pllDiscrFilt = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_dllDiscr   = mxGetPr(result.dllDiscr);
    result.p_dllDiscrFilt   = mxGetPr(result.dllDiscrFilt);
    result.p_pllDiscr   = mxGetPr(result.pllDiscr);
    result.p_pllDiscrFilt   = mxGetPr(result.pllDiscrFilt);

    result.remCarrPhase = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.remCodePhase = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_remCarrPhase   = mxGetPr(result.remCarrPhase);
    result.p_remCodePhase   = mxGetPr(result.remCodePhase);

    result.VSMIndex = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.VSMValue = mxCreateDoubleMatrix(1, codePeriods, 0);
    result.p_VSMIndex = mxGetPr(result.VSMIndex);
    result.p_VSMValue = mxGetPr(result.VSMValue);

    result.PRN = mxCreateDoubleScalar(0);
    result.p_PRN = mxGetPr(result.PRN);

    return result;
}

Result* chn_calc_DP(Channel channel, DataFileInfo fileInfo, Setting setting, Result *result){
    //Result chnlResult = initResult(setting);
    Result *chnlResult = result;
    Result resultIndexer;
    //mexPrintf("%f\n", setting.codeLength);

    resultIndexer.p_absoluteSample  = chnlResult->p_absoluteSample;
    resultIndexer.p_codeFreq        = chnlResult->p_codeFreq      ;
    resultIndexer.p_carrFreq        = chnlResult->p_carrFreq      ;
    resultIndexer.p_I_P             = chnlResult->p_I_P           ;
    resultIndexer.p_I_E             = chnlResult->p_I_E           ;
    resultIndexer.p_I_L             = chnlResult->p_I_L           ;
    resultIndexer.p_Q_P             = chnlResult->p_Q_P           ;
    resultIndexer.p_Q_E             = chnlResult->p_Q_E           ;
    resultIndexer.p_Q_L             = chnlResult->p_Q_L           ;
    resultIndexer.p_dllDiscr        = chnlResult->p_dllDiscr      ;
    resultIndexer.p_dllDiscrFilt    = chnlResult->p_dllDiscrFilt  ;
    resultIndexer.p_pllDiscr        = chnlResult->p_pllDiscr      ;
    resultIndexer.p_pllDiscrFilt    = chnlResult->p_pllDiscrFilt  ;
    resultIndexer.p_remCodePhase    = chnlResult->p_remCodePhase  ;
    resultIndexer.p_remCarrPhase    = chnlResult->p_remCarrPhase  ;
    resultIndexer.p_VSMIndex        = chnlResult->p_VSMIndex      ;
    resultIndexer.p_VSMValue        = chnlResult->p_VSMValue      ;
    
    if(channel.PRN != 0){
        int codePeriods = (int)round(setting.msToProcess / 1000 / setting.intTime);
        double earlyLateSpc = setting.dllCorrelatorSpacing;
        double PDIcode = setting.intTime;

        double tau1code = INT_MIN, tau2code = INT_MIN;
        calcLoopCoef(&tau1code, &tau2code, setting.dllNoiseBandwidth, setting.dllDampingRatio, 1);


        double PDIcarr = setting.intTime;

        double tau1carr = INT_MIN, tau2carr = INT_MIN;
        calcLoopCoef(&tau1carr, &tau2carr, setting.pllNoiseBandwidth, setting.pllDampingRatio, 0.25);

        int dataApdptCoeff = (int)setting.fileType == 1 ? 1 : 2;

        *chnlResult->p_PRN = channel.PRN;

        int idx = (int)(dataApdptCoeff * (setting.skipNumberOfBytes + channel.codePhase - 1));

        double *caCode = (double *)mxMalloc(sizeof(double) * (1023 + 2));
        generateCAcode(channel.PRN, caCode + 1);
        caCode[0] = caCode[1023];
        caCode[1024] = caCode[1];

        double codeFreq = setting.codeFreqBasis;
        double remCodePhase = 0;
        double carrFreq = channel.acquiredFreq;
        double carrFreqBasis = channel.acquiredFreq;
        double remCarrPhase = 0;
        double oldCodeNco = 0;
        double oldCodeError = 0;
        double oldCarrError = 0;
        double oldCarrNco = 0;
        double vsmCnt = 0;

        double codePhaseStep = 0;
        int blksize = 0;

        // if(channel.PRN == 17)
        //     mexPrintf("idx: %d / codePhase: %f\n", idx, channel.codePhase);

        for(int loopCnt = 0; loopCnt < codePeriods; loopCnt++){

            void* dataPtr = setting.dataType[3] != '1' ?
                            ((signed char*)fileInfo.dataPtr) + idx :
                            ((int *)fileInfo.dataPtr) + idx;
            // if(channel.PRN == 17 && loopCnt < 50)
            //     mexPrintf("idx: %d\n", idx);

            *chnlResult->p_absoluteSample++ = (int)(idx / dataApdptCoeff);

            double I_E = 0;
            double Q_E = 0;
            double I_P = 0;
            double Q_P = 0;
            double I_L = 0;
            double Q_L = 0;

            codePhaseStep = codeFreq / setting.samplingFreq;
            blksize = (int)ceil((setting.codeLength - remCodePhase) / codePhaseStep);
            // if(channel.PRN == 17 && loopCnt < 50)
            //     mexPrintf("blksize: (%f 0 %f) / %f = %f\n",setting.codeLength, remCodePhase, codePhaseStep, (setting.codeLength - remCodePhase) / codePhaseStep);
            
            *(resultIndexer.p_remCodePhase++) = remCodePhase;
            *(resultIndexer.p_remCarrPhase++) = remCarrPhase;

            if(idx + blksize >= fileInfo.Length){
                result->ThreadResult = MemoryOut_Data;
                result->ErrorLoopIdx = loopCnt;
                break;
            }            
            // if(channel.PRN == 17 && loopCnt < 50)
            //     mexPrintf("blksize: %d %s\n", dataApdptCoeff,setting.dataType);

            for(int i = 0; i < blksize; i++){
                if(dataApdptCoeff == 1){
                    int IFEN = 0;
                    
                    if(setting.dataType[3] != '1'){
                        IFEN = (int)*(((signed char*)dataPtr)++);
                    }
                    else{
                        IFEN = (int)*(((int*)dataPtr)++);
                    }
                    double earlyCode = caCode[(int)ceil((i * codePhaseStep + remCodePhase - earlyLateSpc))];
                    double lateCode = caCode[(int)ceil((i * codePhaseStep + remCodePhase + earlyLateSpc))];
                    double promptCode = caCode[(int)ceil((i * codePhaseStep + remCodePhase))];
                    // if(channel.PRN == 17 && loopCnt == 0 && blksize < 100){
                    //     mexPrintf("coef[%d]: %d\n", i, IFEN);
                    // }
                    
                    double y = (double)((carrFreq * 2.0 * M_PI) * (double)i / setting.samplingFreq) + remCarrPhase;
                    double iBaseBandSignal = sin(y) * IFEN;
                    double qBaseBandSignal = cos(y) * IFEN;

                    I_E += earlyCode * iBaseBandSignal;
                    Q_E += earlyCode * qBaseBandSignal;
                    I_P += promptCode * iBaseBandSignal;
                    Q_P += promptCode * qBaseBandSignal;
                    I_L += lateCode * iBaseBandSignal;
                    Q_L += lateCode * qBaseBandSignal;

                }
                else{
                    int real = 0;
                    int imag = 0;
                    if(setting.dataType[3] != '1'){
                        real = (int)*(((signed char*)dataPtr)++);
                        imag = (int)*(((signed char*)dataPtr)++);
                    }
                    else{
                        real = (int)*(((int*)dataPtr)++);
                        imag = (int)*(((int*)dataPtr)++);
                    }
                    double earlyCode = caCode[(int)ceil((i * codePhaseStep + remCodePhase - earlyLateSpc))];
                    double lateCode = caCode[(int)ceil((i * codePhaseStep + remCodePhase + earlyLateSpc))];
                    double promptCode = caCode[(int)ceil((i * codePhaseStep + remCodePhase))];
                    

                    double y = (double)((carrFreq * 2.0 * M_PI) * (double)i / setting.samplingFreq) + remCarrPhase;
                    double iBaseBandSignal = real * cos(y) + sin(y) * imag;
                    double qBaseBandSignal = cos(y) * imag - sin(y) * real;

                    // if(channel.PRN == 17 && loopCnt == 1 && i < 100){
                    //     mexPrintf("real[%d]: %d\n", i, real);
                    //     mexPrintf("imag[%d]: %d\n", i, imag);                        
                    //     //mexPrintf("iBaseBandSignal[%d]: %f\n", i, iBaseBandSignal);
                    //     //mexPrintf("qBaseBandSignal[%d]: %f\n", i, qBaseBandSignal);
                    // }
                    I_E += earlyCode * iBaseBandSignal;
                    Q_E += earlyCode * qBaseBandSignal;
                    I_P += promptCode * iBaseBandSignal;
                    Q_P += promptCode * qBaseBandSignal;
                    I_L += lateCode * iBaseBandSignal;
                    Q_L += lateCode * qBaseBandSignal;
                }
            }

            
            remCodePhase = ((blksize-1)*codePhaseStep+remCodePhase) + codePhaseStep - setting.codeLength; 
            remCarrPhase = remainder((double)((carrFreq * 2.0 * M_PI) * (double)blksize / setting.samplingFreq) + remCarrPhase, 2 * M_PI);
            if(remCarrPhase < 0)
                remCarrPhase += 2 * M_PI;
            
            double carrError = atan(Q_P / I_P) / (2.0 * M_PI);
            double carrNco = oldCarrNco + (tau2carr / tau1carr) * (carrError - oldCarrError) + carrError * (PDIcarr / tau1carr);
            oldCarrError = carrError;
            oldCarrNco = carrNco;

            *(resultIndexer.p_carrFreq++) = carrFreq;
            carrFreq = carrFreqBasis + carrNco;
            
            double codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            double codeNco = oldCodeNco + (tau2code/tau1code) * (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;

            *(resultIndexer.p_codeFreq++) = codeFreq;
            codeFreq = setting.codeFreqBasis - codeNco;
            
            *(resultIndexer.p_dllDiscr++)       = codeError;
            *(resultIndexer.p_dllDiscrFilt++)   = codeNco;
            *(resultIndexer.p_pllDiscr++)       = carrError;
            *(resultIndexer.p_pllDiscrFilt++)   = carrNco;
            
            *(resultIndexer.p_I_E++) = I_E;
            *(resultIndexer.p_I_P++) = I_P;
            *(resultIndexer.p_I_L++) = I_L;
            *(resultIndexer.p_Q_E++) = Q_E;
            *(resultIndexer.p_Q_P++) = Q_P;
            *(resultIndexer.p_Q_L++) = Q_L;
            
            idx += (int)dataApdptCoeff*blksize;
        }
        mxFree(caCode);
    }     
    return chnlResult;
}

void calcLoopCoef(double *tau1, double *tau2, double LBW, double zeta, double k){
    double Wn = LBW * 8 * zeta / (4 * zeta * zeta + 1);
    *tau1 = k / (Wn * Wn);
    *tau2 = 2.0 * zeta / Wn;    
}

void calcLoopCoefCarr(double *pf1, double *pf2, double *pf3, Setting setting){
    double a3 = 1.1;
    double b3 = 2.4;
    double Wn = setting.pllNoiseBandwidth / 0.7845;
    *pf1 = b3 * Wn;
    *pf2 = a3 * Wn * Wn * setting.intTime;
    *pf3 = Wn * Wn * Wn * setting.intTime * setting.intTime;
}

void generateCAcode(double PRN, double* output){
    double g2s[51] = 
    {  5,   6,   7,   8,  17,  18, 139, 140, 141, 251,   //...
       252, 254, 255, 256, 257, 258, 469, 470, 471, 472, //...
       473, 474, 509, 512, 513, 514, 515, 516, 859, 860, //...
       861, 862,
       //... end of shifts for GPS satellites 
       //... Shifts for the ground GPS transmitter are not included
       //... Shifts for EGNOS and WAAS satellites (true_PRN = PRN + 87)
                 145, 175,  52,  21, 237, 235, 886, 657, //...
       634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, //...
       386};

    double g2shift = g2s[(int)PRN - 1];
    int g2[1023];
    int reg_g1[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int reg_g2[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    int saveBit_g1 = 0;
    int saveBit_g2 = 0;
    for(int i = 0; i < 1023; i++)
    {
        output[i] = (double)reg_g1[9];
        saveBit_g1 = (int)(reg_g1[2] * reg_g1[9]);
        g2[i] = reg_g2[9];
        saveBit_g2 = (int)(reg_g2[1] * reg_g2[2] * reg_g2[5] * reg_g2[7] * reg_g2[8] * reg_g2[9]);
        
        for(int j = 8; j >= 0; j--){
            reg_g1[j + 1] = reg_g1[j];
            reg_g2[j + 1] = reg_g2[j];
        }
        reg_g1[0] = saveBit_g1;
        reg_g2[0] = saveBit_g2;
    }
    for(int i = 0; i < 1023; i++){
        int idx = (1023 - (int)g2shift + i) % 1023;
        output[i] *= -g2[idx];
    }
}

void gen_legendre_sequence()
{
    int i;
    for (i=0;i<10223;i++)
        legendre[i]=1;

    for(i=0;i<10224;i++)
        legendre[(i*i)%10223]=-1;
    legendre[0]=1;
}

void dispose(Channel channel[], Setting setting){
    mxFree(setting.fileName);
    mxFree(setting.dataType);
}