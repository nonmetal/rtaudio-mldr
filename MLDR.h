#include <complex>

#include <algorithm>

class MLDR {
  int cnt = 0;
  int ref_cnt = 1;

private:
	int    nShift;     
	int    fftSize;     
	int    nFreq;      
	int    nOverlap;    
	int	   nChannel;   

  const double min_pos = 1e-32;
  const double eps = 1e-16;
  int power_method_iter; // 10->1
		
	double gamma;    // forgetting factor for inverse matrix
  double gammaPre, gammaPost;
  int initFrame;
  double alpha;    // smoothing factor for phiY
	double delta_MLDR;    // diagonal loading for inverse matrix
	double epsi_MLDR;     // floor value for lambdaY

  double delta_SVE;
  double theta_SVE;
  double epsi_SVE;

	double *phiY;
	double **W, **XPsi, **WNumer;
	double ***numer, ***Rx, ***PsiRx, ***Psi, ***PsiTmp;

  double* YPred, *denom, *WDenom; // nFreq(complex)
  double* lambdaY, * lambdaYInv;//nFreq(double)

  double* normRxTmp;
  double*** Rn;
  double*** RxTmp, *** RnTmp, ***XFrameTmp,*** RxTmp2, *** RnTmp2,***eigTmp ; //Nfreq*Nch*Nch(complex)
  double**maxEigVec, ** st,**eigVecTmp; //Nch*Nfreq(complex)
  double*LambdaX,*LambdaN, *PhiYHat,*tmpRe,*tmpIm; //Nfreq(real)

  //MLDR v2
  double*** Phi, *** PhiRxPhi, *** PhiHat, *** PhiHatTmp, *** PhiRxPhiTmp2;//nFreq*Nch*NCh(complex)
  double* XPhiX, * YHat;//nFreq (complex)
  double* PsiY, * PsiYHat; //nFreq(double)
  double** XPhiXTmp, ** WDenomTmp, ** PhiRxPhiTmp; //nFreq * ch (complex)

  //SVE v7
  double*** RsTmp;
  int* RsIdx;

  double *SHat;
  double **steerVector;

  int frame;
public:
	inline MLDR(int fftSize_, int nChannel_, 
    // See https://github.com/kooBH/MLDR/issues/1 for parameters history
    /** 20201203 **/
//    double gammaPre_ = 0.3,
//    double gammaPost_ = 0.995,
//    int initFrame_ = 0,
//    double alpha_MLDR_ = 0.2,
//    double delta_MLDR_ = 1,
//    double epsi_MLDR_ = 1e-3,
//    double delta_SVE_ = 1,
//    double epsi_SVE_ = 1e-3,
//    double theta_SVE_ = 0.995,
//    double power_method_iter_ = 1

//    /** 20210304 **/
    double gammaPre_=0.995,
    double gammaPost_=0.995,
    int initFrame_=0,
    double alpha_MLDR_=0.2,
    double delta_MLDR_=1e-2,
    double epsi_MLDR_=1e-6, 
    double delta_SVE_=0.0, 
    double epsi_SVE_=1e-6, 
    double theta_SVE_=0.995,
    double power_method_iter_=1
    ); 
	inline ~MLDR();  
  inline void Process(double** X);
  inline void Process(double** X,int target_channels);

  // initial
  inline void SVE_routine(double** X, double* SHat, double** steerVector, int freq);
  // update 20210303
  inline void SVE_routine_v7(double** X, double* SHat, double** steerVector, int freq);
  // initial
  inline void MLDR_routine(double** X, double* SHat, double** steerVector,int freq);
  // reimpelemnted
  inline void MLDR_routine_v2(double** X, double* SHat, double** steerVector, int freq);

  inline void Clear();
};

inline MLDR::MLDR(int fftSize_, int nChannel_,
  double gammaPre_, double gammaPost_, int initFrame_,
  double alpha_MLDR_, double delta_MLDR_, double epsi_MLDR_,
  double delta_SVE_,double epsi_SVE_,double theta_SVE_,double power_method_iter_)
{
  int channel, channel2, sample, freq;

  fftSize = fftSize_;
  nFreq = fftSize / 2 + 1;
  nChannel = nChannel_;

  gamma = gammaPost_;
  gammaPre = gammaPre_;
  gammaPost = gammaPost_;
  initFrame = initFrame_;
  alpha = alpha_MLDR_;
  delta_MLDR = delta_MLDR_;
  epsi_MLDR = epsi_MLDR_;

  delta_SVE = delta_SVE_;
  epsi_SVE = epsi_SVE_;
  theta_SVE = theta_SVE_;

  power_method_iter = power_method_iter_;

  frame = 0;

  YPred = new double[nFreq * 2];
  memset(YPred, 0, sizeof(double) * nFreq * 2);
  lambdaY = new double[nFreq];
  memset(lambdaY, 0, sizeof(double) * nFreq);
  lambdaYInv = new double[nFreq];
  memset(lambdaYInv, 0, sizeof(double) * nFreq);

  phiY = new double[nFreq];
  memset(phiY, 0, (nFreq) * sizeof(double));

  PhiYHat = new double[nFreq];
  memset(PhiYHat, 0, sizeof(double) * nFreq);
 
  W = new double* [nChannel];
  for (channel = 0; channel < nChannel; channel++) {
    W[channel] = new double[nFreq * 2];
    memset(W[channel], 0, sizeof(double) * (nFreq * 2));
  }

  WNumer = new double* [nFreq*2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    WNumer[freq] = new double[nChannel];
    memset(WNumer[freq], 0, sizeof(double) * nChannel);
  }
  XPsi = new double* [nFreq *2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    XPsi[freq] = new double[nChannel];
    memset(XPsi[freq], 0, sizeof(double) * nChannel);
  }

  PsiRx = new double** [nFreq *2];
  for (freq = 0; freq < nFreq*2; freq++) {
    PsiRx[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      PsiRx[freq][channel] = new double[nChannel];
      memset(PsiRx[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  numer = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    numer[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      numer[freq][channel] = new double[nChannel];
      memset(numer[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  PsiTmp = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    PsiTmp[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      PsiTmp[freq][channel] = new double[nChannel];
      memset(PsiTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }

  normRxTmp = new double[nFreq];

  // NOTE 
  // Rx = repmat(deltaSt*eye(nSensor),[1 1 nFreq]);  
  // Rn = repmat(deltaSt*eye(nSensor),[1 1 nFreq]); 
  //================================================
  Rx = new double** [nFreq*2];
  for (freq = 0; freq < nFreq*2; freq++) {
    Rx[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      Rx[freq][channel] = new double[nChannel];
      memset(Rx[freq][channel], 0, (nChannel) * sizeof(double));
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        if (channel == channel2) {
          if(freq%2==0)
            Rx[freq][channel][channel2] = delta_SVE;
        }
      }
    }
  }
  Rn = new double** [nFreq*2];
  for (freq = 0; freq < nFreq*2; freq++) {
    Rn[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      Rn[freq][channel] = new double[nChannel];
      memset(Rn[freq][channel], 0, (nChannel) * sizeof(double));
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        if (channel == channel2) {
          if (freq % 2 == 0)
            Rn[freq][channel][channel2] = delta_SVE;
        }
      }
    }
  }


	denom = new double[nFreq];
  memset(denom, 0, sizeof(double) * nFreq);
  WDenom = new double[nFreq];
  memset(WDenom, 0, sizeof(double) * nFreq);

	Psi = new double**[fftSize + 2];	

	for (freq = 0; freq < fftSize + 2; freq++)
	{
		Psi[freq] = new double*[nChannel];

		for (channel = 0; channel < nChannel; channel++){
			Psi[freq][channel] = new double[nChannel];
			memset(Psi[freq][channel], 0, (nChannel) * sizeof(double));
			
			if (freq % 2 == 0)
			{
				Psi[freq][channel][channel] = delta_MLDR;
			}
		}
	}

  RxTmp = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    RxTmp[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      RxTmp[freq][channel] = new double[nChannel];
      memset(RxTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  RnTmp = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    RnTmp[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      RnTmp[freq][channel] = new double[nChannel];
      memset(RnTmp[freq][channel], 0, sizeof(double)* nChannel);
    }
  }
  // maxEigVec = ones(nSensor, nFreq)./nSensor;
  maxEigVec = new double*[nFreq * 2];
  for (freq = 0; freq < nFreq; freq++) {
    maxEigVec[freq+ freq] = new double[nChannel];
    maxEigVec[freq+ freq+1] = new double[nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      maxEigVec[freq+ freq][channel] = 1.0/nChannel;
      maxEigVec[freq+ freq+1][channel] = 0.0;
    }
  }
  eigVecTmp = new double* [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    eigVecTmp[freq] = new double[nChannel];
  }
  
  st = new double* [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    st[freq] = new double[nChannel];
    memset(st[freq], 0, sizeof(double) * nChannel);
  }

  XFrameTmp = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    XFrameTmp[freq] = new double* [nChannel];
  printf("");
    for (channel = 0; channel < nChannel; channel++) {
      XFrameTmp[freq][channel] = new double[nChannel];
      memset(XFrameTmp[freq][channel], 0, sizeof(double));
    }
  }
  RxTmp2 = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    RxTmp2[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      RxTmp2[freq][channel] = new double[nChannel];
      memset(RxTmp2[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  RnTmp2 = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    RnTmp2[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      RnTmp2[freq][channel] = new double[nChannel];
      memset(RnTmp2[freq][channel], 0, sizeof(double) * nChannel);
    }
  printf("");
  }
  eigTmp = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    eigTmp[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      eigTmp[freq][channel] = new double[nChannel];
      memset(eigTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }

  LambdaX = new double[nFreq];
  memset(LambdaX, 0, sizeof(double)* nFreq);

  LambdaN = new double[nFreq];
  memset(LambdaN, 0, sizeof(double)* nFreq);

  tmpRe = new double[nFreq];
  tmpIm = new double[nFreq];

  SHat = new double [fftSize+2];
  memset(SHat,0,sizeof(double)*(fftSize+2));

  steerVector = new double*[nChannel];
  for(channel=0;channel <nChannel;channel++){
    steerVector[channel] = new double[nFreq *2];
    for(freq=0;freq<nFreq;freq++){
      steerVector[channel][freq+freq]=1.0/nChannel;
      steerVector[channel][freq+freq+1]=0.0;
    }
  }
  //MLDR v2
  //double *** Phi,*** PhiRxPhi,*** PhiHat, *** PhiHatTmp;//nFreq*Nch*NCh(complex)
  PhiRxPhiTmp2 = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    PhiRxPhiTmp2[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      PhiRxPhiTmp2[freq][channel] = new double[nChannel];
      memset(PhiRxPhiTmp2[freq][channel], 0, sizeof(double) * nChannel);
    }
  }

  Phi = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    Phi[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      Phi[freq][channel] = new double[nChannel];
      memset(Phi[freq][channel], 0, sizeof(double) * nChannel);
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        if (channel == channel2) {
          if (freq % 2 == 0)
            Phi[freq][channel][channel2] = delta_MLDR;
        }
      }
    }
  }
  PhiRxPhi = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    PhiRxPhi[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      PhiRxPhi[freq][channel] = new double[nChannel];
      memset(PhiRxPhi[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  PhiHat = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    PhiHat[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      PhiHat[freq][channel] = new double[nChannel];
      memset(PhiHat[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  PhiHatTmp = new double** [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    PhiHatTmp[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      PhiHatTmp[freq][channel] = new double[nChannel];
      memset(PhiHatTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  // nFreq * ch
  PhiRxPhiTmp = new double* [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    PhiRxPhiTmp[freq] = new double[nChannel];
    memset(PhiRxPhiTmp[freq], 0, sizeof(double) * nChannel);
  }
  XPhiXTmp = new double* [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    XPhiXTmp[freq] = new double[nChannel];
    memset(XPhiXTmp[freq], 0, sizeof(double) * nChannel);
  }
  WDenomTmp = new double* [nFreq * 2];
  for (freq = 0; freq < nFreq * 2; freq++) {
    WDenomTmp[freq] = new double[nChannel];
    memset(WDenomTmp[freq], 0, sizeof(double) * nChannel);
  }
  //double* XPhiX, * YHat;//nFreq (complex)
  //double* PsiY, * PsiYHat;; //nFreq

  XPhiX = new double[nFreq * 2];
  memset(XPhiX, 0, sizeof(double) * (nFreq * 2));
  YHat = new double[nFreq * 2];
  memset(YHat, 0, sizeof(double) * (nFreq * 2));

  PsiY = new double[nFreq];
  memset(PsiY, 0, sizeof(double) * (nFreq));
  PsiYHat = new double[nFreq];
  memset(PsiYHat, 0, sizeof(double) * (nFreq));

  // SVE v7
  RsTmp = new double** [nFreq];
  for (freq = 0; freq < nFreq; freq++) {
    RsTmp[freq] = new double* [nChannel];
    for (channel = 0; channel < nChannel; channel++) {
      RsTmp[freq][channel] = new double[nChannel];
      memset(RsTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  RsIdx = new int[nFreq];
  memset(RsIdx, 0, sizeof(int) * nFreq);
}

inline MLDR::~MLDR()
{
	int channel, freq;

  delete[] YPred;
  delete[] lambdaY;
  delete[] lambdaYInv;

  delete[] phiY;

  delete[] denom;
  delete[] WDenom;

  delete[] normRxTmp;
  
	for (freq = 0; freq < fftSize + 2; freq++)
	{
		for (channel = 0; channel < nChannel; channel++)
		{
			delete[] Psi[freq][channel];
		}
		delete[] Psi[freq];
	}
	delete[] Psi;
  delete[] PhiYHat;

	for (channel = 0; channel < nChannel; channel++){
		delete[] W[channel];
	}
	delete[] W;

  for (freq = 0; freq < nFreq*2; freq++) {
    delete[] WNumer[freq];
    delete[] XPsi[freq];
  }
  delete[] WNumer;
  delete[] XPsi;

  for (freq = 0; freq < nFreq*2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      delete[] PsiRx[freq][channel];
      delete[] numer[freq][channel];
      delete[] PsiTmp[freq][channel];
    }
    delete[] PsiRx[freq];
    delete[] numer[freq];
    delete[] PsiTmp[freq];
  }
  delete[] PsiRx;
  delete[] numer;
  delete[] PsiTmp;

  for (freq = 0; freq < nFreq*2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      delete[] RxTmp[freq][channel];
      delete[] RnTmp[freq][channel];
      delete[] XFrameTmp[freq][channel];
      delete[] RxTmp2[freq][channel];
      delete[] RnTmp2[freq][channel];
      delete[] eigTmp[freq][channel];
    }
    delete[] RxTmp[freq];
    delete[] RnTmp[freq];
    delete[] XFrameTmp[freq];
    delete[] RxTmp2[freq];
    delete[] RnTmp2[freq];
    delete[] eigTmp[freq];
  }

  delete[] RxTmp;
  delete[] RnTmp;
  delete[] XFrameTmp;
  delete[] RxTmp2;
  delete[] RnTmp2;
  delete[] eigTmp;

  for (freq = 0; freq < nFreq*2; freq++) {
    delete[] maxEigVec[freq];
    delete[] st[freq];
    delete[] eigVecTmp[freq];
  }
  delete[] maxEigVec;
  delete[] st;
  delete[] eigVecTmp;

  for (freq = 0; freq < nFreq*2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      delete[] Rx[freq][channel];
      delete[] Rn[freq][channel];
    }
    delete[] Rx[freq];
    delete[] Rn[freq];
  }
  delete[] Rx;
  delete[] Rn;

  delete[] LambdaN;
  delete[] LambdaX;

  delete[] tmpRe;
  delete[] tmpIm;

  delete[] SHat;


  for(channel=0;channel<nChannel;channel++)
    delete[] steerVector[channel];
  delete[] steerVector;

  //ver 2
  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      delete[] Phi[freq][channel];
      delete[] PhiRxPhi[freq][channel];
      delete[] PhiHat[freq][channel];
      delete[] PhiHatTmp[freq][channel];
      delete[] PhiRxPhiTmp2[freq][channel];
    }
    delete[] Phi[freq];
    delete[] PhiRxPhi[freq];
    delete[] PhiHat[freq];
    delete[] PhiHatTmp[freq];
    delete[] PhiRxPhiTmp2[freq];
  }
  delete[] Phi;
  delete[] PhiRxPhi;
  delete[] PhiHat;
  delete[] PhiHatTmp;
  delete[] PhiRxPhiTmp2;

  delete[] XPhiX;
  delete[] YHat;
  delete[] PsiY;
  delete[] PsiYHat;

  for (freq = 0; freq < nFreq * 2; freq++) {
    delete[] XPhiXTmp[freq];
    delete[] WDenomTmp[freq];
    delete[] PhiRxPhiTmp[freq];
  }
  delete[] XPhiXTmp;
  delete[] WDenomTmp;
  delete[] PhiRxPhiTmp;

  //ver 7
  for (freq = 0; freq < nFreq; freq++) {
    for (channel = 0; channel < nChannel;  channel++) {
      delete[] RsTmp[freq][channel];
    }
    delete[] RsTmp[freq];
  }
  delete[] RsTmp;
  delete[] RsIdx;
}


inline void MLDR::Process(double** X){
  cnt++;
  if (frame < initFrame)
    gamma = gammaPre;
  else
    gamma = gammaPost;
  int freq;

//#pragma omp parallel for schedule(static,32)
  for (freq = 0; freq < nFreq; freq++){
    int channel, channel2, channel3, re, im;
    re = freq + freq;
    im = freq + freq + 1;

    //[M] YHat = W(:,freq)'*XFrame(:,freq);
    YPred[re] = 0.0;
    YPred[im] = 0.0;
    for (channel = 0; channel < nChannel; channel++){
      YPred[re] += W[channel][re] * X[channel][re] + W[channel][im] * X[channel][im];
      YPred[im] += W[channel][re] * X[channel][im] - W[channel][im] * X[channel][re];
    }
    
    /********** steering vector estimation *********/
   // SVE_routine(X, SHat, steerVector, freq);
    SVE_routine_v7(X, SHat, steerVector, freq);

    /********** MLDR beamforming *********/
    //MLDR_routine(X, SHat, steerVector, freq);
    
    MLDR_routine_v2(X, SHat, steerVector, freq);
    X[0][re] = SHat[re];
    X[0][im] = SHat[im];

  }
}

inline void MLDR::Process(double** X, int target_channels){
	int tmp = nChannel;
	nChannel = target_channels;
	Process(X);
	nChannel = tmp;
}


inline void MLDR::SVE_routine_v7(double** X, double* SHat, double** steerVector, int freq) {
  int channel, channel2, channel3, re, im;
  re = freq + freq;
  im = freq + freq + 1;
  //[M]PhiYHat = abs(YHat). ^ 2;
  PhiYHat[freq] = YPred[re] * YPred[re] + YPred[im] * YPred[im];
  //[M]lambdaY = max(PhiYHat, epsiSt);
  lambdaY[freq] = std::max(PhiYHat[freq], epsi_SVE);
  //[M]lambdaYInv = 1 ./ lambdaY;
  lambdaYInv[freq] = 1.0 / (lambdaY[freq]);


  // XFrameTmp = (XFrame(:,freq)*XFrame(:,freq)')
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      XFrameTmp[re][channel][channel2] = X[channel][re] * X[channel2][re] + X[channel][im] * X[channel2][im];
      XFrameTmp[im][channel][channel2] = -X[channel][re] * X[channel2][im] + X[channel][im] * X[channel2][re];
    }
  }

  //RxTmp = Rx(:,:,freq) * LambdaX(freq) / (LambdaX(freq) + 1) +   (XFrame(:,freq)*XFrame(:,freq)') / (LambdaX(freq) + 1);
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      // real
      RxTmp[re][channel][channel2] =
        Rx[re][channel][channel2] * LambdaX[freq] / (LambdaX[freq] + 1.0)
        + XFrameTmp[re][channel][channel2] / (LambdaX[freq] + 1.0);
      // imag
      RxTmp[im][channel][channel2] =
        Rx[im][channel][channel2] * LambdaX[freq] / (LambdaX[freq] + 1.0)
        + XFrameTmp[im][channel][channel2] / (LambdaX[freq] + 1.0);
    }
  }


  // RnTmp = Rn(:,:,freq) * LambdaN(freq)  / (LambdaN(freq) + lambdaYInv) +  lambdaYInv.*(XFrame(:,freq) * XFrame(:,freq)') / (LambdaN(freq) + lambdaYInv);
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      //real
      RnTmp[re][channel][channel2] = Rn[re][channel][channel2] * LambdaN[freq] / (LambdaN[freq] + lambdaYInv[freq])
        + lambdaYInv[freq] * XFrameTmp[re][channel][channel2] / (LambdaN[freq] + lambdaYInv[freq]);
      //imag
      RnTmp[im][channel][channel2] = Rn[im][channel][channel2] * LambdaN[freq] / (LambdaN[freq] + lambdaYInv[freq])
        + lambdaYInv[freq] * XFrameTmp[im][channel][channel2] / (LambdaN[freq] + lambdaYInv[freq]);
    }
  }
  //  RxTmp = 0.5 * (RxTmp + RxTmp');
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      RxTmp2[re][channel][channel2] = RxTmp[re][channel][channel2] + RxTmp[re][channel2][channel];
      RxTmp2[re][channel][channel2] *= 0.5;
      if (channel != channel2) {
        RxTmp2[im][channel][channel2] = RxTmp[im][channel][channel2] - RxTmp[im][channel2][channel];
        RxTmp2[im][channel][channel2] *= 0.5;
      }
      else {
        RxTmp2[im][channel][channel2] = 0.0;
      }
    }
  }

  //  RnTmp = 0.5 * (RnTmp + RnTmp');
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      if (channel != channel2) {
        RnTmp2[re][channel][channel2] = RnTmp[re][channel][channel2] + RnTmp[re][channel2][channel];
        RnTmp2[re][channel][channel2] *= 0.5;
        RnTmp2[im][channel][channel2] = RnTmp[im][channel][channel2] - RnTmp[im][channel2][channel];
        RnTmp2[im][channel][channel2] *= 0.5;
      }
      else {
        RnTmp2[re][channel][channel2] = RnTmp[re][channel][channel2];
        RnTmp2[im][channel][channel2] = 0.0;
      }
    }
  }
  //norm(RxTmp)
  normRxTmp[freq] = 0.0;
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      normRxTmp[freq]
        += RxTmp2[re][channel][channel2] * RxTmp2[re][channel][channel2]
        + RxTmp2[im][channel][channel2] * RxTmp2[im][channel][channel2];
    }
  }


  //if norm(RxTmp) ~= 0   % to prevent division by zero (20201213)
  if (normRxTmp[freq] != 0) {
    // RsTmp = RxTmp - theta*RnTmp;
    for (channel = 0; channel < nChannel; channel++)
      memset(RsTmp[freq][channel], 0, sizeof(double) * nChannel);
    for (channel = 0; channel < nChannel; channel++) {
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        RsTmp[freq][channel][channel2]
          = RxTmp2[re][channel][channel2] - theta_SVE * RnTmp2[re][channel][channel2];
      }
    }
    // findIndex = find(diag(RsTmp) < 0);
    RsIdx[freq] = 0;
    for (channel = 0; channel < nChannel; channel++) {
      if (RsTmp[freq][channel][channel] < 0) RsIdx[freq]++;
    }
    // if sum(findIndex) == 0
    if (RsIdx[freq] == 0) {
      // for iterPM = 1 : maxIterPM
      // maxEigVec(:,freq) = (RxTmp - theta*RnTmp)*maxEigVec(:,freq);
      // maxEigVec(:, freq) = maxEigVec(:, freq). / norm(maxEigVec(:, freq));
      // end
      for (channel = 0; channel < nChannel; channel++) {
        for (channel2 = 0; channel2 < nChannel; channel2++) {
          eigTmp[re][channel][channel2] = RxTmp2[re][channel][channel2] - theta_SVE * RnTmp2[re][channel][channel2];
          eigTmp[im][channel][channel2] = RxTmp2[im][channel][channel2] - theta_SVE * RnTmp2[im][channel][channel2];
        }
      }
      for (int i = 0; i < power_method_iter; i++) {
        for (channel = 0; channel < nChannel; channel++) {
          tmpRe[freq] = 0.0;
          tmpIm[freq] = 0.0;
          for (channel2 = 0; channel2 < nChannel; channel2++) {
            tmpRe[freq] += eigTmp[re][channel][channel2] * maxEigVec[re][channel2]
              - eigTmp[im][channel][channel2] * maxEigVec[im][channel2];
            tmpIm[freq] += eigTmp[im][channel][channel2] * maxEigVec[re][channel2]
              + eigTmp[re][channel][channel2] * maxEigVec[im][channel2];
          }
          eigVecTmp[re][channel] = tmpRe[freq];
          eigVecTmp[im][channel] = tmpIm[freq];
        }

        for (channel = 0; channel < nChannel; channel++) {
          maxEigVec[re][channel] = eigVecTmp[re][channel];
          maxEigVec[im][channel] = eigVecTmp[im][channel];
        }
        //norm
        tmpRe[freq] = 0.0;
        for (channel = 0; channel < nChannel; channel++)
          tmpRe[freq] += maxEigVec[re][channel] * maxEigVec[re][channel]
          + maxEigVec[im][channel] * maxEigVec[im][channel];
        tmpRe[freq] = sqrt(tmpRe[freq]);

        for (channel = 0; channel < nChannel; channel++) {
          maxEigVec[re][channel] /= (tmpRe[freq] + min_pos);
          maxEigVec[im][channel] /= (tmpRe[freq] + min_pos);
        }
      }
      // st(:,freq) = maxEigVec(:,freq)./maxEigVec(1,freq);
      tmpRe[freq] = maxEigVec[re][0] * maxEigVec[re][0] + maxEigVec[im][0] * maxEigVec[im][0];
      for (channel = 0; channel < nChannel; channel++) {
        steerVector[channel][re] = (maxEigVec[re][channel] * maxEigVec[re][0] + maxEigVec[im][channel] * maxEigVec[im][0]) / (tmpRe[freq] + min_pos);
        steerVector[channel][im] = (maxEigVec[im][channel] * maxEigVec[re][0] - maxEigVec[re][channel] * maxEigVec[im][0]) / (tmpRe[freq] + min_pos);
      }
      // st(:, freq) = st(:, freq).*sqrt(nSensor / norm(st(:, freq)));
      tmpRe[freq] = 0.0;
      for (channel = 0; channel < nChannel; channel++)
        tmpRe[freq] += steerVector[channel][re] * steerVector[channel][re] + steerVector[channel][im] * steerVector[channel][im];
      tmpRe[freq] = sqrt(tmpRe[freq]);
      tmpRe[freq] = sqrt(nChannel / (tmpRe[freq] + min_pos));
      for (channel = 0; channel < nChannel; channel++) {
        steerVector[channel][re] = steerVector[channel][re] * tmpRe[freq];
        steerVector[channel][im] = steerVector[channel][im] * tmpRe[freq];
      }
     //Rn(:, : , freq) = RnTmp;
      for (channel = 0; channel < nChannel; channel++) {
        for (channel2 = 0; channel2 < nChannel; channel2++) {
          Rn[re][channel][channel2] = RnTmp2[re][channel][channel2];
          Rn[im][channel][channel2] = RnTmp2[im][channel][channel2];
        }
      }
      // LambdaN(freq) = gamma * (LambdaN(freq) + lambdaYInv);
      LambdaN[freq] = gamma * (LambdaN[freq] + lambdaYInv[freq]);
    }//if (RsIdx == 0)
     // Rx(:,:,freq) = RxTmp;
    for (channel = 0; channel < nChannel; channel++) {
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        Rx[re][channel][channel2] = RxTmp2[re][channel][channel2];
        Rx[im][channel][channel2] = RxTmp2[im][channel][channel2];
      }
    }
    // LambdaX(freq) = gamma * (LambdaX(freq) + 1);
    LambdaX[freq] = gamma * (LambdaX[freq] + 1.0);
   }//if (normRxTmp[freq] != 0)

}



inline void MLDR::Clear(){
  int channel, channel2, sample, freq;

  memset(phiY, 0, (nFreq) * sizeof(double));
  for (freq = 0; freq < nFreq * 2; freq++) {
    memset(WNumer[freq], 0, sizeof(double) * nChannel);
  }
  for (freq = 0; freq < nFreq * 2; freq++) {
    memset(XPsi[freq], 0, sizeof(double) * nChannel);
  }

  for (freq = 0; freq < nFreq*2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(PsiRx[freq][channel], 0, sizeof(double) * nChannel);
    }
  }

  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(numer[freq][channel], 0, sizeof(double) * nChannel);
    }
  }

  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(PsiTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  for (freq = 0; freq < nFreq*2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(Rx[freq][channel], 0, (nChannel) * sizeof(double));
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        if (channel == channel2) {
          if(freq%2==0)
            Rx[freq][channel][channel2] = delta_SVE;
        }
      }
    }
  }
for (freq = 0; freq < nFreq*2; freq++) {
  for (channel = 0; channel < nChannel; channel++) {
    memset(Rn[freq][channel], 0, (nChannel) * sizeof(double));
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      if (channel == channel2) {
        if (freq % 2 == 0)
          Rn[freq][channel][channel2] = delta_SVE;
      }
    }
  }
}

  memset(denom, 0, sizeof(double) * nFreq);

  memset(WDenom, 0, sizeof(double) * nFreq);
for (freq = 0; freq < fftSize + 2; freq++)
	{
		for (channel = 0; channel < nChannel; channel++){
			memset(Psi[freq][channel], 0, (nChannel) * sizeof(double));
			
			if (freq % 2 == 0)
			{
				Psi[freq][channel][channel] = delta_MLDR;
			}
		}
	}

  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(RxTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(RnTmp[freq][channel], 0, sizeof(double)* nChannel);
    }
  }
  for (freq = 0; freq < nFreq; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      maxEigVec[freq+ freq][channel] = 1.0/nChannel;
      maxEigVec[freq+ freq+1][channel] = 0.0;
    }
  }
 
  for (freq = 0; freq < nFreq * 2; freq++) {
    memset(st[freq], 0, sizeof(double) * nChannel);
  }

  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(XFrameTmp[freq][channel], 0, sizeof(double));
    }
  }
  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(RxTmp2[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(RnTmp2[freq][channel], 0, sizeof(double) * nChannel);
    }
  }
  for (freq = 0; freq < nFreq * 2; freq++) {
    for (channel = 0; channel < nChannel; channel++) {
      memset(eigTmp[freq][channel], 0, sizeof(double) * nChannel);
    }
  }

  memset(LambdaX, 0, sizeof(double)* nFreq);

  memset(LambdaN, 0, sizeof(double)* nFreq);

  memset(SHat,0,sizeof(double)*(fftSize+2));

  for(channel=0;channel <nChannel;channel++){
    for(freq=0;freq<nFreq;freq++){
      steerVector[channel][freq+freq]=1.0/nChannel;
      steerVector[channel][freq+freq+1]=0.0;
    }
  }

}

void MLDR::MLDR_routine_v2(double** X, double* SHat, double** steerVector, int freq) {
  int channel, channel2, channel3, re, im;
  double t_re, t_im;
  re = freq + freq;
  im = re + 1;

  //[M] PhiRxPhi = Phi(:, : , freq) * XFrame(:, freq) 
  //               * (XFrame(:, freq)')*Phi(:,:,freq);

  for (channel = 0; channel < nChannel; channel++) {
    t_re = 0.0;
    t_im = 0.0;
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      t_re += Phi[re][channel][channel2] * X[channel2][re]
        - Phi[im][channel][channel2] * X[channel2][im];
      t_im += Phi[re][channel][channel2] * X[channel2][im]
        + Phi[im][channel][channel2] * X[channel2][re];
    }
    PhiRxPhiTmp[re][channel] = t_re;
    PhiRxPhiTmp[im][channel] = t_im;
  }
    for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      PhiRxPhiTmp2[re][channel][channel2] = PhiRxPhiTmp[re][channel] * X[channel2][re]
        + PhiRxPhiTmp[im][channel] * X[channel2][im];
      PhiRxPhiTmp2[im][channel][channel2] = -PhiRxPhiTmp[re][channel] * X[channel2][im]
        + PhiRxPhiTmp[im][channel] * X[channel2][re];
    }
  }
  
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      t_re = 0.0;
      t_im = 0.0;
      for (channel3 = 0; channel3 < nChannel; channel3++) {
        t_re += PhiRxPhiTmp2[re][channel][channel3] * Phi[re][channel3][channel2]
          - PhiRxPhiTmp2[im][channel][channel3] * Phi[im][channel3][channel2];
        t_im += PhiRxPhiTmp2[re][channel][channel3] * Phi[im][channel3][channel2]
          + PhiRxPhiTmp2[im][channel][channel3] * Phi[re][channel3][channel2];
      }
      PhiRxPhi[re][channel][channel2] = t_re;
      PhiRxPhi[im][channel][channel2] = t_im;
    }
  }

  //[M]   XPhiX = (XFrame(:, freq)')*Phi(:,:,freq)*XFrame(:,freq);  1x6 6x6 6x1 : 1x1

  for (channel = 0; channel < nChannel; channel++) {
    t_re = 0.0;
    t_im = 0.0;
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      t_re += X[channel2][re] * Phi[re][channel2][channel]
        + X[channel2][im] * Phi[im][channel2][channel];
      t_im += X[channel2][re] * Phi[im][channel2][channel]
        - X[channel2][im] * Phi[re][channel2][channel];
    }
    XPhiXTmp[re][channel] = t_re;
    XPhiXTmp[im][channel] = t_im;
  }
  XPhiX[re] = 0.0;
  XPhiX[im] = 0.0;
  for (channel = 0; channel < nChannel; channel++) {
    XPhiX[re] += XPhiXTmp[re][channel] * X[channel][re] - XPhiXTmp[im][channel] * X[channel][im];
    XPhiX[im] += XPhiXTmp[re][channel] * X[channel][im] + XPhiXTmp[im][channel] * X[channel][re];
  }
  //[M]     YHat = W(:, freq)'*XFrame(:,freq); 1x6 6x1 : 1x1
  YHat[re] = 0.0;
  YHat[im] = 0.0;
  for (channel = 0; channel < nChannel; channel++) {
    YHat[re] += W[channel][re] * X[channel][re]
      + W[channel][im] * X[channel][im];
    YHat[im] += +W[channel][re] * X[channel][im]
      - W[channel][im] * X[channel][re];
  }
  //[M]    PsiYHat = alpha * PsiY(freq) + (1 - alpha) * abs(YHat). ^ 2;
  PsiYHat[freq] = alpha * PsiY[freq] + (1 - alpha) * (YHat[re] * YHat[re] + YHat[im] * YHat[im]);

  //[M] lambdaY = max(PsiYHat, epsi);
  lambdaY[freq] = std::max(PsiYHat[freq], epsi_MLDR);
  //[M]   denom = gamma * lambdaY + XPhiX;
  denom[freq] = gamma * lambdaY[freq] + XPhiX[re];

  //[M] denom = real(denom);
  //[M]   PhiHat = (Phi(:, : , freq) - (PhiRxPhi. / denom)). / gamma;
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      PhiHatTmp[re][channel][channel2] = (Phi[re][channel][channel2] - (PhiRxPhi[re][channel][channel2] / (denom[freq] + min_pos))) / gamma;
      PhiHatTmp[im][channel][channel2] = (Phi[im][channel][channel2] - (PhiRxPhi[im][channel][channel2] / (denom[freq] + min_pos))) / gamma;
    }
  }

  //[M] PhiHat = 0.5.*(PhiHat + PhiHat');
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      PhiHat[re][channel][channel2] = PhiHatTmp[re][channel][channel2] + PhiHatTmp[re][channel2][channel];
      PhiHat[re][channel][channel2] *= 0.5;
      PhiHat[im][channel][channel2] = PhiHatTmp[im][channel][channel2] - PhiHatTmp[im][channel2][channel];
      PhiHat[im][channel][channel2] *= 0.5;
    }
  }
  
  //[M]   W(:, freq) = PhiHat * st(:, freq)
  //                   ./ real(st(:, freq)'*PhiHat*st(:,freq));
 
  for (channel = 0; channel < nChannel; channel++) {
    t_re = 0.0;
    t_im = 0.0;
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      t_re += PhiHat[re][channel][channel2] * steerVector[channel2][re]
        - PhiHat[im][channel][channel2] * steerVector[channel2][im];
      t_im += PhiHat[re][channel][channel2] * steerVector[channel2][im]
        + PhiHat[im][channel][channel2] * steerVector[channel2][re];
    }
    WNumer[re][channel] = t_re;
    WNumer[im][channel] = t_im;
  }

  for (channel = 0; channel < nChannel; channel++) {
    t_re = 0.0;
    t_im = 0.0;
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      t_re += steerVector[channel2][re] * PhiHat[re][channel2][channel]
        + steerVector[channel2][im] * PhiHat[im][channel2][channel];
      t_im += steerVector[channel2][re] * PhiHat[im][channel2][channel]
        - steerVector[channel2][im] * PhiHat[re][channel2][channel];
    }
    WDenomTmp[re][channel] = t_re;
    WDenomTmp[im][channel] = t_im;
  }

  WDenom[freq] = 0.0;
  for (channel = 0; channel < nChannel; channel++) {
    WDenom[freq] += WDenomTmp[re][channel] * steerVector[channel][re]
      - WDenomTmp[im][channel] * steerVector[channel][im];
  }

  for (channel = 0; channel < nChannel; channel++) {
    W[channel][re] = (WNumer[re][channel]+min_pos) / (WDenom[freq] + min_pos);
    W[channel][im] =( WNumer[im][channel]+min_pos) / (WDenom[freq] + min_pos);
  }

  //[M]  YFrame(freq) = W(:, freq)'*XFrame(:,freq);
  SHat[re] = 0.0;
  SHat[im] = 0.0;
  for (channel = 0; channel < nChannel; channel++) {
    SHat[re] += W[channel][re] * X[channel][re] + W[channel][im] * X[channel][im];
    SHat[im] += W[channel][re] * X[channel][im] - W[channel][im] * X[channel][re];
  }

  //[M] PsiY(freq) = alpha * PsiY(freq) + (1 - alpha) * abs(YFrame(freq)). ^ 2;
  PsiY[freq] = alpha * PsiY[freq] + (1 - alpha) * (SHat[re] * SHat[re] + SHat[im] * SHat[im]);

  //[M] Phi(:, : , freq) = PhiHat; PhiHat À» ¾È¾²°í ´Ù Phi ÇØµµ µÉµí. 
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      Phi[re][channel][channel2] = PhiHat[re][channel][channel2];
      Phi[im][channel][channel2] = PhiHat[im][channel][channel2];
    }
  }
}


inline void MLDR::MLDR_routine(double** X, double* SHat, double** steerVector, int freq) {
  int channel, channel2, channel3, re, im;

  re = freq + freq;
  im = re + 1;

  phiY[freq] = alpha * phiY[freq] + (1 - alpha) * (YPred[re] * YPred[re] + YPred[im] * YPred[im]);
  lambdaY[freq] = std::max(phiY[freq], epsi_MLDR);

#if !_SVE
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      XFrameTmp[re][channel][channel2] = X[channel][re] * X[channel2][re] + X[channel][im] * X[channel2][im];
      XFrameTmp[im][channel][channel2] = -X[channel][re] * X[channel2][im] + X[channel][im] * X[channel2][re];
    }
  }
#endif



  // numer
  for (channel = 0; channel < nChannel; channel++)
  {
    memset(PsiRx[re][channel], 0, (nChannel) * sizeof(double));
    memset(PsiRx[im][channel], 0, (nChannel) * sizeof(double));
    memset(numer[re][channel], 0, (nChannel) * sizeof(double));
    memset(numer[im][channel], 0, (nChannel) * sizeof(double));

    // Psi * Rx => PsiRx
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      for (channel3 = 0; channel3 < nChannel; channel3++)
      {
        PsiRx[re][channel][channel2] += Psi[re][channel][channel3] * XFrameTmp[re][channel3][channel2] - Psi[im][channel][channel3] * XFrameTmp[im][channel3][channel2];
        PsiRx[im][channel][channel2] += Psi[re][channel][channel3] * XFrameTmp[im][channel3][channel2] + Psi[im][channel][channel3] * XFrameTmp[re][channel3][channel2];

        //PsiRx[re][channel][channel2] += Psi[re][channel][channel3] * XFrameTmp[freq][channel3][channel2+ channel2]
        //  - Psi[im][channel][channel3 + channel3+1] * XFrameTmp[freq][channel3][channel2+channel2+1];
        //PsiRx[im][channel][channel2] += Psi[re][channel][channel3] * XFrameTmp[freq][channel3][channel2+ channel2+1]
        //  + Psi[im][channel][channel3] * XFrameTmp[freq][channel3][channel2+ channel2];

      }
    }

    // PsiRx * Psi => numer
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      for (channel3 = 0; channel3 < nChannel; channel3++)
      {
        numer[re][channel][channel2] += PsiRx[re][channel][channel3] * Psi[re][channel3][channel2] - PsiRx[im][channel][channel3] * Psi[im][channel3][channel2];
        numer[im][channel][channel2] += PsiRx[re][channel][channel3] * Psi[im][channel3][channel2] + PsiRx[im][channel][channel3] * Psi[re][channel3][channel2];
      }
    }
  }

  //[M]   denom = gamma * lambdaY + XPhiX;    
  // denom
  // X' * Psi = XPsi
  memset(XPsi[re], 0, (nChannel) * sizeof(double));
  memset(XPsi[im], 0, (nChannel) * sizeof(double));
  for (channel = 0; channel < nChannel; channel++)
  {
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      XPsi[re][channel] += X[channel2][re] * Psi[re][channel2][channel] + X[channel2][im] * Psi[im][channel2][channel];
      XPsi[im][channel] += X[channel2][re] * Psi[im][channel2][channel] - X[channel2][im] * Psi[re][channel2][channel];
    }
  }

  // XPsi * X = denom (real)
  denom[freq] = 0.0;
  for (channel = 0; channel < nChannel; channel++)
  {
    denom[freq] += XPsi[re][channel] * X[channel][re] - XPsi[im][channel] * X[channel][im];
  }

  denom[freq] += gamma * lambdaY[freq];

  //[M]PsiHat = (Psi(:,:,freq) - (PhiRxPhi./denom))./gamma;
  // Psi = PsiTmp
  for (channel = 0; channel < nChannel; channel++)
  {
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      PsiTmp[re][channel][channel2] = (Psi[re][channel][channel2] - numer[re][channel][channel2] / denom[freq]) / gamma;
      PsiTmp[im][channel][channel2] = (Psi[im][channel][channel2] - numer[im][channel][channel2] / denom[freq]) / gamma;

    }
  }

  //[M]PsiHat = 0.5.*(PsiHat+PsiHat');
  // Psi hermitian symmetry
  for (channel = 0; channel < nChannel; channel++)
  {
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      Psi[re][channel][channel2] = 0.5 * (PsiTmp[re][channel][channel2] + PsiTmp[re][channel2][channel]);

      if (channel != channel2)
      {
        Psi[im][channel][channel2] = 0.5 * (PsiTmp[im][channel][channel2] - PsiTmp[im][channel2][channel]);
      }
      else
      {
        Psi[im][channel][channel2] = 0;
      }

    }
  }

  // Phi * steerVector => WNumer
  memset(WNumer[re], 0, (nChannel) * sizeof(double));
  memset(WNumer[im], 0, (nChannel) * sizeof(double));
  for (channel = 0; channel < nChannel; channel++)
  {
    for (channel2 = 0; channel2 < nChannel; channel2++)
    {
      WNumer[re][channel] += Psi[re][channel][channel2] * steerVector[channel2][re]
        - Psi[im][channel][channel2] * steerVector[channel2][im];
      WNumer[im][channel] += Psi[re][channel][channel2] * steerVector[channel2][im]
        + Psi[im][channel][channel2] * steerVector[channel2][re];
    }
  }

  // steerVector' * WNumer(Phi * steerVector) => WDenom (real)
  WDenom[freq] = 0.0;
  for (channel = 0; channel < nChannel; channel++)
  {
    WDenom[freq] += steerVector[channel][re] * WNumer[re][channel] + steerVector[channel][im] * WNumer[im][channel];

  }


  // W		
  for (channel = 0; channel < nChannel; channel++) {
    W[channel][re] = (WNumer[re][channel] + min_pos) / (WDenom[freq] + min_pos);
    W[channel][im] = (WNumer[im][channel] + min_pos) / (WDenom[freq] + min_pos);
  }


  // W ' X = SHat 

  SHat[re] = 0.0;
  SHat[im] = 0.0;
  for (channel = 0; channel < nChannel; channel++)
  {
    SHat[re] += W[channel][re] * X[channel][re] + W[channel][im] * X[channel][im];
    SHat[im] += W[channel][re] * X[channel][im] - W[channel][im] * X[channel][re];

  }


  X[0][re] = SHat[re];
  X[0][im] = SHat[im];

}


inline void MLDR::SVE_routine(double** X, double* SHat, double** steerVector, int freq) {
  int channel, channel2, channel3, re, im;
  re = freq + freq;
  im = freq + freq + 1;
  //[M]PhiYHat = abs(YHat). ^ 2;
  PhiYHat[freq] = YPred[re] * YPred[re] + YPred[im] * YPred[im];
  //[M]lambdaY = max(PhiYHat, epsiSt);
  lambdaY[freq] = std::max(PhiYHat[freq], epsi_SVE);
  //[M]lambdaYInv = 1 ./ lambdaY;
  lambdaYInv[freq] = 1.0 / (lambdaY[freq]);

  // XFrameTmp = (XFrame(:,freq)*XFrame(:,freq)')
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      XFrameTmp[re][channel][channel2] = X[channel][re] * X[channel2][re] + X[channel][im] * X[channel2][im];
      XFrameTmp[im][channel][channel2] = -X[channel][re] * X[channel2][im] + X[channel][im] * X[channel2][re];
    }
  }

  //RxTmp = Rx(:,:,freq) * LambdaX(freq) / (LambdaX(freq) + 1) +   (XFrame(:,freq)*XFrame(:,freq)') / (LambdaX(freq) + 1);
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      //real
      RxTmp[re][channel][channel2] =
        Rx[re][channel][channel2] * LambdaX[freq] / (LambdaX[freq] + 1.0)
        + XFrameTmp[re][channel][channel2] / (LambdaX[freq] + 1.0);
      //imag
      RxTmp[im][channel][channel2] =
        Rx[im][channel][channel2] * LambdaX[freq] / (LambdaX[freq] + 1.0)
        + XFrameTmp[im][channel][channel2] / (LambdaX[freq] + 1.0);
    }
  }

  //RnTmp = Rn(:,:,freq) * LambdaN(freq)  / (LambdaN(freq) + lambdaYInv) +  (XFrame(:,freq) * XFrame(:,freq)') / (LambdaN(freq) + lambdaYInv);
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      //real
      RnTmp[re][channel][channel2] = Rn[re][channel][channel2] * LambdaN[freq] / (LambdaN[freq] + lambdaYInv[freq])
        + XFrameTmp[re][channel][channel2] / (LambdaN[freq] + lambdaYInv[freq]);
      //imag
      RnTmp[im][channel][channel2] = Rn[im][channel][channel2] * LambdaN[freq] / (LambdaN[freq] + lambdaYInv[freq])
        + XFrameTmp[im][channel][channel2] / (LambdaN[freq] + lambdaYInv[freq]);
      // printf("%e %+ei ", RnTmp[channel][channel2 + channel2], RnTmp[channel][channel2 + channel2 + 1]);
    }
    // printf("\n");
  }


  //  RxTmp = 0.5 * (RxTmp + RxTmp');
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      RxTmp2[re][channel][channel2] = RxTmp[re][channel][channel2] + RxTmp[re][channel2][channel];
      RxTmp2[re][channel][channel2] *= 0.5;
      if (channel != channel2) {
        RxTmp2[im][channel][channel2] = RxTmp[im][channel][channel2] - RxTmp[im][channel2][channel];
        RxTmp2[im][channel][channel2] *= 0.5;
      }
      else {
        RxTmp2[im][channel][channel2] = 0.0;
      }
    }
  }

  //  RnTmp = 0.5 * (RnTmp + RnTmp');
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      if (channel != channel2) {
        RnTmp2[re][channel][channel2] = RnTmp[re][channel][channel2] + RnTmp[re][channel2][channel];
        RnTmp2[re][channel][channel2] *= 0.5;
        RnTmp2[im][channel][channel2] = RnTmp[im][channel][channel2] - RnTmp[im][channel2][channel];
        RnTmp2[im][channel][channel2] *= 0.5;
      }
      else {
        RnTmp2[re][channel][channel2] = RnTmp[re][channel][channel2];
        RnTmp2[im][channel][channel2] = 0.0;
      }
    }
  }
  //norm(RxTmp)
  normRxTmp[freq] = 0.0;
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      normRxTmp[freq]
        += RxTmp[re][channel][channel2] * RxTmp[re][channel][channel2]
        + RxTmp[im][channel][channel2] * RxTmp[im][channel][channel2];
    }
  }

  //if norm(RxTmp) ~= 0   % to prevent division by zero (20201213)
  if (normRxTmp[freq] != 0) {
    // for iterPM = 1 : maxIterPM
    // maxEigVec(:,freq) = (RxTmp - theta*RnTmp)*maxEigVec(:,freq);
    // maxEigVec(:, freq) = maxEigVec(:, freq). / norm(maxEigVec(:, freq));
    // end
    for (channel = 0; channel < nChannel; channel++) {
      for (channel2 = 0; channel2 < nChannel; channel2++) {
        eigTmp[re][channel][channel2] = RxTmp2[re][channel][channel2] - theta_SVE * RnTmp2[re][channel][channel2];
        eigTmp[im][channel][channel2] = RxTmp2[im][channel][channel2] - theta_SVE * RnTmp2[im][channel][channel2];
      }
    }
    for (int i = 0; i < power_method_iter; i++) {
      for (channel = 0; channel < nChannel; channel++) {
        tmpRe[freq] = 0.0;
        tmpIm[freq] = 0.0;
        for (channel2 = 0; channel2 < nChannel; channel2++) {
          tmpRe[freq] += eigTmp[re][channel][channel2] * maxEigVec[re][channel2]
            - eigTmp[im][channel][channel2] * maxEigVec[im][channel2];
          tmpIm[freq] += eigTmp[im][channel][channel2] * maxEigVec[re][channel2]
            + eigTmp[re][channel][channel2] * maxEigVec[im][channel2];
        }
        eigVecTmp[re][channel] = tmpRe[freq];
        eigVecTmp[im][channel] = tmpIm[freq];
      }

      for (channel = 0; channel < nChannel; channel++) {
        maxEigVec[re][channel] = eigVecTmp[re][channel];
        maxEigVec[im][channel] = eigVecTmp[im][channel];
      }
      //norm
      tmpRe[freq] = 0.0;
      for (channel = 0; channel < nChannel; channel++)
        tmpRe[freq] += maxEigVec[re][channel] * maxEigVec[re][channel]
        + maxEigVec[im][channel] * maxEigVec[im][channel];
      tmpRe[freq] = sqrt(tmpRe[freq]);

      for (channel = 0; channel < nChannel; channel++) {
        maxEigVec[re][channel] /= (tmpRe[freq] + min_pos);
        maxEigVec[im][channel] /= (tmpRe[freq] + min_pos);
      }
    }
    // st(:,freq) = maxEigVec(:,freq)./maxEigVec(1,freq);
    tmpRe[freq] = maxEigVec[re][0] * maxEigVec[re][0] + maxEigVec[im][0] * maxEigVec[im][0];
    for (channel = 0; channel < nChannel; channel++) {
      steerVector[channel][re] = (maxEigVec[re][channel] * maxEigVec[re][0] + maxEigVec[im][channel] * maxEigVec[im][0]) / (tmpRe[freq] + min_pos);
      steerVector[channel][im] = (maxEigVec[im][channel] * maxEigVec[re][0] - maxEigVec[re][channel] * maxEigVec[im][0]) / (tmpRe[freq] + min_pos);
    }
    // st(:, freq) = st(:, freq).*sqrt(nSensor / norm(st(:, freq)));
    tmpRe[freq] = 0.0;
    for (channel = 0; channel < nChannel; channel++)
      tmpRe[freq] += steerVector[channel][re] * steerVector[channel][re] + steerVector[channel][im] * steerVector[channel][im];
    tmpRe[freq] = sqrt(tmpRe[freq]);
    tmpRe[freq] = sqrt(nChannel / (tmpRe[freq] + min_pos));
    for (channel = 0; channel < nChannel; channel++) {
      steerVector[channel][re] = steerVector[channel][re] * tmpRe[freq];
      steerVector[channel][im] = steerVector[channel][im] * tmpRe[freq];
    }
  }
  // else
  // st(:, freq) = ones(nSensor, 1) . / nSensor;
  else {
    for (channel = 0; channel < nChannel; channel++) {
      steerVector[channel][re] = 1 / nChannel;
      steerVector[channel][im] = 0.0;
    }
  }


  //Rx(:,:,freq) = RxTmp;
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      Rx[re][channel][channel2] = RxTmp2[re][channel][channel2];
      Rx[im][channel][channel2] = RxTmp2[im][channel][channel2];
    }
  }

  //Rn(:, : , freq) = RnTmp;
  for (channel = 0; channel < nChannel; channel++) {
    for (channel2 = 0; channel2 < nChannel; channel2++) {
      Rn[re][channel][channel2] = RnTmp2[re][channel][channel2];
      Rn[im][channel][channel2] = RnTmp2[im][channel][channel2];
    }
  }
  LambdaN[freq] = gamma * (LambdaN[freq] + lambdaYInv[freq]);
  LambdaX[freq] = gamma * (LambdaX[freq] + 1.0);

}

