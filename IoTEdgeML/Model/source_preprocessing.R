if(!"seewave" %in% rownames(installed.packages())) {
  install.packages("seewave")
}
if(!"mvtnorm" %in% rownames(installed.packages())) {
  install.packages("mvtnorm")
}
if(!"caret" %in% rownames(installed.packages())) {
  install.packages("caret")
}
if(!"mvnmle" %in% rownames(installed.packages())) {
  install.packages("mvnmle")
}

library(seewave)
library(mvtnorm)
library(caret)
library(mvnmle)

## master_raw 의 하나 제품 받아서 feature 생성
get_features_from_master <- function(df_raw){
  ## Fixed parameters
  Fs = 51200
  
  ## 1. 함수 정의 ####
  ## Preprocessing funtions
  ## raw데이터에서 스위치 on/off 시점을 찾는 함수
  get_start_time <- function(x){
    # 전압값을 기준으로 시작시점 추출
    on_off <- rep("", length(x))
    tmp <- paste0(ifelse(x>=4.00, 1,0), collapse="")
    vol_on <- gregexpr("01", tmp)[[1]]+1
    vol_off <- gregexpr("10", tmp)[[1]]
    
    on_off[vol_on] <- paste("ON", c(1,2,3), sep="_")
    on_off[vol_off] <- paste("OFF", c(1,2,3), sep="_")
    
    return(on_off) 
  }
  
  ## raw데이터에서 스위치 on/off 여부를 테스트 시점별로 구분하여 라벨링하는 함수
  get_using_time <- function(x){
    # 전압값을 기준으로 시작시점 추출
    on_off <- rep("0", length(x))
    tmp <- paste0(ifelse(x>=4.00, 1,0), collapse="")
    vol_on <- gregexpr("01", tmp)[[1]]+1
    vol_off <- gregexpr("10", tmp)[[1]]
    
    on_off[vol_on[1]:vol_off[1]] <- "1"
    on_off[vol_on[2]:vol_off[2]] <- "2"
    on_off[vol_on[3]:vol_off[3]] <- "3"
    
    return(on_off) 
  }
  
  ## Feature generation functions
  ## 1. Raw데이터를 받아 dbA 값 계산(해당 값 검증 필요)
  get_dBA <- function(x, ab, Fs){
    fft_noise <- fft(x)
    Fs_approx <- 1:length(x)*Fs/length(x)
    fft_noise <- fft_noise[Fs_approx>=10 &Fs_approx<=20000]
    ifft_noise <- Re(fft(fft_noise, inverse = TRUE)/length(x))
    
    # ab <- A_weighting(fs=Fs, output="ab")
    ifft_noise_weighted <- signal::filter(filt=Arma(ab$b, ab$a), x=ifft_noise)
    x_dB <- 20*log10(sqrt(mean(ifft_noise_weighted^2, na.rm=TRUE)))
    
    return(x_dB)
  }
  
  ## 2. Raw데이터를 받아 dbA 값 계산(해당 값 검증 필요)
  get_dBA_new <- function(x, Fs){
    tmp_wav <- oscillo(x, f=Fs, plot=FALSE)
    tmp_spec <- seewave::meanspec(tmp_wav, f=Fs, dB="A", dBref=1e-6, plot=FALSE)
    # tmp_specprop <- specprop(tmp_spec, plot=FALSE)
    
    dba <- 10*log10(sum(10^(tmp_spec[,2][!is.na(tmp_spec[,2])]/10)))
    
    return(dba)
  }
  
  ## 3. Raw데이터를 받아 Spectrum으로 변환 후, Kurtosis(첨도) 값 계산
  get_Kurtosis <- function(x, Fs){
    tmp_wav <- oscillo(x, f=Fs, plot=FALSE)
    tmp_spec <- seewave::meanspec(tmp_wav, f=Fs, plot=FALSE)
    tmp_specprop <- specprop(tmp_spec, plot=FALSE)
    
    return(tmp_specprop$kurtosis)
  }
  
  ## 4. Raw데이터를 받아 Spectrum으로 변환 후, Skewness(왜도) 값 계산
  get_Skewness <- function(x, Fs){
    tmp_wav <- oscillo(x, f=Fs, plot=FALSE)
    tmp_spec <- seewave::meanspec(tmp_wav, f=Fs, plot=FALSE)
    tmp_specprop <- specprop(tmp_spec, plot=FALSE)
    
    return(tmp_specprop$skewness)
  }
  
  ## 4. Raw데이터를 받아 Spectrum으로 변환 후, sfm(주파수 스펙트럼의 평편한 정도) 값 계산
  get_Sfm <- function(x, Fs){
    tmp_wav <- oscillo(x, f=Fs, plot=FALSE)
    tmp_spec <- seewave::meanspec(tmp_wav, f=Fs, plot=FALSE)
    tmp_specprop <- specprop(tmp_spec, plot=FALSE)
    
    return(tmp_specprop$sfm)
  }
  
  ## 4. Raw데이터를 받아 Spectrum으로 변환 후, sh(주파수 스펙트럼의 entropy) 값 계산
  get_Sh <- function(x, Fs){
    tmp_wav <- oscillo(x, f=Fs, plot=FALSE)
    tmp_spec <- seewave::meanspec(tmp_wav, f=Fs, plot=FALSE)
    tmp_specprop <- specprop(tmp_spec, plot=FALSE)
    
    return(tmp_specprop$sh)
  }
  
  ## 5. A-weighting 함수
  A_weighting <- function(fs, output="zpk", curve="A"){
    #https://github.com/endolith/waveform_analysis/blob/master/waveform_analysis/weighting_filters/ABC_weighting.py#L29
    
    library(signal)
    library(control)
    
    rms_flat_ <- function(x, na.rm=TRUE){
      x <- unlist(x)
      sqrt(mean(x^2, na.rm=na.rm))
    }
    
    relative_degree_ <- function(z, p){
      degree = length(p) - length(z)
      if (degree < 0){
        stop("Improper transfer function.\n Must have at least as many poles as zeros.")
      }else{
        return(degree)
      }
    }
    
    zpkbilinear_ <- function(z, p, k, fs){
      degree = relative_degree_(z, p)
      
      fs2 = 2.0*fs
      
      # Bilinear transform the poles and zeros
      z_z = (fs2 + z) / (fs2 - z)
      p_z = (fs2 + p) / (fs2 - p)
      
      # Any zeros that were at infinity get moved to the Nyquist frequency
      z_z <- c(z_z, rep(-1, degree))
      
      # Compensate for gain change
      k_z = k * Re(prod(fs2 - z) / prod(fs2 - p))
      
      return(list(z=z_z, p=p_z, k=k_z))
    }
    
    ABC_weighting <- function(curve='A'){
      library(signal)
      
      z = c(0, 0, 0, 0)
      p = c(-2*pi*20.598997057568145,
            -2*pi*20.598997057568145,
            -2*pi*12194.21714799801,
            -2*pi*12194.21714799801,
            -2*pi*107.65264864304628,
            -2*pi*737.8622307362899)
      k = 1.0
      
      tmp_tf <- zp2tf(z,p,k)
      b = tmp_tf$num%>%as.vector()
      a = tmp_tf$den%>%as.vector()
      
      
      tmp_fr <- freqs(filt=Arma(b=b, a=a), W=2*pi*1000)
      k <- k/abs(tmp_fr$H)
      
      return(list(z=z, p=p, k=k))
    }
    
    abc_weight <- ABC_weighting(curve=curve)
    
    tmp_bl <- zpkbilinear_(z=abc_weight$z, p=abc_weight$p, k=abc_weight$k, fs=fs)
    if(output=="zpk"){
      return(list(z=tmp_bl$z, p=tmp_bl$p, k=tmp_bl$k))
    }else if(output=="ab"){
      
      tmp_tf <- zp2tf(tmp_bl$z,tmp_bl$p,tmp_bl$k)
      b = tmp_tf$num%>%as.vector()
      a = tmp_tf$den%>%as.vector()
      
      return(list(b=b,a=a))
    }else{
      return(NULL)
    }
  }
  
  ## 2. 기본 전처리 수행 ####
  ### 2.1. 테스트별 On/Off 시점 추출
  ### 2.2. 테스트별 On 시점으로부터 0.2초-0.4초 데이터만 추출.
  ### 2.3. 2,3 테스트만 추출
  tmp_df_raw <- df_raw%>%
    rename(x1 = d, x2 = r, x3=n)%>%
    arrange(x1)%>%
    mutate(Time_Sep = get_start_time(x2))%>%
    mutate(Time_Used = get_using_time(x2))%>%
    filter(Time_Used=="2"|Time_Used=="3")%>%
    group_by(Time_Used)%>%
    mutate(Test_Time = x1-min(x1))%>%
    filter(Test_Time>=0.2)%>%
    filter(Test_Time<=0.4)%>%
    ungroup()
  
  ## 3. Feature 생성 ####
  # Test2, Test3에 대해 각각 Feature 생성 후 평균값을 계산하여 최종 Feature로 사용
  # 앞에서 정의한 feature 생성 함수들을 적용
  ab <- A_weighting(fs=Fs, output="ab")
  x2_by_test <- split(tmp_df_raw$x2, tmp_df_raw$Time_Used)
  df_raw_features <- data.frame(
    my_dba = sapply(x2_by_test,get_dBA, ab, Fs)%>%mean(),
    my_dba_new = sapply(x2_by_test,get_dBA_new, Fs)%>%mean(),
    my_kt = sapply(x2_by_test,get_Kurtosis, Fs)%>%mean(),
    my_sk = sapply(x2_by_test,get_Skewness, Fs)%>%mean(),
    my_sfm = sapply(x2_by_test,get_Sfm, Fs)%>%mean(),
    my_sh = sapply(x2_by_test,get_Sh, Fs)%>%mean()
  )
  
  return(df_raw_features)
}

## master_cep 의 하나 제품 받아서 feature 생성
get_features_from_cep <- function(df_cep, models=NULL){
  ## cep feature generation
  ## 1. 함수 정의 ####
  
  ## 6_2_Model_Batch_cepstrum_prob.R 에서 생성한 모델을 활용하여, 주어진 cepstrum의 확률값 계산
  get_upper_probability_each <- function(md, x){
    x <- x[names(x)%in%md$var[[1]]]
    
    if(length(x)==1){
      res<- pnorm(q=x,
                  mean = md$muhat,
                  sd = md$sigmahat)
      
      return(1-res)
    }else{
      res <- pmvnorm(lower=x,
                     upper=rep(Inf, length(x)),
                     mean = md$muhat[[1]],
                     sigma=md$sigmahat[[1]])
      
      return(res[1])
    }
  }
  
  ## 데이터로부터 단,다변량 정규분포를 추정하여 데이터프레임 형태로 데이터를 반환하는 함수
  mlest_tidy <- function(x, columns){
    library(fitdistrplus)
    if(length(columns)==1){
      y <- unlist(x[,columns])
      md <- fitdist(y,"norm", lower=c(0,0))
      return(data.frame(var=names(x[,columns]), muhat=md$estimate[1], sigmahat = md$estimate[1]))
    }else{
      md <- mlest(x[,columns])
      
      return(data.frame(var=I(list(names(x[,columns]))), muhat=I(list(md$muhat)), sigmahat = I(list(md$sigmahat))))
      
    }
  }
  
  ## Cepstrum 데이터의 peak 구간을 찾아 labeling하는 함수
  cep_peak_spot <- function(x, peak_inter=0.0002){
    peak_spots <- ifelse(x>(0.00832-peak_inter)&x<(0.00832+peak_inter), "1peak",
                         ifelse(x>(0.01666-peak_inter)&x<(0.01666+peak_inter), "2peak",
                                ifelse(x>(0.02500-peak_inter)&x<(0.02500+peak_inter), "3peak",
                                       ifelse(x>(0.03332-peak_inter)&x<(0.03332+peak_inter), "4peak", "else"))))
    
    return(peak_spots)
  }
  
  
  ## 2. 기본 전처리 수행 ####
  ### cep데이터로부터 peak부분을 찾아서 평균값 계산
  df_p2_my_cep_tmp <- df_cep%>%
    mutate(peak_spots = cep_peak_spot(q, peak_inter = 0.0003))%>%
    filter(peak_spots!="else")%>%
    group_by(peak_spots)%>%
    summarise(x23_mean=(mean(x2)+mean(x3))/2)%>%
    ungroup()
  
  
  ## 2. Feature 생성 ####
  ### Feature 계산
  x23_mean <- df_p2_my_cep_tmp$x23_mean
  names(x23_mean) <- df_p2_my_cep_tmp$peak_spots
  df_cep_features <- lapply(models, get_upper_probability_each, x23_mean)%>%unlist()%>%t()
  
  ### Feature 이름 정의
  feature_names <- paste("my_cep_prob",
                         lapply(models, function(x) paste(x$var[[1]], collapse="")%>%
                                  str_replace_all("peak", ""))%>%unlist()%>%str_pad(width=4, pad="0"),
                         sep="")
  df_cep_features <- df_cep_features%>%
    as.data.frame()%>%
    set_names(feature_names)
  
  return(df_cep_features)
}
## get_feature_from_master 와 get_feature_from_cep 의 wrapping 함수
get_features <- function(df_raw, df_cep, models=NULL){
  df_raw_f <- get_features_from_master(df_raw=df_raw)
  df_cep_f <- get_features_from_cep(df_cep=df_cep, models=models)
  
  df_f <- cbind(df_raw_f, df_cep_f)
  return(df_f)
}

## 예측 함수
predict_fault <- function(df, line, model){
  df_pred <- df%>%
    mutate(predicted=predict(model, df%>%mutate(l=line)))
  
  return(df_pred)
}