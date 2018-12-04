from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2 import robjects
import tzlocal
import pandas as pd
import json
from azureml.core.model import Model

def init():
    r('''
        ## 1. Load libraries and custom functions ####
        if(!"optparse" %in% rownames(installed.packages())) {
          install.packages("optparse")
        }
        if(!"tidyverse" %in% rownames(installed.packages())) {
          install.packages("tidyverse")
        }
        if(!"xgboost" %in% rownames(installed.packages())) {
          install.packages("xgboost")
        }
        if(!"control" %in% rownames(installed.packages())) {
          install.packages("control")
        }

        library(optparse)
        library(tidyverse)
        library(xgboost)
        library(control)

        #전역변수 : Bind 된 경로로 변경 필요
        model_path <- "./model/"

        preprocessfile <- paste0(model_path,"source_preprocessing.R")
        source(preprocessfile)

        select <- dplyr::select
        filter <- dplyr::filter
        ''')
    

def run(input_str):

    try:

            # input_str이 파일 이름을 가진 json으로 입력된다            input_json = json.loads(input_str)

            #전역변수 : Bind 된 경로로 변경 필요
            data_path = './data/'

            raw_file = data_path + input_json['raw']
            cep_file = data_path + input_json['cep']

            r.assign('line_name',input_json['line'])
            r.assign('raw_file',raw_file)
            r.assign('cep_file',cep_file)

            r('''
                tmp_raw <- read.csv(raw_file)%>%set_names(c("d","r","n"))
                tmp_cep <- read.csv(cep_file)%>%set_names(c("q","x1","x2","x3"))

                ## 3.1 Cep 이상치 확률계산 모델 불러오기
                cep_model_fname <- paste0(model_path,"cep_models_2018-11-06-06-29-32.RData")
                load(cep_model_fname)

                ## 3.2 불량 예측 모델 불러오기
                fd_model_fname <- paste0(model_path,"fault_detection_model_2018-11-06-07-39-47.RData")
                load(fd_model_fname)

                ## 3.3 라인기준 모델 선택
                models <- lapply(models_list, function(x) x[x$l==line_name,])

                ## 3.4 Feature 생성
                df_f <- get_features(tmp_raw, tmp_cep, models=models)

                ## 3.5 예측
                df_pred <- predict_fault(df = df_f, line = line_name, model = fault_detection_model)
             ''')

            r_result = robjects.r['df_pred'][21]
            Input_json['predicted']= str(r_result)[4:6]
            
            print("Prediction is ", str(r_result)[4:6])

    except Exception as e:
        
        result = str(e)
                
    return [json.dumps(input_json)]