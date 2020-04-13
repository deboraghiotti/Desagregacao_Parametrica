require(moments)
require(corrgram)
require(pear)
require(gridExtra)
require(knitr)
require(tcltk)
source('C:\Users\debor\Desktop\IC-2019-2020\Desagregacao\Desagregacao_Parametrica\Funcoes_abr2020.R')

div_mensais<-function(sH)
{
  qtd_ano_hist = length(sH[,1])/12                        
  serie_hist = matrix(sH$VAZAO, qtd_ano_hist,byrow = TRUE)
  anos = as.character(sH$MES)
  anos = substr(anos, nchar(anos)-4+1, nchar(anos))
  anos = unique(anos)
  anos = sort(anos)                                       
  
  row.names(serie_hist)= anos
  colnames(serie_hist)=c("JAN","FEV","MAR","ABR","MAI","JUN","JUL","AGO","SET","OUT","NOV","DEZ")
  serie_hist=as.data.frame(serie_hist)                    
  
  return(serie_hist)
}

desag_param_info<-function(){
  Info=list()
  for (i in 1:11){
    if (i==1)
      Info[[i]] = cbind(SeriesDadosHist_normalizada[,i],SeriesDadosHist_normalizada[,12],apply(SeriesDadosHist_normalizada[,-1],1,sum),apply(SeriesDadosHist_normalizada,1,sum))
    else if(i<11)
      Info[[i]] = cbind(SeriesDadosHist_normalizada[,i],SeriesDadosHist_normalizada[,i-1],apply(SeriesDadosHist_normalizada[,((i+1):12)],1,sum),apply(SeriesDadosHist_normalizada[,(i:12)],1,sum))
    else
      Info[[i]] = cbind(SeriesDadosHist_normalizada[,i],SeriesDadosHist_normalizada[,i-1],SeriesDadosHist_normalizada[,12],apply(SeriesDadosHist_normalizada[,(i:12)],1,sum))
  }
  return(Info)
}

desag_param_mult<-function(info){
  aux = (serie_sintetica_normalizada)
  
  Desagregado = list()
  if(input_3<=0){
    print("Numero invalido de series a se desagregar")
    break
  }
  for(j in 1:input_3){
    CONT_ERRO <<- 0
    ERRO_AUX = rnorm (10000*2, 0, 1)
    ERRO_AUX = (ERRO_AUX - mean (ERRO_AUX)) / sd (ERRO_AUX)
    aux1 = matrix (ERRO_AUX, nrow = 2)
    ERRO <<- list()
    for(i in 1:10000) {
      ERRO[[i]] <<- aux1[, i]
    }
    
    ################MUDAR PARA OS ANOS FINAIS AO INVES DOS ANOS INICIAIS###########
    df_serie_desagregada_param = data.frame(matrix(NA, nrow = 10000, ncol = 12))
    cont = 1
    inicio =((10000*(j-1))+1)
    fim = 10000*j
    Zinicial = SeriesDadosHist_normalizada[1,12]
    for(i in inicio:fim){
      df_serie_desagregada_param[cont,] = desag_param(info, aux[i,1], Zinicial)
      Zinicial = df_serie_desagregada_param[cont,12]
      cont = cont+1
    }
    Desagregado[[j]] =  df_serie_desagregada_param
  }
  return(Desagregado)
}


#########Função que desegrega
desag_param<-function(info, ano, Zinicial){
  Meses = 0
  Resto = 0
  for(i in 1:11){
    CONT_ERRO = CONT_ERRO+1
    if(i==1){
      ACF_S = acf(info[[i]],lag.max = 1,type = "covariance", plot = FALSE)
      print(ACF_S)
      ACF_S = ACF_S$acf
      Sxx = ACF_S[1,4,4] # Anual com Anual
      Sxx_ = ACF_S[2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[1,4,1],ACF_S[1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[1,2,1],ACF_S[1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[1,1,1],ACF_S[1,3,1]), rbind(ACF_S[1,1,3],ACF_S[1,3,3]))
      Sx_z = ACF_S[2,2,4] ##### Anual com //CONFERIR
      Sxz_cor = Sxx_%*%(solve(Sxx))%*%(Sx_z)
      Syz_cor = t(Syz)+(t(Syx)%*%solve(Sxx)%*%(Sxz_cor-Sxz))
      A = (t(Syx)-Syz_cor%*%solve(Szz)%*%(Sxz_cor))%*%solve(Sxx-Sxz_cor%*%solve(Szz)%*%Sxz_cor)
      C = (Syz_cor-A%*%Sxz_cor)%*%solve(Szz)
      BBt = (Syy)-(A)%*%(Syx)-C%*%t(Syz_cor)
      B = chol(BBt,pivot = TRUE)
      #print(det(BBt))
      Bt = t(B)
      print(A);print(C)
      print("Bt");print(Bt)
      #erro = rbind(rnorm(1,0,1),rnorm(1,0,1))
      erro = ERRO[[CONT_ERRO]]
      print("i"); print(i)
      print("erro"); print(erro)
      Y = A%*%ano+Bt%*%(erro)+C%*%Zinicial #Y = A%*%serie_sint[k]+Bt%*%erro+C%*%Meses
      Meses[i] = Y[1,1]
      Resto[i] = Y[2,1]
      
    }
    else{
      ACF_S = acf(info[[i]],lag.max = 1,type = "covariance", plot = FALSE)
      #print(ACF_S)
      ACF_S = ACF_S$acf
      Sxx = ACF_S[1,4,4] # Anual com Anual
      Sxx_ = ACF_S[2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[1,4,1],ACF_S[1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[1,2,1],ACF_S[1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[1,1,1],ACF_S[1,3,1]), rbind(ACF_S[1,1,3],ACF_S[1,3,3]))
      Sx_z = ACF_S[2,2,4] ##### Anual com //CONFERIR
      print(Sxx);print(Sxx_);print(Syx);print(Syz);print(Sxz);print(Szz);print(Syy);print(Sx_z);
      A = (t(Syx)-t(Syz)%*%solve(Szz)%*%(Sxz))%*%solve(Sxx-Sxz%*%solve(Szz)%*%Sxz)
      print("A")
      print(A)
      C = (t(Syz) - A%*%Sxz)%*%solve(Szz)
      print(C)
      BBt = (Syy)-(A)%*%(Syx)-C%*%(Syz)
      #print("BBt")
      print(BBt)
      #print(det(BBt))
      B = chol(BBt,pivot=TRUE)
      print("B")
      Bt = t(B)
      print(Bt)
      #  erro = rbind(rnorm(1,0,1),rnorm(1,0,1))
      erro = ERRO[[CONT_ERRO]]
      #      print("i");print(i)
      #      print("erro"); print(erro)
      Y = A%*%Resto[i-1]+Bt%*%erro+C%*%Meses[i-1]
      Meses[i] = Y[1,1]
      Resto[i] = Y[2,1]
    }
    #    if (Meses[i] < 0) {
    #      print(paste ("meses", Meses[i]))
    #      print (Zinicial)
    #      print (i)
    #      print("erro"); print(erro)
    #      print("Bt")
    #      print(Bt)
    #      print("")
    #      print("")
    #    }
    
  }
  Meses[12] = Resto[11]
  print(Meses)
  #print(sum(Meses))
  return(Meses)
}

input_1 <- tk_choose.files()#Arquivo de serie historica
input_2 <- tk_choose.files()#Arquivo de vazao sintetica

entrada <- data.frame(read.csv(input_1, sep=";", dec=",")) #Leitura de dados historicos mensais
serie_sintetica <- data.frame(read.csv(input_2, header = F,sep=";", dec=",")) #Leitura de serie sintetica gerada para essa mesma serie

SeriesDadosHist <- div_mensais(entrada)
SeriesDadosHist_normalizada <- apply(log(SeriesDadosHist),2, function(x) (x-mean(x)))

Parametro_Hist <- desag_param_info()
ParametroP_Hist <- desag_param_info_padronizado()

DesagregadoP = desag_param_mult(Parametro_Hist)
for(i in 1:input_3){
  DesagregadoP[[i]] = abs(DesagregadoP[[i]])
}
