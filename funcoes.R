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

# Funcao parametros_historicos: calcula dos parametros historicos necessarios para a desagregacao parametrica
# Parametro: a serie historica normalizada 
# Retorno: parametros historicos para a desagregacao
parametros_historicos<-function(SeriesDadosHist_normalizada){
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

# Funcao autocovariancia 
# Parametro: os parametros historicos calculados pela funcao parametros_historicos
# Retorno: a autocovariancia desses parametros, necessaria para realizar a desagregacao (de janeiro a novembro)
autocovariancia <- function(info){
  ACF_S = list()
  for(i in 1:11){
    ACF_S[[i]]= acf(info[[i]],lag.max = 1,type = "covariance", plot = FALSE)
    ACF_S[[i]]= ACF_S[[i]]$acf
    print(ACF_S[[i]])
  }
  #print(ACF_S)
  return(ACF_S)
}

# Funcao parametro A : Essa funcao calculo o parametro A da formula da desagregacao parametrica
# Parametro: a autocovariancia calculada na funcao autocovariancia
# Retorno: parametro A da desagregacao
parametro_A <- function(ACF_S){
  A = list()
  for(i in 1:11){
    
    if(i == 1){
      
      Sxx = ACF_S[[i]][1,4,4] # Anual com Anual
      Sxx_ = ACF_S[[i]][2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[[i]][1,4,1],ACF_S[[i]][1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[[i]][1,2,1],ACF_S[[i]][1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[[i]][1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[[i]][1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[[i]][1,1,1],ACF_S[[i]][1,3,1]), rbind(ACF_S[[i]][1,1,3],ACF_S[[i]][1,3,3]))
      Sx_z = ACF_S[[i]][2,2,4] ##### Anual com //CONFERIR
      Sxz_cor = Sxx_%*%(solve(Sxx))%*%(Sx_z)
      Syz_cor = t(Syz)+(t(Syx)%*%solve(Sxx)%*%(Sxz_cor-Sxz))
      A[[i]] = (t(Syx)-Syz_cor%*%solve(Szz)%*%(Sxz_cor))%*%solve(Sxx-Sxz_cor%*%solve(Szz)%*%Sxz_cor)
      
    }else{
      
      Sxx = ACF_S[[i]][1,4,4] # Anual com Anual
      Sxx_ = ACF_S[[i]][2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[[i]][1,4,1],ACF_S[[i]][1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[[i]][1,2,1],ACF_S[[i]][1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[[i]][1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[[i]][1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[[i]][1,1,1],ACF_S[[i]][1,3,1]), rbind(ACF_S[[i]][1,1,3],ACF_S[[i]][1,3,3]))
      Sx_z = ACF_S[[i]][2,2,4] ##### Anual com //CONFERIR
      A[[i]] = (t(Syx)-t(Syz)%*%solve(Szz)%*%(Sxz))%*%solve(Sxx-Sxz%*%solve(Szz)%*%Sxz)
      
    }
    
  }
  
  return(A)
}

# Funcao parametro C : Essa funcao calculo o parametro C da formula da desagregacao parametrica
# Parametros: a autocovariancia calculada na funcao autocovariancia, o parametro A calculado pela funcao parametro_A
# Retorno: parametro C da desagregacao
parametro_C <- function(ACF_S,A){
  C = list()
  for(i in 1:11){
    
    if(i == 1){
      
      Sxx = ACF_S[[i]][1,4,4] # Anual com Anual
      Sxx_ = ACF_S[[i]][2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[[i]][1,4,1],ACF_S[[i]][1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[[i]][1,2,1],ACF_S[[i]][1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[[i]][1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[[i]][1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[[i]][1,1,1],ACF_S[[i]][1,3,1]), rbind(ACF_S[[i]][1,1,3],ACF_S[[i]][1,3,3]))
      Sx_z = ACF_S[[i]][2,2,4] ##### Anual com //CONFERIR
      Sxz_cor = Sxx_%*%(solve(Sxx))%*%(Sx_z)
      Syz_cor = t(Syz)+(t(Syx)%*%solve(Sxx)%*%(Sxz_cor-Sxz))
      C[[i]]= (Syz_cor-A[[i]]%*%Sxz_cor)%*%solve(Szz)
      
      
      
    }else{
      
      Sxx = ACF_S[[i]][1,4,4] # Anual com Anual
      Sxx_ = ACF_S[[i]][2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[[i]][1,4,1],ACF_S[[i]][1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[[i]][1,2,1],ACF_S[[i]][1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[[i]][1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[[i]][1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[[i]][1,1,1],ACF_S[[i]][1,3,1]), rbind(ACF_S[[i]][1,1,3],ACF_S[[i]][1,3,3]))
      Sx_z = ACF_S[[i]][2,2,4] ##### Anual com //CONFERIR
      C[[i]] = (t(Syz) - A[[i]]%*%Sxz)%*%solve(Szz)
      
      
    }
    
  }
  
  return(C)
}

# Funcao parametro_Bt : Essa funcao calculo o parametro B da formula da desagregacao parametrica
# Parametros: a autocovariancia calculada na funcao autocovariancia, o parametro A calculado pela funcao parametro_A e o parametro C calculado pela funcao parametro_C
# Retorno: parametro Bt da desagregacao
parametro_Bt <- function(ACF_S,A,C){
  Bt = list()
  for(i in 1:11){
    
    if(i == 1){
      
      Sxx = ACF_S[[i]][1,4,4] # Anual com Anual
      Sxx_ = ACF_S[[i]][2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[[i]][1,4,1],ACF_S[[i]][1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[[i]][1,2,1],ACF_S[[i]][1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[[i]][1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[[i]][1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[[i]][1,1,1],ACF_S[[i]][1,3,1]), rbind(ACF_S[[i]][1,1,3],ACF_S[[i]][1,3,3]))
      Sx_z = ACF_S[[i]][2,2,4] ##### Anual com //CONFERIR
      Sxz_cor = Sxx_%*%(solve(Sxx))%*%(Sx_z)
      Syz_cor = t(Syz)+(t(Syx)%*%solve(Sxx)%*%(Sxz_cor-Sxz))
      BBt = (Syy)-(A[[i]])%*%(Syx)-C[[i]]%*%t(Syz_cor)
      B = chol(BBt,pivot = TRUE)
      Bt[[i]] = t(B)
      
      
      
      
    }else{
      
      Sxx = ACF_S[[i]][1,4,4] # Anual com Anual
      Sxx_ = ACF_S[[i]][2,4,4] # Anual com Anual lag1
      Syx = cbind(ACF_S[[i]][1,4,1],ACF_S[[i]][1,4,3]) # Sazonal com Anual
      Syz = cbind(ACF_S[[i]][1,2,1],ACF_S[[i]][1,2,3]) # Sazonal com Sazonalidade anterior  ##Alterei aqui
      Sxz = ACF_S[[i]][1,2,4] ##### Anual com Sazonalidade anterior //CONFERIR
      Szz = ACF_S[[i]][1,2,2] # Sazonalidade anterior com Sazonalidade anterior
      Syy = cbind(rbind(ACF_S[[i]][1,1,1],ACF_S[[i]][1,3,1]), rbind(ACF_S[[i]][1,1,3],ACF_S[[i]][1,3,3]))
      Sx_z = ACF_S[[i]][2,2,4] ##### Anual com //CONFERIR
      BBt = (Syy)-(A[[i]])%*%(Syx)-C[[i]]%*%(Syz)
      B = chol(BBt,pivot=TRUE)
      Bt[[i]] = t(B)
      
      
    }
    
  }
  
  return(Bt)
}

# Funcao desagregacao_parametrica_ano: realiza a desagregacao de um ano, gerando 12 valores mensais
# Parametros: a vazao anual(ano), o zInicial, A(calculado com a funcao parametro_A),Bt(calculado com a funcao parametro_Bt)e C(calculado com a funcao parametro_C)
# Retorno: o ano desagregado
# OBS: Essa é a funcao desag_param modificada

desagregacao_parametrica_ano <- function(ano,Zinicial,A,Bt,C){
  Meses = 0
  Resto = 0
  
  # Gerando 22 valoreas aleatorios
  ERRO = rnorm (22, 0, 1)
  #normalizando o ERRO
  ERRO = (ERRO - mean (ERRO)) / sd (ERRO)
  #Criando a matriz com 11 linhas e duas colunas
  ERRO = matrix(ERRO,ncol=2)
  
  for(i in 1:11){
    CONT_ERRO = CONT_ERRO+1
    if(i==1){
      #Pegando a linha correspondente a i
      erro = ERRO[i,]
      Y = A[[i]]%*%ano+Bt[[i]]%*%(erro)+C[[i]]%*%Zinicial #Y = A%*%serie_sint[k]+Bt%*%erro+C%*%Meses
      Meses[i] = Y[1,1]
      Resto[i] = Y[2,1]
      
    }
    else{
      erro = ERRO[i,]
      Y = A[[i]]%*%Resto[i-1]+Bt[[i]]%*%erro+C[[i]]%*%Meses[i-1]
      Meses[i] = Y[1,1]
      Resto[i] = Y[2,1]
    }
    
  }
  Meses[12] = Resto[11]
  return(Meses) 
}

# anual: vazao anual sintetica
# mensal: vetor mensal desagregado da vazao anual sintetica (vetor de 12 valores)
ajuste_proporcional <- function(anual,mensal){
  
  soma_mensal = sum(mensal)
  return(sapply(mensal,function(x)(x*anual/soma_mensal)))
  
}

# Funcao de desagregacao_parametrica : Realiaza a desagregacao da serie_sintetica
# Pametro: a serie sintetica normalizada, o parametro historico calculado pela funcao desag_param_info
# Retorno: a serie desagregada
# OBS: Essa é a funcao desag_param -mult modificada com o calculo de A, Bt, C feitos fora da funcao desag_param

desagregacao_parametrica <- function(serie_sintetica_normalizada,serie_historica_normalizada,info){
  
  #nAnos: numero de anos da serie_sintetica
  nAnos = nrow(serie_sintetica_normalizada)
  
  serie_desagregada_normalizada = data.frame(matrix(NA, nrow = nAnos, ncol = 12))
  inicio = 1
  fim = nAnos
    
  #Calculo dos Parametros A, Bt e C
  ACF_S = autocovariancia(info)
  
  A = parametro_A(ACF_S)
  print('Calculo do parametro A')
  
  C = parametro_C(ACF_S,A)
  print('Calculo do paramtro C')
  
  Bt = parametro_Bt(ACF_S,A,C)
  print('Calculo do parametro Bt')
  
  #Zinicial = a vazão de dezembro do ano 1 da serie historica normalizada
  Zinicial = serie_historica_normalizada[1,12]
  
  for(i in inicio:fim){
    
    # Desagregando o ano i da serie_sintetica_normalizada
    serie_desagregada_normalizada[i,] = desagregacao_parametrica_ano(serie_sintetica_normalizada[i,1], Zinicial,A,Bt,C)
    
    # Zinicial é o valor da vazao calculada pela desagregacao do mes de dezembro. Esse valor será utilizado na desagregacao do ano  i + 1
    # Definicao de Z do artigo: contem a ultima temporada do ano anterior. Esse valor só não segue essa formula para o i = 1 
    
    Zinicial = serie_desagregada_normalizada[i,12]
      
  }
  return(serie_desagregada_normalizada)
}


