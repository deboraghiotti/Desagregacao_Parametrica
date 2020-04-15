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

# Funcao desag_param_info: calcula dos parametros historicos necessarios para a desagregacao parametrica
# Parametro: a serie historica normalizada 
# Retorno: parametros historicos para a desagregacao
desag_param_info<-function(SeriesDadosHist_normalizada){
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
# Parametro: os parametros historicos calculados pela funcao desag_param_info
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

# Funcao desagregacao_parametrica: realiza a desagregacao de um ano, gerando 12 valores mensais
# Parametros: a vazao anual(ano), o zInicial, A(calculado com a funcao parametro_A),Bt(calculado com a funcao parametro_Bt)e C(calculado com a funcao parametro_C)
# Retorno: o ano desagregado
# OBS: Essa é a funcao desag_param modificada

desagregacao_parametrica <- function(ano,Zinicial,A,Bt,C){
  Meses = 0
  Resto = 0
  # CONT_ERRO é inicialmente zero na funcao desag_param_mult
  # o for calcula os valores desagregados de jan a nov
  for(i in 1:11){
    CONT_ERRO = CONT_ERRO+1
    if(i==1){
      erro = ERRO[[CONT_ERRO]]
      # print("i"); print(i)
      # print("erro"); print(erro)
      Y = A[[i]]%*%ano+Bt[[i]]%*%(erro)+C[[i]]%*%Zinicial #Y = A%*%serie_sint[k]+Bt%*%erro+C%*%Meses
      Meses[i] = Y[1,1]
      Resto[i] = Y[2,1]
      
    }
    else{
      
      erro = ERRO[[CONT_ERRO]]
      # print("i");print(i)
      # print("erro"); print(erro)
      Y = A[[i]]%*%Resto[i-1]+Bt[[i]]%*%erro+C[[i]]%*%Meses[i-1]
      Meses[i] = Y[1,1]
      Resto[i] = Y[2,1]
    }
    
  }
  Meses[12] = Resto[11]
  #print(Meses)
  #print(sum(Meses))
  return(Meses) 
}

# Primeiro: retirar o input_3 OK
# Segundo: permitir que o ano da serie sintetica seja diferente de 10000

# Funcao de desagregacao_parametrica_mult : Realiaza a desagregacao da serie_sintetica
# Pametro: a serie sintetica normalizada, o parametro historico calculado pela funcao desag_param_info
# Retorno: a serie desagregada
# OBS: Essa é a funcao desag_param -mult modificada com o calculo de A, Bt, C feitos fora da funcao desag_param

desagregacao_parametrica_mult <- function(serie_sintetica_normalizada,info){
  aux = (serie_sintetica_normalizada)
  #nAnos: numero de anos da serie_sintetica
  nAnos = nrow(serie_sintetica_normalizada)
  #OBS: nem sempre serao 10000, é preciso alterar isso
  Desagregado = list()
  CONT_ERRO <<- 0
  # ERRO_AUX : Funcao rnorm ira gerar 20000 valores aleatorios de media zero e variancia igual a 1 
  ERRO_AUX = rnorm (nAnos*2, 0, 1)
  # Normalizando o ERRO_AUX
  ERRO_AUX = (ERRO_AUX - mean (ERRO_AUX)) / sd (ERRO_AUX)
  # aux1: matriz com duas linhas e 10000 colunas, contendo os valores de ERRO_AUX.
  aux1 = matrix (ERRO_AUX, nrow = 2)
    
  #Criando uma lista de 10000 elementos, onde cada elemento possu dois valores aleatorios
  ERRO <<- list()
  for(i in 1:nAnos) {
    ERRO[[i]] <<- aux1[, i]
  }
    
  ################MUDAR PARA OS ANOS FINAIS AO INVES DOS ANOS INICIAIS###########
  df_serie_desagregada_param = data.frame(matrix(NA, nrow = nAnos, ncol = 12))
  cont = 1
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
  Zinicial = SeriesDadosHist_normalizada[1,12]
    
  # for de 1 a 10000, isto é, cada elemento da serie_sintetica de 10000 anos
  # OBS: cont possui o mesmo valor de i no loop
  for(i in inicio:fim){
      
    # A função desagregacao_parametrica retorna uma lista(ou vetor) de 12 elementos. Essa lista será o valor da serie desagregada do ano i.
    # info = os parametros calculados na funcao desag_param_info
    # aux[i,1] = é o valor da vazao da serie sintetica no ano i
    
    print(paste('Desagregacao do valor',i))
    df_serie_desagregada_param[cont,] = desagregacao_parametrica(aux[i,1], Zinicial,A,Bt,C)
    
    # Realizando o ajuste mensal para preservar a aditividade
    df_serie_desagregada_param[cont,] = ajuste_proporcional(aux[i,1],df_serie_desagregada_param[cont,])
    
    #Zinicial é o valor da vazao calculada pela desagregacao do mes de dezembro. Esse valor será utilizado na desagregacao do ano  i + 1
    #Definicao de Z do artigo: contem a ultima temporada do ano anterior. Esse valor só não segue essa formula para o i = 1 
    
    Zinicial = df_serie_desagregada_param[cont,12]
      
    # O valor de cont é incrementado
    cont = cont+1
  }
    
  Desagregado =  df_serie_desagregada_param
  return(Desagregado)
}

desag_param_mult<-function(info){
  aux = (serie_sintetica_normalizada)
  #OBS: input_3 sempre igual a 1
  #OBS: nem sempre serao 10000, é preciso alterar isso
  Desagregado = list()
  if(input_3<=0){
    print("Numero invalido de series a se desagregar")
    break
  }
  for(j in 1:input_3){
    CONT_ERRO <<- 0
    # ERRO_AUX : Funcao rnorm ira gerar 20000 valores aleatorios de media zero e variancia igual a 1 
    ERRO_AUX = rnorm (10000*2, 0, 1)
    # Normalizando o ERRO_AUX
    ERRO_AUX = (ERRO_AUX - mean (ERRO_AUX)) / sd (ERRO_AUX)
    # aux1: matriz com duas linhas e 10000 colunas, contendo os valores de ERRO_AUX.
    aux1 = matrix (ERRO_AUX, nrow = 2)
    
    #Criando uma lista de 10000 elementos, onde cada elemento possu dois valores aleatorios
    ERRO <<- list()
    for(i in 1:10000) {
      ERRO[[i]] <<- aux1[, i]
    }
    
    ################MUDAR PARA OS ANOS FINAIS AO INVES DOS ANOS INICIAIS###########
    df_serie_desagregada_param = data.frame(matrix(NA, nrow = 10000, ncol = 12))
    cont = 1
    # Como o input_ 3 é 1, jogo o inicio = 1 e o fim = 10000
    inicio =((10000*(j-1))+1)
    fim = 10000*j
    
    #Zinicial = a vazão de dezembro do ano 1 da serie historica normalizada
    Zinicial = SeriesDadosHist_normalizada[1,12]
    
    # for de 1 a 10000, isto é, cada elemento da serie_sintetica de 10000 anos
    # OBS: cont possui o mesmo valor de i no loop
    for(i in inicio:fim){
      
      # A função desag_param retorna uma lista(ou vetor) de 12 elementos. Essa lista será o valor da serie desagregada do ano i.
      # info = os parametros calculados na funcao desag_param_info
      # aux[i,1] = é o valor da vazao da serie sintetica no ano i
      df_serie_desagregada_param[cont,] = desag_param(info, aux[i,1], Zinicial)
      
      #Zinicial é o valor da vazao calculada pela desagregacao do mes de dezembro. Esse valor será utilizado na desagregacao do ano  i + 1
      #Definicao de Z do artigo: contem a ultima temporada do ano anterior. Esse valor só não segue essa formula para o i = 1 
      Zinicial = df_serie_desagregada_param[cont,12]
      
      # O valor de cont é incrementado
      cont = cont+1
    }
    
    # Como j = input_3 = 1, logo essa linha não sera necessaria
    Desagregado[[j]] =  df_serie_desagregada_param
  }
  return(Desagregado)
}


#########Função que desegrega]
# A função desag_param retorna uma lista(ou vetor) de 12 elementos. Essa lista será o valor da serie desagregada do ano i.
# info = os parametros calculados na funcao desag_param_info
# aux[i,1] = é o valor da vazao da serie sintetica no ano i
# zinicial = o valor da vazao calculada pela desagregacao do mes de dezembro. Esse valor será utilizado na desagregacao do ano  i + 1
desag_param<-function(info, ano, Zinicial){
  Meses = 0
  Resto = 0
  # CONT_ERRO é inicialmente zero na funcao desag_param_mult
  # o for calcula os valores desagregados de jan a nov
  for(i in 1:11){
    CONT_ERRO = CONT_ERRO+1
    if(i==1){
      # Esse valor é o mesmo para toda vazao. Não depende do valor da vazao, so do Info
      ACF_S = acf(info[[i]],lag.max = 1,type = "covariance", plot = FALSE)
      print(ACF_S)
      ACF_S = ACF_S$acf
      #Esse valor nao depende da vazao
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

desnormalizar <- function(serie,serie_normalizada){
  media = apply(log(serie),2,mean)
  for(i in 1:nrow(serie)){
    serie_normalizada[i,] = serie_normalizada[i,] + media
  }
  
  serie_normalizada = exp(serie_normalizada)
  
  return(serie_normalizada)
}

# anual: vazao anual sintetica
# mensal: vetor mensal desagregado da vazao anual sintetica (vetor de 12 valores)
ajuste_proporcional <- function(anual,mensal){
  
  soma_mensal = sum(mensal)
  return(sapply(mensal,function(x)(x*anual/soma_mensal)))
  
}
