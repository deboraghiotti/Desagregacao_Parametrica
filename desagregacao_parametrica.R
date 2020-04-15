require(moments)
require(corrgram)
require(pear)
require(gridExtra)
require(knitr)
require(tcltk)
source('funcoes.R')

input_1 <- tk_choose.files()#Arquivo de serie historica
input_2 <- tk_choose.files()#Arquivo de vazao sintetica

entrada <- data.frame(read.csv(input_1, sep=";", dec=",")) #Leitura de dados historicos mensais
serie_sintetica <- data.frame(read.csv(input_2, header = F,sep=";", dec=",")) #Leitura de serie sintetica gerada para essa mesma serie

media_serie_sintetica = apply(serie_sintetica,2,mean)
serie_sintetica_normalizada <- apply(log(serie_sintetica), 2, function(x) (x-mean(x)))

SeriesDadosHist <- div_mensais(entrada)
SeriesDadosHist_normalizada <- apply(log(SeriesDadosHist),2, function(x) (x-mean(x)))

Parametro_Hist <- desag_param_info(SeriesDadosHist_normalizada)

# A Funcao desagregacao_parametrica_mult substitui a funcao desag_param_mult
# Nela, os parametros A, B, C são calculados externamente a funcao desag_param.
# A funcao parametro_A calcula o parametro A
# A funcao parametro_Bt calcula o parametro B
# A funcao parametro_C calcula o parametro C
# a funcao autocovariancia calcula a autocovariancia do Parametro_Hist que é necessario para calcular os parametros A.B e C

DesagregadoP <- desagregacao_parametrica_mult(serie_sintetica_normalizada,Parametro_Hist)
# mediaSS = apply(log(serie_sintetica),2,mean)
# Desagregado = apply(DesagregadoP,1,function(x)(x + mediaSS))
# Desagregado = t(exp(Desagregado))

#DesagregadoP_desnormalizado = desnormalizar(serie_sintetica,DesagregadoP)


