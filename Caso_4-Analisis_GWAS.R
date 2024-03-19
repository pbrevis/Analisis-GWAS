# Análisis GWAS con set de datos 
# m = 3093 marcadores SNP
# n = 281 individuos
# x = 10 cromosomas
# 3 variables fenotípicas
#(EarHT: ear height; dpoll: days to pollen or flowering time; EarDia: ear diameter)

# Instalar paquete rMVP (Memory-Efficient, Visualize-Enhanced,
# Parallel-Accelerated GWAS Tool)
# install.packages("rMVP")

# Cargando bibliotecas
library(rMVP)
library(dplyr)

rm(list = ls()) # Limpiando el ambiente global

# Guardando la ubicación del directorio de trabajo,
# una vez que fuese definida manualmente
directorio <- getwd()
head(directorio)


# PREPARACION DE DATOS PARA rMVP

MVP.Data(fileHMP="mdp_genotype.hmp.txt", # Archivo de datos genotípicos
                                         # en formato HapMap
         filePhe="mdp_traits.txt",       # Archivo de datos fenotípicos
                                         # en formato de texto
         sep.hmp="\t",
         sep.phe="\t",
         SNP.effect="Add",
         fileKin=T,
         out="mvp.hmp")

# Importando datos previamente formateados
geno_data <- attach.big.matrix("mvp.hmp.geno.desc")
pheno_data <- read.table("mvp.hmp.phe",head=TRUE)
map_info <- read.table("mvp.hmp.geno.map" , head = TRUE)
kinship <- MVP.K.VanRaden(geno_data, verbose = T)

# Chequeando clases de datos
class(geno_data)
class(pheno_data)
class(map_info)
class(kinship)

# Vista rápida a los datos importados
head(geno_data[,1:10])
dim(geno_data)
head(kinship[,1:10])
dim(kinship)
glimpse(pheno_data)
dim(pheno_data)
glimpse(map_info)
dim(map_info)


# Distribución de marcadores SNP en los 10 cromosomas
freq_table <- map_info %>% group_by(CHROM) %>% summarise(frequency=n()) %>%
  arrange(desc(frequency))
print(freq_table)

# INICIO DEL GWAS

# 1. Establecer PCA como covariante
dir.create("Pca")
setwd("./Pca/")
getwd()

GWAS_PCA_mvp <- MVP(
  phe = pheno_data, 
  geno = geno_data, 
  map = map_info,
  method =  c("MLM"),
  nPC.GLM = 5, 
  nPC.MLM = 5, 
  nPC.FarmCPU = 5)
# MLM = Mixed Linear Model
# nPC = número de "principal components" agregados como "fixed effects"

# Otros modelos estadísticos disponibles:
# GLM = Generalized Linear Model
# FarmCPU = Fixed and random model Circulating Probability Unification


# 2. Establecer parentesco (kinship) como covariante
setwd(directorio)
dir.create("Kinship")
setwd("./Kinship/")
getwd()

GWAS_Kin_mvp <- MVP(
  phe = pheno_data, 
  geno = geno_data, 
  map = map_info,
  method =  c("MLM"),
  nPC.GLM = 0,
  nPC.MLM = 0, 
  nPC.FarmCPU = 0,
  K = kinship)


# 3. Establecer PCA y Kinship como covariantes
setwd(directorio)
dir.create("./PcaKinship" )
setwd("./PcaKinship/")

GWAS_PCA_Kin_mvp <- MVP(
  phe = pheno_data, 
  geno = geno_data, 
  map = map_info,
  method =  c("MLM"),
  nPC.GLM = 5,
  nPC.MLM = 5, 
  nPC.FarmCPU = 5,
  K = kinship)