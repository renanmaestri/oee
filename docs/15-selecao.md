# 4. Seleção e composição fenotípica {-}

Análises de diversidade e/ou composição funcional envolvem calcular medidas que 
consideram conjuntos de espécies descritas tanto por um ou mais atributos quanto sua 
distribuição num conjunto de assembleias. Neste sentido, o parentesco evolutivo entre 
as espécies influenciará tanto a evolução fenotípica da linhagem estudada (evolução de 
nicho), quanto a distribuição das linhagens entre as assembleias (diversificação e 
dispersão). Portanto, análises de determinantes ambientais da diversidade e/ou 
composição funcional podem apresentar inflação do erro tipo I. No script abaixo 
mostramos algumas maneiras de contornar o problema tanto para medidas de 
diversidade (Diniz-Filho et al. 2011) quanto para composição fenotípica (Duarte et al. 
2018). 



``` r
# require(phytools) ou require(ape)
require(picante)
# Carregar dados
comm<-read.table("dados/community.txt",h=T)
dim(comm)
comm
# Carregar variável ambiental
env<-read.table("dados/environment.txt",h=T)
dim(env)
env
# Carregar filogenia (duas alternativas)
tree<-phytools::read.newick("bird_orders.new")
#ou
tree<-ape::read.tree("dados/bird_orders.txt")
phytools::plotTree(tree,fsize=0.4,ftype="i",type="fan",lwd=1)
# Carregar atributos
traits<-read.table("dados/traits.txt",h=T)
dim(traits)
traits
```

Organizar os conjuntos de dados: 



``` r
require(SYNCSA)
org<-organize.syncsa(comm,traits,phylodist=cophenetic(tree),envir=env)
comm<-org$community
phydist<-org$phylodist
traits<-org$traits
env<-org$environmental
```

Quantificar sinal filogenético nos traits: 



``` r
require(phytools)
k.T1<-phylosig(tree,traits$T1,method="K",test=T)
k.T2<-phylosig(tree,traits$T2,method="K",test=T)
```

Calcular autovetores filogenéticos: 



``` r
require(PVR)
pvr_class<- PVRdecomp(tree)
pvr_reg.T1<-PVR(pvr_class,trait=traits$T1,envVar=NULL,method="moran")
sel.vec.brown<-pvr_reg.T1@Selection$Vectors
r2.pvr.T1<-pvr_reg.T1@PVR$R2
res.pvr.T1<-pvr_reg.T1@PVR$Residuals
```

Calcular resíduos da diversidade fenotípica (controlando efeito da filogenia) e analisar 
associação com o ambiente: 



``` r
require(picante)
require(cluster)
# Calcular distâncias fenotípicas:
dis.T1<-as.matrix(daisy(as.data.frame(res.pvr.T1),
metric="euclidean",stand=FALSE))
rownames(dis.T1)<-colnames(comm)
colnames(dis.T1)<-colnames(comm)
# Calcular mfd – mean functional distance usando resíduos filogenéticos de T1:
ses.mfd.T1=ses.mpd(comm,
dis.T1,null.model="independentswap",abundance.weighted=TRUE,runs=999,iteratio
ns=1000)
ses.mfd.T1=as.numeric(ses.mfd.T1[,6])
# Estimar a associação entre ses.mfd.T1 e env:
lm.ses.mfd.T1<-lm(ses.mfd.T1~E,data=env)
summary.lm(lm.ses.mfd.T1)
```

Observação: Nesta análise, apenas o efeito do sinal filogenético do atributo é 
controlado. A distribuição das linhagens no espaço não é afetada. 
Analisar associação entre valores médios do atributo “T1” e variável ambiental “E” 
controlando ou não o efeito da filogenia: 



``` r
source(“cwm.sig.R”)
# Ignorando o efeito da filogenia:
cwmsig.T1.phylo<-cwm.sig(comm=comm,traits=traits,envir=env, formula="T1~E",
PGLS=FALSE,tree=tree,runs = 999)
summary.lm(cwmsig.T1.phylo$model)
# Controlando o efeito da filogenia:
cwmsig.T1.nophylo<-cwm.sig(comm=comm,traits=traits,envir=env,
formula="T1~E", PGLS=TRUE,tree=tree,runs = 999)
summary.lm(cwmsig.T1.nophylo$model)
```


## Exercício

1) Calcule a relação entre diversidade fenotípica e ambiente sem controlar o efeito 
da filogenia e compare os resultados. Reanalise usando T2. 
2) De que forma a função cwm.sig permite controlar o efeito da filogenia sobre a 
relação entre valores médios de atributos nas comunidades e a variável 
ambiental? Qual o argumento chave para controlar o efeito da filogenia? 

## Referências

Diniz-Filho, J. A. F., M. V. Cianciaruso, T. F. Rangel, and L. M. Bini. 2011. 
Eigenvector estimation of phylogenetic and functional diversity. Functional 
Ecology 25:735-744. 
Duarte, L. D. S., V. J. Debastiani, M. B. Carlucci, and J. A. F. Diniz-Filho. 2018. 
Analyzing community-weighted trait means across environmental gradients: 
should phylogeny stay or should it go? Ecology 99:385–398. 

