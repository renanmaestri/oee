# 2. Diversidade filogenética {-}

Existem diversas medidas de diversidade filogenética (PD), e várias delas apresentam 
correspondência com medidas de diversidade funcional. Além disso, algumas são 
correlacionadas com a riqueza de espécies observada nas comunidades. 
Para exercitarmos o cálculo de PD, tomemos como exemplo 50 comunidades de aves 
descritas pela abundância de indivíduos de diferentes ordens. 



``` r
# require(phytools) ou require(ape)
require(picante)
# Carregar dados
comm<-read.table("dados/community.txt",h=T)
dim(comm)
comm
# Carregar filogenia (duas alternativas)
tree<-phytools::read.newick("bird_orders.new")
#ou
tree<-ape::read.tree("dados/bird_orders.txt")
plotTree(tree,fsize=0.4,ftype="i",type="fan",lwd=1)
```

Para verificar a correspondência entre espécies nos dados e na filogenia usar a função 
match.phylo.comm. 



``` r
match.species<-picante::match.phylo.comm(tree,comm)
tree<-match.species$phy
comm<- match.species$comm
```




``` r
pd_comm<-picante::pd(comm, tree)
sespd_comm<-ses.pd(comm, tree, null.model="taxa.labels", runs=999,
iteractions=100)
```

Para calcular o tamanho de efeito padronizado (“standardized effect size” – SES), a 
função ses.pd disponibiliza várias opções de modelos nulos. Antes de rodar a análise, 
avaliar qual modelo mais se adequa ao objetivo da análise. No manual do pacote picante 
(https://cran.r-project.org/web/packages/picante/picante.pdf), os diferentes modelos 
nulos são apresentados em detalhe. 

## Mean Phylogenetic Distance (MPD)




``` r
mpd_comm<-picante::mpd(comm,cophenetic(tree), abundance.weighted = FALSE)
ses_mpd_comm<-ses.mpd(comm, cophenetic(tree), null.model = "taxa.labels",
abundance.weighted = FALSE, runs = 999, iterations = 1000)
```

Verificar a opção de modelo nulo mais adequada (https://cran.r- 
project.org/web/packages/picante/picante.pdf). 

## Mean Nearest Taxon Distance (MNTD)




``` r
mntd_comm<-picante::mntd(comm,cophenetic(tree), abundance.weighted = FALSE)
ses_mntd_comm<-ses.mntd(comm, cophenetic(tree), null.model = "taxa.labels",
abundance.weighted = FALSE, runs = 999, iterations = 1000)
```

Verificar a opção de modelo nulo mais adequada (https://cran.r- 
project.org/web/packages/picante/picante.pdf 

## PD de Rao (Rao’s D)




``` r
raoD_comm<-picante::raoD(comm,tree)$Dkk
```


## Variabilidade filogenética (PSV)




``` r
psv_comm<-picante::psv(comm,tree,compute.var=TRUE,scale.vcv=TRUE)
psr_comm<-picante::psr(comm,tree,compute.var=TRUE,scale.vcv=TRUE)
```

A função psr calcula uma medida de riqueza filogenética ao multiplicar a riqueza 
observada pelo valor de PSV. 

## Exercício

1) Calcule as diferentes medidas de diversidade filogenética, e correlacione os 
valores obtidos para cada métrica de diversidade. O que torna alguns mais 
correlacionados entre si do que outros? Quais são mais correlacionados com o 
número de espécies observado em cada comunidade? 
2) Calcule SES.MPD utilizando diferentes modelos nulos e compare os resultados. 
Por que alguns modelos dão resultados tão diferentes dos outros? 

