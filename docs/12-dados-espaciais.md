# (PART) Módulo II {-}

# 1. Manipulando as bases de dados {-}

Em estudos sobre macroevolução, normalmente as bases de dados envolvem 
informações sobre atributos e árvores filogenéticas. Porém, quando integramos o 
componente espacial em estudos sobre ecologia evolutiva, outras informações são 
também necessárias. 
A primeira e mais básica informação que precisamos organizar é sobre a ocorrência ou 
abundância das espécies nas assembleias, que podem ser comunidades, sítios de 
amostragem, células num grid sobreposto num mapa, etc. Normalmente, essa base de 
dados é uma matriz na qual as linhas são as assembleias, e as colunas são as espécies: 



``` r
# Carregar dados de espécies por assembleias
comm<-read.table("dados/community.txt",h=T)
dim(comm)
comm
```

Em alguns casos, especialmente quando se trabalha com muitas assembleias e/ou 
espécies, costuma-se organizar os dados de assembleias em data frames, nas quais a 
primeira coluna informa a assembleia (sítio), a segunda coluna informa a identidade da 
espécie, e a terceira informa o valor, seja de presença/ausência da espécie no sítio, ou de 
abundância de indivíduos, da seguinte forma: 
| Site | Species | Value |
|---|---|---|
| Com1 | Sp_1 | 2 |
| Com1 | Sp_2 | 10 |
| Com2 | Sp_1 | 5 |
| Com3 | Sp_3 | 1 |
| Com3 | Sp_3 | 8 | 
Neste caso, para podermos seguir nossas análises, precisamos primeiramente converter 
o data frame em matriz, o que pode ser feito da seguinte forma: 



``` r
require(reshape2)
comm.df<-read.table("dados/community_dataframe.txt",h=T)
dim(comm.df)
comm.df
comm<-acast(data=comm.df, formula=site~spp, value.var="abundance")
comm
```

O arquivo contendo a árvore filogenética pode ou não conter apenas as espécies 
ocorrentes nas assembleias (no exemplo, sim), mas em alguns casos a filogenia contém 
espécies adicionais. Em outros casos, nem todas as espécies na matriz de assembleias 
estão presentes na árvore filogenética. 



``` r
require(ape)
require(phytools)
tree<-read.tree("dados/bird_orders.txt")
plotTree(tree,fsize=0.4,ftype="i",type="phylogram",lwd=1)
```

Note que, se a análise envolver atributos das espécies, este(s) deverá(ão) ser definido(s) 
num arquivo no qual as espécies serão as linhas e o(s) atributo(s) a(s) coluna(s): 



``` r
atributos<-read.table("dados/traits.txt",h=T)
atributos
```

Outras duas matrizes de dados importantes são as que descrevem as assembleias a partir 
de variáveis ambientais e de coordenadas geográficas, que poderão ser necessárias 
dependendo da análise que realizarmos. Estas matrizes podem ser organizadas como um 
único data frame, ou separadamente em dois arquivos. Por simplicidade, aqui utilizamos 
apenas um único arquivo: 



``` r
env<- read.table("dados/environment.txt",h=T)
ambiente<- as.data.frame(env[,3])
rownames(ambiente)<-rownames(env)
colnames(ambiente)<- “E”
ambiente
coordenadas<- env[,1:2]
coordenadas
```

Como visto, análises eco-evolutivas que incorporam a dimensão espacial demandam 
muitas matrizes de dados distintas. Isso pode gerar problemas, se as linhas e colunas nas 
matrizes não forem correspondentes entre os diversos objetos (comunidades, atributos, 
filogenia, ambiente, coordenadas espaciais. Uma forma de garantir que os dados estejam 
todos organizados é utilizando a função organize.syncsa: 



``` r
require(SYNCSA)
org<-organize.syncsa(comm=comm, traits = atributos, phylodist = cophenetic(tree),
envir = env, check.comm = TRUE)
```

Caso as comunidades apresentem espécies não presentes na filogenia, ou a filogenia 
apresente espécies não presentes nas comunidades, a função removerá essas espécies e 
informará ao usuário. 
Como podemos notar, a função não organiza os dados de comunidades de acordo com a 
árvore filogenética, mas a partir da matriz de distâncias filogenéticas. Em boa parte das 
análises isso é suficiente. Mas em outros casos os nomes das espécies na matriz de 
comunidades devem corresponder à estrutura da árvore. A função match.phylo.comm 
pode resolver o problema neste caso: 



``` r
match.species<-picante::match.phylo.comm(tree,comm)
tree<-match.species$phy
comm<- match.species$comm
```

De toda forma, sempre vale conferir se os nomes estão todos alinhados antes de rodar as 
análises. 
Com base apenas nos objetos comm, atributos e ambiente, e sem o uso de funções 
específicas, calcule apenas com operações simples, no R (ou no Excel) o valor médio 
dos atributos T1 e T2 nas comunidades, e correlacione cada valor médio com a variável 
ambiental. 

