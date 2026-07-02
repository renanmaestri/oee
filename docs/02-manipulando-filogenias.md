# (PART) Módulo I {-}

# 1. Manipulando filogenias no R - Sinal Filogenético {-}


## O ambiente R

O R pode ser instalado a partir do website: https://www.r-project.org/ 
Nós vamos usar também o RStudio: https://www.rstudio.com/ 
R é o programa estatístico predominante em ecologia e evolução, em larga medida 
porque é flexível e adaptável para várias necessidades, e existe uma comunidade grande 
de indivíduos dispostos a fornecer ajuda. No contexto de análises comparativas, o R tem 
sido a área primária de desenvolvimento de novos métodos. 

## Set working directory

É conveniente criar uma pasta única com todos os arquivos das aulas práticas. 
Nós vamos informar ao R para carregar os dados a partir desta pasta. 
Isso pode ser feito no RStudio com o menu: Session>Set Working Directory>Choose 
Directory. 
Ou com o código setwd(). 

## Scripts

Um script é um bloco de texto que pode ser editado em qualquer editor de texto comum 
(ex: bloco de notas), ou com o editor de texto do próprio R ou RStudio. Você vai 
escrever (ou copiar) os códigos dentro do editor de texto, e depois colar no console do 
R; alternativamente é possível selecionar a parte do código a ser enviada ao console do 
R e pressionar crtl+r ou ctrl+Enter. 
Ao salvar o script você mantém um registro do código utilizado. 
Se você quiser adicionar comentários ao código basta adicionar o sinal # na frente do 
comentário. 

## Funcionamento básico do R

Objetos e funções são os principais componentes do R. Objetos são variáveis, dados e 
resultados. Funções são tipos especiais de objetos que usam um ou mais argumentos 
para fazer alguma coisa dentro do R. Compilados de funções são organizados dentro de 
livrarias (library) criadas por múltiplos autores contribuintes. 
O R funciona através de comandos, que podem ser funções ou operações (+,-,*,etc). 
Um comando retorna um objeto na tela, ou armazena o objeto na memória quando um 
nome é atribuído com o operador <-. 



``` r
2+5
x<-2+5
x
# Os objetos são case-sensitive
x<-1
X<-10
x
X
```


## Tipos de objeto no R

Alguns tipos de objetos são mais comuns no R, e vamos encontrá-los frequentemente: 
vector, factor, matrix, data frame, e list. 



``` r
# Vetores: séries de elementos do mesmo tipo
x<-c(1,2,3,4)
x<-1:5
x
class(x)
length(x)
# A função which é útil para selecionar alguns dados de um vetor
which(x>3)
which(x>=3)
# Fatores possuem níveis
f<-c("Sp1","Sp1","Sp2","Sp3","Sp3")
f
f<-factor(f)
f
table(f)
class(f)
# Os vetores podem conter nomes, o que é bastante útil em análises filogenéticas
names(x)=f
x
# Matrizes
m<-matrix(1:9,3,3) # dados, n-linhas, n-colunas
m
class(m)
m[,2] # a segunda coluna
m[3,] # a terceira linha
# Data frames suportam vários tipos de dados distintos
d<-data.frame(f,x)
d
class(d)
# Lista é uma estrutura geral que pode carregar vários tipos de dados, incluindo árvores
# filogenéticas
l<-list(f,x,m,d)
l
l[[4]] # quarto elemento da lista
# Função ls() lista todos os objetos criados na seção
ls()
?ls # As funções possuem páginas de ajuda descrevendo a função.
```

Algumas operações são úteis para reproduzir códigos. O for loop repete uma 
determinada operação por um número fixo de vezes. Por exemplo, para repetir a mesma 
operação em todas as filogenias de um conjunto de dados. 



``` r
# For loop para calcular a média de atributos, ex: colunas de uma matriz
X<-matrix(1:16,4,4)
X
# vetor vazio para salvar os resultados
medias<-vector(mode="numeric",length(ncol(X)))
medias
# loop
for(i in 1:ncol(X)){
medias[i]<-mean(X[,i])
}
# resultado
medias
```


## Pacotes filogenéticos no R

Existem muitos pacotes lidando com métodos filogenéticos no R. Nós vamos trabalhar 
com alguns pacotes mais importantes neste momento: ape, phytools e geiger. 



``` r
# Instalando pacotes:
install.packages("ape")
install.packages("phytools")
install.packages("geiger")
```

A maioria dos pacotes é constantemente atualizado. É recomendável checar 
frequentemente por atualizações. Através do comando ? seguido do nome do pacote (ex: 
?phytools) é possível acessar o índice com a descrição de todas as funções daquele 
pacote. 

## O objeto “phylo”

Árvores filogenéticas são lidas no R como um objeto list de classe phylo . 



``` r
# Carregar ape
require(ape)
# Construindo uma árvore newick via texto
texto<-"(((((((vaca,
cavalo),baleia),(morcego,(camundongo,humano))),(pardal,tartaruga)),celacanto),nemo),
tubarão);"
# Função read.tree lê o arquivo contendo a árvore
arvore<-read.tree(text=texto)
# Representando a árvore graficamente
plot(arvore,edge.width=2)
# Estrutura da árvore
arvore
str(arvore)
# Decompondo a árvore
class(arvore)
plot(arvore)
tiplabels()
nodelabels()
# ramos
arvore$edge # matriz contendo os índices do nó inicial e final para cada ramo
# O número de linhas em tree$edge é igual o número de ramos; cada ramo começa e
# termina com um par de índices (nó-nó ou nó-terminal)
# terminais
arvore$tip.label
# nós
arvore$Nnode
# Um objeto "phylo" também pode conter o comprimento dos ramos: edge.length
```


## Lendo e manipulando árvores filogenéticas

É aconselhável criar um único diretório contendo a(s) árvore(s) e os dados fenotípicos, e 
informar ao R de qual diretório os arquivos devem ser lidos. 
Árvores filogenéticas podem ser importadas para dentro do ambiente R através de 
funções como read.tree e read.nexus. 



``` r
# Ler árvore contendo 31 espécies do gênero de roedores Neotropicais Akodon
akodon.tree<-read.tree(file="akodon.tree")
akodon.tree
# Existem muitas formas de plotar a árvore
plot(akodon.tree)
# cladograma orientado pra esquerda
plot(akodon.tree, type = "c", use.edge.length = FALSE,direction="l")
# nomes alinhados à esquerda
plot(akodon.tree, type = "c", use.edge.length = FALSE,direction="l",adj=0)
# visualizar partes da árvore com zoom
zoom(akodon.tree, 13:17, subtree = FALSE, cex = 0.8)
par(mfrow=c(1,1))
# mudando cores
plot(akodon.tree, type = "fan", edge.color = "violetred", tip.color = "turquoise1",
cex = 0.5)
# com o pacote phytools
require(phytools)
plotTree(akodon.tree,ftype="i",fsize=0.8,lwd=1)
# filograma arredondado
roundPhylogram(akodon.tree)
# cladograma desenraizado (unrooted)
plot(unroot(akodon.tree),type="unrooted",cex=0.8,use.edge.length=FALSE,lab4ut="axi
al")
# árvore circular (fan tree)
plotTree(akodon.tree,type="fan",fsize=0.8,lwd=1,ftype="i")
# Legenda dos comprimento de ramos
plotTree(akodon.tree,ftype="i",fsize=0.8,lwd=1,mar=c(3.0,0.1,0.1,0.1))
axisPhylo()
```


## Extraindo clados ou espécies da árvore

Para extrair espécies da árvore podemos identificar o conjunto a ser removido. 
Por exemplo, as espécies de Akodon do sul do Brasil: A. paranaensis, A. azarae, A. 
reigi, A. montensis, C. angustidens. 



``` r
akodon.tree$tip.label
sul.especies<-
c("Akodon_paranaensis","Akodon_azarae","Akodon_reigi","Akodon_montensis","Cast
oria_angustidens")
# remover essas espécies da árvore
akodon.drop<-drop.tip(akodon.tree,sul.especies)
akodon.drop
plotTree(akodon.tree,ftype="i",fsize=0.8,lwd=1)
plotTree(akodon.drop,ftype="i",fsize=0.8,lwd=1)
# ou manter apenas essas espécies na árvore
akodon.drop1<-drop.tip(akodon.tree,setdiff(akodon.tree$tip.label,sul.especies))
akodon.drop1
plotTree(akodon.drop1,ftype="i",fsize=0.8,lwd=1)
# extrair um clado da árvore
plotTree(akodon.tree,ftype="i",fsize=0.8,lwd=1)
nodelabels()
clado1<-extract.clade(akodon.tree,node=35)
plotTree(clado1,ftype="i",fsize=0.8,lwd=1)
# remover o clado
akodon.noclado1<-drop.tip(akodon.tree,clado1$tip.label)
akodon.noclado1
plotTree(akodon.noclado1,ftype="i",fsize=0.8,lwd=1)
# função interativa do phytools para remover clados
akodon.col<-collapseTree(akodon.tree)
plotTree(akodon.col,ftype="i",fsize=0.8,lwd=1)
```


## Resolvendo politomias

Várias funções implementadas no R para métodos comparativos vão lidar com 
filogenias contendo politomias. No entanto, nem todas as funções vão funcionar; 
algumas exigem árvores sem politomias. 
É possível resolver aleatoriamente as politomias com a função multi2di do ape. 



``` r
ex<-read.tree(text="((A,B,C),(D,E));")
plot(ex)
is.binary(ex)
# resolver ao acaso
ex1<-multi2di(ex)
plot(ex1)
is.binary(ex1)
```

Existem outros métodos que permitem gerar múltiplas árvores com as politomias 
resolvidas ao acaso. Por exemplo: http://wsmartins.net/sunplin/ - Rangel, T.F. et al. 
(2015) Phylogenetic uncertainty revisited: implications for ecological analyses. 
Evolution, 69, 1301–1312. 
Lembre-se que ao resolver apenas 1 vez ao acaso não existe avaliação de incerteza 
filogenética. 

## Múltiplas árvores

Múltiplas filogenias podem ser úteis se queremos incorporar incerteza filogenética, por 
exemplo, replicando as mesmas análises em múltiplas árvores. 
Árvores múltiplas pertencem ao objeto de classe multiPhylo, que é uma lista de 
filogenias da classe phylo. Muitas funções podem prontamente ser aplicadas tanto para 
phylo como para multiPhylo. 



``` r
akodon.trees<-c(akodon.tree,akodon.drop,akodon.drop1,clado1,akodon.noclado1)
akodon.trees
class(akodon.trees)
print(akodon.trees,details=TRUE)
str(akodon.trees[1])
# Extraindo clados de múltiplas árvores
trees<-pbtree(n=50,nsim=100)
trees
plot(trees[[1]])
tree_cortado<- vector(mode = "list", length(trees))
for(i in 1:length(trees)){
tree_cortado[[i]]<-keep.tip(trees[[i]], tip=trees[[1]]$tip.label[1:12])
}
tree_cortado
plot(tree_cortado[[1]])
```


## Sinal Filogenético

Sinal filogenético é meramente um padrão que mostra se espécies mais próximas têm 
valores de atributo mais similares do que espécies mais distantes. É comum encontrar 
menções de que as regressões foram corrigidas pela filogenia porque um sinal 
filogenético foi encontrado nas variáveis brutas. No entanto, o uso de regressão 
filogenética deve estar associado aos resíduos do modelo de interesse apresentarem 
sinal filogenético (ver Revell 2010 Methods in Ecology and Evolution). 
Para calcular o sinal filogenético, nós vamos usar alguns dados de borboletas do artigo 
de Swanson et al 2016 Proc R Soc B. 



``` r
require(ape)
# Carregar dados
dados<-read.table("dados/butterfly-data.txt",h=T,row.names = 1)
dados
attach(dados)
wing_length<-setNames(wing_length,rownames(dados))
temp<-setNames(temp,rownames(dados))
eye_width<-setNames(eye_width,rownames(dados))
# Carregar árvore
tree<-read.tree("dados/butterfly-tree.txt")
plot(tree)
axisPhylo()
```

Checar correspondência entre espécies na filogenia e nos dados. Todas as espécies 
precisam estar tanto na árvore como nos dados, e também tem que existir 
correspondência exata nos nomes. 



``` r
require(geiger)
obj<-name.check(tree,dados)
obj
# 59 espécies na árvore não contêm dados
```

Nós podemos remover estas espécies da árvore para seguir com as análises. 



``` r
# remover espécies da árvore
tree.pruned<-drop.tip(tree,obj$tree_not_data)
name.check(tree.pruned,dados)
plot(tree.pruned)
```

Vários métodos existem para calcular sinal filogenético. Os dois métodos mais comuns 
são a estatística K (Blomberg et al 2003) e o lambda de Pagel (1999). 
K 
K=1 indica que espécies próximas são tão parecidas quanto prevê o modelo Browniano 
K<1 indica sinal filogenético menor do que esperado pelo modelo Browniano 
K>1 indica sinal filogenético maior do que esperado pelo modelo Browniano 
λ 
λ=0 indica ausência de sinal filogenético (sinal esperado em uma filogenia estrela) 
λ=1 indica sinal filogenético correspondente ao modelo Browniano 



``` r
require(phytools)
# K
phylosig(tree.pruned,wing_length,method="K")
# Pagel's λ
phylosig(tree.pruned,wing_length,method = "lambda")
# Visualizando o atributo na filogenia
plotTree.barplot(tree.pruned,wing_length,
args.plotTree=list(fsize=0.5,ftype="i"))
par(mfrow=c(1,1))
```


## Exercício 1 - Árvore filogenética no R

Carregue uma árvore filogenética para o seu grupo de estudo e observe a estrutura dela. 
Represente a árvore graficamente de diferentes formas. Se você tiver variáveis 
fenotípicas para o seu grupo de estudo, tente carregar no R e testar o sinal filogenético. 

