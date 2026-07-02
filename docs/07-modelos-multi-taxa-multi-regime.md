# 6. Modelos multi-taxa e multi-regime {-}

Modelos mais complexos que permitem que os parâmetros σ 2 , α, e θ (taxa, seleção, 
ótimo) variem podem ser implementados com o pacote OUwie. 
Descrição dos modelos: 
1. BM1 - Modelo Browniano com um único valor de taxa evolutiva (rate) para toda a 
árvore; 
2. BMS - Modelo Browniano com diferentes taxas evolutivas de acordo com os estados 
na árvore; 
3. OU1 - Modelo Ornstein-Uhlenbeck com um único ótimo para todas as espécies; 
4. OUM - Modelo OU com diferentes ótimos de acordo com os estados, mas um único 
alpha e sigma2 para todos os regimes seletivos; 
4.1. OUMV - Modelo OU com diferentes médias de estados (ótimos), diferentes 
sigma2, mesmo alpha; 
4.2. OUMA - Modelo OU com diferentes médias de estados (ótimos), diferentes alpha, 
mesmo sigma2; 
4.3. OUMVA - Modelo OU com múltiplos theta, alpha e sigma2. 
Em resumo, é possível criar hipóteses onde os ótimos, as taxas ou os alpha variem, em 
qualquer combinação. No entanto, modelos mais complexos vão estimar muitos 
parâmetros e podem gerar over-fitting se o número de parâmetros for parecido com o 
número de observações (espécies). Modelos mais complexos (com muitos parâmetros) 
também podem ser favorecidos artificialmente na comparação entre modelos (ver 
Cooper et al. 2016). 
O OUwie precisa de um dataframe com três colunas: nomes das espécies, o código do 
regime, e o atributo contínuo de interesse. 



``` r
require(ape)
require(geiger)
require(phytools)
require(OUwie)
# Carregar dados
dados<-read.table("dados/sigmodontinae-atributos.txt",h=T,row.names=1)
head(dados)
diet<-as.factor(dados$Diet) # A dieta é um código (Onivoro-1; Herbivoro-2;
Insetivoro-3)
names(diet)=rownames(dados)
size<-dados$Size
names(size)=rownames(dados)
```

O OUwie requer uma filogenia onde os nós internos estejam rotulados com o regime 
seletivo ancestral, o que fizemos na reconstrução ancestral de atributos discretos. 



``` r
# Carregar filogenia
tree<-read.nexus("dados/sigmodontinae-tree.tre")
plotTree(tree,fsize=0.4,ftype="i",type="fan",lwd=1)
str(tree)
# Match Species
match.species<-treedata(tree,dados)
tree<-match.species$phy
str(tree)
plotTree(tree,fsize=0.5,ftype="i",type="fan",lwd=1)
# Reconstruir evolução da dieta
diet.anc<-ace(diet,tree,type="discrete",model="ER")
diet.anc
```

O ace retorna a probabilidade de cada estado em cada nó. Precisamos de uma árvore 
única com o "melhor" estado para cada nó (alternativamente, podemos usar o 
make.simmap como fizemos anteriormente, e ajustar um modelo para cada N histórias). 
Nós podemos encontrar o estado com maior probabilidade em cada nó, e atribuir esse 
valor aos nós de uma árvore em formato phylo. 



``` r
# Primeiro, corrigimos a ordem númerica no nome dos nós para n+1
node.states<-diet.anc$lik.anc
rownames(node.states) <- seq(1:nrow(node.states)) + length(tree$tip.label)
node.states
# Extraímos o estado mais provável de cada nó
best<-apply(node.states,1,which.max)
best
# Incluimos esses valores na árvore
tree$node.label<-best
tree
```

Ajustando modelos com o OUwie 



``` r
# Brownian motion 1 taxa (equivalente ao modelo BM do fitContinuous)
BM1<-OUwie(tree,dados,model="BM1")
BM1
```

Apha=NA pois este é um modelo Browniano. Mesmo valor de sigma.sq = taxa única. 
Optima é o valor estimado na raiz. 



``` r
# Brownian motion múltiplas taxas
BMS<-OUwie(tree,dados,model="BMS")
BMS
```

Taxa de evolução de tamanho para insetívoros é mais rápida que para outros grupos. 
Mas ainda não sabemos se esse é o modelo evolutivo com maior suporte. 



``` r
# OU 1 otimo (equivalente ao modelo OU do fitContinuous)
OU1<-OUwie(tree,dados,model="OU1")
OU1
```

OU com um ótimo, único alpha e único sigma2. 



``` r
# OU múltiplos otimos
OUM<-OUwie(tree,dados,model="OUM")
OUM
```

Tamanho 'ótimo' para cada dieta é estimado, porém o alpha neste caso é muito próximo 
de zero. 



``` r
# OU múltiplos otimos, multiplos sigma2
OUMV<-OUwie(tree,dados,model="OUMV")
OUMV
# OU múltiplos otimos, multiplos alpha
OUMA<-OUwie(tree,dados,model="OUMA")
OUMA
# OU multiplos otimos, alpha e sigma2
OUMVA<-OUwie(tree,dados,model="OUMVA")
OUMVA
# Comparando modelos
aic.scores<-
setNames(c(BM1$AICc,BMS$AICc,OU1$AICc,OUM$AICc,OUMV$AICc,OUMA$A
ICc,OUMVA$AICc),
c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"))
aic.scores
aicw(aic.scores) # Delta AICc e AICc weights
```

Livre. 

