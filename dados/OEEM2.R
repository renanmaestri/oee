#ATIVIDADE 1

# Carregar dados de espécies por assembleias
comm<-read.table("community.txt",h=T)
dim(comm)
comm

# Carregar árvore
require(ape)
require(phytools)
tree<-read.tree("bird_orders.txt")
plotTree(tree,fsize=0.4,ftype="i",type="phylogram",lwd=1)

#Carregar atributos
atributos<-read.table("traits.txt",h=T)
atributos

# Carregar ambiente
env<- read.table("environment.txt",h=T)
ambiente<-as.data.frame(env[,3])
rownames(ambiente)<-rownames(env)
colnames(ambiente)<-"E"
ambiente

# Carregar coordenadas geográficas
coordenadas<- env[,1:2]
coordenadas

# Organiza dados
require(SYNCSA)
org<-organize.syncsa(comm=comm, traits = atributos, phylodist = cophenetic(tree),
                     envir = env, check.comm = TRUE)

match.species<-picante::match.phylo.comm(tree,comm)
tree<-match.species$phy
comm<- match.species$comm

# ATIVIDADE 2  - DIVERSIDADE FILOGENÉTICA

# Faith's PD
pd_comm<-picante::pd(comm, tree)
sespd_comm<-ses.pd(comm, tree, null.model="taxa.labels", runs=999,
                   iteractions=100)

# MPD
mpd_comm<-picante::mpd(comm,cophenetic(tree), abundance.weighted = F)

ses_mpd_comm<-ses.mpd(comm, cophenetic(tree), null.model = "taxa.labels",
                      abundance.weighted = FALSE, runs = 999, iterations = 1000)

mntd_comm<-picante::mntd(comm,cophenetic(tree), abundance.weighted = FALSE)
ses_mntd_comm<-ses.mntd(comm, cophenetic(tree), null.model = "taxa.labels",
                        abundance.weighted = FALSE, runs = 999, iterations = 1000)


#PD de Rao (Rao’s D)
raoD_comm<-picante::raoD(comm,tree)$Dkk

#PSV
psv_comm<-picante::psv(comm,tree,compute.var=TRUE,scale.vcv=TRUE)
psr_comm<-picante::psr(comm,tree,compute.var=TRUE,scale.vcv=TRUE)


plot(mpd_comm,pd_comm[,2])
plot(pd_comm[,1],pd_comm[,2])

# ATIVIDADE 3 - FILOBETADIVERSIDADE

require(picante)

phylosor_comm<-phylosor(comm,tree)
phylosor_comm

# Unifrac para dados de presença/ausência:

require(picante)
unifrac_comm<-unifrac(comm, tree)
unifrac_comm

# Ou

require(GUniFrac)
unifrac_presab_comm<-GUniFrac(comm, tree)$unifracs
unifrac_presab_comm[,,"d_UW"]
unifrac_presab_comm

# Unifrac para dados de abundância:

require(GUniFrac)
unifrac_abund_comm<-GUniFrac(comm, tree)$unifracs
unifrac_abund_comm[,,"d_1"]
unifrac_abund_comm

source("beta.pd.decompo.R")
beta_pd=beta.pd.decompo(comm, tree, type="both",output.dist=F, random=F)
beta_pd$betadiv

require(picante)
comdist_comm<-comdist(comm, cophenetic(tree), abundance.weighted = FALSE)
                      

require(picante)
comdistnt_comm<-comdistnt(comm, cophenetic(tree), abundance.weighted = FALSE, exclude.conspecifics = FALSE)

require(picante)
# Rao D
raoD_beta<-picante::raoD(comm, tree)$Dkl

# Rao H
raoH_beta<- picante::raoD(comm,tree)$H

require(devtools)
install_github("vanderleidebastiani/PCPS")
require(PCPS)
pcps_comm<-PCPS::pcps(comm,cophenetic(tree),method="bray",squareroot = TRUE)
#Matriz P:
P<-pcps_comm$P

#Matriz DP de dissimilaridades pareadas entre comunidades:
require(vegan)
DP<- sqrt(vegdist(P,method="bray"))

# Autovalores dos PCPS:
AV<-pcps_comm$values

# Autovetores dos PCPS (scores das comunidades):
Vec<- pcps_comm$vectors

# Scores das espécies no espaço dos vetores de sítio:
Sp<- summary(pcps_comm)$scores$scores.species

# Correlações entre as abundâncias das espécies ponderadas pela filogenia e os vetores dos sítios:
Cor_sp<- pcps_comm$correlations

# Usando ADONIS - Matriz DP (fuzzy weighting):
p.sig_comm<-matrix.p.sig(comm, phylodist=cophenetic(tree), envir=as.data.frame(env$E), formula= p.dist~env$E, FUN = FUN.ADONIS2.margin, method.p="bray", sqrt.p=TRUE, runs = 99, checkdata=FALSE)
p.sig_comm

# Usando vetores de PCPS:
pcps.sig_comm <- pcps.sig(comm, phylodist=cophenetic(tree), envir = as.data.frame(env$E), FUN = FUN.GLM, method = "bray", squareroot=TRUE, formula = pcps.1~env$E, choices = 1, runs = 99,checkdata=FALSE)
pcps.sig_comm

# ATIVIDADE 4

require(SYNCSA)
org<-organize.syncsa(comm,atributos,phylodist=cophenetic(tree),envir=env)
comm<-org$community
phydist<-org$phylodist
traits<-org$traits
env<-org$environmental

require(phytools)
k.T1<-phylosig(tree,traits$T1,method="K",test=T)
k.T2<-phylosig(tree,traits$T2,method="K",test=T)

require(PVR)
pvr_class<- PVRdecomp(tree)
pvr_reg.T1<-PVR(pvr_class,trait=traits$T1,envVar=NULL,method="moran")
sel.vec.brown<-pvr_reg.T1@Selection$Vectors
r2.pvr.T1<-pvr_reg.T1@PVR$R2
res.pvr.T1<-pvr_reg.T1@PVR$Residuals

require(picante)
require(cluster)

# Calcular distâncias fenotípicas:
dis.T1<-as.matrix(daisy(as.data.frame(res.pvr.T1), metric="euclidean",stand=FALSE))
rownames(dis.T1)<-colnames(comm)
colnames(dis.T1)<-colnames(comm)
# Calcular mfd – mean functional distance usando resíduos filogenéticos de T1:
ses.mfd.T1=ses.mpd(comm, dis.T1,null.model="independentswap",abundance.weighted=TRUE,runs=999,iterations=1000)
ses.mfd.T1=as.numeric(ses.mfd.T1[,6])

# Estimar a associação entre ses.mfd.T1 e env:
lm.ses.mfd.T1<-lm(ses.mfd.T1~E,data=env)
summary.lm(lm.ses.mfd.T1)

source("cwm.sig.R")

# Ignorando o efeito da filogenia:
cwmsig.T1.phylo<-cwm.sig(comm=comm,traits=traits,envir=env, formula="T1~E", PGLS=FALSE,tree=tree,runs = 999)
summary.lm(cwmsig.T1.phylo$model)
# Controlando o efeito da filogenia:
cwmsig.T1.nophylo<-cwm.sig(comm=comm,traits=traits,envir=env, formula="T1~E", PGLS=TRUE,tree=tree,runs = 999)
summary.lm(cwmsig.T1.nophylo$model)

# ΑTIVIDADE 5

require(ade4)

tree<-("(((Struthioniformes:21.8,Tinamiformes:21.8):4.1,((Craciformes:21.6,Galliformes:21.6):1.3,Anseriformes:22.9):3)A:2.1,(Turniciformes:27,(Piciformes:26.3,((Galbuliformes:24.4,((Bucerotiformes:20.8,Upupiformes:20.8):2.6,(Trogoniformes:22.1,Coraciiformes:22.1):1.3):1)B:0.6,(Coliiformes:24.5,(Cuculiformes:23.7,(Psittaciformes:23.1,(((Apodiformes:21.3,Trochiliformes:21.3):0.6,(Musophagiformes:20.4,Strigiformes:20.4):1.5):0.6,((Columbiformes:20.8,(Gruiformes:20.1,Ciconiiformes:20.1):0.7):0.8,Passeriformes:21.6):0.9):0.6):0.6):0.8):0.5):1.3):0.7):1);")

tree.phylog<-newick2phylog(tree,add.tools = T)


ori<-originality(tree.phylog,method=5)
# method = 1 = Vane-Wright et al.'s (1991) node-counting index 
# 2 = May's (1990) branch-counting index
# 3 = Nixon and Wheeler's (1991) unweighted index, based on the sum of units in binary values 
# 4 = Nixon and Wheeler's (1991) weighted index 
# 5 = QE-based index 
# 6 = Isaac et al. (2007) ED index 
# 7 = Redding et al. (2006) Equal-split index

# Ou

require (picante)
es.spp<-evol.distinct(tree,type="equal.splits",scale = FALSE,use.branch.lengths = TRUE)[,2]
# type = equal.splits (opção 7 da função anterior) ou fair.proportion (opção 6 da função anterior).


require(Herodotools)
require(PCPS)


#Idade das assembleias (Van Dijk et al. 2021)

require (Herodotools)

# Carregar objeto descrevendo sítios por região biogeográfica

biogeo<-read.table("biogeo.txt",h=T)
class(biogeo)

# Carregar objeto descrevendo os nós da filogenia por região biogeográfica, a partir de uma reconstrução ancestral, via BioGeoBears ou alguma outra ferramenta de reconstrução de caracteres ancestrais (e.g. função ace do pacote ape):

ancestral.area<-read.table("node_biogeo.txt",h=T)

# Calcula a idade das assembleias baseada na chegada do clado:
Idade<- calc_age_arrival(comm, tree, ancestral.area, biogeo)


# Definir número máximo de clusters a serem analisados (argumento max.n.clust da função calc_evoregions):
pcps_comm<-PCPS::pcps(comm,cophenetic(tree),method="bray",squareroot = TRUE)
Val<-pcps_comm$values
# Quatro eixos de pcps  apresentam autovalores > 5% da variação total em P.
Vec<- as.data.frame(pcps_comm$vectors)
dim(comm)
dim(Vec)

class(Vec)

n.clust3<-find_max_nclust(x=Vec, threshold= 4,max.nclust=3)
n.clust4<-find_max_nclust(x=Vec, threshold= 4,max.nclust=4)
n.clust5<-find_max_nclust(x=Vec, threshold= 4,max.nclust=5)
n.clust6<-find_max_nclust(x=Vec, threshold= 4,max.nclust=6)
n.clust7<-find_max_nclust(x=Vec, threshold= 4,max.nclust=7)
n.clust8<-find_max_nclust(x=Vec, threshold= 4,max.nclust=8)
n.clust9<-find_max_nclust(x=Vec, threshold= 4,max.nclust=9)

# Estimar evoregions:
evoreg<-calc_evoregions(comm, tree, max.n.clust = 8)


biogeo<-read.table("biogeo.txt",h=T)

# Carregar objeto descrevendo os nós da filogenia por região biogeográfica, a partir de uma reconstrução ancestral, via BioGeoBears ou alguma outra ferramenta de reconstrução de caracteres ancestrais (e.g. função ace do pacote ape):
ancestral.area<-read.table("node_biogeo.txt",h=T)

require(Herodotools)
calc_insitu_metrics(comm, tree, ancestral.area, biogeo)

devtools::install_github("GabrielNakamura/Herodotools",force = TRUE)

