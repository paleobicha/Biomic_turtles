# -----------------------------------------------------------------------
# Project: Sebas Tortugas
# File name: 
# Last updated: 2024-01-17
# Author: Sara Gamboa
# Email: saragamb@ucm.es
# Repository: https://github.com/paleobicha/Biomic_turtles
# -----------------------------------------------------------------------
###Cargamos las librerías
library(ape)
library(picante)
library(dplyr)
library(geiger)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)

###Cargamos nuestros archivos (filogenia, datasets, etc.)
#La filogenia, y ploteamos para ver que todo va bien
phy <- read.nexus("turtletree.tre")
plot.phylo(phy, show.tip.label = T, font=1, cex=0.3)
###Ahora anadimos nuestra informacion de BSI. 
###Lo hacemos con una tabla tabulada donde aparezca el nombre del taxon, tal y como aparece en el arbol.
bsi<-read.csv("TBSI.csv" ,header=TRUE, sep=";")
biomas <- read.csv("TBiomas.csv", header = T, sep=";")
head(bsi)
head(biomas)

###Calculamos la diferencia evolutiva (Evolutionary Distinctiveness) de los taxones de nuestro arbol
###Luego hacemos la inversa para sacar DR
ED_tips <- evol.distinct(phy, type = "equal.splits", scale = FALSE, use.branch.lengths = TRUE)
ED <- ED_tips$w
DR <- 1/ED

###Damos los nombres de nuestros taxones a la tabla con los DR
names(DR) <- ED_tips$Species
head(DR)
hist(DR)

###Ahora asignamos a cada especie su DR

bsi$DR <- 0
head(bsi)
dim(bsi)

i <- 1
for (i in 1:length(bsi$tip.label)) {
  bsi$DR[i] <- DR[which(names(DR)==bsi$tip.label[i])]
}

plot(bsi[,2:3])
bsi <- bsi[!(bsi$bsi %in% 0),]
bsi$bsi <- as.factor(bsi$bsi)
unique(bsi$bsi)

###Vamos a hacer unos boxplots de nuestros resultados. Siguiendo otros trabajos como Gamboa et al., 2022, he coloreado los especialistas de rojo,
## ls generalistas moderados de amarillo y los generalistas extremos (BSI>=5) de azul.
bsi_colour <- c("#FF5959", "#FFAD5A", "#FFAD5A", "#FFAD5A", "#4F9DA6")
boxplot(bsi$DR~bsi$bsi, log="y", axes=T, ylim=c(0.01,1))

p <- ggplot(bsi, aes(x=bsi, y=DR, group=bsi)) + 
  geom_boxplot(fill=bsi_colour, color = bsi_colour, alpha = 0.8, outlier.colour = "darkgrey", outlier.shape = 20, outlier.alpha = 0.5) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color=bsi_colour, fill=bsi_colour, show.legend=TRUE)

###Sobre un boxplot más basico he creado dos alternativas, una normal y otra logaritmica para que puedas elegir la que mas te guste
p + scale_x_discrete(name = "BSI") +
  scale_y_continuous(name = "Diversification rate", trans = "log10") + 
  theme(panel.background = element_rect(fill = "white", colour = "darkgrey", linewidth = 0.8), 
        panel.grid.major.y = element_line(colour = "darkgrey", linewidth  = 0.3), panel.grid.major.x = element_line(colour = "NA")) +
  scale_y_sqrt(name = "Diversification rate \n (log scaled)")

p + scale_x_discrete(name = "BSI") +
  scale_y_continuous(name = "Diversification rate") + 
  theme(panel.background = element_rect(fill = "white", colour = "darkgrey", linewidth = 0.8), 
        panel.grid.major.y = element_line(colour = "darkgrey", linewidth  = 0.3), panel.grid.major.x = element_line(colour = "NA"))

####Vamos a calcular unas varaibles estadísticas para hacer plots más bonitos
means <- aggregate(DR ~  bsi, bsi, mean)
max <-aggregate(DR ~  bsi, bsi, max)
min <-aggregate(DR ~  bsi, bsi, min)
sd <-aggregate(DR ~  bsi, bsi, sd)

###Creamos la funcion Summary para sacar los valores de la muestra 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###Vamos a sacar una cajita para plotear y guardar el plot más facilmente
###Esta funcion (quartz) es especifica para Mac, si tienes PC debes buscar como se llama, porque no lo se
###Las medidas se dan en pulgadas, así que basta con poner el nº de cm y dividir por 2.54
quartz(width = 8/2.54, height = 8/2.54)

pd <- position_dodge(0.1)
tgc <- summarySE(bsi, measurevar = "DR", groupvars = "bsi")


####Vamos a hacer un plot muy sencillo y visual.
##En esta versión los bigotes del boxplot representan el error estandar de la muestra
ggplot(tgc, aes(x=bsi, y=DR)) + 
  geom_errorbar(aes(ymin=DR-se, ymax=DR+se, color=factor(bsi)), width=.1, position=pd) +
  scale_color_manual("BSI", breaks=c(1:5),values=bsi_colour) +
  geom_line(position=pd) +
  geom_point(position=pd, colour=bsi_colour, cex=1.8) +
  scale_fill_manual("Biome Specialization Index", breaks=c(1:5),values=bsi_colour) +
  xlab("Biome Specialization Index") + ylab("Diversification rate")+
  theme(panel.background = element_rect(fill = "white", size = 0.8), 
        panel.grid.major.y = element_line(colour = "darkgray", size = 0.3), panel.grid.major.x = element_line(colour = "NA"),
        axis.line = element_line(colour = "dimgrey"), panel.border = element_blank(),
        text = element_text(size=8))

###En este segundo caso los bigotes representan el intervalo de confianza al 95%
ggplot(tgc, aes(x=bsi, y=DR)) + 
  geom_errorbar(aes(ymin=DR-ci, ymax=DR+ci, color=factor(bsi)), width=.1, position=pd) +
  scale_color_manual("BSI", breaks=c(1:10),values=bsi_colour) +
  geom_line(position=pd) +
  geom_point(position=pd, colour=bsi_colour, cex=1.8) +
  scale_fill_manual("Biome Specialization Index", breaks=c(1:10),values=bsi_colour)+
  xlab("Biome Specialization Index") + ylab("Diversification rate")+
  theme(panel.background = element_rect(fill = "white", size = 0.8), 
        panel.grid.major.y = element_line(colour = "darkgray", size = 0.3), panel.grid.major.x = element_line(colour = "NA"),
        axis.line = element_line(colour = "dimgrey"), panel.border = element_blank(),
        text = element_text(size=8))
###Como ves los resultados son muy similares, ya que son los mismos datos, pero el error estándar suele magnificar las diferencias entre los grupos
##Mejorando la visualización. Elije el que prefirais.


######
##Ahora vamos a los análisis por biomas
biomes <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII")
col_biomes <- c("#813BB0", "#84B43C", "#1C542D", "#EEB003", "#A43749", "#A13678", "#2E8A3C", "#CB9162", "#2D866A", "#308091")

bsi_or <- bsi[order(bsi$tip.label),]
par(mfrow=c(4, 2), mar=c(2,4,0.5,0.5))

for(b in biomes){
  #quartz(height = 5.6/2.54, width = 9/2.54)
  presentes <- biomas[,which(colnames(biomas)==b)]==1
  bsi_bioma <- bsi$bsi[presentes]
  DR_bioma <- bsi$DR[presentes]
  boxplot(DR_bioma~bsi_bioma,col=col_biomes[which(biomes==b)])
}

dim(biomas)
especialistas <- biomas[rowSums(biomas[,2:11])==1,]
names(especialistas) <- c("Species", "I", "II", "II/III", "III", "IV", "V", "VI", "VII", "VIII", "IX")
head(especialistas)
dim(especialistas)
m1 <- col(especialistas[-1]) * especialistas[-1]
i1 <- m1 != 0 
especialistas[-1][i1] <- rep(colnames(m1), each = nrow(m1))[i1]
especialistas

colapse <- function(x)
{
  y <- paste0(names(x)[x!=0],collapse = '-')
  if(y=='') {y=0}
  return(y)
}
especialistas$bio <- apply(especialistas[,-1],1,colapse )

especialistas$DR <- 0
for (i in 1:length(especialistas$Species)) {
  especialistas$DR[i] <- DR[which(names(DR)==especialistas$Species[i])]
}

stn_data <-data.frame(especialistas$Species, especialistas$DR, as.factor(especialistas$bio))
names(stn_data) <- c("Species", "DR", "bio")
###Graficos por biomas
tgc_s <- summarySE(stn_data, measurevar = "DR", groupvars = "bio")
tgc_s$bio <- as.factor(tgc_s$bio)
pd <- position_dodge(0.1) # move them .05 to the left and right
col_biomes <- c("#813BB0", "#84B43C", "#1C542D", "#EEB003", "#A43749", "#A13678", "#2E8A3C")

quartz(width = 8/2.54, height = 8/2.54)


ggplot(tgc_s, aes(x=bio, y=DR, group=bio)) + 
  geom_errorbar(aes(ymin=DR-se, ymax=DR+se, color=factor(bio)), width=.1, position=pd) +
  scale_color_manual("Bioma", breaks=c(1:7),values=col_biomes) +
  geom_line(position=pd) +
  geom_point(position=pd, colour=col_biomes, cex=1.8) +
  scale_fill_manual("Bioma",values=col_biomes) +
  xlab("Biome") + ylab("Diversification rate") +
  theme(panel.background = element_rect(fill = "white", linewidth = 0.8), 
        panel.grid.major.y = element_line(colour = "darkgray", linewidth = 0.3), panel.grid.major.x = element_line(colour = "NA"),
        axis.line = element_line(colour = "dimgrey"), panel.border = element_blank(),
        text = element_text(size=8))


ggplot(tgc_s, aes(x=bio, y=DR)) + 
  geom_errorbar(aes(ymin=DR-ci, ymax=DR+ci, color=factor(bio)), width=.1, position=pd) +
  scale_color_manual("Bioma", breaks=c(1:10),values=col_biomes) +
  geom_line(position=pd) +
  geom_point(position=pd, colour=col_biomes, cex=1.8) +
  scale_fill_manual("Bioma", breaks=c(1:10),values=col_biomes)+
  xlab("Biome") + ylab("Diversification rate")+
  theme(panel.background = element_rect(fill = "white", size = 0.8), 
        panel.grid.major.y = element_line(colour = "darkgray", size = 0.3), panel.grid.major.x = element_line(colour = "NA"),
        axis.line = element_line(colour = "dimgrey"), panel.border = element_blank(),
        text = element_text(size=8))


#####Análisis estadísticos

##PGLS
bsi$bsi <- as.numeric(bsi$bsi)
pglsModel <- gls(DR ~ bsi, correlation = corBrownian(phy = phy),
                 data = bsi, method = "REML")
pglsModel <- gls(DR ~ bsi, correlation = corPagel(1, phy = phy, fixed = FALSE),
                 data = bsi, method = "REML")
summary(pglsModel)


##ANOVA
#Extraemos y preparamos los datos
dim(especialistas)
stat_bsi <- rep(0, 149)
names(stat_bsi) <- stn_data$Species

for (i in 1:149) {
  stat_bsi[i] <- as.character(stn_data$bio[which(stn_data$Species==names(stat_bsi)[i])])
}
stat_bsi <- as.factor(stat_bsi)

stat_DR <- rep(0, 149)
names(stat_DR) <- stn_data$Species
for (i in 1:149) {
  if ((names(stat_DR)[i] %in% bsi$tip.label)==F) {
    stat_DR[i] <- 0
  } else {stat_DR[i] <- stn_data$DR[which(stn_data$Species==names(stat_DR)[i])]}
}

for (i in 1:149) {
  stat_DR[i] <- stn_data$DR[which(stn_data$Species==names(stat_DR)[i])]
}

summary(stat_bsi)

anova_m <- aov.phylo(stat_DR ~ stat_bsi, phy, nsim = 1000, 
                 test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"))
print(attributes(anova_m)$summary)

ggdensity(stat_bsi, 
          main = "Density plot DR",
          xlab = "DR")


phy$tip.label

#######DR TREE#########

or <- match(names(DR),phy$tip.label)
ladd_tree <- ladderize(phy)
anc <- ace(DR,ladd_tree, type = "continuous",method="ML")$ace
ancestry <- fastAnc(ladd_tree, DR)
acr <- contMap(ladd_tree, log(DR), show.tip.label=F, 
               fsize=c(0.2,0.7), outline=F, method="fastAnc", plot=F, ftype ="off", sig = 2,
               anc.states=anc, type="fan",rotate.tree=270, open.angle=180,direction="rightwards",legend = F, part= 0.5) 
plotBranchbyTrait(ladd_tree, log(DR_ordenado), mode="tips", palette=palette, legend=F, type="fan", show.tip.label=F,open.angle=180,rotate=180, y.lim=100, edge.width=1.6, direction="rightwards")
col_scale <- rev(brewer.pal(10, "RdYlBu"))
col_scale <- col_scale[1:10]
palette <- colorRampPalette(col_scale)
scale_cols <- palette(21)
acr$cols[] <- palette(1001)

biomes <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X")
col_biomes <- c("#813BB0", "#84B43C", "#1C542D", "#EEB003", "#A43749", "#A13678", "#2E8A3C", "#CB9162", "#2D866A", "#308091")

quartz(height = 29.7/2.54,width = 21/2.54) ##DINA4 normal
quartz(height = 21/2.54,width = 29.7/2.54) ##DINA4 apaisado
quartz(height = 10/2.54,width = 10/2.54)
quartz(height = 20/2.54,width = 20/2.54)
par(oma=c(4,2,4,2))
par(oma=c(2,2,2,2))
par(mar=c(8,8,8,8))

plotBranchbyTrait(ladd_tree, log(DR_ordenado), mode="tips", palette=palette, legend=F, type="fan", show.tip.label=F,open.angle=180,rotate=180, y.lim=100, edge.width=1.6, direction="rightwards")

offset <- 1.5

for(b in biomes){
  trans <- rep(0.1, length(ladd_tree$tip.label))
  presentes <- biomas[,which(colnames(biomas)==b)]==1
  presentes <- biomas$Species[presentes]
  trans[ladd_tree$tip.label %in% presentes] <- 1
  col <- col_biomes[biomes %in% b]
  tiplabels(pch=19, thermo=NULL, offset=offset, cex=0.25, col=alpha(col, trans))
  offset <- offset+1.8
}

which(names(DR_ordenado)=="Graptemys_pearlensis_RCT242")

DR_mod <- DR_ordenado[-c(1,2,5,33,82,86,123)]
length(DR_mod)
tips <- c("AA_Alligator","AA_Ggallus","Aldabrachelys_gigantea_HBS118599","Chelodina_timorensis_AMNH167209",
          "Elseya_albagula_AMR123067","Elseya_lavarackorum_321", "Graptemys_pearlensis_RCT242")
phy_mod <- drop.tip(phy, tips)

ladd_mod <- ladderize(phy_mod)
plotBranchbyTrait(ladd_mod, log(DR_mod), mode="tips", palette=palette, legend=F, type="fan", show.tip.label=F,open.angle=180,rotate=180, y.lim=100, edge.width=1.6, direction="rightwards")

offset <- 1.5

for(b in biomes){
  trans <- rep(0.1, length(ladd_tree$tip.label))
  presentes <- biomas[,which(colnames(biomas)==b)]==1
  presentes <- biomas$Species[presentes]
  trans[ladd_tree$tip.label %in% presentes] <- 1
  col <- col_biomes[biomes %in% b]
  tiplabels(pch=19, thermo=NULL, offset=offset, cex=0.35, col=alpha(col, trans))
  offset <- offset+1.8
}

pdf(file="DRtree_Testu.pdf", width=200, height=200)
dev.off()
