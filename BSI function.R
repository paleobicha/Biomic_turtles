
#en el paréntesis debe ir el nombre del archivo de texto, que debe estar guardado en la carpeta datos BSI MC-primero cargar la función y luego ejecutar la familia.
BSI("Testudinidae")

BSI<-function(n){
name<-n	
archivo<-paste(name,".txt", sep="")
setwd("~/Desktop/datos BSI MC")
a<-read.table(file=archivo, header=FALSE)
bsidea<-rowSums(a)
maxbsi<-max(bsidea)
sp<-nrow(a)
sumbiom<-colSums(a)
sim=10000
mc<-array(NA, dim=c(sp,10,sim))
m<-array(NA, dim=dim(a))
for (i in 1:sim){
for (k in 1:10){
	m[,k]<-sample(a[,k])
	}
	mc[,,i]<-m
	}
		
	
bsimc<-matrix(data=NA, nrow=sp, ncol=sim)
bsi<-1
for(j in 1:sim)
	{
	bsi<-(rowSums(mc[,,j]))
	bsimc[,j]<-(bsi)
	}
		
bsimc2<-matrix(data=NA, nrow=sp, ncol=sim)
for(i in 1:(sp*sim)){
	if(bsimc[i]==0){
		bsimc2[i]<-NA
		
		}
	else(bsimc2[i]=bsimc[i])
	
	}	
	BSI<-matrix(data=NA, nrow=1, ncol=10)
	for(b in 1:maxbsi){
	BSI[b]=sum(bsidea==b)
	}
	BSI[is.na(BSI)] <- 0		

##tabla BSI		
table<-matrix(data=NA, nrow=10, ncol=6)
table2<-matrix(data=NA, nrow=10, ncol=2)
for(k in 1:10){
sumx<-1
sm<-1
for(i in 1:sim){
	sumx[i]=sum(bsimc[,i]==0)
	sm[i]=100*(sum(bsimc[,i]==k))/(sp-sumx[i])	
	}

obs=100*BSI[k]/sp
obs2<-format(obs, digits=2)


if(max(sm)==0){pvalue=1}
else{if(obs==0){pvalue=sum(sm==obs)/sim}
	
	else{if (sum(sm>obs)/sum(sm<=obs)>1){
		pvalue=sum(sm<obs)/sum(sm>=obs)
		}
		else{pvalue=sum(sm>obs)/sum(sm<=obs)}}}

	
	
pvalue2<-format(pvalue, digits=3)
table[k,1]<-k
table[k,2]<-format(100*BSI[k]/sp, digits=3)
table[k,3]<-format(mean(sm), digits=2)
table[k,4]<-format(sd(sm), digits=2)
max2<-format(max(sm), digits=2)
min2<-format(min(sm), digits=2)
table[k,5]<-paste(min2,"-",max2)
table[k,6]<-pvalue2
table2[k,1]<-k
table2[k,2]<-mean(sm)*sp/100
}

show(table)
write.table(table, file=paste(name," BSI.txt", sep=""), quote=FALSE,sep="	",  row.names=FALSE, col.names=FALSE)	


hist(bsidea, breaks=c(0,1,2,3,4,5,6,7,8,9,10),xaxp=c(0,10,10), labels=FALSE, xlab="BSI", ylab="% species", ylim=c(0,sp/2),main=paste(name), right=TRUE,las=1,axes=F)
axis(2,c(0,sp/4,sp/2),labels=c(0,25,50),las=1)
axis(1,c(0.5:9.5),labels=c(1:10))
points(table2[,1]-0.5, table2[,2], pch=20, cex=2)
legend("topright",c(paste("mean =",format(mean(bsidea),digits=3)),paste("sd =",format(sd(bsidea), digits=3)),paste("n =",nrow(a))),bty="n",cex=1.1,yjust=1)
dev.copy2pdf(device=quartz,file=paste(name,".pdf",sep=""))
	
		
###TABLA BSI=1
tablabsi1<-matrix(data=NA, nrow=10, ncol=8)
tablabsi1[,1]<-c("I","II","II/III","III","IV","V","VI","VII","VIII","IX")
tablabsi1[,2]<-sumbiom##sp en cada bioma
bioma<-1
bsis<-1
biomebsi<-1
col5<-1
for(v in 1:10){
	bioma=mc[,v,]*bsimc2
	nbsi<-matrix(data=NA, nrow=sim, ncol=1)
	biomebsi=a[,v]*bsidea
		for(z in 1:1){
		tablabsi1[v,3]<-sum(biomebsi==z)##especialistas en cada bioma
		
		
		nbsi<-colSums(bioma==z, na.rm=TRUE)/(sumbiom[v])*100##% de especialistas dentro de las sp de ese bioma de las simulaciones
		
					
		}
	
	
	percentage<-(sum(biomebsi==z)/(sumbiom[v]))*100##% de especialistas dentro de las sp de ese bioma OBSERVADOS para los cálculos
	
	if(sumbiom[v]==0){pvalue=1
		pvalue2<-format(pvalue, digits=3)
		tablabsi1[v,4]<-0##% de especialistas dentro de las sp de ese bioma OBSERVADOS##
		tablabsi1[v,5]<-0
		tablabsi1[v,6]<-0
		min<-0
		max<-0
		tablabsi1[v,7]<-paste(min,"-",max)
		tablabsi1[v,8]<-pvalue2	
		
		
		
		}
		else{
			tablabsi1[v,4]<-format(((sum(biomebsi==z))/(sumbiom[v])*100), digits=4)##% de especialistas dentro de las sp de ese bioma OBSERVADOS##
					
			if (sum(nbsi>percentage)/sum(nbsi<=percentage)>1){
			pvalue=sum(nbsi<percentage)/sum(nbsi>=percentage)
			}
			else{pvalue=sum(nbsi>percentage)/sum(nbsi<=percentage)}
			pvalue2<-format(pvalue, digits=3)
	
		
		
		tablabsi1[v,5]<-format(mean(nbsi), digits=3)
		tablabsi1[v,6]<-format(sd(nbsi), digits=2)
		min<-format(min(nbsi), digits=3)
		max<-format(max(nbsi), digits=3)
		tablabsi1[v,7]<-paste(min,"-",max)
		tablabsi1[v,8]<-pvalue2	
			}
		
			
		}
show(tablabsi1)		
write.table(tablabsi1, file=paste(name," BSI=1.txt", sep=""), quote=FALSE,sep="	",  row.names=FALSE, col.names=FALSE)		
		
		
		}
		


