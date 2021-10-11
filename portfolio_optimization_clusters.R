rm(list=ls())
graphics.off()
library(ggplot2)
library(DEoptim)
#first need many assets returns which we can get from quantmod
require(quantmod)

#csv for all s&p500 tickers is available at https://datahub.io/core/s-and-p-500-companies
temp = readr::read_csv('constituents_csv.csv')
(tickers = as.array(temp$Symbol))
#some tickers arent available so remove them from the array
tickers = tickers[tickers!="BRK.B"]
tickers = tickers[tickers!="BF.B"]
tickers = tickers[tickers!='CARR']
tickers = tickers[tickers!='CTVA']
tickers = tickers[tickers!='DOW']
tickers = tickers[tickers!='OGN']
tickers = tickers[tickers!='OTIS']


#the next line will take a minute to get all the data
myStocks = lapply(tickers, function(x) {getSymbols(x, 
                                                   from = "2016/12/31",to = "2019/12/31",
                                                   periodicity = "daily",
                                                   auto.assign=FALSE)} )
names(myStocks) = tickers
head(myStocks$TSLA)
tail(myStocks$TSLA)
head(myStocks$AMD)
head(myStocks$INTC)

#we have a list of stock prices and many vars but really only want 
# adjusted close.
adjusted.cls = lapply(myStocks, Ad) #Ad is an extraction function from quantmod (extracts adjusted price)
adjusted.cls = do.call(merge,adjusted.cls)
#View(adjusted.cls)
dim(adjusted.cls)
is.matrix(adjusted.cls)
#now need to get returns for all
returns500 = apply(adjusted.cls, 2,Delt)
returns500 = returns500[2:dim(returns500)[1],] #remove NA row created
summary(returns500[,1:5]*252)
summary(returns500[,450:455]*252)
#need to remove NAs
dim(returns500[complete.cases(returns500),]) #would remove 500ish scenarios
dim(returns500[,colSums(is.na(returns500))==0]) #would remove 5 assets
#opt to remove 5 assets with missing values from the data
returns500 = returns500[,colSums(is.na(returns500))==0]

#now we can look for an optimal portfolio within the S&P500
#Guy Yollins optimization function for DE:
# Note: L is the loss threshold
# more info on Omega here http://www.performance-measurement.org/Winton2003.pdf
optOmega = function(x,rets,L){ #make sure returns data units match L units
  ret = rets%*%x
  obj = -(mean(pmax(ret-L,0))/mean(pmax(L-ret,0)))
  pen = 1000*(1-sum(x))^2
  return(obj+pen)
}

cvarOpt = function(rmat,alpha=.05,rmin = 0,wmin = 0, wmax = 1,weight.sum=1){
  rmat = as.matrix(rmat)
  require(Rglpk) #works extrememly similar to lp_solve
  n = ncol(rmat) #number of assets
  s = nrow(rmat) #number of scenarios
  averet = colMeans(rmat)
  #build A matrix
  Amat = rbind(cbind(rbind(1,averet),matrix(data = 0,nrow = 2,ncol = s+1)),
               cbind(rmat,diag(s),1))
  objL = c(rep(0,n),rep(-1/(alpha*s),s),-1) #objective
  bvec = c(weight.sum,rmin,rep(0,s)) #rhs vector
  dir.vec = c('==','>=',rep('>=',s)) #direction vector
  bounds = list(lower = list(ind = 1:n,val=rep(wmin,n)),
                upper = list(ind=1:n,val=rep(wmax,n))) #only for rglpk
  
  
  res = Rglpk_solve_LP(obj = objL,mat=Amat,dir = dir.vec,rhs = bvec,
                       types = rep('C',length(objL)),max=T,bounds=bounds)
  w = as.numeric(res$solution[1:n])
  return(list(w=w,status=res$status,obj_val=res$optimum))
}

sharpeMaxDE = function(x,rets,riskfreerate){ #maximizes sharpe ratio (minimzes -returns/risk)
  rets = as.matrix(rets)
  z=-(mean(rets%*%x)-riskfreerate)/(var(rets%*%x)^.5) #assumes constant risk free rate
  pen = 1000*(1-sum(x))^2 #weights constraint
  return(z+pen)
}

lower = rep(0,dim(returns500)[2]) #now have 498 assets
upper = rep(1,dim(returns500)[2])

#generate random portfolios to visualize the efficient frontier
###############################################################################################
#start function
makeRandPort = function(rets,nAssets,nPorts,Opts = FALSE){#inputs:returns data,number of assets,num of portfolios to generate,whether or not to optimize each random portfolio
  nr = nrow(rets)
  nc = ncol(rets)
  if(nAssets>nc) print('CANNOT HAVE MORE ASSETS THAN COLUMNS IN RETURNS DATA')
  w = runif(nAssets*nPorts)
  w = matrix(w,ncol = nAssets)
  w=t(apply(w, 1, function(x){x/sum(x)}))#make weights sum to one
  #print(w)
  results = cbind.data.frame('PortNum'=NA,'MU'=NA,'SD'=NA,'Type'=NA)
  
  if(Opts==TRUE){
    w.all = vector(mode = 'list',length = 5*nPorts)
    for (i in 1:nPorts) {
      indexes = sample(1:nc,nAssets,replace = FALSE)
      p = rets[,indexes]
      #rand ports
      results[i,] =c(i,mean(p%*%w[i,]),var(p%*%w[i,])^.5,'Rand.')
      w.all[[i]]= cbind.data.frame('weights'=w[i,],'Stock'=colnames(p)) #make list to store all weights and stocks
      #cvar opt
      wcvar = cvarOpt(p)
      wcvar = wcvar$w
      results[i+nPorts,]=c(i,mean(p%*%wcvar),var(p%*%wcvar)^.5,'cVaR')
      w.all[[i+nPorts]]= cbind.data.frame('weights'=wcvar,'Stock'=colnames(p))
      #omega opt using DEoptim and L = 0
      wOmega = DEoptim(optOmega,lower = rep(0,nAssets) ,upper = rep(1,nAssets),rets = p, L=0, 
                       control = list(NP=10*(nAssets),itermax=12000,trace = FALSE,
                                      reltol=.0001,steptol=250,parallelType=1,
                                      F=.2,CR=.8))
      wOmega = wOmega$optim$bestmem
      results[i+2*nPorts,]=c(i,mean(p%*%wOmega),var(p%*%wOmega)^.5,'Omega')
      w.all[[i+2*nPorts]]= cbind.data.frame('weights'=wOmega,'Stock'=colnames(p))
      #Sharpe Opt assuming risk free rate of .02
      wSharpe = DEoptim(sharpeMaxDE,lower = rep(0,nAssets),upper=rep(1,nAssets),rets = p,riskfreerate=.02, 
                        control = list(NP=10*(nAssets),itermax=12000,trace = FALSE,
                                       reltol=.0001,steptol=250,parallelType=1,
                                       F=.2,CR=.8))
      wSharpe = wSharpe$optim$bestmem
      results[i+3*nPorts,] = c(i,mean(p%*%wSharpe),var(p%*%wSharpe)^.5,'Sharpe')
      w.all[[i+3*nPorts]]= cbind.data.frame('weights'=wSharpe,'Stock'=colnames(p))
      #even portfolio
      wEven = rep(1/nAssets,nAssets)
      results[i+4*nPorts,]=c(i,mean(p%*%wEven),var(p%*%wEven)^.5,'EqualW')
      w.all[[i+4*nPorts]]=cbind.data.frame('weights'=wEven,'Stock'=colnames(p))
      
      print('PROGRESS:')
      print(100*i/nPorts)
    }
  }
  
  if(Opts == FALSE){
    w.all = vector(mode = 'list',length = nPorts)
    for (i in 1:nPorts) {
      indexes = sample(1:nc,nAssets,replace = FALSE)
      p = rets[,indexes]
      results[i,] =c(i,mean(p%*%w[i,]),var(p%*%w[i,])^.5,'Rand.')
      w.all[[i]] = cbind.data.frame('weights'=w[i,],'Stock'=colnames(p)) #make list to store all weights and stocks
      
      print('PROGRESS:')
      print(100*i/nPorts)
    }
  }
  
  
  return(list(res=results,WandS = w.all))#return mu / sd dataframe and list of dataframes containing weights and stocks
}
#end function
##########################################################################################
ports = makeRandPort(rets=252*returns500,nAssets=6,nPorts=10,Opts = T)
ports$res$MU = as.numeric(ports$res$MU)
ports$res$SD = as.numeric(ports$res$SD)
ports$res$PortNum = as.integer(ports$res$PortNum)

ggplot(data = ports$res,aes(x=SD,y=MU,color = factor(Type)))+geom_point()+
  theme_classic()+
  scale_x_continuous(breaks = round(seq(0, max(ports$res$MU), by = 0.1),1)) +
  scale_y_continuous(breaks = round(seq(0, max(ports$res$SD), by = 0.1),1))

#show how each opt method affects the mean variance relative to equal weights
mv = ports$res # (didnt want to keep accessing the ports list) mv is mean varaince dataframe
centers = subset(mv,Type == 'EqualW')
legs = subset(mv,Type != 'EqualW')
#make cluster variables data
numClust = max(mv$PortNum)
mv.clusters = data.frame()
for (i in 1:numClust) {
  mini = data.frame(PortNum = c(rep(centers$PortNum[i], (length(unique(mv$Type))-1)), legs$PortNum[legs$PortNum==i]),
                    mu = c(rep(centers$MU[i], (length(unique(mv$Type))-1)), legs$MU[legs$PortNum==i]),
                    sd = c(rep(centers$SD[i], (length(unique(mv$Type))-1)), legs$SD[legs$PortNum==i]),
                    Type = c(rep(centers$Type[i], (length(unique(mv$Type))-1)), legs$Type[legs$PortNum==i]))
  mini$Clust = rep(1:(length(unique(mv$Type))-1),2) + (i-1)*(length(unique(mv$Type))-1)
  #print(mini)
  mv.clusters = rbind.data.frame(mv.clusters,mini)
}
rm(centers)
rm(legs)

ggplot(data = mv.clusters,aes(x=sd,y=mu))+geom_point(size=3)+
  geom_line(group = mv.clusters$Clust,aes(color=factor(PortNum),alpha=.7),size=1.5)+
  geom_text(aes(label=Type),hjust=0, vjust=0)+theme_classic()+
  scale_x_continuous(breaks = round(seq(0, max(mv.clusters$sd), by = 0.1),1)) +
  scale_y_continuous(breaks = round(seq(0, max(mv.clusters$mu), by = 0.1),1))
#Notice optimized portfolios' points are generally pushed up and left towards efficient frontier.
#Now lets visualize the frontier using only random portfolios