# Code for the Bayesian estimation of Main effects

# Specify the no of MCMC iterations as "ite"
ite=50000;

# Specify the no of MCMC iterations to burn-in (remove from the final calculation) as "br"
br=10000;

# Read the phenotype file 
dat=read.table("pheno.txt",header=T)

# Read the marker information 
M=read.table("geno_ord_2331_epi.txt",header=T,row.names=1)

# read the map file 
map=read.table("map_ord_2331_epi.txt",header=T)

# Define the no of markers
mark=dim(M)[1]

#Define the no of lines with the phenotype information
line=dim(M)[2]

#Define the marker information as a matrix
M=as.matrix(M)

# Define the phenotype as a vector called 'y' 
y=as.vector(dat$phe)

#Create vector 'A' to store the mean of the phenotype
A=vector()
#Define a vector 'va' to store residual
va=vector()

t=1
#Assign the phenotype mean to 'A'
A[1]=mean(y)

#Matrices to store the current and proposed values of marker effects respectively
B=matrix(0,nrow=mark,ncol=2) #marker x 2
vb=matrix(0,nrow=mark,ncol=2) #marker x 2

#Set first values as 0
B[,1]=rep(0,mark)
#initial values of marker effects
vb[,1]=0.5*rep(1,mark)
va[1]=0.5
#Define vectors to store the final values
fA=vector()
fva=vector()

#Define a matrix 'res' to store the marker effects.
res=matrix(0,nrow=mark,ncol=1)

#Start the MCMC chain
for (s in 1:ite ) 
	{
	#update A(mean)
	mA=mean(y-t(B[,t])%*%M)
	vA=va[t]/line
	A[t+1]=mA+sqrt(vA)*rnorm(1)

	#update B (marker effect)
	B[,t+1]=B[,t]
	Bmxj=t(B[,t+1])%*%M

	for (j in 1:mark ) 
		{	
		Bmxj=Bmxj-B[j,t+1]*M[j,] # to keep k not j
		e=sum(M[j,]*(y-A[t+1]-Bmxj))

		d1=sum(M[j,]^2) #Sig^2_ij

		mB=vb[j,t]*e/(vb[j,t]*d1+va[t])

		vB=vb[j,t]*va[t]/(vb[j,t]*d1+va[t])

		B[j,t+1]=mB+sqrt(vB)*rnorm(1)
		Bmxj=Bmxj+B[j,t+1]*M[j,]
		}
# update va the residual error
	c2=rchisq(1,line)
	e2=sum((y-A[t+1]-t(B[,t+1])%*%M)^2)
	va[t+1]=e2/c2
#update vb
	vb[,t+1]=(B[,t+1]^2)/rchisq(mark,1)
#store data
	fA[s]=A[t+1]
	fva[s]=va[t+1]
#Save the marker effects after the burn-in period 
	if(s>=br)
		{
	res=res+B[,t+1]
		}
	
#next round initial values	
	A[t]=A[t+1]
	B[,t]=B[,t+1]
	va[t]=va[t+1]
	vb[,t]=vb[,t+1]

	}

#Vector m_eff contains the final results (as the average of all MCMC chains after the burn-in period)
m_eff=res/(ite-br)


# Run the MCMC chain for 5 times and store the marker effects in a file called "all_run_main.txt"
############################################################



