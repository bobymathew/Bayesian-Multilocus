# Code for the Bayesian estimation of three-way epistasis

# Specify the no of MCMC iterations as "ite"
ite=50000;

# Specify the no of MCMC iterations to burn-in (remove from the final calculation) as "br"
br=10000;

#the number of two-way pseudo markers 
ns=1000

#the number of three-way pseudo markers to be considered for the three-way epistasis
nth=1000


# Read the phenotype file 
dat=read.table("pheno.txt",header=T)

# Read the marker information 
M=read.table("geno_ord_2331_epi.txt",header=T,row.names=1)

# Define the no of markers
mark=dim(M)[1]

#Define the no of lines with the phenotype information
line=dim(M)[2]

#Define the marker information as a matrix
M=as.matrix(M)

# Define the phenotype as a vector called 'y' 
y=as.vector(dat$phe)


marker=t(M)
n=dim(marker)[1]
p=dim(marker)[2]


#Read the file with the marker effects for the 2331 markers based on 5 MCMC chains. Each column contains the marker effects(2331 markers) corresponding to a MCMC chain
eff=read.table("all_run_main.txt",header=T)

#Read the file with the epistatic effects based 1000 pseudo markers based on 5 MCMC chains. Each column contains the epistatic effects(1000 markers) corresponding to a MCMC chain
epi=read.table("all_run_epi.txt",header=T)


#Read the most important two-way interacting marker coordinates
sure_two=read.table("sure_sel_markers_two_way.txt",header=T)

#Read the most important two-way interacting marker coordinates found using program Three_way_SIS.R"
sel_three=read.table("three_way_sel_cordinates.txt",header=T)

# Calculate the mean of marker effects from 5 MCMC chains
m_eff=rowMeans(eff)
# Calculate the mean of two-way interaction effects from 5 MCMC chains
m_epi=rowMeans(epi)

# Calculate residual error(E1i)
err=y-(t(M)%*%m_eff)

#Create the pseudo marker based on the selected indices using the program 'Two_way_SIS.R'
newM=matrix(0,nrow=n,ncol=ns)

for (i in 1:ns ) 
	{

	newM[,i]=marker[,sel[i,1]]*marker[,sel[i,2]]

	}

#Calculate epistasis residuals (E2i)
err_th=err-(newM%*%m_epi)



#Create the pseudo marker for the Three-way epistasis search based on the selected indices using the program 'Three_way_SIS.R'
newM_TH=matrix(0,nrow=line,ncol=nth)

for (i in 1:nth ) 
	{

	newM_TH[,i]=(marker[,sel_th[i,1]]*marker[,sel_th[i,2]])*marker[,sel_th[i,3]]

	}

# Take the transpose of the pseudo marker 
eM=t(newM_TH) # should be p x n
mark=dim(eM)[1]
line=dim(eM)[2]

#Define the residuals(after correcting main and two-way epistasis) as the phenotype for the epistasis search
y=as.vector(err_th)

#The following steps are the same as in  "Bayesian_main_effec.R"

#Create vector 'A' to store the mean of the phenotype
A=vector()
#Define a vector 'va' to store residual
va=vector()

t=1
#Assign the phenotype mean to 'A'
A[1]=mean(y)

#Matrices to store the current and proposed values of marker effects respectively
B=matrix(0,nrow=mark,ncol=2) #marker x 2 order 
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

for (s in 1:ite ) 
	{
	#update A(mean)
	mA=mean(y-t(B[,t])%*%eM)
	vA=va[t]/line
	A[t+1]=mA+sqrt(vA)*rnorm(1)

	#update B (marker effect)
	B[,t+1]=B[,t]
	Bmxj=t(B[,t+1])%*%eM

	for (j in 1:mark ) 
		{	
		Bmxj=Bmxj-B[j,t+1]*eM[j,] # to keep k not j
		e=sum(eM[j,]*(y-A[t+1]-Bmxj))

		d1=sum(eM[j,]^2) #Sig^2_ij

		mB=vb[j,t]*e/(vb[j,t]*d1+va[t])

		vB=vb[j,t]*va[t]/(vb[j,t]*d1+va[t])


		B[j,t+1]=mB+sqrt(vB)*rnorm(1)
		Bmxj=Bmxj+B[j,t+1]*eM[j,]
		}
# update va the residual error
	c2=rchisq(1,line)
	e2=sum((y-A[t+1]-t(B[,t+1])%*%eM)^2)
	va[t+1]=e2/c2
#update vb
	vb[,t+1]=(B[,t+1]^2)/rchisq(mark,1)
#store data
	fA[s]=A[t+1]
	fva[s]=va[t+1]
	if(s>=th)
		{
	res=res+B[,t+1]
		}
#Save the marker effects after the burn-in period	
	A[t]=A[t+1]
	B[,t]=B[,t+1]
	va[t]=va[t+1]
	vb[,t]=vb[,t+1]
	}
#Vector m_epi_eff contains the final results(as the average of all MCMC chains after the burn-in period)
m_epi_eff=res/(ite-th)

############################################################
