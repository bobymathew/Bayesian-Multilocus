#Sure Independence Screening (SIS) algorithm for Three-Way Epistasis Effects

#Define the no of three-way pseudo markers to be saved
nth=1000

#the number of two-way pseudo markers 
ns=1000

# Read the phenotype file 
dat=read.table("pheno.txt",header=T)

# Read the marker information 
M=read.table("geno_ord_2331_epi.txt",header=T)


# Define the no of markers
mark=dim(M)[1]

#Define the no of lines with the phenotype information
line=dim(M)[2]



#Define the marker information as a matrix
M=as.matrix(M)

# Define the phenotype as a vector called 'y' 
y=as.vector(dat$phe)


#Read the file with the marker effects for the 2331 markers based on 5 MCMC chains. Each column contains the marker effects(2331 markers) corresponding to a MCMC chain
eff=read.table("all_run_main.txt",header=T)

#Read the file with the epistatic effects based 1000 pseudo markers based on 5 MCMC chains. Each column contains the epistatic effects(1000 markers) corresponding to a MCMC chain
epi=read.table("all_run_epi.txt",header=T)

#Read the most important 1000 two-way interacting marker coordinates
sure=read.table("sure_sel_markers_two_way.txt",header=T)


# Calculate the mean of marker effects from 5 MCMC chains
m_eff=rowMeans(eff)
# Calculate the mean of two-way interaction effects from 5 MCMC chains
m_epi=rowMeans(epi)

# Calculate residual error(E1i)
err=y-(t(M)%*%m_eff)


marker=t(M)
n=dim(marker)[1]
p=dim(marker)[2]

#Create the pseudo marker to correct for the two-way interaction effects.
newM=matrix(0,nrow=n,ncol=ns)

for (i in 1:ns ) 
	{

	newM[,i]=marker[,sel[i,1]]*marker[,sel[i,2]]

	}

#Calculate epistasis residuals (E2i)
err_th=err-(newM%*%m_epi)




#three way sure#########################################################33




#Create a matrix to store the indices of the most three-way correlated markers(1st and 2nd and 3rd column) and the corresponding marginal correlation(4th column of the matrix)
sel_th=matrix(0,nrow=nth,ncol=4)

	for (i in 1:p ) 
	{
		
		  for (j in (i+1):p )
			{	 
				for (l in (j+1):p ) 
					{	

					if(l>p) { break; }
		 		else 
					{

					epiX= (marker[,i]*marker[,j])*marker[,l]

					Nepix=(epiX-mean(epiX))/sd(epiX)
					# compute marginal correlation 
					val = abs(Nepix%*%err_th)
					if (is.nan(val)) {val=0}
					if (val >=sel_th[1,4])
						{
						sel_th[1,4]=val

						sel_th[1,1:3]=c(i,j,l)

						sel_th=sel_th[ order(sel_th[,4],decreasing=F), ]
						}
		  			}
				 }
			}
	}
sel_th=sel_th[ order(sel_th[,3],decreasing=F), ]

sel_th=sel_th[ order(sel_th[,2],decreasing=F), ]
sel_th=sel_th[ order(sel_th[,1],decreasing=F), ]

#Create  pseudo markers for the three-way epistasis search
newM_TH=matrix(0,nrow=n,ncol=nth)

for (i in 1:nth ) 
	{

	newM_TH[,i]=(marker[,sel_th[i,1]]*marker[,sel_th[i,2]])*marker[,sel_th[i,3]]

	}
#Store the selected marker indices in a file called "three_way_sel_cordinates.txt"
############################################################



