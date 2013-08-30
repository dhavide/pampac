        COMMON /CONST/ANU,DT
        common /period/ disc_per	
        common /section/ secxtol,secftol,xsec,pmax,nnwt,nfix,npsfix
        common /fixed/ ifix
        common /gmres_pars/ gmrestol,maxdim,maxit
c  set time step
	dt=0.01d0
c  set discrete period
        disc_per=5
c  set Poincare intersection data
        secxtol=1d-11
	secftol=1d-11
c  NOTE xsec depends on the initial data
	xsec=-4d-2
	pmax=10d0
	nnwt=5
	nfix=1
	ns_dum=2**nn/6
	npsfix=2*(ns_dum+1)+(ns_dum+1)**2+1
c  set fixed forcing component NOTE depends on truncation
        ifix=3291
c  set GMRES parameters
	maxdim=40
	maxit=1
	gmrestol=1d-5

