# BACI power but combining multiple watersheds. 
# The ns_t and ns_c are now vectors with the number of streams in each watershed. Length of the vectors must be the same
#    for example ns_c=(2,3,2) and ns_t=(1,1,1)
# indicates three watersheds with 2/1, 3/1, and 2/1 streams in control vs treatment


# BACI Power function
# 2018-08-09 CJS removed using lmer() and did the computations directly using matrices
#            CJS did not expand the individual observations since quadrats would be averaged in the analysis.
# 2016-06-01 CJS Converted output to a data.frame rather than a list
# 2014-02-26 CJS Added control () to the lmer() call to avoid problems with no sub-samples
# 2013-02-04 CJS Added code to compute the one-sided power for the specified BACI contrast as well
# 2013-12-16 CJS Modified code to deal with changes from new release of lmer()

# This function computes the power for a BACI design with multiple treatment
# and control sites, multiple years measured before and after treatment is applied, and
# sub-sampling at each year-site combination.
#
# The information we need is:
#     Alpha level (usually .05 or .10)
#     Variance components (these were obtained from the an analysis of previous data)
#        sdSite      -> site to site STANDARD DEVIATION
#        sdYear      -> year to year STANDARD DEVIATION
#        sdSiteYear  -> site-year interaction STANDARD DEVIATION
#        sdResid     -> within site sub-sampling STANDARD DEVIATION
#     Sub-sampling sample sizes, i.e. the number of sub-samples taken at each site-period combinations.
#       We allow for different sample sizes in the treatment-time combinations
#        but all sites in that combination must have the same number of samples.
#        It is possible to generalize this; contact me for details.
#        n_TA   -> number of subsamples in Treatment-After  combination
#        n_TB   -> number of subsamples in Treatment-Before combination
#        n_CA   -> number of subsamples in Control-After    combination
#        n_CB   -> number of subsamples in Control-Before   combination
#     Number of sites
#        We allow for a different number of sites for Treatment and Control areaa
#        ns_T     -> number of treatment sites (i.e. downstram of the project)
#        ns_C     -> number of control sites (i.e. upstream of the project)
#     Number of years of monitoring before/after impact
#        ny_B     -> number of years monitoring before project starts
#        ny_A     -> number of years monitoring after  project starts
#     Marginal Means
#        These are used to form the BACI contrast of interest
#        mu_TA, mu_TB, mu_CA, mu_CB (i.e mean of Treatment-After,
#        Treatment-Before, Control-After, Control_before)
#        These are chosen based on the size of impact that may be biologically important
 
# This code was originally created by 
#    Tony Booth
#    Department of Ichthyology and Fisheries Science, 
#    Rhodes University, 
#    PO Box 94, 
#    Grahamstown 6140 SOUTH AFRICA
#    t.booth@ru.ac.za

# The computations are based on the 
#    Stroup, W. W. (1999)
#    Mixed model procedures to assess power, precision, and sample size in the design of experiments.
# paper where "dummy" data is generated and "analyzed" and the resulting F-statistics etc
# provide information needed to compute the power


#-----------------------------------------------

baci.power2 <- function(n_TA,n_TB,n_CA,n_CB,nsw_T,nsw_C,ny_B,ny_A,mu_TA,mu_TB,mu_CA,mu_CB,
                      sdYear, sdSite, sdSiteYear, sdResid, sdWS.BACI=0, alpha=0.05){
   if(length(nsw_T) != length(nsw_C))stop("Number of watersheds (length of ns_T and ns_C) must be equal")
# This computes the power of the baci-fry example. Arguments are defined above
  
   ns_T = sum(nsw_T)
   ns_C = sum(nsw_C)
  
# First generate a "dummy" dataset of the appropriate size with the response variable
# set to the mean as needed. Arguments of the function were defined above.      
  
#  Total number ob observations
   n <- ns_T*ny_A*n_TA+   # treatment sites x years after  x traps/site/year
        ns_T*ny_B*n_TB+   # treatment sites x years before x traps/site/year
        ns_C*ny_A*n_CA+   # control   sites x years after  x traps/site/year
        ns_C*ny_B*n_CB    # control   sites x years before x traps/site/year
   
  testdata <- data.frame(Site=rep(1:(ns_T+ns_C), each=ny_B+ny_A),
                          Year= 1:(ny_B+ny_A),
                          CI  = c(rep("T", ns_T*(ny_B+ny_A)),rep("C", each=ns_C*(ny_B+ny_A))),
                          BA  = c( rep("B",ny_B),rep("A",ny_A)),
                          mu  = c(rep(c(rep(mu_TB, ny_B), rep(mu_TA,ny_A)),ns_T), rep(c(rep(mu_CB, ny_B), rep(mu_CA,ny_A)),ns_C)),
                          nQ  = c(rep(c(rep(n_TB, ny_B),  rep(n_TA,ny_A)),ns_T),  rep(c(rep(n_CB, ny_B),  rep(n_CA, ny_A)),ns_C)),
                          WS  = c(rep( rep(1:length(nsw_T), nsw_T), each=ny_B+ny_A), rep( rep(1:length(nsw_C), nsw_C), each=ny_B+ny_A)),
                          stringsAsFactors=FALSE)


# Now to fit the data and extract the various bits of interest
# We don't actually need the results, but we want the various matrices
  nL <- ns_T+ns_C     # number of sites for the control and treatment
  nT <- ny_A + ny_B   # number of times before and after
  nLT <- nL*nT

#  Get the design matrices for fixed and random effects

   X <- model.matrix(~-1+BA:CI, data=testdata)
   
   Location    <- model.matrix(~0+as.factor(Site), data=testdata)
   Time        <- model.matrix(~0+as.factor(Year), data=testdata)
   LocationTime<- model.matrix(~0+as.factor(Site):as.factor(Year), data=testdata)
   WS.BACI     <- model.matrix(~0+as.factor(WS):as.factor(BA):as.factor(CI), data=testdata)


   # Add the variance components together to generate the variance-covariance matrix
   V <- diag(sdResid^2/testdata$nQ,nrow(testdata),nrow(testdata)) + Time%*%t(Time)*sdYear^2 + Location%*%t(Location)*sdSite^2 + 
       LocationTime%*%t(LocationTime)*sdSiteYear^2 +
       WS.BACI %*% t(WS.BACI)*sdWS.BACI^2
   
# Get fixed effects and fixed effects varcovar matrix
# because the design is balanced, it not necessary to use weighted least squares
  VI <- solve(V)
  beta<- solve(t(X) %*% VI %*% X)   %*%  t(X) %*% VI %*% testdata$mu

# the contrast vector for the BACI effect
  K <- c(1,-1,-1,1)
  baci=mu_TA-mu_TB-mu_CA+mu_CB

#  calculate the non-centrality parameter, and then the power
  ncp <- as.numeric(t(t(K)%*%beta)%*%solve(t(K) %*%   solve(t(X) %*% VI %*% X)%*% K)  %*% (t(K) %*% beta))
  
  dfdenom <- (ny_B-1 + ny_A-1)*(ns_T-1 + ns_C-1) +  # site*time(BACI) term
             (ny_B-1 + ny_A-1)*1 +                  # BA*site(CI) term
             (ns_T-1 + ns_C-1)*1                    # CI*time(BA) term
  if(dfdenom==0){ # simple baci design with
  	 dfdenom = n-1-3;
  }           

  Fcrit <- qf(1-alpha, 1, dfdenom)
  ts.power <- 1 - pf(Fcrit, 1, dfdenom,ncp)

#  Compute the one-sided power, i.e. to detect the change in one direction only
#  Because I don't know which direction of interest, the one-sided power is 
#  computed in both directions.
#
  Tcrit <- qt(1-alpha,dfdenom)
  os.power.a <- 1-pt(Tcrit,dfdenom,sqrt(ncp))
  os.power.b <- pt(-Tcrit,dfdenom,sqrt(ncp))

  return(data.frame(alpha=alpha, 
     sdSite=sdSite, sdYear=sdYear, sdSiteYear=sdSiteYear, sdResid=sdResid,sdWS.BACI=sdWS.BACI,
     n_TA=n_TA, n_TB=n_TB, n_CA=n_CA,n_CB=n_CB,
     ns_T=ns_T, ns_C=ns_C, ny_B=ny_B, ny_A=ny_A, 
     mu_TA=mu_TA, mu_TB=mu_TB, mu_CA=mu_CA, mu_CB=mu_CB, 
     baci=baci,
     dfdenom=dfdenom, ncp=ncp, Fcrit=Fcrit, ts.power=ts.power,
     Tcrit=Tcrit, os.power.a=os.power.a, os.power.b=os.power.b))
}

