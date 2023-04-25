# Volume & Surface Area Calculations for PB2017
# Arina Favilla
# Written on 30-Apr-2022
# Edited on 
  # 1-May-2022 - Work in progress
  #2-May-2022 - Finalized code

# Load Morph data ---------------------------------------------------------

dat <- read.csv("DiseaseAtSea_AllMorphs_PB2017.csv")
n=length(unique(dat$TOPPID))

# Calculate section Lengths from cumulative lengths
dat$Tail_Ank <- dat$Length_Ank
dat$Ank_Pel <- dat$Length_Pel - dat$Length_Ank
dat$Pel_Umb <- dat$Length_Umb - dat$Length_Pel
dat$Umb_Mid <- dat$Length_Mid - dat$Length_Umb
dat$Mid_Ster <- dat$Length_Ster - dat$Length_Mid 
dat$Ster_Ax <- dat$Length_Ax - dat$Length_Ster 
dat$Ax_Neck <- dat$Length_Neck - dat$Length_Ax  
dat$Neck_Ears <- dat$Length_Ears - dat$Length_Neck 
dat$Ears_Nose <- dat$Curve - dat$Length_Ears 

# Set up dataframes ---------------------------------------------------------

SA_CircTC <- as.data.frame(matrix(nrow=n*2,ncol=16))
colnames(SA_CircTC)<-c("TOPPID","SealID","Procedure","Mass","Standard","Curve","Method","Tail_Ank","Ank_Pel","Pel_Umb","Umb_Mid","Mid_Ster","Ster_Ax","Ax_Neck","Neck_Ears","Ears_Nose")
SA_CircTC$TOPPID <- rep(dat$TOPPID)
SA_CircTC$SealID <- rep(dat$SealID)
SA_CircTC$Procedure <- rep(dat$Procedure)
SA_CircTC$Mass <- rep(dat$Mass)
SA_CircTC$Standard <- rep(dat$Standard)
SA_CircTC$Curve <- rep(dat$Curve)
SA_CircTC$Method <- rep("CircTC",times=n*2)

SA_EllipTC <- SA_CircTC
SA_EllipTC$Method <- rep("EllipTC",times=n*2)

V_CircTC <- SA_CircTC
V_EllipTC <- SA_EllipTC


# Functions for straight lengths of cones ---------------------------------

# Determine straight lengths (i.e., height) of each cone
sLCone = function(r,s){
  sqrt(s^2-r^2)
}
sLTruncCone = function(R,r,s){
  sqrt(s^2-(R-r)^2)
}
sLECone = function(w,h,s){
  mean(c(sqrt(s^2-w^2),sqrt(s^2-h^2)))
}
SLETruncCone = function(W,H,w,h,s){
  mean(c(sqrt(s^2-(W-w)^2),sqrt(s^2-(H-h)^2)))
}


# Functions for V & SA of Circular Cones & Cylinder ----------------------------------

# Volume and Lateral SA of full cone
VolCone = function(r,L){
  (1/3)*pi*L*r^2
}
LatSACone = function(r,s){
  pi*r*s
}

# Volume and Lateral SA of truncated cone
VolTruncCone = function(R,r,L){
  (1/3)*pi*L*(R^2+R*r+r^2)
}
LatSATruncCone = function(R,r,L){
  pi*(R+r)*sqrt((R-r)^2+L^2)
}

# Volume and Lateral SA of cylinder
VolCyl = function(r,s){
  pi*s*r^2
}
LatSACyl = function(r,s){
  2*pi*r*s
}


# Functions for V & SA of Elliptical Cones --------------------------------

library(pracma) # surface area requires complete elliptic integral of 2nd kind

# Volume and Lateral SA of elliptical full cone
VolECone = function(w,h,L){
  (1/3)*pi*L*w*h
}
LatSAECone = function(w,h,L){
  ke = ellipke(sqrt((1-(h^2/w^2))/(1+(h^2/L^2))))
  2*w*sqrt(h^2+L^2)*ke$e
}

# Volume and Lateral SA of elliptical truncated cone
VolETruncCone = function(W,H,w,h,L){
  (1/6)*pi*L*(2*w*h+W*h+H*w+2*W*H)
}
LatSAETruncCone = function(W,H,w,h,L){ # Assumes W > H 
  ke = ellipke(sqrt((1-((H^2)/(W^2)))/(1+((H^2)*(W-w)^2)/((L^2)*(W^2)))))
  (2*W*sqrt(H^2+((L*W)/(W-w))^2)-2*w*sqrt(((H^2)*(w^2))/(W^2)+((L*w)/(W-w))^2))*ke$e
}

# Volume and Lateral SA of elliptical cylinder
VolECyl = function(w,h,s){
  pi*w*h*s
}
LatSAECyl = function(w,h,s){
  ke = ellipke(sqrt(1-(h/w)^2))
  4*w*s*ke$e
}


# Calculate SA & V --------------------------------------------------------

Slant = dat[,39:47]
Girth = dat[,15:22]
Radii = Girth/(2*pi)
Height = dat[,23:30]
Width = dat[,31:38]


# Circular truncated cone for sections: Ank_Pel, Pel_Umb, Umb_Mid, Mid_Ster, Ster_Ax, Ax_Neck, Neck_Ears
# Circular cone for sections: Tail_Ank, Ears_Nose
for(j in 1:16){
  # Cone for Tail_Ank
  V_CircTC[j,1+7] = VolCone(Radii[j,1],sLCone(Radii[j,1],Slant[j,1])) 
  SA_CircTC[j,1+7] = LatSACone(Radii[j,1],Slant[j,1])
  # Truncated Cone Ank_Pel, Pel_Umb, Umb_Mid, Mid_Ster, Ster_Ax, Ax_Neck, Neck_Ears
  for(i in 1:7){
    if(Girth[j,i+1]>Girth[j,i]){ # checks which base has a larger radius
      R = Radii[j,i+1]
      r = Radii[j,i]
    }else {
      R = Radii[j,i]
      r = Radii[j,i+1]
    }
    V_CircTC[j,i+8] = VolTruncCone(R,r,sLTruncCone(r,R,Slant[j,i+1])) 
    SA_CircTC[j,i+8] = LatSATruncCone(R,r,sLTruncCone(r,R,Slant[j,i+1]))
  }
  # Cone for Ears_Nose
  V_CircTC[j,9+7] = VolCone(Radii[j,8],sLCone(Radii[j,8],Slant[j,9])) 
  SA_CircTC[j,9+7] = LatSACone(Radii[j,8],Slant[j,9])
}

  
# Elliptical truncated cone for sections: Ank_Pel, Pel_Umb, Umb_Mid, Mid_Ster, Ster_Ax, Ax_Neck, Neck_Ears
# Elliptical cone for sections: Tail_Ank, Ears_Nose
for(j in 1:16){ # ran this separately so that I could change j when it gets stuck -- best solution for now
  # Cone for Tail_Ank
  if(Height[j,1]>Width[j,1]){
    D = Height[j,1]
    d = Width[j,1]
  }else{
    D = Width[j,1]
    d = Height[j,1]
  }
  R = D/2
  r = d/2
  V_EllipTC[j,1+7] = VolECone(R,r,sLECone(R,r,Slant[j,1]))
  SA_EllipTC[j,1+7] = LatSAECone(R,r,sLECone(R,r,Slant[j,1])) # unsolvable for 2017004Rec(j=8,i=1), 2017005Rec(j=10,i=1), 2017006Dep (j=11,i=1)
}

for(j in 1:16){
  for(i in 1:7){
    if(Girth[j,i+1]>Girth[j,i]){ # checks which base has a larger radius
      if(Height[j,i+1]>Width[j,i+1]){ # checks which slump is larger for larger base
        D1 = Height[j,i+1]
        d1 = Width[j,i+1]
      }else{
        D1 = Width[j,i+1]
        d1 = Height[j,i+1]
      }
      if(Height[j,i]>Width[j,i]){ # checks which slump is larger for smaller base
        D2 = Height[j,i]
        d2 = Width[j,i]
      }else{
        D2 = Width[j,i]
        d2 = Height[j,i]
      }
    }else{
      if(Height[j,i]>Width[j,i]){ # checks which slump is larger for larger base
        D1 = Height[j,i]
        d1 = Width[j,i]
      }else{
        D1 = Width[j,i]
        d1 = Height[j,i]
      }
      if(Height[j,i+1]>Width[j,i+1]){
        D2 = Height[j,i+1]
        d2 = Width[j,i+1]
      }else{
        D2 = Width[j,i+1]
        d2 = Height[j,i+1]
      }
    }
    R1 = D1/2
    r1 = d1/2
    R2 = D2/2
    r2 = r1/2
    V_EllipTC[j,i+8] = VolETruncCone(R1,r1,R2,r2,SLETruncCone(R1,r1,R2,r2,Slant[j,i+1]))
    SA_EllipTC[j,i+8] = LatSAETruncCone(R1,r1,R2,r2,Slant[j,i+1])
    
  }
  # Cone for Ears_Nose
  if(Height[j,8]>Width[j,8]){
    D = Height[j,8]
    d = Width[j,8]
  }else{
    D = Width[j,8]
    d = Height[j,8]
  }
  R = D/2
  r = d/2
  V_EllipTC[j,9+7] = VolECone(R,r,sLECone(R,r,Slant[j,9]))
  SA_EllipTC[j,9+7] = LatSAECone(R,r,sLECone(R,r,Slant[j,9]))
}


# Combine results ---------------------------------------------------------

SA <- rbind(SA_CircTC, SA_EllipTC)
SA <- SA[order(SA$TOPPID),]

V <- rbind(V_CircTC, V_EllipTC)
V <- V[order(V$TOPPID),]

for(i in 1:nrow(SA)){
  SA$Ank_Ears[i] <- sum(SA[i,9:15])
  SA$Ank_Nose[i] <- sum(SA[i,9:16])
  SA$Tail_Nose[i] <- sum(SA[i,8:16])
}

for(i in 1:nrow(V)){
  V$Ank_Ears[i] <- sum(V[i,9:15])
  V$Ank_Nose[i] <- sum(V[i,9:16])
  V$Tail_Nose[i] <- sum(V[i,8:16])
}

write.csv(SA,"PB2017_SA.csv",row.names = FALSE)
write.csv(V,"PB2017_V.csv",row.names = FALSE)
