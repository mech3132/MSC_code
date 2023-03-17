
# Function to get Info on Edges
#---------------------------------
getTreeInfo=function(tr2)
{
  NH=node.depth.edgelength(tr2)
  treeH=max(NH)
  edgeInfo=cbind(tr2$edge,treeH-NH[tr2$edge[,1]],treeH-NH[tr2$edge[,2]])
  row.names(edgeInfo)=as.character(1:dim(edgeInfo)[1])
  return(edgeInfo)
}

# Function to get all branches and their descendant node at a given slice
#------------------------------------------------------------------------
GetBranchNode=function(slice,edgeInfo) {
  a=(edgeInfo[,3]>=slice)*(edgeInfo[,4]<=slice)
  return(edgeInfo[a==1,])
}

# Functions to fill the branch*sites matrices 
# (sites can be environmental sites or hosts for microbial communities)
#----------------------------------------------------------------------
vectorPres=function(node,tree)
{
  pres=rep(0,length(tree$tip.label))
  names(pres)=tree$tip.label
  pres[clade.members(node,tree,tip.label=T)]=1 
  return(pres)
}

vectorPresBigMat=function(nodes,tree,mat)
{
  for (i in nodes)
  {
    #print(100*match(i,nodes)/length(nodes))
    dd=clade.members(i,tree,tip.label=T)
    mat[as.character(i),dd]=1 
  }  
  
  return(mat)
}

# Function 'GetBranchOcc': computes a Branch*Sites matrix and directly save it in a file (not in Renv)
# Branch*sites matrices are also frequently named 'OTU tables' in the field of microbiology
#--------------------------------------------------------------------------------------------------------

#   slice: the age of the desired slice
#   tree: the community phylogenetic tree (must be ultrametric)
#   sitesp: the site * species matrice. Species names must match the tip names in the phylogenetic tree
#            In microbiology, sitesp is equivalent to the OTU table determining the distribution of unique 16S sequences across samples (100% similarity OTUs) 
#   pathtoSaveBranchOcc: directory where Branch*sites matrices are stored
#   bigmatrix=F: if site*species is a VERY big matrix, it is highly recommended to set bigmatrix=T (use of the bigmemory package)

GetBranchOcc=function(slice,tree,sitesp,pathtoSaveBranchOcc,bigmatrix=F)
{
  edgeInfo=getTreeInfo(tree)
  
  if (slice==0)
  {OTUmat=sitesp}
  else 
  {
    brCorres=GetBranchNode(slice=slice,edgeInfo=edgeInfo)
    NodeDes=brCorres[,2]    
    #print("getting OTUS host matrix")
    OTUnames=as.character(brCorres[,2])
    if (bigmatrix==F) {branchPAtotal=matrix(0,ncol=length(tree$tip.label),nrow=length(NodeDes),dimnames=list(OTUnames,tree$tip.label))}
    else if (bigmatrix==T) {branchPAtotal=filebacked.big.matrix(ncol=length(tree$tip.label),nrow=length(NodeDes),dimnames=list(OTUnames,tree$tip.label),backingfile = paste("branchPAtotal",slice,sep=""),backingpath =paste(pathtoSaveBranchOcc), descriptorfile = paste("branchPAtotalsauv",slice,sep=""))}   
      
      mat1=vectorPresBigMat(node=NodeDes,tree=tree,mat=branchPAtotal)
      mat1=mat1[,colnames(sitesp)]
      mat2=as.matrix(mat1)
      mat2=t(mat2)
      occM=as.matrix(sitesp)
      OTUmat=(occM)%*%(mat2)
  } 
      save(OTUmat,file=paste(pathtoSaveBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
}    


# Function 'GetBranchOcc': computes a Branch*Sites matrix and directly save it in a file (not in Renv)
# Branch*sites matrices are also frequently named 'OTU tables' in the field of microbiology
#--------------------------------------------------------------------------------------------------------

#   slice: the age of the desired slice (reverse from above; larger # is deeper from tip)
#   tree: the community phylogenetic tree (NON-ultrametric)
#   sitesp: the site * species matrice. Species names must match the tip names in the phylogenetic tree
#            In microbiology, sitesp is equivalent to the OTU table determining the distribution of unique 16S sequences across samples (100% similarity OTUs) 
#   pathtoSaveBranchOcc: directory where Branch*sites matrices are stored
#   bigmatrix=F: if site*species is a VERY big matrix, it is highly recommended to set bigmatrix=T (use of the bigmemory package)

GetBranchOcc_fromtips=function(slice,tree,sitesp,pathtoSaveBranchOcc,bigmatrix=F){
  # edgeInfo=getTreeInfo(tree)
  # bigmatrix=F
  # intree = bacttree_raw
  # sitesp = t(otu)
  # slice=slices[3]
  # intree = testtree
  # sitesp = t(testotu)
  # slice=slices[6]
  # 
  if (slice==0)
  {OTUmat=sitesp}
  else 
  {
    collapsed <- collapse_tree_at_resolution(tree, resolution = slice, rename_collapsed_nodes = TRUE, shorten = FALSE)
    # plot(intree)
    # plot(collapsed$tree)
    if ( length(collapsed$tree$tip.label)==1 ) {
      OTUmat = data.frame("A"=rowSums(sitesp))
    } else {
      if ( length(collapsed$collapsed_nodes)>1 ) {
        # oldNamesCollapse <- lapply(get_tips_for_mrcas(intree, collapsed$collapsed_nodes), function(x) {intree$tip.label[seq(min(x), max(x))]})
        # collapseIntoNames <- intree$tip.label[collapsed$farthest_tips] # names from ORIGINAL tree
        # names(oldNamesCollapse) <- collapseIntoNames
        newNames <- collapsed$tree$tip.label
        
        # brCorres=GetBranchNode(slice=slice,edgeInfo=edgeInfo)
        # NodeDes=brCorres[,2]    
        #print("getting OTUS host matrix")
        # OTUnames=as.character(brCorres[,2])
        if (bigmatrix==F) {branchPAtotal=matrix(0,nrow=length(tree$tip.label),ncol=length(newNames),  dimnames=list(tree$tip.label, newNames))}
        else if (bigmatrix==T) {branchPAtotal=filebacked.big.matrix(nrow=length(tree$tip.label),ncol=length(newNames),  dimnames=list(tree$tip.label, newNames),backingfile = paste("branchPAtotal",slice,sep=""),backingpath =paste(pathtoSaveBranchOcc), descriptorfile = paste("branchPAtotalsauv",slice,sep=""))}   
        nodesToCollapse <- collapsed$collapsed_nodes+length(tree$tip.label)
        names(nodesToCollapse) <- tree$tip.label[collapsed$farthest_tips]
        # Fill matrix
        mat1 = getNewMat(nodesToCollapse = nodesToCollapse,intree = tree, mat=branchPAtotal)
        # mat1=vectorPresBigMat(node=NodeDes,tree=tree,mat=branchPAtotal)
        mat1=mat1[colnames(sitesp),]
        mat2=as.matrix(mat1)
        # mat2=t(mat2)
        occM=as.matrix(sitesp)
        OTUmat=(occM)%*%(mat2)
        # sum(OTUmat)
        # sum(occM)
      } else {
        OTUmat=sitesp
      }
    }
    }
    # plot(collapsed$tree)
    # plot(intree)
    
  # View(OTUmat)
  save(OTUmat,file=paste(pathtoSaveBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
  # return(list(mat=OTUmat, tree=collapsed$tree))
}    

# Function to fill in new matrix given old and new tips
getNewMat <- function(nodesToCollapse,intree,mat=branchPAtotal){
  for ( n in colnames(mat)) {
    if ( n %in% names(nodesToCollapse) ) {
      currTips <- clade.members(nodesToCollapse[[n]], phy=intree, tip.labels = TRUE)
      mat[currTips, n] <- 1
    } else {
      mat[n,n] <- 1
    }
  }
  return(mat)
}

# Function to compute raw Jaccard and Sorensen Beta diversities 
#----------------------------------------------------------
getBeta=function(mat,ab=F)
{
  if (ab==T)
  {   
    h=bray.part(data.frame(mat))
    bctu=as.matrix(h[[1]])
    bcne=as.matrix(h[[2]])
    bc=as.matrix(h[[3]])
    res=abind(bctu,bcne,bc,along=0)
    dimnames(res)[[1]]=c("bctu","bcne","bc")
  }
  if (ab==F)
  {
    h=beta.pair(data.frame(mat), index.family="jaccard")
    hh=beta.pair(data.frame(mat), index.family="sorensen")
    jtu=as.matrix(h[[1]])
    jne=as.matrix(h[[2]])
    jac=as.matrix(h[[3]])
    stu=as.matrix(hh[[1]])
    sne=as.matrix(hh[[2]])
    sor=as.matrix(hh[[3]])  
    res=abind(jtu,jne,jac,stu,sne,sor,along=0)
    dimnames(res)[[1]]=c("jtu","jne","jac","stu","sne","sor")  
  }
  return(res)
}

# Function 'GetBetaDiv' to compute Beta-diversities from site*species matrices and to directly save the matrix of beta-diversities 
#---------------------------------------------------------------------------------------------------------------------------------

# INPUT VARIABLES
#   slice: the age of the desired slice
#   pathtoGetBranchOcc: the input file where the branch*sites matrix is saved (result of the function GetBranchOcc)
#   pathtoSaveBeta: the output file where the matrices of beta diversities are saved (several beta-diversity metrics are used)
#   bigmatrix=F: if you used bigmatrix=T to create the branch*sites matrix with 'GetBranchOcc' (OTU table), use bigmatrix=T again.

# OUTPUT
# save Beta diversity matrices as a 3D array in 'pathtoSaveBeta' 
# Array = array of sites*sites*beta diversity metrics 

#beta diversity metrics 

#"jtu" : True Turnover component of Jaccard (Presence/Absence)
#"jne" : Nestedness component of Jaccard (Presence/Absence)
#"jac" : Jaccard (Presence/Absence)
#"stu" : True Turnover component of Sorensen (Presence/Absence)
#"sne" : Nestedness component of Sorensen (Presence/Absence)
#"sor" : Sorensen (Presence/Absence)
#"bctu" : True Turnover component of Bray-Curtis (Abundance version of Sorensen)
#"bcne" : Nestedness component of Bray-Curtis (Abundance version of Sorensen)
#"bc" : Bray-Curtis (Abundance version of Sorensen)

GetBetaDiv=function(slice,pathtoGetBranchOcc,pathtoSaveBeta)
{
  # slice=0.44
  # OTUmat = OTUmat
  load(paste(pathtoGetBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
  OTUmatPA=OTUmat
  OTUmatPA[OTUmatPA>0]=1
  Be1=getBeta(OTUmatPA,ab=F)  
  Be=getBeta(OTUmat,ab=T)
  Betaa=abind(Be1,Be,along=1)
  save(Betaa,file=paste(pathtoSaveBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
}   

GetBetaDiv_bcjac=function(slice,doBeta,pathtoGetBranchOcc,pathtoSaveBeta)
{
  # slice=slices[length(slices)]
  # OTUmat = OTUmat
  load(paste(pathtoGetBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
  OTUmatPA=OTUmat
  OTUmatPA[OTUmatPA>0]=1
  print(paste0("Getting beta div ",slice))
  if ( "bc" %in% doBeta ) {
    bc=as.matrix(vegdist(OTUmat, method = "bray"))
  } else if ( "jac" %in% doBeta ) {
    jac=as.matrix(vegdist(OTUmat, method="jaccard"))
  }
  if ( length(doBeta)==2 ) {
    Betaa=abind(bc=bc, jac=jac, along=0)
  } else if ( doBeta[1]=="jac") {
    Betaa=abind(jac=jac, along=0)
  } else if ( doBeta[1]=="bc") {
    Betaa=abind(bc=bc, along=0)
  }
  save(Betaa,file=paste(pathtoSaveBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
}   


WriteOTU=function(slice,pathtoGetBranchOcc,pathtoSaveOTU)
  {
    # slice=0.44
    # OTUmat = OTUmat
    load(paste(pathtoGetBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))
    OTUmatPA=OTUmat
    OTUmatPA[OTUmatPA>0]=1
    # Be1=getBeta(OTUmatPA,ab=F)  
    # Be=getBeta(OTUmat,ab=T)
    # Betaa=abind(Be1,Be,along=1)
    write.table(t(OTUmat) %>% as.data.frame() %>% rownames_to_column(var="#OTUID"),file=paste(pathtoSaveOTU,"OTUTable_counts",slice,".txt",sep=""), row.names=FALSE, quote=FALSE, sep="\t")
    write.table(t(OTUmatPA) %>% as.data.frame() %>% rownames_to_column(var="#OTUID"),file=paste(pathtoSaveOTU,"OTUTable_PA",slice,".txt",sep=""), row.names=FALSE, quote=FALSE, sep="\t")
}   


GetBetaDiv_brayonly=function(slice,pathtoGetBranchOcc,pathtoSaveBeta)
{
  # slice=slices[5]
  # pathtoGetBranchOcc="05b_BDTT/fullSet_raw/branchSitesMatrices/"
  load(paste(pathtoGetBranchOcc,"Branch_Site_matrix_SliceNo",slice,".rdata",sep=""))

  Betaa <- as.matrix(vegdist(OTUmat, method="bray"))
  # Be1=getBeta(OTUmatPA,ab=F)  
  # Be=getBeta(OTUmat,ab=T)
  # Betaa=abind(Be1,Be,along=1)
  save(Betaa,file=paste(pathtoSaveBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
}   

# FUNCTION 'GetCorrelations': computes correlations between beta-diversities and environmental distances at the desired slice
#----------------------------------------------------------------------------------------------------------------------------

# INPUT VARIABLES
#   slice: the age of the desired slice
#   indice: the betadiv metric you want:
#		"jtu": True Turnover component of Jaccard (Presence/Absence)
#		"jne": Nestedness component of Jaccard (Presence/Absence)
#		"jac": Jaccard (Presence/Absence)
#		"stu": True Turnover component of Sorensen (Presence/Absence)
#		"sne": Nestedness component of Sorensen (Presence/Absence)
#		"sor": Sorensen (Presence/Absence)
#		"bctu": True Turnover component of Bray-Curtis (Abundance version of Sorensen)
#		"bcne": Nestedness component of Bray-Curtis (Abundance version of Sorensen)
#		"bc": Bray-Curtis (Abundance version of Sorensen)
#   pathtoGetBeta : the file where Beta-Diversity matrices were previously stored (output of the 'GetBetaDiv' function)
#   EnvDist: matrix of environmental distances. (e.g. geographic or climatic distances between sites, dietary or phylogenetic distances between hosts, etc).
#   TypeofMantel: type of correlation used to run the Mantel test (either "Spearman" or "Pearson")
#   nperm: number of permutations to perform to compute a p-value

# OUTPUT
# A 2*n matrix with R2 coefficients and their associated p-values


GetCorrelations=function(slice,indice="sor",pathtoGetBeta="",EnvDist,TypeofMantel="Spearman",nperm=1000)
{
  print(paste0("Getting corr ",slice))
  # slice=slices[5]
  load(file=paste(pathtoGetBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
  #get same sites
  site=intersect(dimnames(Betaa)[[2]],colnames(EnvDist))
  # Remove sites that are na
  
  naEnv=names(which(is.na(EnvDist[,1]))) ## MYC
  naBeta=names(which(is.na(rowSums(Betaa[indice,,])))) ## MYC
  site <- site[!site%in% c(naEnv, naBeta)]
  Betaa=Betaa[indice,site,site]
  EnvDist=EnvDist[site,site]
  test <- tryCatch(MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm,mrank=T), error = function(e) e) ### MYC
  
  
  if ( class(test)[1]=="simpleError" ) {
    return(c(R2=NA, pval=NA))
  } else {
    if (TypeofMantel=="Spearman") {  
      multiMant_SE=test
      # multiMant_SE=MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm,mrank=T)
      }
    if (TypeofMantel=="Pearson") {  multiMant_SE=MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm)}
    
  }
    return(multiMant_SE$r.squared)
}

# ## not finished
# GetCorrelations_allDat=function(slice,indice="bc",pathtoGetBeta="",EnvDistList,TypeofMantel="Spearman",nperm=1000)
# {
#   print(paste0("Getting corr ",slice))
#   # slice=slices[5]
#   load(file=paste(pathtoGetBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
#   Bdist <- as.matrix(Betaa[indice,,])
#   # Filter out nas from Bdist
#   naBeta=names(which(is.na(rowSums(Bdist)))) ## MYC
#   
#   #get same sites
#   allnames <- colnames(Bdist)[!colnames(Bdist)%in%naBeta]
#   for ( p in names(EnvDistList) ) {
#     tempEnv <- as.matrix(EnvDistList[[p]])
#     removeNames <- names(which(is.na(tempEnv[,1])))
#     intersectNames <- intersect(colnames(BDist), colnames(tempEnv))
#     intersectNames <- intersectNames[!intersectNames%in%removeNames]
#     allnames <- c(allnames, intersectNames)
#   }
#   filtSites <- names(which(table(allnames)==(length(EnvDistList)+1)))
# 
#   # naEnv=names(which(is.na(EnvDist[,1]))) ## MYC
#   filtSites <- filtSites[!filtSites%in% c(naEnv, naBeta)]
#   Bdist_filt=as.dist(Bdist[filtSites,filtSites])
#   
#   for ( p in names(EnvDistList) ) {
#     tempEnv <- as.matrix(EnvDistList[[p]])
#     assign(paste0(p,"_filt"), as.dist(tempEnv[filtSites,filtSites]))
#   }
#   frml <- paste0("Bdist_filt ~ ",paste0(paste0(names(EnvDistList),"_filt"),collapse=" + "))
#   test <- tryCatch(MRM(formula(frml),nperm = nperm,mrank=T), error = function(e) e) ### MYC
#   
#   if ( class(test)[1]=="simpleError" ) {
#     return(c(R2=NA, pval=NA))
#   } else {
#     if (TypeofMantel=="Spearman") {  
#       multiMant_SE=test
#       # multiMant_SE=MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm,mrank=T)
#     }
#     if (TypeofMantel=="Pearson") {  multiMant_SE=MRM(formula(frml),nperm = nperm,mrank=T)}
#   }
#   #####return(multiMant_SE$r.squared)
# }


#### MY OWN 
# 
# force.ultrametric<-function(tree,method=c("nnls","extend")){
#   method<-method[1]
#   if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
#                                      rooted=TRUE,trace=0)
#   else if(method=="extend"){
#     h<-diag(vcv(tree))
#     d<-max(h)-h
#     ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
#                y=tree$edge[,2])
#     tree$edge.length[ii]<-tree$edge.length[ii]+d
#   } else 
#     cat("method not recognized: returning input tree\n\n")
#   tree
# }
# # 
# runBDTT <- function(slices, intree, sitesp, output, setName, metrics=c("bc","bctu","bcne","sor","jac"), predictors=c("phyloDist","geoDist","ecoDist","timeDist","climDist"), sliceMethod=c("bdtt","chen")) {
#   if ( sliceMethod == "bdtt" ) {branchFUNC=GetBranchOcc}
#   if ( sliceMethod == "chen" ) {branchFUNC=GetBranchOcc_fromtips}
#   if ( !sliceMethod %in% c("bdtt","chen")) { 
#     print("Choose bdtt or chen for slice method")
#     stop()
#   }
#   
#   if ( ! file.exists(paste(output,"Beta/BetaDiv_BetaDivSliceNo",slices[1],".rdata",sep=""))  ) {
#     dir.create(output)
#     dir.create(paste0(output, "branchSitesMatrices/"))
#     dir.create(paste0(output, "Beta/"))
#     print("Slicin' and dicin'; making some betas")
#     lapply(slices,branchFUNC,tree=intree,sitesp=sitesp,pathtoSaveBranchOcc=paste0(output,"branchSitesMatrices/"),bigmatrix=T)
#     lapply(slices,GetBetaDiv,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices/"),pathtoSaveBeta=paste0(output,"Beta/"))
#   } else {
#     print("Already ran Branchocc and Betadiv")
#   }
#   
#   if ( file.exists(paste(output,"Correlations/allDat_",setName,".RData", sep="")) ) {
#     print("Loading correlations")
#     load(paste(output,"Correlations/allDat_",setName,".RData", sep=""))
#   } else {
#     dir.create(paste0(output,"Correlations"))
#     print("Conducting Correlations")
#     allDat <- data.frame()
#     for ( m in metrics ) {
#       for ( p in predictors) {
#         if ( ! file.exists(paste0(output,"Correlations/Cors",m,p,".RData")) ) {
#           Cors=sapply(slices,GetCorrelations,indice=m, EnvDist=get(p),pathtoGetBeta=paste0(output,"Beta/"),nperm=1000)
#           save(Cors, file=paste0(output,"Correlations/Cors",m,p,".RData"))
#         } else {
#           load(paste0(output,"Correlations/Cors",m,p,".RData"))
#         }
#         colnames(Cors)=slices
#         
#         newDat <- t(Cors) %>% as.data.frame() %>% rownames_to_column(var="slice")%>% 
#           mutate(type=p, metric=m, slice=as.numeric(slice))
#         if (nrow(allDat)==0) {
#           allDat <- newDat
#         } else {
#           allDat <- rbind(allDat, newDat)
#         }
#       }
#       
#     }
#     save(allDat, file=paste(output,"Correlations/allDat_",setName,".RData", sep=""))
#   }
#   return(allDat)
# }

# 
# runCor_useExitingBDTT <- function(slices, sites, subsetName, output, setName, metrics=c("bc","bctu","bcne","sor","jac"), predictors=c("phyloDist","geoDist","ecoDist","timeDist","climDist","timeDist")) {
#   if ( file.exists(paste(output,"Correlations_",subsetName,"/allDat_",setName,".RData", sep="")) ) {
#     print("Loading correlations")
#     load(paste(output,"Correlations_",subsetName,"/allDat_",setName,".RData", sep=""))
#   } else {
#     dir.create(paste0(output,"Correlations_",subsetName))
#     print("Conducting Correlations")
#     allDat <- data.frame()
#     for ( m in metrics ) {
#       for ( p in predictors) {
#         if ( ! file.exists(paste0(output,"Correlations_",subsetName,"/Cors",m,p,".RData")) ) {
#           filtEnv <- get(p)[sites,sites]
#           Cors=sapply(slices,GetCorrelations,indice=m, EnvDist=filtEnv,pathtoGetBeta=paste0(output,"Beta/"),nperm=1000)
#           save(Cors, file=paste0(output,"Correlations_",subsetName,"/Cors",m,p,".RData"))
#         } else {
#           load(paste0(output,"Correlations_",subsetName,"/Cors",m,p,".RData"))
#         }
#         colnames(Cors)=slices
#         
#         newDat <- t(Cors) %>% as.data.frame() %>% rownames_to_column(var="slice")%>% 
#           mutate(type=p, metric=m, slice=as.numeric(slice))
#         if (nrow(allDat)==0) {
#           allDat <- newDat
#         } else {
#           allDat <- rbind(allDat, newDat)
#         }
#       }
#       
#     }
#     save(allDat, file=paste(output,"Correlations_",subsetName,"/allDat_",setName,".RData", sep=""))
#   }
#   return(allDat)
# }


#### FINAL VERSIONS BELOW #####

getDM <- function(output,slice,indice,sites) {
  pathtoGetBeta=paste0(output,"Beta/")
  load(file=paste(pathtoGetBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
  #get same sites
  bextract=Betaa[indice,sites,sites]
  return(bextract)
}

runBDTTOnly <- function(slices, intree, sitesp, output, sliceMethod=c("bdtt","chen"), doBeta=c()) {
  if ( sliceMethod == "bdtt" ) {branchFUNC=GetBranchOcc}
  if ( sliceMethod == "chen" ) {branchFUNC=GetBranchOcc_fromtips}
  if ( !sliceMethod %in% c("bdtt","chen")) { 
    print("Choose bdtt or chen for slice method")
    stop()
  }
  
  if ( ! file.exists(paste(output,"branchSitesMatrices/","Branch_Site_matrix_SliceNo",slices[length(slices)],".rdata",sep=""))  ) {
    dir.create(output)
    dir.create(paste0(output, "branchSitesMatrices/"))
    dir.create(paste0(output, "Beta/"))
    
    print("Slicin' and dicin'")
    lapply(slices,branchFUNC,tree=intree,sitesp=sitesp,pathtoSaveBranchOcc=paste0(output,"branchSitesMatrices/"),bigmatrix=T)
  } else {
    print("Already ran Branchocc")
  }
  if (length(doBeta)>0  ){
    print("Makin' some betas")
    # lapply(slices,GetBetaDiv,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices/"),pathtoSaveBeta=paste0(output,"Beta/"))
    lapply(slices,GetBetaDiv_bcjac,doBeta=doBeta,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices/"),pathtoSaveBeta=paste0(output,"Beta/"))
  } else {
    dir.create(paste0(output, "OTU_tables/"))
    print("Writin' OTUs") 
    
    if ( ! file.exists(paste(output,"OTU_tables/","OTUTable_counts",slices[length(slices)],".txt",sep=""))) {
      lapply(slices, WriteOTU,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices/"), pathtoSaveOTU=paste0(output,"OTU_tables/"))
    } else {
      print("Tables already printed")
    }
    
  }
  
  print("Finished saving betas")
}


runBDTTOnly_withhighlow <- function(slices, intree, sitesp, output, sliceMethod=c("bdtt","chen"), doBeta=c(), levels=c("Full","High","Low")) {
  if ( sliceMethod == "bdtt" ) {branchFUNC=GetBranchOcc}
  if ( sliceMethod == "chen" ) {branchFUNC=GetBranchOcc_fromtips}
  if ( !sliceMethod %in% c("bdtt","chen")) { 
    print("Choose bdtt or chen for slice method")
    stop()
  }
  dir.create(output)
  
  #### Get high, low abund ###
  sitesp_withzero <- sitesp
  sitesp_withzero[sitesp_withzero==0] <- NA
  sitespRanked <- t(apply(sitesp_withzero, 1, function(x) rank(x, na.last = "keep")/length(x[!is.na(x)])))
  sitespRanked[is.na(sitespRanked)] <- 0
  # dim(sitespRanked)
  for ( otuType in levels) {
    # otuType="High"
    otuTemp <- NULL
    if ( otuType == "Full") {
      otuTemp<- sitesp
      # dim(sitesp)
    } else if (otuType == "Low") {
      sitesp_slice <- sitesp
      sitesp_slice[sitespRanked>0.5] <- 0
      spKeep <- which(colSums(sitesp_slice)>0)
      sitesp_slice <- sitesp_slice[,spKeep]
      colnames(sitesp_slice) <- names(spKeep)
      otuTemp <- sitesp_slice
      # dim(otuTemp)
    } else if (otuType == "High") {
      sitesp_slice <- sitesp
      sitesp_slice[sitespRanked<=0.5] <- 0
      spKeep <- which(colSums(sitesp_slice)>0)
      sitesp_slice <- sitesp_slice[,spKeep]
      colnames(sitesp_slice) <- names(spKeep)
      otuTemp <- sitesp_slice
      # dim(otuTemp)
    }
    # dim(otuTemp)
    bacttree_filt <- keep.tip(intree, colnames(otuTemp))
    # print(slices)
  ######### BRANCH MATRICES
    if ( ! file.exists(paste(output,"branchSitesMatrices_/",otuType,"/Branch_Site_matrix_SliceNo",slices[length(slices)],".rdata",sep=""))  ) {
      dir.create(paste0(output, "branchSitesMatrices_",otuType,"/"))
      
      print("Slicin' and dicin'")
      lapply(slices,branchFUNC,tree=bacttree_filt,sitesp=otuTemp,pathtoSaveBranchOcc=paste0(output,"branchSitesMatrices_",otuType,"/"),bigmatrix=T)
    } else {
      print("Already ran Branchocc")
    }
    
    ########## BETAS
    
    if (length(doBeta)>0  ){
      print("Makin' some betas")
      dir.create(paste0(output, "Beta_",otuType,"/"))
      
      # lapply(slices,GetBetaDiv,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices/"),pathtoSaveBeta=paste0(output,"Beta/"))
      lapply(slices,GetBetaDiv_bcjac,doBeta=doBeta,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices_",otuType,"/"),pathtoSaveBeta=paste0(output,"Beta_",otuType,"/"))
    } else {
      dir.create(paste0(output, "OTU_tables_",otuType,"/"))
      print("Writin' OTUs") 
      
      if ( ! file.exists(paste(output,"OTU_tables_",otuType,"/","OTUTable_counts",slices[length(slices)],".txt",sep=""))) {
        lapply(slices, WriteOTU,pathtoGetBranchOcc=paste0(output,"branchSitesMatrices_",otuType,"/"), pathtoSaveOTU=paste0(output,"OTU_tables_",otuType,"/"))
      } else {
        print("Tables already printed")
      }
      
    }
  }

  
  print("Finished saving betas")
}





runCorOnly <- function(slices, sites, subsetName, output, metrics=c("bc","jac"), predictors=c("phyloDist","geoDist","ecoDist","timeDist","climDist","timeDist")
                       , phyloDist=NULL, geoDist=NULL, ecoDist=NULL, timeDist=NULL, climDist=NULL) {

  
  if ( file.exists(paste(output,"Correlations_",subsetName,"/allDat.RData", sep="")) ) {
    print("Loading correlations")
    load(paste(output,"Correlations_",subsetName,"/allDat.RData", sep=""))
  } else {
    # Check intersection of all site names
    dir.create(paste0(output,"Correlations_",subsetName))
    print("Conducting Correlations")
    allDat <- data.frame()
    for ( m in metrics ) {
      # m="jac"
      for ( p in predictors) {
        # p="phyloDist"
        if ( ! file.exists(paste0(output,"Correlations_",subsetName,"/Cors",m,p,".RData")) ) {
          allEnv <- as.matrix(get(p))
          filtsites <- intersect(colnames(allEnv), sites) 
          # Warn if sites are not in both
          if ( length(filtsites) != length(sites) ) {
            print("WARNING: missing sites in data that are not in site list or loaded betas")
          } 
          filtEnv <- allEnv[filtsites,filtsites]
          Cors=sapply(slices,GetCorrelations,indice=m, EnvDist=filtEnv,pathtoGetBeta=paste0(output,"Beta/"),nperm=1000)
          save(Cors, file=paste0(output,"Correlations_",subsetName,"/Cors",m,p,".RData"))
        } else {
          load(paste0(output,"Correlations_",subsetName,"/Cors",m,p,".RData"))
        }
        colnames(Cors)=slices
        
        newDat <- t(Cors) %>% as.data.frame() %>% rownames_to_column(var="slice")%>% 
          mutate(type=p, metric=m, slice=as.numeric(slice))
        if (nrow(allDat)==0) {
          allDat <- newDat
        } else {
          allDat <- rbind(allDat, newDat)
        }
      }
      
    }
    save(allDat, file=paste(output,"Correlations_",subsetName,"/allDat.RData", sep=""))
  }
  return(allDat)
}


runCorOnly_withPrevSlice <- function(slices, sites, subsetName, output, metrics=c("bc","jac"), levels=c("Full","High","Low"), predictors=c("phyloDist","geoDist","ecoDist","timeDist","climDist","timeDist","studyDist")
                       , phyloDist=NULL, geoDist=NULL, ecoDist=NULL, timeDist=NULL, climDist=NULL, studyDist=NULL, shuffle=FALSE) {
  
  if ( shuffle ) {
    randomOrder <- sample(1:length(sites))
  }
  if ( file.exists(paste(output,"Correlations_",subsetName,"/allDat.RData", sep="")) ) {
    print("Loading correlations")
    load(paste(output,"Correlations_",subsetName,"/allDat.RData", sep=""))
  } else {
    # Check intersection of all site names
    dir.create(paste0(output,"Correlations_",subsetName))
    print("Conducting Correlations")
    allDat <- data.frame()
    for ( m in metrics ) {
      # m="bc"
      for ( p in predictors) {
        # p="timeDist"
        print(c("DOING",p))
        allEnv <- as.matrix(get(p))
        filtsites <- intersect(colnames(allEnv), sites) 
        
        # Warn if sites are not in both
        if ( length(filtsites) != length(sites) ) {
          print("WARNING: missing sites in data that are not in site list or loaded betas")
        } 
        filtEnv <- allEnv[filtsites,filtsites]
        # Shuffle sites of required
        if ( shuffle ) {
          print("SHUFFLING")
          randomOrder <- randomOrder[!randomOrder>length(filtsites)]
          rownames(filtEnv) <- filtsites[randomOrder]
          colnames(filtEnv) <- filtsites[randomOrder]
        }
        ## CHECK FULL
        for ( l in levels ) {
          if ( ! file.exists(paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_Full.RData")) ) {
            Cors=sapply(slices,GetCorrelations,indice=m, EnvDist=filtEnv,pathtoGetBeta=paste0(output,"Beta_",l,"/"),nperm=1000)
            save(Cors, file=paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_",l,".RData"))
          } else {
            load(paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_",l,".RData"))
          }
        }
        # 
        # ## CHECK HIGH
        # if ( ! file.exists(paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_High.RData")) ) {
        #   Cors_high=sapply(slices,GetCorrelations,indice=m, EnvDist=filtEnv,pathtoGetBeta=paste0(output,"Beta_High/"),nperm=1000)
        #   save(Cors_high, file=paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_High.RData"))
        # } else {
        #   load(paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_High.RData"))
        # }
        # ### CHECK LOW
        # if ( ! file.exists(paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_Low.RData")) ) {
        #   Cors_low=sapply(slices,GetCorrelations,indice=m, EnvDist=filtEnv,pathtoGetBeta=paste0(output,"Beta_Low/"),nperm=1000)
        #   save(Cors_low, file=paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_Low.RData"))
        # } else {
        #   load(paste0(output,"Correlations_",subsetName,"/Cors",m,p,"_Low.RData"))
        # }
        
        # } else {
        #   print(paste0("loading ",m," ",p))
        # }
        # colnames(Cors_full)=slices
        # colnames(Cors_low)=slices
        # colnames(Cors_high)=slices
        # 
        # newDat_full <- t(Cors_full) %>% as.data.frame() %>% rownames_to_column(var="slice")%>% 
        #   mutate(type=p, metric=m, slice=as.numeric(slice), group="full")
        # newDat_high <- t(Cors_high) %>% as.data.frame() %>% rownames_to_column(var="slice")%>% 
        #   mutate(type=p, metric=m, slice=as.numeric(slice), group="high")
        # newDat_low <- t(Cors_low) %>% as.data.frame() %>% rownames_to_column(var="slice")%>% 
        #   mutate(type=p, metric=m, slice=as.numeric(slice), group="low")
        # newDat <- rbind(newDat_full, newDat_high, newDat_low)
        if (nrow(allDat)==0) {
          allDat <- Cors
        } else {
          allDat <- rbind(allDat, Cors)
        }
      }
      
    }
    save(allDat, file=paste(output,"Correlations_",subsetName,"/allDat.RData", sep=""))
  }
  return(allDat)
}


runCorOnly_withPrevSlice_allpredictors <- function(slices, sites, subsetName, output, metrics=c("bc","jac"), predictorNamedList, nperm=1000, TypeofMantel="Spearman",levels=c("Full","High","Low")) {
  
  if ( file.exists(paste(output,"Correlations_",subsetName,"/allDat.RData", sep="")) ) {
    print("Loading correlations")
    load(paste(output,"Correlations_",subsetName,"/allDat.RData", sep=""))
  } else {
    # Check intersection of all site names
    dir.create(paste0(output,"Correlations_",subsetName))
    print("Conducting Correlations")
    allDat <- data.frame()
    for ( m in metrics ) {
      # m="bc"
      for ( slice in slices ) {
        # slice=slices[20]
        print(slice)
        for (l in levels) {
          # l="Full"
          if ( !file.exists(paste0(output,"Correlations_",subsetName,"/Cors_",m,l,slice,".RData")) ) {
            
            pathtoGetBeta=paste0(output,"Beta_",l,"/")
            # print(pathtoGetBeta)
            ######## Insert correlation~~~~~~~~~~~~~~~~~~~~~~~~~
            load(file=paste(pathtoGetBeta,"BetaDiv_BetaDivSliceNo",slice,".rdata",sep=""))
            #get same sites
            # REMOVE BETAs
            nonNAnames <- names(which(!is.na(rowSums(Betaa[m,,]))))
            allColNames <- c(nonNAnames, sites)
            for ( p in names(predictorNamedList)) {
              # p="geoDist"
              tempMat <- as.matrix(predictorNamedList[[p]])
              keepNames <- names(which(apply(tempMat, 1, function(x) any(!is.na(x)))))
              # keepNames <- names(which(!is.na(rowSums(tempMat))))
              # print(length(keepNames))
              allColNames <- c(allColNames, keepNames)
            }
            filtsites <- names(which(table(allColNames)==length(predictorNamedList)+2))
            # find common sites between all of these
            # Get function for MRM
            for ( p in names(predictorNamedList)) {
              # p="timeDist"
              # temp <- as.matrix(predictorNamedList[[p]])[filtsites,filtsites]
              assign(paste0(p,"_filt"), as.dist(as.matrix(predictorNamedList[[p]])[filtsites,filtsites]))
            }
            # Warn if sites are not in both
            if ( length(filtsites) != length(sites) ) {
              print("WARNING: missing sites in data that are not in site list or loaded betas")
            } 
            
            # Betaa1=Betaa[indice,filtsites,filtsites]
            # print(c("betaDist",dim(Betaa1)))
            betaDist=as.dist(Betaa[m,filtsites,filtsites])
            # Get environmental formula
            EnvFrml <-paste0("betaDist ~ ", paste(paste0(names(predictorNamedList),"_filt"), collapse = " + "))
            # nperm=10
            test <- tryCatch(MRM(formula(EnvFrml),nperm = nperm,mrank=T), error = function(e) e) ### MYC
            
            if ( class(test)[1]=="simpleError" ) {
              corrresults <- data.frame(predictor=NA, Coef=NA, pval_coef=NA, R2=NA, pval=NA, type=l, metric=m, slice=slice)
            } else {
              if (TypeofMantel=="Spearman") {  
                multiMant_SE=test
                # multiMant_SE=MRM(as.dist(Betaa)~as.dist(EnvDist),nperm = nperm,mrank=T)
              }
              if (TypeofMantel=="Pearson") {  multiMant_SE=MRM(frml,nperm = nperm,mrank=T)}
              corrresults <- as.data.frame(multiMant_SE$coef) %>% rownames_to_column(var="predictor") %>%
                rename(Coef="betaDist", pval_coef=pval) %>% cbind(as.data.frame(t(multiMant_SE$r.squared))) %>%
                mutate(type=l, metric=m, slice=slice) %>% filter(predictor!="Int")
            }
            
            save(corrresults, file=paste0(output,"Correlations_",subsetName,"/Cors_",m,l,slice,".RData"))
            
          } else {
            load(paste0(output,"Correlations_",subsetName,"/Cors_",m,l,slice,".RData"))
          }
          
          
          allDat <- rbind(allDat, corrresults)
          ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        }
        
      }
      
    }
    save(allDat, file=paste(output,"Correlations_",subsetName,"/allDat.RData", sep=""))
  }
  return(allDat)
}

