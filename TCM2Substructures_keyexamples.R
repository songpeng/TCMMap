# Get the PubChem Substructures of 800 compounds from TCM.
# Songpeng Zu
# 2015-04-07

#-- Load package
require(ChemmineR)
require(xlsx)
#-- Load TCM and drugs data.
tcmcid <- read.xlsx("812种化学成分和靶标集合.xlsx",1,colClasses="numeric",
                    colIndex=1)

#-- Function to extract the substructures.
compounds <- getIds(unlist(tcmcid$CID))
cid(compounds) <- sdfid(compounds)
fpset3 <- fp2bit(compounds)

#-- Write the substructures.
# The matrix is saved as column-wise.
tcmdata <- matrix(fpset3@fpma,nrow=length(tcmcid$CID),byrow=FALSE)
rownames(tcmdata) <- unlist(dimnames(fpset3@fpma)[1])
