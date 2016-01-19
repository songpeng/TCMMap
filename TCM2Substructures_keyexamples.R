# Get the PubChem Substructures of 800 compounds from TCM.
# Songpeng Zu
# 2015-04-07

#-- Load package
require(ChemmineR)
require(xlsx)

#-- Load TCM data.
tcmcid <- read.xlsx("812种化学成分和靶标集合.xlsx",1,colClasses="numeric",
                    colIndex=1)

#-- Function to extract the substructures.
writecomsub <- function(cids,writefilenm){
    compounds <- getIds(cids)
    cid(compounds) <- sdfid(compounds)
    fpset3 <- fp2bit(compounds)
    data <- matrix(fpset3@fpma,nrow=length(cids),byrow=FALSE)
    rownames(data) <- unlist(dimnames(fpset3@fpma)[1])
    write.table(data,file=writefilenm,quote=FALSE,sep="\t",row.names = TRUE,
                col.names = FALSE)
}

#-- Get the subs of compounds from TCM by mannually.
compounds <- getIds(unlist(tcmcid$CID))
cid(compounds) <- sdfid(compounds)
fpset3 <- fp2bit(compounds)

#-- Write the substructures.
# The matrix is saved as column-wise.
tcmdata <- matrix(fpset3@fpma,nrow=length(tcmcid$CID),byrow=FALSE)
rownames(tcmdata) <- unlist(dimnames(fpset3@fpma)[1])
write.table(tcmdata,file="tcm2subs800_pubchem.txt",quote=FALSE, sep="\t",
            row.names = TRUE,col.names = FALSE)

#-- Load Drugs data.
drugcid <- read.table("ANdrugsfromBaiMing.txt",colClasses="numeric",header=FALSE)
writecomsub(drugcid$V1,"ANdrugs2sub.txt")
