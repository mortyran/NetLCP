\name{binaryRegulation}
\alias{binaryRegulation}
\title{Binary Regulatory in Local Heterogeneous network}
\usage{binaryRegulation(transcriptomeList, regulationType, allRegulation)}
\arguments{
  \item{transcriptomeList}{a vector consist of transcriptome including circRNA, lncRNA, miRNA ,mRNA or KEGG/Reactome/Wikipathway pathway IDs.}
  \item{regulationType}{character string naming the type of regulation which can be located including "circRNA-miRNA", "lncRNA-miRNA", "lncRNA-mRNA", "miRNA-mRNA", "miRNA-pathway", "mRNA-pathway".}
  \item{allRegulation}{a logical value indicates locating all related regulations from the storage or internal regulations between input elements.}
}
\value{
The dataframe of regulatory data with the colnames of Node1, Node2, Datasoure and RegulationType.
}
\description{
Locate the experimental binary regulatory data among biological elements in regional Heterogeneous network which is integrated from a wide range of bioinformatics databases.
}
\examples{
binaryRegulation(transcriptomeList = c("2309", "3838", "MIMAT0000255"), regulationType = "miRNA-mRNA",  allRegulation = FALSE)
}
