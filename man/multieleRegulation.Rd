\name{multieleRegulation}
\alias{multieleRegulation}
\title{Multielements Regulatory Data Extracting}
\usage{multieleRegulation(elementList, regulationType, allRegulation)}
\arguments{
  \item{elementList}{a vector consist of elements including circRNA, lncRNA, miRNA, mRNA or pathway IDs.}
  \item{regulationType}{character string naming the type of regulation which can be located including "circRNA-miRNA-mRNA", "lncRNA-miRNA-mRNA", "lncRNA-miRNA-pathway" and "circRNA-miRNA-pathway".}
  \item{allRegulation}{a logical value indicates locating all related regulations from the storage or internal regulations between input elements.}
}
\value{
The dataframe of regulatory data with the colnames of Node1, Node2, Datasoure and RegulationType.
}
\description{
Locate the experimental Multielement regulatory data among biological elements in regional Heterogeneous network which is integrated from a wide range of bioinformatics databases.
}
\examples{
multieleRegulation(elementList = c("2309", "3838", "ENSG00000259366", "MIMAT0000255"), regulationType = "lncRNA-miRNA-mRNA",  allRegulation = FALSE)
}
