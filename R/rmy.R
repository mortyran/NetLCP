.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    crayon::green(
      "NetLCP,
More information can be found at https://github.com/rmyhandsome/NetLCP
If you use NetLCP in you publication, please cite this publication:
NetLCP: A bioinformatics tool for prioritizing regulations among diverse biological regulatory layers network with genetic variant ‘switches’ detection
Authors: MingYu Ran (rmyhandsome@163.com)
Maintainer: MingYu Ran."
    )
    ,
    crayon::red(
      "
Please read the tutorial in https://mortyran.github.io/NetLCP/ for data preparation before using NetLCP."
    )
  )
}
