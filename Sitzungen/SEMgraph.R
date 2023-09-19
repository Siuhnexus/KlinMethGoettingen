library(lavaan)
library(DiagrammeR)

SEMgraph = function(lavaanFit, showPar=TRUE, showSignif=TRUE, showCov=TRUE, showInt=TRUE, showVar=FALSE, filterPhant=TRUE) {
  generateEdges = function(perows, showPar=T, showSignif=T, startDir=F, endDir=T) {
    if(endDir == T) {
      endDir = startDir
    }
    hasStartPort = startDir != F
    hasEndPort = endDir != F
    makeBold = showPar & showSignif & !is.na(perows$stars) & perows$stars != ""
    hasLabel = perows$label != ""
    if(length(hasLabel) == 0) hasLabel = F
    showEstimate = showPar & !is.na(perows$stars)
    
    paste0(
      perows$rhs, ifelse(hasStartPort, paste0(":", startDir), ""), "->",
      perows$lhs, ifelse(hasEndPort, paste0(":", endDir), ""),
      "[label=",
      ifelse(makeBold, "< <b>", "' "),
      ifelse(hasLabel, paste0(perows$label, ifelse(showEstimate, " (", "")), ""),
      ifelse(showEstimate, paste0(
        format(round(perows$est, 3), nsmall=3, trim=T),
        perows$stars
      ), ""),
      ifelse(showEstimate & hasLabel, ")", ""),
      ifelse(makeBold, "</b>>", "'"), "]",
      ";",
      collapse="", recycle0 = T)
  }
  
  pe = parameterEstimates(lavaanFit)
  fd = pe[pe$op == "=~",]
  pe = pe[!pe$op == "=~" & (!is.na(pe$z) | pe$op == "~"),]
  # Filter phantom latents
  if(filterPhant) {
    for(f in unique(fd$lhs)) {
      if(nrow(fd[fd$lhs == f,]) == 1) {
        fd = fd[fd$lhs != f,]
      }
    }
  }
  pe = rbind(fd, pe)
  
  pe$stars = ifelse(pe$pvalue < 0.001, "***", ifelse(pe$pvalue < 0.01, "**", ifelse(pe$pvalue < 0.05, "*", "")))
  
  nodeNames = unique(c(pe$lhs, pe$rhs))
  nodeNames = nodeNames[nodeNames != ""]
  
  isFactor = rep(0, length(nodeNames))
  for(i in 1:length(nodeNames)) {
    isFactor[i] = nrow(pe[pe$op == "=~" & pe$lhs == nodeNames[i],]) > 0
  }
  
  addIntercept = nrow(pe[pe$op == "~1",]) > 0 && showInt
  
  regressions = pe[pe$op == "~",]
  factors = pe[pe$op == "=~",]
  covariances = pe[pe$op == "~~" & pe$lhs != pe$rhs,]
  variances = pe[pe$op == "~~" & pe$lhs == pe$rhs,]
  
  vizStr = paste0(
    "digraph{node[penwidth=2];",
    ifelse(addIntercept, "subgraph{rank=source;node[shape=triangle];1;}", ""),
    "node[shape=box]",
    paste0(nodeNames[!isFactor], ";", collapse="", recycle0=T),
    "node[shape=circle];",
    paste0(nodeNames[isFactor == T], ";", collapse="", recycle0=T),
    "edge[minlen=2,constraint=false];",
    ifelse(addIntercept, paste0("1->", pe[pe$op == "~1",]$lhs, ";", collapse="", recycle0 = T), ""),
    "edge[constraint=true];",
    generateEdges(regressions, showPar=showPar, showSignif=showSignif),
    "edge[dir=back];",
    generateEdges(factors, showPar=showPar, showSignif=showSignif),
    "edge[dir=both,constraint=false];",
    ifelse(showCov, generateEdges(covariances, showPar=showPar, showSignif=showSignif), ""),
    ifelse(showVar, generateEdges(variances, showPar=showPar, showSignif=showSignif, startDir="s"), ""),
    "}"
  )
  grViz(vizStr)
}
