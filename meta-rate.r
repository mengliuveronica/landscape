Break = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "\"meta_rate\", An R program for interrater reliability in L2 Meta-Analyses.
 Copyright (C) 2019-present  Reza Norouzian, rnorouzian@gmail.com\n"

message(Break, notice, Break)


Break = "\n*********************************************************************************\n"

cite <- "To cite this package use:\n\nNorouzian,R.(2021). Interrater reliability in second language meta-analyses:\nThe case of categorical moderators. Studies in Second Language Acquisition, 43,896-915."

cat(Break, cite, Break)

#=============================================================================================================================

trim <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

#===============================================================================================================================


d.prepos <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, stder = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post = NA, control = NA, outcome = NA, time = NA, ...) 
{
}  

#========================================================================================

rm.allrowNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE])
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}

#===============================================================================================================================

rm.allcolNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[, colSums(is.na(i) | i == "") != nrow(i), drop = FALSE])
    
  } else { X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE] }
}

#===============================================================================================================================

rm.colrowNA <- function(X){
  
  r <- rm.allrowNA(X)
  rm.allcolNA(r)  
  
}                                      

#================================================================================================================================

drop.col <- function(dat, vec){
  
  vec <- trimws(vec)
  names(dat) <- trimws(names(dat))
  
  f <- function(dat, vec) {
    i1 <- !names(dat) %in% vec
    setNames(dat[i1], names(dat)[i1])
  }
  
  if(inherits(dat, "list")) { lapply(dat, f, vec = vec)
  } else { f(dat, vec) }
}               

#================================================================================================================================              

full.clean <- function(X, omit, all = TRUE, omit.auto.suffix = TRUE)
{
  
  X <- rm.colrowNA(X)
  
  X <- if(inherits(X, "list") & omit.auto.suffix){ lapply(X, function(x) trim(setNames(x, sub("\\.\\d+$", "", names(x))))) 
    
  } else if(inherits(X, "data.frame") & omit.auto.suffix) { trim(setNames(X, sub("\\.\\d+$", "", names(X)))) } else { X }
  
  if(all){ X } else { 
    
    drop.col(X, vec = omit)
  }
}                

#===============================================================================================================================

detail2 <- function(X, useNA = "ifany"){
  
  nr <- nrow(X)
  nc <- ncol(X)
  tab <- table(row(X), unlist(X), useNA = useNA)
  pj <- colSums(tab)/(nr * nc)
  pjk <- (colSums(tab^2) - nr * nc * pj)/(nr * nc * (nc - 1) * pj)
  K <- (pjk - pj)/(1 - pj)
  h <- names(K)
  h[is.na(h)] <- "NA"
  setNames(K, h)
}         

#===============================================================================================================================

detail <- function(X, useNA = "ifany") {
  X <- as.matrix(X)
  tab <- table(row(X), unlist(X), useNA = useNA)
  w <- diag(ncol(tab))
  rosum <- rowSums(tab)
  obs_oc <- tab * (t(w %*% t(tab)) - 1)
  obs_c <- colSums(obs_oc)
  max_oc <- tab * (rosum - 1)
  max_c <- colSums(max_oc)
  SA <- obs_c / max_c
  h <- names(SA)
  h[is.na(h)] <- "NA"
  setNames(SA, h)
}                                                       

#===============================================================================================================================                                                          

set.margin <- function() 
{
  par(mgp = c(1.5, 0.14, 0), mar = c(2.5, 2.6, 1.8, .5), 
      tck = -0.02)
}                                                         

#===============================================================================================================================                                                          


splot <- function(y, main, lwd = 5, lend = 2, show.sa = FALSE, digits = 3, cex.sa = .9){
  
  ll <- length(y)
  
  x <- seq_len(ll)
  
  plot(x, y, type = "h", main = main, xlim = c(.95, 1.02*max(x)), ylim = 0:1,
       ylab = "SA%", xaxt = "n", xlab = "Category", lend = lend, lwd = lwd,
       col = colorRampPalette(c(4, 2))(ll), font.lab = 2, 
       panel.first = abline(h = 0, col = 8), las = 1, cex.axis = .9)
  
  if(show.sa) text(x[y != 0]-.015, .4, round(y[y != 0], digits), pos = 2, xpd = NA, srt = 90, font = 2, cex = cex.sa)
  
  axis(1, at = x, labels = names(y))
}


#===============================================================================================================================

irr <- int <- function (X, nsim = 1e3, useNA = "ifany", level = .95, digits = 6, raw = TRUE) 
{
  
  if(!inherits(X, c("data.frame", "matrix", "table"))) stop("Ratings must be 'data.frame', 'matrix', and if not raw, a 'table'.", call. = FALSE)
  
  if(raw) X <- table(row(X), unlist(X), useNA = useNA)
  
  X2 <- X * (X - 1)
  sumcol <- colSums(X)
  sumrow <- rowSums(X)
  nc <- ncol(X)
  nr <- nrow(X)
  tot <- sum(X)
  pij <- X2/(sumrow * (sumrow - 1))
  pi <- rowSums(pij)
  p <- mean(pi)
  pj <- sumcol/tot
  pj2 <- pj^2
  pe <- sum(pj2)
  KAPPA <- (p - pe)/(1 - pe)
  s <- (nc * p - 1)/(nc - 1)
  pi.v.boot <- replicate(nsim, pi.boot <- sample(pi, size = nr, replace = TRUE))
  p.boot <- colMeans(pi.v.boot)
  s.boot <- sapply(seq_len(nsim), function(i) (nc * p.boot[i] - 1)/(nc - 1))
  
  p <- (1 - level) / 2
  s.boot.ci <- quantile(s.boot, probs = c(p, 1-p), na.rm = TRUE)
  
  return(round(c(Fleiss_KAPPA = KAPPA, 
                 Sindex = s, 
                 lower = s.boot.ci[[1]], 
                 upper = s.boot.ci[[2]], 
                 conf.level = level), digits))
}                                      


#===============================================================================================================================                   


int2 <- function(X, level = .95, useNA = "ifany", nsim = 1e3, digits = 4, raw = TRUE){ 
  
  X <- table(row(X), unlist(X), useNA = useNA)
  
  agree.mat <- as.matrix(X) 
  n <- nrow(agree.mat) # number of studies or groups within studies
  q <- ncol(agree.mat) # number of categories
  f <- 0               # population correction 
  
  weights.mat <- diag(q)
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  ac1 <- (pa-pe)/(1-pe)
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec == 0) # this replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  
  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-ac1) * (pe.ivec-pe)/(1-pe)
  
  var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - ac1)^2)
  stderr <- sqrt(var.ac1)
  p.value <- 2*(1-pt(ac1/stderr,n-1))
  
  lower <- ac1 - stderr*qt(1-(1-level)/2,n-1)
  upper <- min(1,ac1 + stderr*qt(1-(1-level)/2,n-1))
  
  return(round(c(AC = ac1, lower = lower, upper = upper, conf.level = level), digits))
}                   

#===============================================================================================================================

is.constant <- function(x) length(unique(x)) == 1L 

#===============================================================================================================================                   

drop.inner.list <- function(L, what, omit.auto.suffix = TRUE) {
  
  if(omit.auto.suffix) L <- lapply(L, function(x) setNames(x, sub("\\.\\d+$", "", names(x))))
  
  L[!names(L) %in% what]
}


#===============================================================================================================================

meta_rate <- function(..., sub.name = "group.name", nsim = 1e3, level = .95,
                      useNA = "ifany", type = c("s", "ac"), na.rm = FALSE, 
                      digits = 3, common = FALSE, all = TRUE, drop = NULL,
                      plot = TRUE, lwd = 5, lend = 1, show.sa = TRUE, 
                      sub.level = NULL, study.level = NULL, file.name = NULL,
                      reset = TRUE, rev.page = FALSE, cex.sa = .9)
{
  
  r <- list(...) 
  
  type <- trimws(type)
  type <- match.arg(type)
  
  if(!all(sapply(r, inherits, c("data.frame", "matrix")))) stop("Coding sheet(s) must be 'Excel CSV' files, 'data.frame' or 'matrix'.", call. = FALSE)
  
  n.df <- length(r)
  
  r <- lapply(r, as.data.frame)
  
  ar <- formalArgs(d.prepos)[-c(2, 22)]
  
  r <- full.clean(r, ar, all)
  
  check <- all(sapply(r, function(i) "study.name" %in% names(i)))
  
  if(!check) stop("Add a new column named 'study.name' to the coding sheet(s).", call. = FALSE)
  
  r <- lapply(r, function(x) do.call(rbind, c(split(x, x$study.name), make.row.names = FALSE)))
  
  drop <- trimws(drop)              
  drop <- drop[!drop %in% "study.name"]
  
  if(length(drop) != 0) r <- drop.col(r, drop)   
  
  r <- unname(r)
  
  sub.name <- trimws(sub.name)
  if(n.df == 1) tbl <- table(names(r[[1]])[!names(r[[1]]) %in% c("study.name", sub.name)])
  
  com.names <- if(n.df >= 2) { 
    
    ok <- is.constant(sapply(r, nrow))
    
    if(!ok) stop("The coding sheets don't have the same number of rows.", call. = FALSE)
    
      vec <- names(unlist(r, recursive = FALSE))
      unique(vec[duplicated(vec)])
    
  } else { 
    
      names(which(tbl >= 2))
  }
  
  dot.names <- if(all) com.names else com.names[!com.names %in% ar]
  
  if(length(dot.names) == 0) stop("No 2 raters detected OR no two moderators names match.", call. = FALSE)
  
  if(n.df >= 2) { 
    
    r <- do.call(cbind, r)
    
    tbl <- table(names(r)[!names(r) %in% c("study.name", sub.name)]) 
    
  } else { r <- r[[1]]
  
  }
  
  n.coder <- tbl[tbl >= 2]
  
  i1 <- colnames(r) != 'study.name'
  st.level <- names(which(sapply(split.default(r[i1], names(r)[i1]), function(x) 
   base::all(!colSums(!aggregate(.~ study.name, transform(x, study.name = r$study.name), FUN = is.constant)[-1])))))
  
  st.level <- st.level[st.level %in% dot.names]
  
  exclude <- trimws(sub.level)
  
  st.level <- st.level[!st.level %in% c(exclude,"study.name", sub.name)]
  
  L <- split.default(r[names(r) %in% dot.names], names(r)[names(r) %in% dot.names])
  
  if(length(st.level) != 0) L[st.level] <- lapply(L[st.level], function(x) x[ave(seq_along(x[[1]]), r$study.name, FUN = seq_along) == 1, ]) 
  
  L <- drop.inner.list(L, c("study.name", sub.name))
  
  if(na.rm) L <- lapply(L, na.omit)
  
  f <- if(type == "s") int else int2
  out <- lapply(L, f, nsim = nsim, level = level, digits = digits, useNA = useNA, raw = TRUE)
  
  A <- lapply(L, detail, useNA = useNA)
  
  study.level <- sapply(seq_along(out), function(i) names(out)[[i]] %in% st.level)
  
  d <- data.frame(out)
  
  d[] <- lapply(d, as.list)
  
  if(plot){
    
    n <- length(L)
    
    if(reset){
      graphics.off()
      org.par <- par(no.readonly = TRUE)
      on.exit(par(org.par))
    }
    dev <- if(!rev.page) n2mfrow(n) else rev(n2mfrow(n))
    if(n > 1L) { par(mfrow = dev) ; set.margin() }
    
    invisible(mapply(splot, y = A, main = names(A), lwd = lwd, lend = lend, show.sa = show.sa, digits = digits, cex.sa = cex.sa))
  }
  
  res <- data.frame(t(rbind(d, row.comprd = sapply(L, nrow), min.cat = sapply(A, function(i) if(any(i < 1)) names(i)[which.min(i)] else "--"), 
                            n.coder = n.coder, study.level = ifelse(study.level, "Yes", "No"))))
  
  output <- data.frame(lapply(res, unlist))
  
  if(common) output <- output[output$n.coder == max(output$n.coder),]
  
  file.name <- trimws(file.name)
  
  if(length(file.name) != 0){
    nm <- paste0(file.name, ".csv")
    ur <- try(write.csv(output, nm), silent = TRUE)
    if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
    message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
  }
  
  return(output)
}


#================================================================================================================================================================

irr.diag <- function(X, useNA = "ifany"){
  
  a <- detail2(X, useNA = useNA)
  b <- detail(X, useNA = useNA)
  
  round(data.frame(Fleiss_KAPPA_cat. = a, SA = b), 3)
}

#==========================================================================================================================================

find.irr <- function(X, what, sub.name = "group.name"){
  
  if(!inherits(X, "data.frame")) stop("Data must be an Excel CSV file or a 'data.frame'.", call. = FALSE)
  
  X <- full.clean(X)
  
  s <- as.list(substitute(what))  
  
  res <- Filter(NROW, X[rowSums(X[grep(as.character(s[[2]]), names(X))] == s[[3]], na.rm = TRUE) > 0,][c("study.name", sub.name)])
  
  if(length(res) == 0) NULL else res
}        

#===========================# Datasets # ===================================================================================== 

#table1 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/t1.csv", row.names = 1)
#table2 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/t2.csv", row.names = 1)          
#table3 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/t3.csv", row.names = 1)
#table5 <- read.csv('https://raw.githubusercontent.com/hkil/m/master/t5.csv', row.names = 1)
#c1 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/c1.csv")
#c2 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/c2.csv")
#c3 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/c3.csv")
#c4 <- read.csv("https://raw.githubusercontent.com/hkil/m/master/c4.csv")           
