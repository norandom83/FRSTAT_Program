# Size and power of the two between-cointegration tests.
# Run this file from the folder in which results should be saved.

############################################################
## Settings
############################################################

seed = 123456
seed_MC = 654321
set.seed(seed)

titer = 1000                         # Number of Monte Carlo repetitions
nobslist = c(200,400,800)
nbasis = 20                         # Number of nonconstant orthonormal Fourier directions
kx = 2
alphalist = c(.01,.05,.10)
reportalpha = .05

R_MC = 2000                         # auxiliary null draws per sample for feasible critical values 
ngrid = 499                         # number of Brownian subintervals
maglist = c(0,.30)                   # 0: size; .30: power
powerdir = 1                        # additional trend along phi_1 only

design = data.frame(name = c("IID", "Serial", "Serial + endogeneity"),
                    rho = c(0,.3,.3), eta = c(0,0,.3))

outdir = "results_between_short"
print_every = 50

############################################################
## Functions
############################################################

center = function(x) sweep(x,2,colMeans(x),"-")
csum = function(x) apply(x,2,cumsum)

parzen = function(x) {
  x = abs(x)
  ans = numeric(length(x))
  ind = x <= .5
  ans[ind] = 1-6*x[ind]^2+6*x[ind]^3
  ind = x > .5 & x < 1
  ans[ind] = 2*(1-x[ind])^3
  ans
}

psdroot = function(x) {
  x = (x+t(x))/2
  ee = eigen(x,symmetric = TRUE)
  # Truncate negative numerical eigenvalues at zero; no ridge is added.
  sweep(ee$vectors,2,sqrt(pmax(ee$values,0)),"*") %*% t(ee$vectors)
}

dgp = function(nobs,nbasis,rho,eta,mag,direction) {
  # x_0 = 0. Stationary AR states are initialized from their exact joint law.
  # Rows are times and columns are orthonormal Fourier coefficients.
  eps = matrix(rnorm((nobs+1)*2),nobs+1,2)
  z = matrix(rnorm((nobs+1)*nbasis),nobs+1,nbasis)
  z[,1:2] = eta*eps+sqrt(1-eta^2)*z[,1:2]
  Sigma_x = matrix(c(1,.3,.3,1),2,2)
  xinnovation = eps %*% chol(Sigma_x)

  dx = matrix(0,nobs+1,2)
  ecoef = matrix(0,nobs+1,nbasis)
  dx[1,] = xinnovation[1,]
  ecoef[1,] = z[1,]
  for (tt in 2:(nobs+1)) {
    dx[tt,] = rho*dx[tt-1,]+sqrt(1-rho^2)*xinnovation[tt,]
    ecoef[tt,] = rho*ecoef[tt-1,]+sqrt(1-rho^2)*z[tt,]
  }
  dx = dx[-1,,drop = FALSE]
  ecoef = ecoef[-1,,drop = FALSE]
  x = csum(dx)
  u = sweep(ecoef,2,seq_len(nbasis),"/")
  q = cumsum(rnorm(nobs))
  if (mag > 0) u[,direction] = u[,direction]+mag*q

  beta = matrix(0,nbasis,2)
  beta[c(1,3),1] = c(1,.5)
  beta[c(2,4),2] = c(.8,-.4)
  beta0 = c(.2,-.1,rep(0,nbasis-2))
  y = sweep(x %*% t(beta)+u,2,beta0,"+")
  list(x = x,dx = dx,y = y)
}

sample_inputs = function(x,y,dx,h) {
  nobs = nrow(x)
  p = ncol(x)
  j = ncol(y)
  xc = center(x)
  yc = center(y)
  # Ordinary least squares with an intercept, not fully modified residuals.
  Bpre = qr.solve(xc,yc)
  upre = yc-xc %*% Bpre
  Vhat = crossprod(upre)/nobs
  eV = eigen(Vhat,symmetric = TRUE)
  vV = eV$vectors[,1]
  varV = eV$values[1]
  supre = csum(upre)
  Kscaled = crossprod(supre)/nobs^2
  statK = eigen(Kscaled,symmetric = TRUE,only.values = TRUE)$values[1]
  statV = drop(crossprod(vV,Kscaled %*% vV))/varV

  # Paper convention: lag j uses an earlier dx and a later residual.
  # dx is not sample-demeaned; all lags use divisor nobs, not nobs-j.
  Z = cbind(dx,upre)
  Omega = crossprod(Z)/nobs
  Omegaplus = crossprod(upre,dx)/nobs
  if (h >= 1) for (lag in seq_len(h)) {
    weight = parzen(lag/h)
    if (weight == 0) next
    early = seq_len(nobs-lag)
    late = early+lag
    G = crossprod(Z[late,,drop = FALSE],Z[early,,drop = FALSE])/nobs
    Omega = Omega+weight*(G+t(G))
    Omegaplus = Omegaplus+weight*G[p+seq_len(j),seq_len(p),drop = FALSE]
  }
  list(stat = c(K = statK,V = statV),Omega = (Omega+t(Omega))/2,
       Omegaplus = Omegaplus,vV = vV,varV = varV,
       residuals = upre,Vhat = Vhat,Kscaled = Kscaled)
}

primitive_paths = function(dW,kx) {
  # Independent standard Brownian increments; all ordinary sums use left endpoints.
  n = nrow(dW)
  dt = 1/n
  r = (0:(n-1))/n
  p = seq_len(kx)
  W = rbind(rep(0,ncol(dW)),csum(dW))
  Wleft = W[seq_len(n),,drop = FALSE]
  Wbridge = Wleft-outer(r,W[n+1,])
  dWbridge = sweep(dW,2,W[n+1,]/n,"-")
  Xc = center(Wleft[,p,drop = FALSE])
  J = crossprod(Xc)*dt
  Q = rbind(rep(0,kx),csum(Xc)[seq_len(n-1),,drop = FALSE])*dt
  H = Q %*% solve(J)
  A = crossprod(dWbridge,Xc)
  F = cbind(Wbridge-H %*% t(A),-H)
  list(F = F,K = crossprod(F)*dt)
}

make_bank = function(R_MC,ngrid,kx,nbasis) {
  # Cache Gram matrices instead of all Brownian paths.
  # This is the same joint Brownian simulation as in the Supplement, written
  # after whitening x. Sample covariance inputs are still reestimated EVERY time.
  bank = array(0,c(nbasis+2*kx,nbasis+2*kx,R_MC))
  for (ei in seq_len(R_MC)) {
    dW = matrix(rnorm(ngrid*(kx+nbasis)),ngrid,kx+nbasis)/sqrt(ngrid)
    bank[,,ei] = primitive_paths(dW,kx)$K
  }
  bank
}

null_loading = function(inputs,kx) {
  Om = inputs$Omega
  ix = seq_len(kx)
  iu = (kx+1):nrow(Om)
  Lx = t(chol(Om[ix,ix,drop = FALSE]))
  C = t(forwardsolve(Lx,t(Om[iu,ix,drop = FALSE])))
  S = psdroot(Om[iu,iu,drop = FALSE]-tcrossprod(C))
  D = t(forwardsolve(Lx,t(inputs$Omegaplus)))
  # Joint W_U = C*B_x + S*B_e; the one-sided term contributes -D*H.
  cbind(C,S,D)
}

simulate_null = function(inputs,bank,kx) {
  loading = null_loading(inputs,kx)
  loading_t = t(loading)
  av = drop(loading_t %*% inputs$vV)
  vals = matrix(NA_real_,dim(bank)[3],2,dimnames = list(NULL,c("K","V")))
  for (ei in seq_len(nrow(vals))) {
    K0 = bank[,,ei]
    KR = loading %*% K0 %*% loading_t
    vals[ei,1] = eigen(KR,symmetric = TRUE,only.values = TRUE)$values[1]
    vals[ei,2] = drop(crossprod(av,K0 %*% av))/inputs$varV
  }
  vals
}

############################################################
## LaTeX tables
############################################################

write_tables = function(ans,settings,outdir) {
  a = ans[abs(ans$alpha-settings$reportalpha) < 1e-12,]
  stamp = paste0("Rejection frequencies in percent at the ",
                 100*settings$reportalpha,"\\% nominal level.")
  size_note = paste0("Each cell uses ",
    format(settings$titer,big.mark = ",",scientific = FALSE,trim = TRUE),
    " independent data replications. Critical values are reestimated for each sample using ",
    format(settings$R_MC,big.mark = ",",scientific = FALSE,trim = TRUE),
    " auxiliary draws and the Parzen kernel with $h$ chosen as the nearest integer to $T^{1/4}$. Parentheses give Monte Carlo standard errors in percentage points, conditional on the common auxiliary draws.")
  entries = function(row) sprintf("%.1f (%.1f) & %.1f (%.1f)",
                                  row$rateK,row$mcseK,row$rateV,row$mcseV)
  table_lines = function(rows,caption,label,note) {
    lines = c("\\begin{table}[h!]", "\\centering",
      "\\renewcommand{\\arraystretch}{0.7}",
      paste0("\\caption{",caption," ",stamp,"}"),
      paste0("\\label{",label,"}"), "\\begin{tabular}{lccccc}","\\toprule",
      "Design & $\\rho$ & $\\eta$ & $T$ & $\\widehat{\\mathcal T}_K$ & $\\widehat{\\mathcal T}_V$ \\\\",
      "\\midrule")
    for (ii in seq_len(nrow(rows))) {
      row = rows[ii,]
      lines = c(lines,paste0(row$design," & ",row$rho," & ",row$eta," & ",
                             row$T," & ",entries(row)," \\\\"))
    }
    c(lines,"\\bottomrule","\\end{tabular}",
      paste0("\\par\\smallskip\\begin{minipage}{\\textwidth}\\footnotesize ",
             note,"\\end{minipage}"),"\\end{table}")
  }
  lines = table_lines(a[a$mag == 0,],
    "Empirical size of the between-cointegration tests.","tab:bc_size",size_note)
  writeLines(lines,file.path(outdir,"table_size.tex"))

  power_caption = paste0("Empirical power of the between-cointegration tests with $\\delta=",
                         settings$maglist[settings$maglist > 0],"$.")
  power_note = paste0("The Monte Carlo settings are as in Table~\\ref{tab:bc_size}. Power is the rejection frequency under $U_t=U_t^S+",
    settings$maglist[settings$maglist > 0],"q_t\\phi_",settings$powerdir,"$ and is not size-adjusted.")
  lines = table_lines(a[a$mag > 0,],power_caption,"tab:bc_power",power_note)
  writeLines(lines,file.path(outdir,"table_power.tex"))
  writeLines(c("\\documentclass[11pt]{article}","\\usepackage[margin=1in]{geometry}",
               "\\usepackage{amsmath,booktabs}","\\begin{document}",
               "\\input{table_size.tex}","\\clearpage","\\input{table_power.tex}",
               "\\end{document}"),file.path(outdir,"tables.tex"))
}

############################################################
## Simulation
############################################################

dir.create(outdir,recursive = TRUE,showWarnings = FALSE)
settings = list(seed = seed,seed_MC = seed_MC,titer = titer,nobslist = nobslist,nbasis = nbasis,kx = kx,
                alphalist = alphalist,reportalpha = reportalpha,R_MC = R_MC,
                ngrid = ngrid,maglist = maglist,powerdir = powerdir,design = design)

# Reuse independent auxiliary draws across samples (common random numbers).
# This reduces runtime without replacing fitted inputs with population values.
set.seed(seed_MC)
bank = make_bank(R_MC,ngrid,kx,nbasis)
set.seed(seed)

alternatives = rbind(data.frame(mag = 0,direction = 0),
                     expand.grid(mag = maglist[-1],direction = powerdir))
result = data.frame()
raw = list()
cell = 0
for (iset in seq_len(nrow(design))) {
  rho = design$rho[iset]
  eta = design$eta[iset]
  for (model in seq_along(nobslist)) {
    nobs = nobslist[model]
    h = max(1,as.integer(round(nobs^(1/4))))
    for (im in seq_len(nrow(alternatives))) {
      mag = alternatives$mag[im]
      direction = alternatives$direction[im]
      cell = cell+1
      testK = testV = criticalK = criticalV = matrix(NA_real_,titer,length(alphalist))
      statistics = pvalues = matrix(NA_real_,titer,2,dimnames = list(NULL,c("K","V")))

      for (iter in seq_len(titer)) {
        dat = dgp(nobs,nbasis,rho,eta,mag,direction)
        fit = sample_inputs(dat$x,dat$y,dat$dx,h)
        draws = simulate_null(fit,bank,kx)
        cvK = as.numeric(quantile(draws[,1],1-alphalist,type = 1,names = FALSE))
        cvV = as.numeric(quantile(draws[,2],1-alphalist,type = 1,names = FALSE))
        statistics[iter,] = fit$stat
        criticalK[iter,] = cvK
        criticalV[iter,] = cvV
        testK[iter,] = fit$stat[1] > cvK
        testV[iter,] = fit$stat[2] > cvV
        pvalues[iter,] = c((1+sum(draws[,1] >= fit$stat[1]))/(R_MC+1),
                           (1+sum(draws[,2] >= fit$stat[2]))/(R_MC+1))
        if (iter %% print_every == 0 || iter == titer) {
          ia = match(reportalpha,alphalist)
          cat(design$name[iset],"| T =",nobs,"| h =",h,
              "| delta =",mag,"| direction =",direction,
              "| iter =",iter,"| level =",100*reportalpha,"percent",
              "| K =",round(100*mean(testK[seq_len(iter),ia]),1),
              "| V =",round(100*mean(testV[seq_len(iter),ia]),1),"percent\n")
        }
      }

      rateK = colMeans(testK)
      rateV = colMeans(testV)
      result = rbind(result,data.frame(design = design$name[iset],rho = rho,eta = eta,
        T = nobs,h = h,mag = mag,direction = direction,alpha = alphalist,
        rateK = 100*rateK,rateV = 100*rateV,
        mcseK = 100*sqrt(rateK*(1-rateK)/titer),
        mcseV = 100*sqrt(rateV*(1-rateV)/titer)))
      raw[[cell]] = list(design = design$name[iset],T = nobs,mag = mag,direction = direction,
                         statistics = statistics,pvalues = pvalues,
                         criticalK = criticalK,criticalV = criticalV,
                         rejectK = testK,rejectV = testV)
      # Save every completed cell, including the RNG state and session details.
      saveRDS(list(settings = settings,result = result,raw = raw,
                   rng = .Random.seed,session = sessionInfo()),file.path(outdir,"results.rds"))
      write.csv(result,file.path(outdir,"rejection_rates.csv"),row.names = FALSE)
    }
  }
}
write_tables(result,settings,outdir)
capture.output(sessionInfo(),file = file.path(outdir,"sessionInfo.txt"))
print(result[result$alpha == reportalpha,],row.names = FALSE)
