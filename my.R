# A: nxn, real symmetric matrix, e.g. covariance matrix
# W: nxn, precomputed set of eigenvectors (eigenvectors are columns)
# return_extra: boolean, whether to return  (Lambdas, R, S) or just return Lambdas;
estimate_eigenvalues = function(A, W, return_extra = TRUE) {
  n = dim(A)[1];
  R = diag(1, n) - t(W) %*% W;
  S = t(W) %*% A %*% W;
  L = diag(S) / (1- diag(R));
  if (return_extra) {
    return(list(L = L, R = R, S = S));
  } else {
    return(L);
  }
}

ogita_aishima_step = function(A, W) {
  n = dim(A)[1];
  temp = estimate_eigenvalues(A, W, TRUE);
  L = temp$L;
  R = temp$R;
  S = temp$S;
  D = diag(L);
  delta = 2 * (norm(S-D, type = "F") + norm(A, type = "F") * norm(R, type = "F"));
  E = matrix(, nrow = n, ncol=n);
  for(i in 1:n) {
    for (j in 1:n) {
      if (abs(L[i] - L[j]) > delta) {
        E[i,j] = (S[i,j] + L[j]*R[i,j])/(L[j] - L[i]);
      } else {
        E[i,j] = R[i,j]/2;
      }
    }
  }
  return(W + W %*% E);
}

ogita_aishima = function(A, W, max_iter_count = NA, tol = 1e-6, sort_by_eigenvalues = FALSE) {
  itrCount =0;
  while (TRUE) {
    itrCount = itrCount + 1;
    W2 = ogita_aishima_step(A, W);
    if (!is.na(max_iter_count) && itrCount == max_iter_count) {
      break;
    }
    epsilon = norm(W2-W, type = "F");
    if (epsilon < tol) {
      break;
    }
    W = W2;
  }
  if (sort_by_eigenvalues) {
    lmbda = estimate_eigenvalues(A, W, FALSE);
    indices_orders = order(lmbda, decreasing = TRUE);
    W = W[, indices_orders];
  }
  return(W);
}


ipca = function(X, W = NA, tol = 1e-6, max_iter_count = NA) {
  T = dim(X)[1];
  n = dim(X)[2];
  if (is.na(W)) {
    pca = prcomp(X, retx = TRUE);
    W = as.matrix(pca$rotation);
  }
  V = cov(X);
  W = ogita_aishima(V, W, max_iter_count, tol);
  
  means = as.vector(colMeans(X));
  X_centered = X - matrix(rep(means, T), ncol = length(means), byrow = TRUE);
  
  P = X %*% W;
  return(list(P = P, W = W));
}



if (!require("arrow")) {
  install.packages("arrow");
  library(arrow);
}


first_year = 2008; last_year = 2020;

# read the equity_indices
df = read_parquet('/Users/lyf/projects/MQF/STAT-970/jupyter/xpca/data/equity_indices.parquet')

df = df[df$Datetime >= "2008-01-01" & df$Datetime <= "2020-12-31 23:59:59",];

Z = prcomp(df[,1:9], retx = TRUE)

Z_periods_ipca = c();
W = NA;
for (year in first_year:last_year) {
  startTime = sprintf("%s-01-01", year);
  endTime = sprintf("%s-12-31 23:59:59", year);
  subdf = df[df$Datetime >= startTime & df$Datetime <= endTime,];
  result = ipca(as.matrix(subdf[,1:9]), W);
  Z_period = result$P;
  Z_periods_ipca = rbind(Z_periods_ipca, Z_period);
  W = result$W;
}

# par(mfrow=c(3,3))
for (i in 1:9) {
  plot(Z$x[, i], Z_periods_ipca[, i])
}

