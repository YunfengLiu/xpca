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
    return(W);
  }
  return(W);
}


ipca = function(X, W = NA, tol = 1e-6, max_iter_count = NA) {
  T = dim(X)[1];
  n = dim(X)[2];
  if (!is.matrix(W)) {
    pca = prcomp(X, retx = TRUE);
    W = as.matrix(pca$rotation);
  } else {
    W = as.matrix(W);
  }
  V = cov(X);
  W = ogita_aishima(V, W, max_iter_count, tol);

  means = as.vector(colMeans(X));
  X_centered = X - matrix(rep(means, T), ncol = length(means), byrow = TRUE);

  P = X_centered %*% W;
  return(list(P = P, W = W));
}

# alpha: scaler
# X: Txn data frame
ewmcov = function(alpha, X) {
  X = as.matrix(X);
  T = dim(X)[1];
  n = dim(X)[2];
  m = list();
  x_1 = matrix(X[1,], n);
  m[[1]] = x_1;
  S = list();
  S_1 = matrix(rep(0, n*n), n, n);
  S[[1]] = S_1;
  for (t in 2:T) {
    x_t = matrix(X[t,], n);
    m_t = (1-alpha) * x_t + alpha * m[[t-1]];
    m[[t]] = m_t;
    x_t_centered = x_t - m_t;
    single_cov = (x_t_centered %*% t(x_t_centered));
    S[[t]] = (1-alpha) * single_cov + alpha * S[[t-1]];
  }
  return(list(means = m, covs = S));
}

sorted_eig = function(A) {
  cov_init = cov(A);
  eigens_init = eigen(cov_init);
  ev_init = eigens_init$vectors;
  ev_vals = eigens_init$values;
  indices_orders = order(ev_vals, decreasing = TRUE);
  W = ev_init[, indices_orders];
  return(W);
}

ewmpca = function(X, alpha, W_initial, tol = 1e-6, max_iter_count=NA) {
  X = as.matrix(X);
  T = dim(X)[1];
  n = dim(X)[2];
  m = list();
  x_1 = matrix(X[1,], n);
  m[[1]] = x_1;
  S = list();
  S_1 = matrix(rep(0, n*n), n, n);
  S[[1]] = S_1;
  Z = matrix(, T, n);
  Z[1,] = rep(0, n);
  W = W_initial;
  for (t in 2:T) {
    x_t = matrix(X[t,], n);
    m_t = (1-alpha) * x_t + alpha * m[[t-1]];
    m[[t]] = m_t;
    x_t_centered = x_t - m_t;
    single_cov = (x_t_centered %*% t(x_t_centered));
    S[[t]] = (1-alpha) * single_cov + alpha * S[[t-1]];
    W = ogita_aishima(S[[t]], W, max_iter_count, tol, sort_by_eigenvalues = TRUE);
    Z[t,] = t(x_t_centered) %*% W;
  }
  return(Z); # EWM PCAs
}

