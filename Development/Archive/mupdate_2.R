mu_mat.copy <- mu_mat

for (l in 1:n){

  # the l-th row of mu_mat, alpha_mat stacked n times
  mu_stack <- matrix(mu_mat.copy[l, ], n, p, T)
  alpha_stack <- matrix(alpha_mat[l, ], n, p, T)

  # the element-wise product of X_mat, mu_stack, and alpha stack;
  # the i,j entry is x_i,j * mu_l,j * alpha_l,j
  X_mu_alpha <- X_mat * mu_stack * alpha_stack

  # the k-th column is the rowSums of X_mu_alpha minus the k-th column of
  # X_mu_alpha (accounts for m \neq k in summation)
  X_mu_alpha_k <- matrix(rowSums(X_mu_alpha), n, p) - X_mu_alpha

  # the k-th column is y minus the k-th column of X_mu_alpha_k
  y_k <- matrix(y, n, p) - X_mu_alpha_k

  # the k-th column is d_:,l * x_:,k * y_k_:,k
  d_x_y <- D[ , l] * X_mat * y_k

  # the update of the l-th row of mu
  mu_mat[l, ] <- S_sq[l, ] / sigmasq * colSums(d_x_y)

}
