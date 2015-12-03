(n k) = 2^n / sqrt(1/2*n*pi) * exp(-(k-(n/2))^2 / n/2)

// P(x+X successes in n+k trials | x successes in n trials)

P(n, x, k, X) := ((n+1)*binomial(k, X)*binomial(n, x))
              / ((k+n+1)*binomial(k+n, X+x));

// Shitty: when n is high and k is low, the approximation is really bad
binomial_approx(n, k) := 2^n / sqrt(1/2*n*%pi) * %e^(-(k-(n/2))^2/(n/2));
P_approx(n, x, k, X) := ((n+1)*binomial_approx(k, X)*binomial_approx(n, x))
                  / ((k+n+1)*binomial_approx(k+n, X+x));

// With beta
binomial_beta(n, k) := 1 / ((n+1) * beta(n-k+1, k+1));
P_beta(n, x, k, X) := ((n+1)*binomial_beta(k, X)*binomial_beta(n, x))
                  / ((k+n+1)*binomial_beta(k+n, X+x));

                   beta(- X - x + n + k + 1, X + x + 1)
P_beta = ---------------------------------------------------------
         (k + 1) beta(- x + n + 1, x + 1) beta(- X + k + 1, X + 1)

b = 0.9