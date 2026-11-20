# The largest totient ratio with prescribed 2-adic valuation

Code for the paper "The largest totient ratio with prescribed 2-adic valuation" (Y. Shmalo).

For s >= 0 let M(s) be the largest value of n/phi(n) over the integers n with v_2(phi(n)) = s. The paper shows that M(s) is twice the optimum of a knapsack problem over the odd primes, and that, with L = log s, M(s)/e^gamma = L + log L - C + (log L - kappa)/L + O((log L)^2/L^2).

The paper source and PDF are in `paper/`.

`check.py` uses only the Python standard library. It does the following:

- It solves the knapsack problem exactly for s <= 12 and prints M(s) with an extremal integer n (Table 1 of the paper).
- It checks that n/phi(n) <= M(v_2(phi(n))) for every 2 <= n <= 200000.
- It checks the values against exhaustive search over subsets for s <= 6.
- It computes the constants C, D and kappa to 40 digits.
- It computes M(s) at s = 16, 32, ..., 512 by dynamic programming, together with the three- and four-term approximations (the remark at the end of Section 5).

The output of `python check.py` is in `check_output.txt`. It runs in a few seconds.
