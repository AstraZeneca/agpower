
# lam.type = "base"
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "blind", lam.type = "base")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "blind", lam.type = "base")

# lam.type = "pool"
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "blind", lam.type = "pool")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "blind", lam.type = "pool")

# lam.type = "base"
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "blind", lam.type = "base")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "blind", lam.type = "base")

# lam.type = "pool"
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = NULL, ar = 0.5, frailty.type = "blind", lam.type = "pool")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = c(1000,1100.1), RR = 0.7, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = NULL, pow = 0.8, ar = 0.5, frailty.type = "blind", lam.type = "pool")


# lam.type = "base"
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "blind", lam.type = "base")
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "blind", lam.type = "base")
power.lwyy.test(N = NULL, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "unbli", lam.type = "base")
power.lwyy.test(N = NULL, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "blind", lam.type = "base")

# lam.type = "pool"
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = NULL, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "blind", lam.type = "pool")
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = 1000, RR = NULL, thta = 1, L = NULL, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "blind", lam.type = "pool")
power.lwyy.test(N = NULL, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "unbli", lam.type = "pool")
power.lwyy.test(N = NULL, RR = NULL, thta = 1, L = 1000, tau = 0.9, lam = 1.1, alp = 0.05, pow = c(0.8, 0.9, 0.99), ar = 0.5, frailty.type = "blind", lam.type = "pool")



