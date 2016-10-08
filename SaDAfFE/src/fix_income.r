bondvalue = function(c, t, r, p){
# Compute bv = bond values (current prices) corresponding
# to all values of yield to maturity in the
# input vector r
#
# INPUT
#  c = coupon payment (semiannual)
#  t = time to maturity (in years)
#  r = vector of yieds to maturity(semiannual rates)
#  p = par value
  bv = c / r + (p - c / r) * (1 + r) ^ (-2 * t)
  bv
}
