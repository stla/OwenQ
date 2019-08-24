## ----setup, include=FALSE------------------------------------------------
library(OwenQ)
knitr::opts_chunk$set(collapse=TRUE)

## ------------------------------------------------------------------------
ptOwen(q=1, nu=3, delta=2)
pt(q=1, df=3, ncp=2)

## ------------------------------------------------------------------------
p1 <- pt(q=80, df=4, ncp=70)
p2 <- ptOwen(q=80, nu=4, delta=70)
wolfram <- 0.54742763380700947685
p1 - wolfram
p2 - wolfram

## ----ptOwen_fails--------------------------------------------------------
ptOwen(q=50, nu=3500, delta=50)
ptOwen(q=50, nu=3600, delta=50)
ptOwen(q=50, nu=3650, delta=50)
ptOwen(q=50, nu=3660, delta=50)
ptOwen(q=50, nu=3670, delta=50)
ptOwen(q=50, nu=3680, delta=50)

## ----pt_boost------------------------------------------------------------
OwenQ:::pt_boost(q=50, nu=3500, delta=50)
OwenQ:::pt_boost(q=50, nu=3600, delta=50)
OwenQ:::pt_boost(q=50, nu=3650, delta=50)
OwenQ:::pt_boost(q=50, nu=3660, delta=50)
OwenQ:::pt_boost(q=50, nu=3670, delta=50)
OwenQ:::pt_boost(q=50, nu=3680, delta=50)

