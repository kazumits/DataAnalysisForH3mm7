# Functional toys
# compose: f . g
"%.%" <- function(f,g) function(x) f(g(x))
# apply function: f $ x
"%$%" <- function(f,x) f(x)
# partial application:: function -> function
papp <- function(f,...) function(x) f(x,...)
# `dot(anyfunc(.))` generates `function(.) anyfunc(.)`
dot <- function(f) eval(substitute(function(.) f))
# col-wise function -> row-wise function
flipwise <- function(f) t %.% f %.% t

negate <- function(x) -x

# DEMO:
# flipwise %$% papp(scale,scale=FALSE,center=TRUE) %.% dot(.*10) %$% matrix(1:4,2)

# compose with dots: imcomplete, fails when #chains >= 3
dot_ <- function(f,env) eval(function(.) f,list(f=f)) # SE
"%dot%" <- function(f,g) {
  print(match.call())
  f_ <- eval(substitute(function(.) f))
  g_ <- eval(substitute(function(.) g))
  #browser()
  function(x) f_(g_(x))
}

