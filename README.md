# orthoMTL

## Description

OrthoMTL is an R package to solve penalized learning tasks, assuming orthogonality between tasks support. It also allows to perform feature selection and sparsity-inducing schemes.

Source extending the method presented in '
[On Learning Matrices with Orthogonal Columns or Disjoint Supports' (Vervier et al., 2014)] (http://link.springer.com/chapter/10.1007%2F978-3-662-44845-8_18)

The main addition is the possibility to work on survival data.

## Installation

To install OrthoMTL from R, type:

```{r}
cred <- git2r::cred_user_pass(rstudioapi::askForPassword("username"), rstudioapi::askForPassword("Password"))
devtools::install_git("http://gitlabce.apps.dit-prdocp.novartis.net/VERVIKE1/orthomtl.git", credentials = cred,lib = 'YOURPATH')

library(orthoMTL,lib.loc = 'YOURPATH')
```

Check the [vignettes](vignettes/) for examples.
