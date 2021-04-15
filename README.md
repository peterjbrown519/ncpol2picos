# ncpol2picos
Python package to relax noncommutative polynomial optimization problems to picos SDPs

The package was created to fix some bugs in Peter Wittek's ncpol2sdpa. In particular ncpol2sdpa no longer supports complex moment matrices. I was personally unable to fix these bugs and so instead rewrote the package (in a suboptimal way) but the main design ideas of ncpol2sdpa are still there and some code has been borrowed. The package differs from ncpol2sdpa in that it handles the noncommutative polynomial algebra natively rather than relying on sympy. At the moment this means that it is much slower than ncpol2sdpa and we are looking to improve this. 

For almost all applications however, I would recommend a user to use ncpol2sdpa if they can. The package ncpol2sdpa has been tested much more than this one, it is faster and has greater functionality! If for whatever reason you can't use ncpol2sdpa then you are welcome to try this package. 
