### Troubleshooting MOFA2

Rflomics uses the R package [MOFA2](https://www.bioconductor.org/packages/release/bioc/html/MOFA2.html). 
This package depends on a python script and this
can lead to several issues. The first step is to consult the MOFA2 troubleshooting
[FAQ](https://biofam.github.io/MOFA2/troubleshooting.html) on their website.

If none of the proposed solutions are resolving your issues, you can try some
additional steps and verifications. 

1. Install [Python](https://www.python.org/downloads/) or make sure you have it installed.

You can check you have python installed by running `which python` in a shell terminal. 
On windows, replace with `where python`. You should be able to see the path(s) were
Python binaries are installed. You can have multiple python version installed at once.
You have to make sure the one used by R is this correct one. 

2. Install a python environment (conda, miniconda, anaconda, ...)
3. Install pip command (if not already done). For windows users, you can use [this guide](https://phoenixnap.com/kb/install-pip-windows).
4. Install mofapy 2: in a shell terminal, run `pip install mofapy2`.
5. Set the value of the RETICULATE_PYTHON environment variable to a Python 
binary and restart your session:
* either manually indicate it in your .Rprofile (hidden file) 
* or use the following command in a terminal, in the .Rprofile folder:

```
echo "Sys.setenv(RETICULATE_PYTHON = \"path_to_python_bin\")" >> .Rprofile
```
6. Install MOFA2 using `BiocManager::install("MOFA2")` in R. 

You can check your configuration using `reticulate::py_config()` and 
`reticulate::py_list_packages()` in R. 
Your Python binary path and mofapy2 (not mofapy!) path should appear 
after step 5. If it's not the case, something went wrong. 
