.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "We strongly recommend that prior to using this package, you familiarise ", 
    "yourself with the underlying concepts of dBBMMs and dBGB.\n",
    "Based on the characteristics of your data, your data may be ",
    "inadequate for your model choice.\n",
    "Please read the package vignette and familiarise yourself with the ",
    "different options before deciding to run a dBBMM or dBGB.\n",
    "See the references listed in the GitHub README at ",
    "https://github.com/SimonDedman/movegroup ",
    "for further background on these methods and use ?movegroup for guidance.\n"
  )
}