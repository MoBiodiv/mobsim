## Test environments

## R CMD check results

Status: 1 WARNING, 4 NOTEs

See
  'C:/FelixMay/MoBiodiv/mobsim.Rcheck/00check.log'
for details.


checking for hidden files and directories ... NOTE
Found the following hidden files and directories:
  revdep/.cache.rds
These were most likely included in error. See section 'Package
structure' in the 'Writing R Extensions' manual.

checking top-level files ... NOTE
Non-standard file/directory found at top level:
  'revdep'

checking dependencies in R code ... NOTE
Missing or unexported object: 'sads::Svolkov'

checking compiled code ... NOTE
File 'mobsim/libs/x64/mobsim.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

It is good practice to register native routines and to disable symbol
search.

See 'Writing portable packages' in the 'Writing R Extensions' manual.
 WARNING
'qpdf' is needed for checks on size reduction of PDFs
R CMD check results
0 errors | 0 warnings | 4 notes

R CMD check succeeded

## Downstream dependencies
