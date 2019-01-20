# The following excerpt of testData_t5 is used in PCMBaseCpp for a
# minibenchmark of the likelihood calculation for the MGPM model:

library(MGPMSimulations)
library(data.table)
idx <- c(289,
         449,
         833,
         929,
         1281,
         1377,
         1409,
         1857,
         1897,
         1921)
testData_t5_miniBenchmark <- testData_t5[idx]
testData_t5_miniBenchmark[, id:=idx]
usethis::use_data(testData_t5_miniBenchmark, overwrite = TRUE)
