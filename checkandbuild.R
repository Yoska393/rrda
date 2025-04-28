library(devtools)
library(attachment)

# when i make a change in code, it goes to this check and build
# when i make change in description -> run the function document
load_all()
document()
att_amend_desc()

check()

# when i am sure that its done
build()

# for the install
#install.packages("~/Documents/R/HayatoINRAE/rrda_0.0.0.9000.tar.gz", repos = NULL, type = "source")
install.packages("~/Documents/R/HayatoINRAE/rrda_0.1.1.tar.gz", repos = NULL, type = "source")

library(rrda)

