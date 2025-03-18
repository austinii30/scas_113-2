
source("code.R")

# 3.(a)
beg <- Sys.time()
rN.1(1000000)
end <- Sys.time()
cat("rN.1, time: ", beg-end, "\n")


beg <- Sys.time()
rN.2(1000000)
end <- Sys.time()
cat("rN.2, time: ", beg-end, "\n")

# 3.(b)

