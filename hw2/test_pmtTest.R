###########################################################
# Filename: test_pmtTest.R
# Course  : Unit test for the permutation test.
# Author  : Potumas Liu 
# Date    : 2025/03/24
###########################################################

source("code.R")


random_numbers <- runif(1000)

# results from the function
funcTest <- pmtTest(random_numbers)

# results from plane code
random_numbers <- random_numbers[1:(length(random_numbers) - length(random_numbers) %% 3)]
k <- 3
data_matrix <- matrix(random_numbers, ncol = k, byrow = TRUE)
rank_matrix <- t(apply(data_matrix, 1, function(x) order(x)))
# 將排名結果轉換為字符格式，便於統計頻率
rank_strings <- apply(rank_matrix, 1, paste, collapse = "")
rank_table <- table(rank_strings)
print(rank_table)
expected_counts <- rep(nrow(data_matrix) / factorial(k), factorial(k))
chi_test <- chisq.test(rank_table, p = rep(1/factorial(k), factorial(k)))


print(funcTest == chi_test$p.value)
