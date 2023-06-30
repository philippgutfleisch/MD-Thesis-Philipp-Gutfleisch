data <- gsva_matrix

char_columns <- sapply(data, is.character)             # Identify character columns
data_chars_as_num <- data                              # Replicate data
data_chars_as_num[ , char_columns] <- as.data.frame(   # Recode characters as numeric
        apply(data_chars_as_num[ , char_columns], 2, as.numeric))
sapply(data_chars_as_num, class) 

gsva_matrix <- data_chars_as_num
