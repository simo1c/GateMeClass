# scGateMe
scGateMe is an R package for the classification of Flow Cytometry data

Example:

gates <- data.frame(read_excel("gates.xlsx"))
colnames(gates) <- gsub("[.]", "-", colnames(gates))
gates

system.time(
  res <- scGateMe(m,
                  gates, 
                  refine = F,
                  verbose = T,
                  narrow_gate_table = T, 
                  seed = 1,
                  marker_seq_eval = F)
)
