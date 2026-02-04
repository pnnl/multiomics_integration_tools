do.call(rbind,
  lapply(readRDS("Models/MOFA/splits.RDS")$splits,
         function(x) {x$in_id})
) %>%
  write.table("Models/MultiMLP/splits.txt", quote = F, row.names = F, col.names = F)
