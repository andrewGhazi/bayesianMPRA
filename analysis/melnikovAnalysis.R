# Well let's try to use the data from Melnikov 2012

melnFiles = list.files('/mnt/bigData2/andrew/MPRA/Melnikov/',
                       pattern = '.txt')

meln = melnFiles %>% map_df(~read_tsv(paste0('/mnt/bigData2/andrew/MPRA/Melnikov/', .x),
                                      col_types = cols(ID = col_character(),
                                                       Sequence = col_character(),
                                                       Tags = col_character(),
                                                       Counts = col_character())) %>% mutate(file = .x))
