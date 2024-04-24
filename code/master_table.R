#Master table
# This code produces a master table saved as file='Rdata/master_table.Rda'

KO.table <- read.table(paste(locn,"1_VIRGO/8.A.kegg.ortholog.txt", sep=""), 
                       header=F, row.names=1, sep="\t", stringsAsFactors=F, quote='')
KO.table$row_names <- row.names(KO.table)

path.table <- read.table(paste(locn,'1_VIRGO/8.C.kegg.pathway.copy.txt', sep=""), 
                         sep="\t", header=T, row.names=1, fill=TRUE)
path.table$row_names <- row.names(path.table)

tax.table$row_names <- row.names(tax.table)

#Make table of K# with taxon
KO.tax <- merge(KO.table, tax.table, by.x = "row_names", by.y = "row_names", all = T)


# ref.table should have 4 variables - ko, pathway, map, modules
master_table <- merge(KO.table, path.table, by.x = "row_names", by.y = "row_names", all = T)
master_table <- master_table[!duplicated(master_table$V2),]

# need to make module column
master_table$module <- master_table$map_number

##################### Assigning module column ########################
#this code fills in module values to the module column
master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00010","map00020","map00030","map00040","map00051","map00052","map00053","map00500",
                                   "map00520","map00620","map00630","map00640","map00650","map00660","map00562"), 
                               "Carbohydrate Metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00190","map00195","map00196","map00710",
                                   "map00720","map00680","map00910","map00920"),
                               "Energy metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00061","map00062","map00071","map00073",
                                   "map00100","map00120","map00121","map00140",
                                   "map00561","map00564","map00565","map00600",
                                   "map00590","map00591","map00592","map01040"),
                               "Lipid metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00230","map00240"),"Nucleotide metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00250","map00260","map00270","map00280",
                                   "map00290","map00300","map00310","map00220",
                                   "map00330","map00340","map00350","map00360",
                                   "map00380","map00400"),
                               "Amino acid metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00410","map00430","map00440","map00450",
                                   "map00460","map00470","map00480"),
                               "Metabolism of other amino acids")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00510","map00513","map00512","map00515",
                                   "map00514","map00532","map00534","map00533",
                                   "map00531","map00563","map00601","map00603",
                                   "map00604","map00511","map00540","map00542",
                                   "map00541","map00550","map00552","map00571",
                                   "map00572","map00543"),
                               "Glycan biosynthesis and metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00730","map00740","map00750","map00760",
                                   "map00770","map00780","map00785","map00790",
                                   "map00670","map00830","map00860","map00130"),
                               "Metabolism of cofactors and vitamins")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00900","map00902","map00909","map00904",
                                   "map00906","map00905","map00981","map00908",
                                   "map00903","map00907","map01052","map00522",
                                   "map01051","map01059","map01056","map01057",
                                   "map00253","map00523","map01054","map01053","map01055"),
                               "Metabolism of terpenoids and polyketides")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00940","map00945","map00941","map00944",
                                   "map00942","map00943","map00946","map00901",
                                   "map00403","map00950","map00960","map00996",
                                   "map00232","map00965","map00966","map00402",
                                   "map00311","map00332","map00261","map00331",
                                   "map00521","map00524","map00525","map00401",
                                   "map00404","map00405","map00333","map00254",
                                   "map00998","map00999","map00997"),
                               "Biosynthesis of other secondary metabolites")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map00362","map00627","map00364","map00625",
                                   "map00361","map00623","map00622","map00633",
                                   "map00642","map00643","map00791","map00930",
                                   "map00363","map00621","map00626","map00624",
                                   "map00365","map00984","map00980","map00982","map00983"),
                               "Xenobiotics biodegradation and metabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map03020","map03022","map03040"),
                               "Transcription")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map03010","map00970","map03013","map03015","map03008"),
                               "Translation")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map03060","map04141","map04130","map04120",
                                   "map04122","map03050","map03018"),
                               "Folding, sorting and degradation")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map03030","map03410","map03420","map03430",
                                   "map03440","map03450","map03460"),
                               "Replication and repair")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map03082","map03083"),
                               "Chromosome")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map03230","map03240","map03250","map03260",
                                   "map03262","map03261","map03263","map03264",
                                   "map03265","map03266","map03268","map03267","map03269"),
                               "Information processing in viruses")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map02010","map02060","map03070"),
                               "Membrane transport")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map02020","map04010","map04013","map04016",
                                   "map04011","map04012","map04014","map04015",
                                   "map04310","map04330","map04340","map04341",
                                   "map04350","map04390","map04391","map04392",
                                   "map04370","map04371","map04630","map04064",
                                   "map04668","map04066","map04068","map04020",
                                   "map04070","map04072","map04071","map04024",
                                   "map04022","map04151","map04152","map04150","map04075"),
                               "Signal transduction")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04080","map04060","map04061","map04512","map04514"),
                               "Signaling molecules and interaction")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04144","map04145","map04142","map04146",
                                   "map04140","map04138","map04136","map04137","map04139"),
                               "Transport and catabolism")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04110","map04111","map04112","map04113",
                                   "map04114","map04210","map04214","map04215",
                                   "map04216","map04217","map04115","map04218"),
                               "Cell growth and death")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04510","map04520","map04530","map04540","map04550"),
                               "Cellular community - eukaryotes")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map02024","map05111","map02025","map02026"),
                               "Cellular community - prokaryotes")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map02030","map02040","map04814","map04810"),
                               "Cell motility")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04640","map04610","map04611","map04613",
                                   "map04620","map04624","map04621","map04622",
                                   "map04623","map04625","map04650","map04612",
                                   "map04660","map04658","map04659","map04657","map04662",
                                   "map04664","map04666","map04670","map04672","map04062"),
                               "Immune system")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04911","map04910","map04922","map04923",
                                   "map04920","map03320","map04929","map04912",
                                   "map04913","map04915","map04914","map04917",
                                   "map04921","map04926","map04935","map04918",
                                   "map04919","map04928","map04916","map04924",
                                   "map04614","map04925","map04927"),
                               "Endocrine system")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04260","map04261","map04270"),
                               "Circulatory system")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04970","map04971","map04972","map04976","map04973",
                                   "map04974","map04975","map04979","map04977","map04978"),
                               "Digestive system")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04962","map04960","map04961","map04964","map04966"),
                               "Excretory system")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04724","map04727","map04725","map04728","map04726",
                                   "map04720","map04730","map04723","map04721","map04722"),
                               "Nervous system")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04744","map04745","map04740","map04742","map04750"),
                               "Sensory system")

master_table$module <- replace(master_table$module, master_table$module %in% c("map04320","map04360","map04361","map04380"),
                               "Development and regeneration")

master_table$module <- replace(master_table$module, master_table$module %in% c("map04211","map04212","map04213"),
                               "Aging")

master_table$module <- replace(master_table$module, master_table$module %in% c("map04710","map04713","map04711","map04712","map04714","map04626"),
                               "Environmental adaptation")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05200","map05202","map05206","map05205",
                                   "map05204","map05207","map05208","map05203",
                                   "map05230","map05231","map05235"),
                               "Cancer: overview")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05210","map05212","map05225","map05226",
                                   "map05214","map05216","map05221","map05220",
                                   "map05217","map05218","map05211","map05219",
                                   "map05215","map05213","map05224","map05222","map05223"),
                               "Cancer: specific types")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05166","map05170","map05161","map05160",
                                   "map05171","map05164","map05162","map05168",
                                   "map05163","map05167","map05169","map05165"),
                               "Infectious disease: viral")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05110","map05120","map05130","map05132",
                                   "map05131","map05135","map05133","map05134",
                                   "map05150","map05152","map05100"),
                               "Infectious disease: bacterial")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05146","map05144","map05145","map05140",
                                   "map05142","map05143"),
                               "Infectious disease: parasitic")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05310","map05322","map05323","map05320",
                                   "map05321","map05330","map05332","map05340"),
                               "Immune disease")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05010","map05012","map05014","map05016",
                                   "map05017","map05020","map05022"),
                               "Neurodegenerative disease")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05030","map05031","map05032","map05033","map05034"),
                               "Substance dependence")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map05417","map05418","map05410","map05412",
                                   "map05414","map05415","map05416"),
                               "Cardiovascular disease")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map04930","map04940","map04950","map04936",
                                   "map04932","map04931","map04933","map04934"),
                               "Endocrine and metabolic disease")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map01501","map01502","map01503"),
                               "Drug resistance: antimicrobial")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map01521","map01524","map01523","map01522"),
                               "Drug resistance: antineoplastic")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07011","map07012","map07013","map07021",
                                   "map07019","map07020","map07014","map07023",
                                   "map07026","map07044","map07053"),
                               "Chronology: Antiinfectives")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07040","map07041","map07042","map07043", "map07045"),
                               "Chronology: Antineoplastics")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07032","map07030","map07033","map07015",
                                   "map07039","map07028","map07029","map07031",
                                   "map07027","map07056","map07057"),
                               "Chronology: Nervous system agents")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07055","map07016","map07017","map07018","map07037",
                                   "map07038","map07046","map07047","map07048","map07049",
                                   "map07050","map07051","map07052","map07054"),
                               "Chronology: Other drugs")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07220","map07215","map07214","map07213","map07212",
                                   "map07227","map07211","map07228","map07224","map07229"),
                               "Target-based classification: G protein-coupled receptors")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07225","map07226","map07223","map07222"),
                               "Target-based classification: Nuclear receptors")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07221","map07230","map07036",
                                   "map07231","map07232","map07235"),
                               "Target-based classification: Ion channels")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07233","map07234"),
                               "Target-based classification: Transporters")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07216","map07219","map07024","map07217","map07218"),
                               "Target-based classification: Enzymes")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07025","map07034","map07035"),
                               "Structure-based classification")

master_table$module <- replace(master_table$module, master_table$module %in% 
                                 c("map07110","map07112","map07114","map07117"),
                               "Skeleton-based classification")


#################### Re-assigning pathways ###########################
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K03367", "K03739","K03740", "K14188" ,"K14198", "K01401", "K14205"),
                                "Cationic antimicrobial peptide (CAMP) resistance")
master_table$pathway <- replace(master_table$pathway,master_table$V2 %in% c("K00688", "K00693"),"Starch and sucrose metabolism")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K03841"),"Glycolysis / Gluconeogenesis")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K11262"),"Fatty acid biosynthesis")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K02991"),"Ribosome")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01897"),"Fatty acid metabolism")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K07359"),"AMPK signaling pathway")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01759"),"Pyruvate metabolism")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K13288"),"RNA degradatio")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01580"),"Alanine aspartate and glutamate metabolism")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K03364"),"Cell Cycle")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K01078"),"Thiamine metabolism")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K04079"),"Protein processing in endoplasmic reticulum")
master_table$pathway <- replace(master_table$pathway, master_table$V2 %in% c("K00863"),"Fructose and mannose metabolism")


#################### Re-assigning modules ###########################
master_table <- replace(master_table, is.na(master_table), "Unassigned")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K03367", "K03739","K03740", "K14188" ,"K14198", "K01401", "K14205"),
                               "Drug resistance: antimicrobial")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K00688", "K03841", "K00693", "K01759", "K00863"),"Carbohydrate Metabolism")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K11262", "K01897"),"Lipid metabolism")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K02991"),"Translation")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K13288", "K04079"),"Folding, sorting and degradation")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K01580"),"Amino acid metabolism")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K07359"),"Signal transduction")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K01078"),"Metabolism of cofactors and vitamins")
master_table$module <- replace(master_table$module, master_table$V2 %in% c("K03364"),"Cell growth and death")


# save the table as a table
# then we can call it easy
save(master_table, file='Rdata/master_table.Rda')

