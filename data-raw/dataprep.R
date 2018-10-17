live <- read.csv('C:/Users/kowalewski/Dropbox (UFL)/R files/Final_All_Taxa_BetaDiv.csv', head=TRUE)
dead <- read.csv('C:/Users/kowalewski/Dropbox (UFL)/R files/Final_Dead.csv', head=TRUE)
refs <- read.csv('C:/Users/kowalewski/Dropbox (UFL)/R files/Refs_local.csv')
taxa <- read.csv('C:/Users/kowalewski/Dropbox (UFL)/R files/Refs.csv')
rownames(live) <- paste('locality',live[,1])
rownames(dead) <- paste('locality',live[,1])
live <- live[,-1]
dead <- dead[,-1]
keep <- which(refs$Habitat %in% c(4, 5))
live <- live[keep,]
dead <- dead[keep,]
refs <- refs[keep,]
refs <- droplevels(refs)
keept <- which((colSums(live) + colSums(dead))>0)
live <- as.matrix(live[,keept])
dead <- as.matrix(dead[,keept])
taxa <- taxa[keept,]
Habitat <- as.character(refs$Habitat)
Habitat[Habitat=='4'] <- 'nearshore'
Habitat[Habitat=='5'] <- 'offshore'
Habitat <- as.factor(Habitat)
Fp <- as.character(taxa$Preservation)
Fp[Fp %in% c(1,2)] <- 'fragile'
Fp[Fp %in% c(3,4)] <- 'durable'
Fp <- as.factor(Fp)
colnames(live) <- paste('species',1:ncol(live))
colnames(dead) <- paste('species',1:ncol(live))
FidData <- list(live=as.matrix(live), dead=as.matrix(dead), habitat=Habitat, fossiltype=Fp)
devtools::use_data(FidData)

