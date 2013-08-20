#players <- read.csv("player_data.csv")
#elections <- read.csv("election_data.csv")
HOFdat <- read.csv("HOFvotingdata.csv", stringsAsFactors=FALSE)
HOFdat$P <- as.numeric(HOFdat$position == "P")
#we only consider elections after 1967
election_dat <- subset(HOFdat, Year >= 1967)
#disregard players that get no votes?
player_dat <- subset(election_dat, Votes != 0)
#format years on ballot
player_dat$YoB <- as.integer(gsub("[a-z]+", "", player_dat$YoB))
#grab columns of interest and sort them by player name and years on ballot
election_vars <- c("Year", "Votes", "YoB", "NumBallots")
stats <- c("WAR", "P")
standard <- c("WAR") #variables that need to be standardized
datAll <- player_dat[c("Name", election_vars, stats)]
datAll <- plyr::arrange(datAll, Name, YoB)
#break up variables into their inherit levels of info
election_dat <- datAll[c("Name", election_vars)]
election_dat$prop <- with(election_dat, Votes/NumBallots)
player_dat <- unique(datAll[c("Name", stats)])
#standardize the covariates
for (i in standard) player_dat[,i] <- scale(player_dat[,i])
#rescale time to [0, 1]
election_dat$YoB <- (election_dat$YoB-1)/14
n.players <- dim(player_dat)[1]
n.predictors <- dim(player_dat)[2]