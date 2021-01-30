# cviceni 03
rm(list = ls())
library(tidyverse)
source("Support/my_NLRM.R") # nacteni my_NLRM funkce do Environment
source("Support/norm_rnd.R") # nacteni norm_rnd funkce do Environment
source("Support/gamm_rnd_Koop.R") # nacteni gamm_rnd_Koop funkce do Environment
# Nacteni dat
data <- read_delim("capm2_data.csv",
                    ";", escape_double = FALSE,
                    locale = locale(date_format = "%YM%m",
                                    decimal_mark = ",", grouping_mark = "."),
                    trim_ws = TRUE)

# Vykreslseni MSFT
ggplot(data = data, aes(x = obs, y = MSFT)) +
    geom_line(color = "darkblue", size = 0.75) +
    theme_bw()


# Priprava promennych
y <- data$GM - data$RKFREE
X <- cbind(rep(1,nrow(data)),(data$MKT - data$RKFREE))

# Nastaveni apriornich hyperparametru
beta_0 <- c(0,1)
cov_beta_0 <- diag(c(0.05^2,0.5^2))
s2_0 <- 0.2^2
h_0 <- 1/s2_0
nu_0 <- 10
V_0 <- cov_beta_0*(nu_0-2)/nu_0*h_0

# Odhad NLRM s prirozene konjugovanou apriorni hustotou
res_GM <- my_NLRM(y,X,beta_0,V_0,h_0,nu_0)

# 1. Prezentace vysledku

cat('\nParametr    Prior   Prior std.    Posterior   Posterior std.\n',
    '============================================================\n',
    'Alpha:     ',res_GM$beta_0[1],'     ',res_GM$b0_std[1],'    ',res_GM$beta_1[1],'    ',res_GM$b1_std[1],'\n',
    'Beta:      ',res_GM$beta_0[2],'     ',res_GM$b0_std[2],'    ',res_GM$beta_1[2],'    ',res_GM$b1_std[2],'\n',
    'h:         ',res_GM$h_0,'     ',res_GM$h0_std,'    ',res_GM$h_1,'    ',res_GM$h1_std)
# 2. Test hypotezy beta = 1
y <- (data$GM - data$RKFREE) - (data$MKT - data$RKFREE)
X <- rep(1,nrow(data))
beta0_res <- beta_0[1]
V0_res <- V_0[1,1]
res_GM_rest <- my_NLRM(y,X,beta0_res,V0_res,h_0,nu_0)

# Prezentace vysledku pro omezeny model beta = 1
cat('\n\nOmezeny model beta = 1',
    '\nParametr    Prior   Prior std.    Posterior   Posterior std.\n',
    '============================================================\n',
    'Alpha:     ',res_GM_rest$beta_0,'     ',res_GM_rest$b0_std,'    ',res_GM_rest$beta_1,'    ',res_GM_rest$b1_std,'\n',
    'h:         ',res_GM_rest$h_0,'     ',res_GM_rest$h0_std,'    ',res_GM_rest$h_1,'    ',res_GM_rest$h1_std)

log_BF = res_GM_rest$log_ML - res_GM$log_ML;
BF = exp(log_BF);

cat('\n\nBayesuv faktor porovnavajici omezeny a neomezeny model\n',
    'Pr = ',BF)

# 3. Hypoteza beta > 1
MC <- 100000
beta_sim <- matrix(0,nrow = 2,ncol = MC)

for (i in 1:MC) {
    h_sim <- gamm_rnd_Koop(res_GM$h_1,res_GM$nu_1,1)
    beta_sim[,i] <- norm_rnd(h_sim^(-1)*res_GM$V_1) + res_GM$beta_1
}

# Vypocet pravdepodobnosti beta > 1
pr_beta <- sum(beta_sim[2,] > 1)/MC
cat('\n\nPravdepodobnost beta > 1',
    '\nPr = ',pr_beta)

zscore <- (1 - res_GM$beta_1[2])/res_GM$b1_std[2]
pr_beta_analyticky = 1 - pt(zscore,res_GM$nu_1)
cat('\nPr = ',pr_beta_analyticky,' (analyticky)')
