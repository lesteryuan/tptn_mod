## 11.6.2019: Production version of TP-chl model
## 12.17.2019Cleaned and commented
## 1.6.2021: Testing adjustment for dtp
## 1.21.2021: Limit only model for TP
ntumodel <- function(df1, runmod = T) {

    require(rstan)

    nchains <- 3    # select number of chains

    ## omit records that are missing data
    incvec <- ! is.na(df1$ptl.result) & ! is.na(df1$turb.result) &
         ! is.na(df1$chl) & ! is.na(df1$index.site.depth)
    cat("N omitted due to missing records:", sum(!incvec), "\n")
    df1 <- df1[incvec,]

    ## select index sites
    df1 <- subset(df1, sample.type == "MICX")

    ## drop below detection limit turb and tp
    norig <- nrow(df1)
    incvec <- df1$turb.result > 0.01
    df1 <- df1[incvec,]
    incvec <- df1$ptl.result > 1
    df1 <- df1[incvec,]
    incvec <- df1$chl >0
    df1 <- df1[incvec,]
    cat("N dropped for detection limit:", norig - nrow(df1), "\n")

    ## drop samples with chl out of Mo range
    incvec <- df1$chl < 108
    df1 <- df1[incvec,]
    incvec <- df1$chl > 1
    df1 <- df1[incvec,]

    ## scale chl and tp
    chlmn <- mean(log(df1$chl))
    tpsc <- exp(mean(log(df1$ptl.result)))
    chlsc <- exp(chlmn)
    df1$chl.sc <- df1$chl/chlsc
    df1$tp.sc <- df1$ptl.result/tpsc

    ##drop HI, only dealing with conterminous US
    incvec <- df1$st.nla2012 == "HI"
    df1 <- df1[!incvec,]

    ## reset factor levels for state
    df1$state <- factor(df1$st.nla2012)
    df1$statenum <- as.numeric(df1$state)

    ## reset L3 ecoregion codes
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)

    datstan <- list(n = nrow(df1),
                    neco = max(df1$econum), econum = df1$econum,
                    tp = log(df1$tp.sc),
                    chl = log(df1$chl.sc))

    print(str(datstan))

    modstan <- '
        data {
            int n;                 // number of samples
            int neco;              // number of L3 ecoregions
            int econum[n];         // ecoregion assigment
            vector[n] tp;          // TP
            vector[n] chl;         // Scaled Chl
        }
        parameters {
            real muk;  // exponents on chl and u in models
            real mud[2];
            real<lower = 0> sigd[2];// SD of ecoregion- or depth-specific coef
            vector[neco] etad1;
            vector[n] etad1a;

            real<lower = 0> sigtp;    // measurement error of tp

        }
        transformed parameters {
            vector[neco] d1;
            vector[n] d1a;

            d1 = mud[1] + sigd[1]*etad1;
            d1a = d1[econum] + sigd[2]*etad1a;
        }
        model {
            vector[n] tp_mn;

            muk ~ normal(1,1);    // muk should be somewhere around 1
            mud ~ normal(0,4);

            sigd[1] ~ cauchy(0,3);
            sigd[2] ~ cauchy(0,3);

            etad1 ~ normal(0,1);
            etad1a ~ normal(0,1);

            sigtp ~ normal(0.1, 0.002);

           for (i in 1:n) {
               tp_mn[i] = log_sum_exp(d1a[i], mud[2] + muk*chl[i]);
           }

            tp ~ student_t(4,tp_mn,sigtp);
        }
    '

    if (runmod) {

        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 2400, chains = nchains,
                    warmup = 600, thin = 3)
        return(fit)
    }

}

## runmod variable set to T to run simulation and set to F to
##  run post processing.
fitout <- ntumodel(dat.merge.all, runmod = T)
## post processing
varout.p.limnat <- extract(fitout, pars = c("muk", "mud", "sigd", "d1", "d1a"))







