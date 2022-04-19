## Mo Chl-TP model
## Commented and cleaned 4.18.2022

tss.explore <- function(df1,  runmod = T) {

    df1$tss <- as.numeric(as.character(df1$tss))
    df1$dtp <- as.numeric(as.character(df1$dtp))
    df1$nvss <- as.numeric(as.character(df1$nvss))
    df1$tp <- as.numeric(as.character(df1$tp))
    df1$chl <- as.numeric(as.character(df1$chl))
    df1$srp <- as.numeric(as.character(df1$srp))

    date0 <- strptime(paste(df1$month, df1$day, "2004", sep = "-"),
                      format = "%m-%d-%Y")
    df1$yday <- date0$yday

    ## set up month of sampling as a factor
    df1$yday.q <- df1$month
    df1$yday.q <- factor(df1$yday.q)

    ## drop missing data
    incvec <- ! is.na(df1$vss) & ! is.na(df1$chl) & ! is.na(df1$tp) &
        ! is.na(df1$nvss)
    df1 <- df1[incvec,]

    df1$lake <- factor(df1$lake)
    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    ## drop one big outlier
    incvec <- log(df1$tp - df1$dtp) < 0
    df1 <- df1[!incvec, ]

    ## center all data
    varlist<- c("tss", "chl", "tp")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    save(mn.val, file = "mn.val.mo.rda")

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$dtp <- df1$dtp/mn.val["tp"]
    df1$srp <- df1$srp/mn.val["tp"]

    ## stan model
    modstan <- '
        data {
            int n;               // number of samples
            int nlake;           // number of lakes
            int lakenum[n];      // identifier for lake
            int nseas;           // number of months
            int seasnum[n];      // identifier for number of seasons
            vector[n] tp;        // total P
            vector[n] dtp;       // dissolved P
            vector[n] chl;       // Chl-a
            vector[n] tss;       // total suspended solids
        }
        parameters {
            real muu;                // grand mean of non phytoplankton sediment
            real<lower = 0> sigu;
            vector[n] etau;

            real mub;                // grand mean of TSS per unit of Chl
            vector[nseas] etab;
            real<lower = 0> sigb;

            real k[3];               // model exponents

            vector[2] mud;             // grand mean of TP model coefficients
            real<lower = 0> sigd[2];
            vector[nseas] etad1;
            vector[nlake] etad2;

            real<lower = 0> sigtss;     // measurement error for TSS
            real<lower = 0> sigtp;      // measurement error for TP

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] tss_mn;
            vector[nseas] d1;          // season specific coefficients
            vector[nlake] d2;
            vector[nseas] b;

            b = mub + etab*sigb;
            u = muu + etau*sigu;

            d1 = mud[1] + etad1*sigd[1];
            d2 = mud[2] + etad2*sigd[2];

            for (i in 1:n) {
               tss_mn[i] = log_sum_exp(b[seasnum[i]] + k[1]*chl[i], u[i]);

               // configuration of the TP model specified here.
               // mud[1] as the coefficient for chl term indicates that
               // only one coefficient is estimate for all P-Chl relationships
               // d2[seasnum[i]] for the non-phytoplankton P term indicates
               // that d2 for different lakes are estimated
               tp_mn[i] = log_sum_exp(mud[1]+k[2]*chl[i],
                                      d2[lakenum[i]]+k[3]*u[i]);

            }
        }
        model {
            // non-informative priors
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            mub ~ normal(0,3);
            etab ~ normal(0,1);
            sigb ~ cauchy(0,3);

            mud ~ normal(0,3);
            sigd ~ cauchy(0,3);
            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);

            k[1] ~ normal(1,1);     // exponents should be close to 1 so set
            k[2] ~ normal(1,1);     // a wide prior around 1
            k[3] ~ normal(1,1);

            sigtp ~ cauchy(0,3);
            sigtss ~ normal(0.1, 0.02);   // prior for tss sampling error

            // assume tss and tp observations are t-distributed with
            // 4 degrees of freedom to make the regressions robust.
            tss ~ student_t(4,tss_mn, sigtss);
            tp ~ student_t(4,tp_mn, sigtp);
        }
    '
    ## function to compute RMS error
    rmsout <- function(x,y) sqrt(sum((x-y)^2, na.rm = T)/
                                     min(sum(! is.na(x)), sum(! is.na(y))))

    ## function to compute predicted TP from model
    gettp <- function(df, varout.loc) {

        k <- apply(varout.loc$k, 2, mean)
        mud <- apply(varout.loc$mud, 2, mean)
        d1 <- apply(varout.loc$d1, 2, mean)
        d2 <- apply(varout.loc$d2, 2, mean)
        mub <- mean(varout.loc$mub)
        b <- apply(varout.loc$b, 2, mean)
        u <- exp(apply(varout.loc$u, 2, mean))

        ## change formula below to match model specification
        tppred <- exp(mud[1])*df$chl^k[2] + df$dtp +
            exp(d2[df$lakenum])*u^k[3]

        return(list(tppred = tppred))
    }
    extractvars <- c("k", "u", "d1", "d2", "sigd", "sigtp", "sigtss","mud",
                     "mub", "sigb", "b")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
        options(mc.cores = nchains)

        datstan <- list(n = nrow(df1),
                        nlake = max(df1$lakenum),lakenum = df1$lakenum,
                        nseas = max(df1$seasnum), seasnum = df1$seasnum,
                        nvss = df1$nvss,
                        tp = log(df1$tp - df1$dtp),
                        dtp = df1$dtp,
                        vss = df1$vss,
                        chl = log(df1$chl),
                        tss = log(df1$tss))
        print(str(datstan))

        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1800, chains = nchains,
                    warmup = 600, thin= 2,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))

        varout <- extract(fit, pars = extractvars)

        tp.pred <- gettp(df =df1, varout.loc = varout)[[1]]
        dev.new()
        plot(log(tp.pred), log(df1$tp), xlab = "Predicted centered TP",
             ylab = "Observed centered TP")
        abline(0,1)
        cat("RMS:", rmsout(log(tp.pred), log(df1$tp)), "\n")
        return(varout)
    }

}

varout.mo.test <- tss.explore(moi3.all, runmod = T)
