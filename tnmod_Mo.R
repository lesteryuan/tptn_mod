## 4.9.2022: R script for Mo TN model

tnmod <- function(df1, runmod = T) {

    print(nrow(df1))

    nperiod <- 6
    df1$yday.q <- df1$month
    df1$yday.q <- factor(df1$yday.q)
    print(table(df1$yday.q))

    ## drop missings
    incvec <- ! is.na(df1$don) & ! is.na(df1$chl) & ! is.na(df1$din) &
        ! is.na(df1$tn) & ! is.na(df1$doc) & ! is.na(df1$vss)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    ## center variables
    varlist<- c("tn", "chl", "doc", "tss")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    save(mn.val, file = "mn.val.tnmo.rda")

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$din <- df1$din/mn.val["tn"]
    df1$don <- df1$don/mn.val["tn"]
    df1$dtn <- df1$dtn/mn.val["tn"]
    df1$vss <- df1$vss/mn.val["tss"]

    ## drop pn measurements near zero or negative
    pn <- df1$tn - df1$din - df1$don
    incvec <- pn < 1e-10
    print(sum(incvec))

    df1 <- df1[!incvec,]
    print(nrow(df1))

    modstan <- '
        data {
            int n;
            int nlake;
            int lakenum[n];
            int nseas;
            int seasnum[n];
            vector[n] tn;                // particulate N
            vector[n] chl;               // Chl
            vector[n] tss;               // VSS
        }
        parameters {
            real k[3];
            vector[2] mud;
            real<lower = 0> sigd[2];
            vector[nlake] etad1;
            vector[nseas] etad2;

            real<lower = 0> sigtn;

            real muu;
            real<lower =0> sigu;
            vector[n] etau;

            real mub;
            real<lower = 0> sigb;
            vector[nseas] etab;
            real<lower = 0> sigtss;

        }
        transformed parameters {
            vector[n] tn_mn;
            vector[n] u;
            vector[n] tss_mn;

            vector[nseas] b;
            vector[nlake] d1;
            vector[nseas] d2;

            b = mub + sigb*etab;

            d1 = mud[1] + etad1*sigd[1];
            d2 = mud[2] + etad2*sigd[2];

            u = muu + sigu*etau;

            for (i in 1:n) {
                tss_mn[i] = log_sum_exp(b[seasnum[i]] + k[3]*chl[i], u[i]);

                tn_mn[i] = log_sum_exp(mud[1] + k[1]*chl[i],
                        d2[seasnum[i]] + k[2]*u[i]);
            }
        }
        model {
            mud ~ normal(0,3);
            sigd ~ cauchy(0,3);
            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);

            k[3] ~ normal(0.85,0.01);
            k[1] ~ normal(1,1);
            k[2] ~ normal(1,1);
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);
            sigtss ~ normal(0.1,0.002);

            mub ~ normal(0,3);
            etab ~ normal(0,1);
            sigb ~ cauchy(0,3);

            sigtn ~ cauchy(0,3);

            tss ~ student_t(4,tss_mn, sigtss);
            tn ~ normal(tn_mn, sigtn);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2, na.rm = T)/
                                     min(sum(! is.na(x)), sum(!is.na(y))))
    gettn <- function(df, varout) {
        k <- apply(varout$k,2,mean)
        mud <- apply(varout$mud, 2, mean)
        mub <- mean(varout$mub)
        d1 <- apply(varout$d1, 2, mean)
        d2 <- apply(varout$d2, 2, mean)

        u <- exp(apply(varout$u, 2, mean))
        tnpred <- exp(mud[1])*df$chl^k[1] +
            exp(d2[df$lakenum])*u^k[2]
        return(tnpred)
    }
    extractvars <- c("mud", "k", "u", "mub", "d1", "d2", "sigd")

    nchains <- 3
    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)

        options(mc.cores = nchains)

        datstan <- list(n = nrow(df1),
                        nlake = max(df1$lakenum),lakenum = df1$lakenum,
                        nseas = max(df1$seasnum), seasnum = df1$seasnum,
                        tn = log(df1$tn-df1$dtn),
                        chl = log(df1$chl),
                        tss = log(df1$vss))
        print(str(datstan))

        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1800, chains = nchains,
                    warmup = 600, thin= 1,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))
        save(fit, file = "fitout.rda")
        varout <- extract(fit, pars = extractvars)

        tnpred <- gettn(df1, varout)
        dev.new()
        plot(log(tnpred), log(df1$tn-df1$dtn))
        abline(0,1)
        cat("Internal RMS:", rmsout(log(tnpred), log(df1$tn-df1$dtn)), "\n")

        return(varout)
    }

}
varout.test<- tnmod(moi3.all, runmod = T)
