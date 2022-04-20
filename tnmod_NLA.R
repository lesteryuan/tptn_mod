## TN model using NLA data

tn.model <- function(df1, runmod = F) {
    require(rstan)
    nchains <- 3    # number of chains


    ## omit records that are missing data
    incvec <- ! is.na(df1$ntl.result) & ! is.na(df1$no3no2.result) &
        ! is.na(df1$chl) & ! is.na(df1$doc.result) & !is.na(df1$us.l3code)

    cat("N omitted due to missing records:", sum(!incvec), "\n")
    df1 <- df1[incvec,]
    print(nrow(df1))

    ## drop chl = 0
    incvec <- df1$chl > 0
    df1 <- df1[incvec,]

    ## compute TN-DIN and drop values that are <= 0
    df1$tkn <- df1$ntl.result - df1$no3no2.result
    incvec <- df1$tkn <= 0
    cat("TKN <= 0:", sum(incvec), "\n")
    df1 <- df1[!incvec,]

    ## select index sites
    incvec <- df1$sample.type == "MICX"
    df1 <- df1[incvec,]

    ## reset state and L3 ecoregion factors
    df1$state <- factor(df1$st.nla2012)
    df1$statenum <- as.numeric(df1$state)
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)

    ## drop samples with chl outside of Mo range
    incvec <- df1$chl < 108
    df1 <- df1[incvec,]
    incvec <- df1$chl > 1
    df1 <- df1[incvec,]

    ## center chl and doc
    ## use same chl scale as TP model
    load("chlsc.rda")
    df1$chl.sc <- df1$chl/chlsc

    docsc <- exp(mean(log(df1$doc.result)))
    df1$doc.sc <- df1$doc.result/docsc

    ## center tn and nox
    tnsc <- exp(mean(log(df1$ntl.result)))
    df1$tn.sc <- df1$ntl.result/tnsc
    df1$nox.sc <- df1$no3no2.result/tnsc

    save(tnsc, file = "tnsc.rda")

    ## drop 5 outliers in tn doc relationship
    incvec <- (log(df1$tn.sc - df1$nox.sc) - log(df1$doc.sc)) > 3 |
        (log(df1$tn.sc - df1$nox.sc) - log(df1$doc.sc)) < -2
    print(sum(incvec))

    df1 <- df1[!incvec,]

    tnchldat <- df1[, c("chl", "chl.sc", "tkn", "doc.sc", "doc.result",
                        "statenum", "state", "ntl.result", "tn.sc",
                        "no3no2.result", "nox.sc", "index.lat.dd",
                        "index.lon.dd", "econum", "us.l3code")]
    save(tnchldat, docsc,tnsc, file = "tnchldat.rda")

    datstan <- list(n = nrow(df1),
                    neco = max(df1$econum),econum = df1$econum,
                    tn = log(df1$tn.sc - df1$nox.sc),
                    nox = log(df1$nox.sc),
                    doc = log(df1$doc.sc),
                    chl = log(df1$chl.sc))

    print(str(datstan))

    modstan <- '
        data {
            int n;                // number of samples
            int neco;            // number of ecoregions
            int econum[n];       // ecoregion of each sample
            vector[n] tn;        // scaled TN
            vector[n] chl;       // scaled chl
            vector[n] nox;       // scaled NOx
            vector[n] doc;       // scaled DOC
        }
        parameters {
            real muk;                // mean value of exponent on chl

            real mud[2];              // mean value of model coefficients
            real<lower = 0> sigd[2]; // SD of model coefficients among ecoregions
            vector[neco] etad2;
            vector[n] etad2a;

            real<lower = 0> sigtn;  // measurement error of tn
            real muu;
            vector[n] etau;
            real<lower = 0> sigu;

        }
        transformed parameters {

            vector[neco] d2;
            vector[n] d2a;
            vector[n] u;

            d2 = mud[2] + sigd[1]*etad2;

            d2a = d2[econum] + sigd[2]*etad2a;

            u = muu + etau*sigu;
        }
        model {
            matrix[n,3] temp;
            vector[n] tnmean;

            mud ~ normal(0,4);
            muk ~ normal(1,1);

            sigd ~ cauchy(0,4);
            etad2 ~ normal(0,1);
            etad2a ~ normal(0,1);

            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            sigtn ~ normal(0.1, 0.002);

           temp[,1] = mud[1] + muk*chl;
           temp[,2] = d2a + doc;
           temp[,3] = u;
           for (i in 1:n) tnmean[i] = log_sum_exp(temp[i,]);

            tn ~ student_t(4,tnmean, sigtn);
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

## save extracted variables to varout to post-process
fitout <- tn.model(dat.merge.all, runmod = T)

varout.n.limnat <- extract(fitout, pars = c("u","muk", "mud",  "d2", "d2a",
                                       "sigd", "u", "muu"))


