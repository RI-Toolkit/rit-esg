#' esg_var_simulator
#'
#' @description
#' The discrete-time economic scenario generator simulates the trajectories of
#' 11 Australian economic and financial variables (in brackets are the `$names`
#' of the output dataframes):
#'
#' (1) 3-month zero-coupon bond yields (`$zcp3m_yield`),
#'
#' (2) 10-year zero-coupon bond spread (`$zcp10y_spread`),
#'
#' (3) NSW home value index (`$home_index`),
#'
#' (4) NSW home rental yields (`$rental_yield`),
#'
#' (5) Australia GDP (`$GDP`),
#'
#' (6) Australia CPI (`$CPI`),
#'
#' (7) S&P/ASX200 closing price (`$ASX200`),
#'
#' (8) Australian dollar trade-weighted index (`$AUD`),
#'
#' (9) Australia mortgage rate (`$mortgage_rate`),
#'
#' (10) NSW unemployment rate (`$unemployment_rate`),
#'
#' (11) Stochastic discount factors (`$discount_factors`).
#'
#' Simulations are based on a Vector Autoregression model, refer to section
#' `details` for explanations of the mathematical model. This function uses
#' the package `zoo` to convert the frequnency units. Period-by-period summary
#' statistics can be obtained from \code{esg_summary}.
#'
#' @details
#' Factors (1)-(8) (in rates) were fitted using a Vector Autoregressive model
#' (VAR), factors (9)-(10) were respectively expressed as a fixed margin over
#' factors (1)-(2) due to strong correlations, while factor (11) is derived
#' from the fitted VAR model with arbitrage-free assumptions. Further details
#' on model parameter estimation and forecasts can be found in note (a) below.
#'
#' Vector Autoregression (VAR) is a regression of a time series where the ouput
#' depends linearly on the past values of itself, and the past values of other
#' variables, up to some specfied order:
#'
#' \deqn{\mathbf{z}_t = \mathbf{\mu} + \Phi_1 \mathbf{z}_{t-1} + \Phi_2 \mathbf{z}_{t-2} + \cdots + \Phi_p \mathbf{z}_{t-p} + \mathbf{\epsilon},}
#'
#' where
#'
#' * \eqn{\mathbf{z}_{t}} is the vector of economic variables,
#'
#' * \eqn{\mathbf{\mu}} is the vector of intercepts,
#'
#' * \eqn{\Phi_{i}, i=1,\cdots,p} are coefficient matrices of size \eqn{n \times n}
#' with \eqn{n} being the number of economic variables and \eqn{p} the lags.
#'
#' * \eqn{\mathbf{\epsilon}} is a vector of white noises.
#'
#' The stochastic discount factor is defined as:
#'
#' \deqn{\mathbf{s}_{t+1} = \exp(-\mathbf{e}_1^\top \mathbf{z}_t - \frac{1}{2} \lambda_t^\top \lambda_t - \lambda_t^\top \mathbf{\epsilon}_{t+1}),}
#'
#' where \eqn{\mathbf{e}_1^\top \mathbf{z}_t} and \eqn{\mathbf{\epsilon}_t}
#' respectively denote the 3-month zero-coupon bond rates and white noises from
#' the fitted VAR model, and \eqn{\mathbf{\lambda}_t} is the market price of
#' risk process which is assumed to be affine over factors (1)-(8).
#'
#'**Notes:**
#'
#'(a) Procedure for parameter estimation:
#'
#'* Transformed factors (3),(5)-(8) to continuously compounded growth rates,
#'tested for correlation, causality, and stationarity.
#'
#'* Found the optimal lag order for Vector Autoregression using AIC, SIC, HQC.
#'
#'* Fitted the VAR model using ordinary least squares. This was followed by
#'evaluation.
#'
#'* Associated the stochastic discount factors with VAR factors and nominal bond
#'prices, estimated the market price of risk parameters by minimising the sum
#'of squared error.
#'
#'     Procedure for forecasts: Generated factors (1)-(8) using Vector
#' Autoregression formula. From the generated paths and noises, derived
#' stochastic discount factors.
#'
#'     Detailed R codes for parameter estimation can be found in the economic tutorial/economic.
#'
#' (b) Large values of percentage change can appear if the original forecasts
#' are near-zero, or if the Gaussian noise is large, though with low probabilities.
#' This happens especially for interest rates in the first few periods due to
#' historical-low rates in 2021.
#'
#'(c) Negative values for rate factors i.e., factors (1)(2)(4)(9)(10), are
#'theoretically allowed as Vector Autoregression models assume that the noise
#'follow a Gaussian distribution. Index factors, i.e., factors (3)(5)-(8), on
#'the other hand, are all positive.
#'
#'(d) Choosing between discrete- and continuous-time models:
#'
#'* The outputs are different for the two simulators, users should choose the
#'model based on their objectives.
#'
#'* The base time step for discrete-time model is one quarter, whereas there is
#'no such as a base for continuous-time. For time steps smaller than one quarter,
#'discrete-time model will interpolate the quarterly statistics, whereas the
#'continuous-time model simply generates random noises for each specific time
#'step. Consequently, for large time steps, the executing time for continuous-time
#'models are shorter than the dicrete-time model.
#'
#' @references
#'
#'Daniel H Alai, Hua Chen, Daniel Cho, Katja Hanewald, and Michael Sherris. Developing equity release markets: Risk analysis for reverse mortgages and home reversions. _North American Actuarial Journal_, 18(1):217–241, 2014.
#'
#'Andrew Ang and Monika Piazzesi. A no-arbitrage vector autoregression of term structure dynamics with macroeconomic and latent variables. _Journal of Monetary economics_, 50(4):745–787, 2003.
#'
#'Andrew Ang, Monika Piazzesi, and Min Wei. What does the yield curve tell us about gdp growth? _Journal of econometrics_, 131(1-2):359–403, 2006.
#'
#' @param num_years integer denoting number of years to forecast from 01-01-2021,
#' default 5 years
#' @param num_paths integer denoting number of simulations to make for each variable,
#' default 10 paths'year', 'quarter' or 'month' denoting the simulation
#' frequency, default 'quarter'. The base simulation time step is one quarter,
#' linear interpolation will be used if the required frequency is higher,
#' whereas arithmetic average will be used if the frequency is lower.
#' @param perc_change set TRUE for outputs to be expressed as period-by-period
#' percent change, default FALSE. The reference level, i.e., the original values
#' in the first output period, will be appended above the percentage changes for
#' each variable and each trajectory. See note (b) in section `details`.
#' @param return_sdf set TRUE to return the stochastic discount factors,
#' default TRUE
#' @param seed Specify the seed for simulations, no default
#'
#' @return A list of 10 dataframes containing simulated trajectories from
#' 01-01-2021 of the 10 variables, or a list of 11 dataframes including the
#' simulated stochastic discount factors if `return_sdf` is set TRUE. Rows are
#' the trajectories (e.g., `trajectory_1`), columns are the time steps (e.g.,
#' `2021-01-01`). See note (c) for explanations on the negativity of output values.
#' @export esg_var_simulator
#'
#' @examples # simulate 10 years of data
#'
#' sim <- esg_var_simulator(num_years = 10, num_paths = 10, frequency = 'year')
#'
#' # suppose we wish to look at the 3 months zero coupon bonds
#' sim$zcp3m_yield
#'
#' # if we wanted a specific trajectory, say 3
#' sim$zcp3m_yield[3,]
esg_var_simulator = function (num_years = 5, num_paths = 10, frequency = "quarter", perc_change = FALSE, return_sdf = TRUE, seed = NULL) {

    ################
    # error messages
    is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (num_years <= 0 | num_paths <= 0 | !is.wholenumber(num_years)
        | !is.wholenumber(num_paths)) {
        stop("Number of years and paths to simulate must be positivie integers. ")

    } else if (!frequency %in% c("year", "quarter", "month")) {
        stop ("Frequency must be one of 'year', 'quarter', and 'month'. ")

    } else if (!is.logical(perc_change) | !is.logical(return_sdf)) {
        stop ("perc_change and return_sdf must be logical. ")

    }

    ##########################################################
    # VAR(2) calibrated coefficients (for stationary series) #
    ##########################################################

    # variable names
    var_names = c("zcp3m_yield", "zcp10y_spread", "home_index", "rental_yield", "GDP", "CPI", "ASX200", "AUD")
    sim_var_names = c(var_names, "mortgage_rate", "unemployment_rate")

    VAR = var_model()
    {
        intercept = VAR$intercept
        coef = VAR$coef
        covres = VAR$covres
        init_stat_2020q4 = VAR$init_stat_2020q4
        init_stat_2021q1 = VAR$init_stat_2021q1
        init_orig = VAR$init_orig
        mortgage_rate = VAR$mortgage_rate
        unemployment_rate = VAR$unemployment_rate
        lambda_0 = VAR$lambda_0
        lambda_1 = VAR$lambda_1
        init_lambdat = VAR$init_lambdat
        init_st = VAR$init_st
    }

    # intercept in VAR(1)
    # coefficient matrix in VAR(1)
    # residual covariance matrix in VAR(1)
    colnames(coef) = var_names
    colnames(covres) = var_names

    coef = as.matrix(t(coef))
    coef1 = coef
    coef1[1,1] = coef1[1,1] + 1
    coef1[4,4] = coef1[4,4] + 1
    coef2 = matrix(0,nrow = length(var_names), ncol = length(var_names))
    coef2[,1] = -coef[,1]
    coef2[,4] = -coef[,4]

    ##################
    # initialisation #
    ##################

    # starting quarter
    init_qtr = as.Date("2021-01-01")
    num_pred = 4 * num_years
    time_index = seq(from = init_qtr, length.out = num_pred + 1, by = "quarter")
    path_index = paste("trajectory_", 1:num_paths, sep = "")
    progression = floor(num_paths / 5)

    ############################
    # step-by-step simulations #
    ############################

    # white noise
    set.seed(seed)
    noise = matrix(data = stats::rnorm(length(intercept) * num_pred * num_paths, 0, 1),
                   nrow = length(intercept))
    noise = lapply(seq(from = 1, to = num_paths * num_pred, by = num_pred),
                   function (x) {noise[, x:(x+num_pred-1)]})
    noise = lapply(noise,
                   function(x) {row.names(x) = var_names;
                                colnames(x) = as.character(time_index[-1]);
                                return (x)})
    names(noise) = path_index

    # whole path (inputs/outputs are both stationary)
    var_path = function (num_pred, noise_index) {
        path = as.data.frame(matrix(NA, nrow = num_pred, ncol = length(var_names)))
        row.names(path) = time_index[-1]
        colnames(path) = var_names

        # simulate for num_pred steps
        new_init = init_stat_2021q1 # z_{t-1}
        old_init = init_stat_2020q4 # z_{t-2}

        for (i in 1:num_pred) {
            e = as.vector(noise[[noise_index]][,i])
            zt = intercept + coef1 %*% new_init + coef2 %*% old_init + as.matrix(chol(covres)) %*% e
            path[i,] = zt
            old_init = new_init # z_{t-2} <- z_{t-1}
            new_init = zt # z_{t-1} <- z_t

        }
        return (path)
    }

    ########################################
    # simulation for the stationary series #
    ########################################

    prog_ind = 1; cat("Progress: 0% \n")
    var_sim_stationary = function (num_pred, num_paths) {

        # loops thru the series (separate lists)
        v_path = replicate(n = num_paths,
                            expr = {data.frame(matrix(NA, nrow = num_pred, ncol = length(intercept)))},
                            simplify = F)
        v_path = lapply(1:num_paths,
                        function (x) {
                            if (x == progression * prog_ind && prog_ind < 5) {
                                cat(paste(20*prog_ind, "%\n", sep = ""))
                                prog_ind <<- prog_ind + 1
                            }
                            var_path(num_pred, x)

                        })
        return (v_path)
    }
    stat = var_sim_stationary(num_pred, num_paths)
    stat = lapply(stat, function (x) {cbind(x[,1:2]/100, x[,3:8])})

    ################################################
    # convert forecast variables -> original units #
    ################################################

    index2grow_inv = function (x, init) {
        # reverse index to growth rate: home_value, gdp, cpi, asx200, aud
        Reduce (function (init, x) {init * exp(x)}, c(init, x), accumulate = TRUE)
    }


    # reorganise the results
    sim = replicate(n = length(var_names),
                    expr = {data.frame(matrix(NA, nrow = num_pred, ncol = num_paths))},
                    simplify = F)
    sim = lapply(1:8, function (y) { lapply(1:num_paths, function (x) {stat[[x]][,y]}) })
    sim = lapply(sim, function (x) {as.data.frame(x)})
    sim[[1]] = rbind(init_orig[1]/100, as.data.frame(sim[[1]])) # zcp3m
    sim[[2]] = rbind(init_orig[2]/100, as.data.frame(sim[[2]])) # zcp10y
    sim[[3]] = apply(sim[[3]], 2, function (x) {index2grow_inv(x, init_orig[3])}) # home_index
    sim[[4]] = rbind(init_orig[4], as.data.frame(sim[[4]])) # rental
    sim[[5]] = apply(sim[[5]], 2, function (x) {index2grow_inv(x, init_orig[5])}) # GDP
    sim[[6]] = apply(sim[[6]], 2, function (x) {index2grow_inv(x, init_orig[6])}) # CPI
    sim[[7]] = apply(sim[[7]], 2, function (x) {index2grow_inv(x, init_orig[7])}) # ASX200
    sim[[8]] = apply(sim[[8]], 2, function (x) {index2grow_inv(x, init_orig[8])}) # AUD
    sim[[9]] = sim[[1]] + mortgage_rate # mortage_rate
    sim[[10]] = sim[[2]] + unemployment_rate # unemployment_rate
    sim = lapply(sim, function(x){row.names(x) = time_index; colnames(x) = path_index; return (x)})
    names(sim) = sim_var_names

    ###############################
    # stochastic discount factors #
    ###############################

    if (isTRUE(return_sdf)) {
        # find lambda_t's for different trajectories
        lambda_t = replicate(n = num_paths,
                                expr = {matrix(NA, ncol = num_pred+1, nrow = 8)},
                                simplify = F)
        lambda_t = lapply(1:num_paths,
                          function (x) {lambda_t[[x]] = sapply(1:num_pred,
                          function (y) {lambda_t[[x]][,y] = lambda_0 + lambda_1 %*%
                                                            t(as.matrix(stat[[x]][y,])) })})
        lambda_t = lapply(lambda_t,
                          function (x) {x = cbind(init_lambdat, as.data.frame(x));
                                        row.names(x) = var_names; colnames(x) = as.character(time_index);
                                        x = x[,-ncol(x)]; return (x)})

        # find s_t for different trajectories
        st = as.data.frame(matrix(NA, nrow = num_pred, ncol = num_paths))
        st_expn = function (time,path) {
            # finds s(t+1)
            exp(- stat[[path]][time,1] - 1/2 * sum(lambda_t[[path]][,time]^2)
                - sum(lambda_t[[path]][,time] * noise[[path]][,time]))
        }
        st = sapply(1:num_paths,
                    function (x) {sapply(1:num_pred,
                    function (y) {st[y,x] = st_expn(y,x)}, simplify = TRUE)},
                    simplify = TRUE)
        st = rbind(init_st,st)

        st[2,] = ifelse(st[2,] > 1.2,1.2,ifelse(st[2,] < 0.8, 0.8, st[2,])) # trim the irregular values: historical data: interest rate < 0
        row.names(st) = as.character(time_index)
        colnames(st) = path_index

        ########
        # output
        sim[[length(sim_var_names) + 1]] = st
        names(sim)[length(sim_var_names) + 1] = "discount_factors"
    }

    #################
    # Adj frequency #
    #################

    output = list()
    if (frequency == "month") {
        time_index_month = seq(from = init_qtr, length.out = num_years * 12 + 1, by = "month")

        qtr2month = function (x) {
            # transforms quarterly data to monthly data
            qtr_data = zoo::zoo (x, time_index)
            month_data = zoo::zoo (NA, time_index_month)
            data = merge (qtr_data, month_data)
            data$month_data = zoo::na.approx(data$qtr_data, rule=12)
            return (as.vector(data$month_data))
        }
        output = lapply(sim,
                        function (x) apply(x, 2, qtr2month))
        output = lapply(output,
                        function(x) {row.names(x) = as.character(time_index_month);
                                     return (x[-nrow(x), ])})

    } else if (frequency == "quarter") {
        output = lapply(sim, function(x) {x = x[-nrow(x), ]}) # remove the last row (1 Jan)

    } else if (frequency == "year") {
        time_index_year = seq(from = init_qtr, length.out = num_years, by = "year")
        output = lapply(sim,
                        function (x) apply(x[-nrow(x), ], 2,
                        function (y) {colMeans(matrix(y, nrow=4))} ))
        output = lapply(output,
                        function(x) {row.names(x) = as.character(time_index_year); x})
    }
    output = lapply(output, function(x){x = t(as.data.frame(x))})

    #############
    # Adj units #
    #############

    if (isTRUE(perc_change)) {
        ref_level = lapply(output, function (x) {x = as.data.frame(x[,1]); colnames(x) = paste("ref_level", time_index[1]);x})
        output = lapply(output, function (x) {(x[,-1] - x[,-ncol(x)]) / x[,-ncol(x)]})
        output = lapply(1:length(output), function (x) {cbind(ref_level[[x]],output[[x]])}) # include the reference level in outputs
        names(output) = names(sim)
    }

    cat("100% \n")
    return (output)
}
