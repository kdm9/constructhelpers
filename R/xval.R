#' csh.x.validation
#'
#' @param n.reps            Number of replicate runs/mcmc chains
#' @param K                 Vector of layers/ancestral pops. (e.g. 1:6)
#' @param train.prop        Proportion of SNPs to use for training vs test. typically 0.8-0.9
#' @param freqs             Per-locality allele frequencies, see conStruct::conStruct()
#' @param data.partitions   Pre-computed data partitions, see conStruct::x.validation()
#' @param geoDist           Geographic distance matrix, see conStruct::conStruct()
#' @param coords            X/Y or Long/Lat coordiate matrix, see conStruct::conStruct()
#' @param prefix            File prefix for any output files, see conStruct::x.validation()
#' @param n.iter            Number of MCMC interations
#' @param make.figs         Save plots to disk
#' @param save.files        Save intermediates to disk
#'
#' @return tibble with one row per MCMC run (so length(K) * n.reps * 2)
#' @export
csh.x.validation <- function(freqs = NULL, data.partitions = NULL, geoDist, coords, prefix, K, n.reps=8,  n.iter=10000, 
                             train.prop = 0.9, make.figs = FALSE, save.files = FALSE,
                              ...) {
    #######################################
    #  Check args (from check.xval.call)  #
    #######################################
    # these two are needed to ensure that the check....() functions work.
    n.nodes=2
    parallel=T
    args = as.list(environment())
    conStruct:::check.for.files(args)
    conStruct:::check.genetic.data.arg(args)
    args$spatial <- TRUE
    conStruct:::check.geoDist.arg(args)
    conStruct:::check.coords.arg(args)

    #####################
    #  Make partitions  #
    #####################
    if (is.null(data.partitions)) {
        data.partitions = conStruct:::make.data.partitions(n.reps, freqs, train.prop)
    }
    conStruct:::check.data.partitions.arg(args = as.list(environment()))

    ###############
    #  Main loop  #
    ###############
    x.val = foreach::foreach(rep.no=1:n.reps, .combine=bind_rows) %:%
        foreach::foreach(k=K, .combine=bind_rows) %:%
        foreach::foreach(mdl=c("sp", "nsp"), .combine=bind_rows) %dopar% {
            
            # Prep data structures
            spatial = mdl == "sp"
            dat = data.partitions[[rep.no]]
    
            # Run construct mcmc chain
            cs = conStruct:::xval.conStruct(
                spatial = spatial,
                K = k, 
                data = dat$training, 
                geoDist = geoDist,
                coords = coords, 
                prefix = paste0(prefix, "_", mdl, "_", "rep", rep.no, "K", k), 
                n.iter = n.iter,
                make.figs = make.figs,
                save.files = save.files)

            # calculate fit
            fit = unlist(conStruct:::fit.to.test(dat$testing, cs[[1]]))
            mean.fit = mean(fit, na.rm=T)

            # re-create data block (downstream stuff needs it)
            dblock = conStruct:::xval.make.data.block(k, dat$training, coords, spatial, geoDist)
            dblock = conStruct:::unstandardize.distances(dblock)

            cat(paste0("mdl=", mdl, " K=",k, " rep=", rep.no, " Done\n"))
            tibble_row(mdl=mdl, K=k, rep=rep.no,
                       construct.res=cs, mdl.fit=list(fit), mean.fit=mean.fit,
                       data.part=list(dat), data.block = dblock)
    }
    return(x.val)
}
