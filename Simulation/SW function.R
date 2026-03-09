createS <- function(n, p,
                    topology = "identity",
                    dataset = FALSE,
                    precision = FALSE,
                    nonzero = 0.25,
                    m = 1L,
                    banded.n = 2L,
                    invwishart = FALSE,
                    nu = p + 1,
                    Plist) {
        ##############################################################################
        # - Simulate one or more random symmetric square (or datasets) matrices from
        #   various models.
        # - n          > A vector of sample sizes
        # - p          > An integer giving the dimension. p should be greater than
        #                or equal to 2.
        # - topology   > character. Specify the topology to simulate data from.
        #                See details.
        # - dataset    > logical. Should dataset instead of its sample covariance be
        #                returned? Default is FALSE.
        # - precision  > logical. Should the constructed precision matrix be returned?
        # - nonzero    > numeric of length 1 giving the value of the non-zero entries
        #                for some topologies
        # - m          > integer. The number of conditionally independent subgraphs.
        #                I.e. the number of blocks.
        # - banded.n   > interger. The number of off-diagonal bands used if
        #                topology is "banded". Use as paremeter if topology is
        #                "Watt-Strogatz".
        # - invwishart > logical. If TRUE the constructed precision matrix is
        #                used as the scale matrix of an inverse Wishart distribution
        #                and class covariance matrices are drawn from this
        #                distribution.
        # - nu         > The degrees of freedom in the inverse wishart distribution.
        #                A large nu implies high class homogeneity.
        #                Must be greater than p + 1.
        # - Plist      > A list of user-supplied precision matrices. Should be the
        #                same length as n.
        #
        # NOTES:
        # - Returns a list of matrices if length(n) > 1. The output is simplified if
        #   n has length 1 where only the matrix is returned
        # EXAMPLE:
        # - G1 <- createS(n=10, p=100, topology="small-world", m=5, banded.n=3, precision=T)
        # - Theta1 <- (G1!=0)*1
        ##############################################################################
        
        if (missing(p) && !missing(Plist)) {
                p <- nrow(Plist[[1]])
        }
        stopifnot(p > 1)
        stopifnot(m >= 1)
        G <- length(n)
        
        if (dataset && precision) {
                stop("dataset and precision cannot be TRUE at the same time.")
        }
        
        if (invwishart && missing(nu)) {
                stop("argument 'nu' is missing. Supply the degrees of freedom 'nu' for ",
                     "the inverse Wishart distribution.")
        }
        
        topology <- match.arg(topology,
                              c("identity", "star", "clique", "complete",
                                "chain", "banded", "Barabassi", "small-world",
                                "scale-free", "Watts-Strogatz", "random-graph",
                                "Erdos-Renyi"))
        
        # Construct the precision matrix "constructor"
        if (topology == "identity") {
                
                submat <- function(p) diag(p)
                
        } else if (topology == "star") {
                
                submat <- function(p) {
                        subP <- diag(p)
                        subP[1, seq(2, p)] <- subP[seq(2, p), 1] <- 1/seq(2, p)
                        return(subP)
                }
                
        } else if (topology == "chain") {
                
                submat <- function(p) {
                        s <- seq_len(p - 1)
                        subP <- diag(p)
                        subP[cbind(s, s + 1)] <- subP[cbind(s + 1, s)] <- nonzero
                        return(subP)
                }
                
        } else if (topology == "clique" || topology == "complete") {
                
                submat <- function(p) {
                        subP <- diag(p)
                        subP[lower.tri(subP)] <- subP[upper.tri(subP)]  <- nonzero
                        return(subP)
                }
                
        } else if (topology == "banded") {
                
                submat <- function(p) {
                        if (banded.n > p) {
                                stop("The number of bands cannot exceed the dimension of each block")
                        }
                        subP <- diag(p)
                        for (j in seq(1, banded.n)) {
                                s <- seq_len(p - j)
                                subP[cbind(s, s + j)] <- subP[cbind(s + j, s)] <- 1/(j + 1)
                        }
                        return(subP)
                }
                
        } else if (topology == "Barabassi" || topology == "scale-free") {
                
                submat <- function(p) {
                        G <- barabasi.game(p, power = 1, directed = FALSE)
                        adj <- get.adjacency(G, sparse = FALSE)
                        return(diag(p) + nonzero*adj)
                }
                
        } else if (topology == "Watts-Strogatz" || topology == "small-world") {
                
                submat <- function(p) {
                        G <- watts.strogatz.game(1, p, banded.n, 0.05)
                        adj <- get.adjacency(G, sparse = FALSE)
                        return(diag(p) + nonzero*adj)
                }
                
        } else if (topology == "Erdos-Renyi" || topology == "random-graph") {
                
                submat <- function(p) {
                        G <- erdos.renyi.game(p, 1/p)
                        adj <- get.adjacency(G, sparse = FALSE)
                        return(diag(p) + nonzero*adj)
                }
                
        } else {
                
                stop("requested topology not implemented yet.")
                
        }
        
        # Construct the block split
        blocks <- split(seq_len(p), ceiling(m*seq_len(p)/p))
        
        # Fill in blocks to construct full precisions
        P <- diag(p)
        for (b in blocks) {
                P[b, b] <- submat(length(b))
        }
        if (rcond(P) < sqrt(.Machine$double.eps)) {  # Check condition number
                warning("The generated precision matrix has a very high condition number ",
                        "and the generated data might be unreliable.")
        }
        S <- solve(P)
        
        # Construct names
        n.letters <- which.max(p <= 26^(1:3))
        x <- expand.grid(rep(list(LETTERS), n.letters))
        nms <- do.call(paste0, x)
        
        # Construct list to fill and iterate over all classes
        ans <- vector("list", G)
        names(ans) <- paste0("class", seq_len(G))
        for (g in seq_len(G)) {
                
                if (!missing(Plist)) {
                        stopifnot(length(Plist) == length(n))
                        stopifnot(nrow(Plist[[g]]) == ncol(Plist[[g]]))
                        stopifnot(nrow(Plist[[g]]) == p)
                        
                        Sg <- solve(Plist[[g]])
                } else if (invwishart) {
                        stopifnot(nu - p - 1 > 0)
                        Sg <- drop(.armaRInvWishart(n = 1, psi = (nu - p - 1)*S, nu = nu))
                } else {
                        Sg <- S
                }
                
                if (precision) {
                        
                        if (invwishart) {
                                ans[[g]] <- solve(Sg)
                        } else {
                                ans[[g]] <- P
                        }
                        
                } else {
                        
                        ans[[g]] <- rmvnormal(n = n[g], mu = rep(0, p), sigma = Sg)
                        
                        if (!dataset) {
                                ans[[g]] <- covML(ans[[g]])
                        }
                        
                }
                
                if (p <= 17576) {  # Only give names for "small" dimensions
                        colnames(ans[[g]]) <- nms[1:p]
                        if (!dataset) {
                                rownames(ans[[g]]) <- nms[1:p]
                        }
                }
        }
        
        if (G == 1) {  # Simplify output if ns is length 1
                ans <- ans[[1]]
        }
        
        return(ans)
}
