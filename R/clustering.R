#' Calculates a self-organizing map
#'
#' @param data A numerical matrix or data.frame.
#' @param sizemultiplier A size multiplier.
#' @param num_epochs How long the training should last.
#' @param train Row indices of training data.
#' @param seed A random seed (default NULL).
#' @param initrand Whether initialization should be random (not recommended).
#' @param intoroidal If a toroidal map should be created.
#'
#' @return A som object.
#'
#' @export
calcSOM <- function(data, sizemultiplier = 1, num_epochs = 200, train = NULL, seed = NULL, initrand = FALSE, intoroidal = FALSE) {
    if (num_epochs < length(train)/25) {
        lower_bound = floor(length(train)/25)
        warning(paste(num_epochs, "epochs is low for", length(train), "training genes, it was set to",
            lower_bound))
        num_epochs = lower_bound
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(train)) {
        train = 1:min(nrow(datamat), 5000)
    }
    data.train = data[train, ]

    init.map = mapinit(data.train, sizemultiplier, intoroidal)
    # Handle output of mapinit

    maxr = min(0.5 * init.map$h, init.map$w)
    test.som <- kohonen::som(data = data.train, grid = somgrid(init.map$w, init.map$h, "hexagonal"), rlen = num_epochs,
        radius = c(maxr, 1), init = init.map$initgrid, toroidal = intoroidal)

    par(mfrow = c(2, 2))
    plot(test.som, type = "changes", main = "t")
    plot(test.som, type = "count")
    plot(test.som, type = "dist.neighbours")
    plot(test.som, main = "t")
    par(mfrow = c(1, 1))

    return(test.som)
}

#' Initializes a SOM map with two PCA components
#'
#' @param datamat A numerical matrix or data.frame.
#' @param prefactor A factor.
#' @param toroidal If toroidal.
#'
#' @return A list with intialization parameters.
mapinit <- function(datamat, prefactor = 1, toroidal = FALSE) {
    munits <- 5 * prefactor * sqrt(length(datamat[, 1]))
    pcar <- prcomp(datamat)
    ev <- pcar$sdev^2
    if (ev[1] > 5 * ev[2]) {
        h <- round(sqrt(munits * 5))
    } else {
        h <- round(sqrt(munits * ev[1]/ev[2]))
    }
    w <- round(munits/h)
    initcodes <- array(0, dim = c(h, w, length(datamat[1, ])))
    gridshape <- array(0, dim = c(h * w, length(datamat[1, ])))
    for (x in 1:w) {
        for (y in 1:h) {
            if (toroidal) {
                initcodes[y, x, ] = sin((pi * x)/w) * pcar$rotation[, 2] + sin((pi * y)/h) * pcar$rotation[,
                  1]
            } else {
                initcodes[y, x, ] = (x/w) * pcar$rotation[, 2] + (y/h) * pcar$rotation[, 1]
            }
            gridshape[(y - 1) * w + x, ] = initcodes[y, x, ]
        }
    }
    list(h = h, w = w, initgrid = gridshape)
}
