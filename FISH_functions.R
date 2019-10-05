## ---- EDA.histograms
EDA_histograms = function(var='', group='',dat=NULL, var.lookup) {
    pretty.name = as.character(var.lookup$pretty.name[var.lookup$Abbreviation==var])
    if (is.numeric(dat[,var])) {
        grid.arrange(
            ggplot(dat) + geom_histogram(aes(x=!!sym(var))) +
            scale_x_continuous(pretty.name),
            ggplot(dat) + geom_histogram(aes(x=!!sym(var))) +
            scale_x_log10(pretty.name),
            nrow=1
        )
    } else {
        grid.arrange(
            ggplot(dat) + geom_bar(aes(x=!!sym(var), fill=!!sym(group))) +
            scale_x_discrete(pretty.name),
            nrow=1
        )
    }
}
## ----end

## ---- EDA.density
EDA_density = function(var='', group='', dat=NULL, var.lookup) {
    pretty.name = as.character(var.lookup$pretty.name[var.lookup$Abbreviation==var])
    if (is.numeric(dat[,var])) {
        grid.arrange(
            ggplot(dat) + geom_density(aes(x=!!sym(var), fill=!!sym(group)), alpha=0.5) +
            scale_x_continuous(pretty.name),
            ggplot(dat) + geom_density(aes(x=!!sym(var), fill=!!sym(group)), alpha=0.5) +
            scale_x_log10(pretty.name),
            nrow=1
        )
    } else {
    }
}
## ----end

## ---- assignMonotone
assignMonotone = function(data,formula) {
    MF=model.frame(formula, data)
    dataClasses=attr(terms(MF),'dataClasses')[-1]
    MF=MF %>% mutate_if(is.factor, as.numeric)
                                        #MM=model.matrix(formula,data)
    VAR_MONOTONE <- cor(MF[, 1], MF[, -1], method = 'spearman') / abs(cor(MF[, 1], MF[, -1], method = 'spearman'))
    VAR_MONOTONE[dataClasses=='factor'] = 0
    VAR_MONOTONE[is.na(VAR_MONOTONE)]=0
    as.vector(VAR_MONOTONE)
}
## ----end



## ---- plotABT
plot.abt <-
function (x, i.var = 1, npts, trim.grid = FALSE, return.grid = FALSE, 
    center = TRUE, cvalue, cfun = mean, add.zero = FALSE, xlab, 
    ylab, zlab = "Partial Dependence", clip = TRUE, type = c("link", 
        "response")[1], ptype = c("contour", "persp", "lines", 
        "trellis")[1], pcol = "vary", fpcol = heat.colors(50), 
    d = 6, theta = 40, n.trees, add.legend = TRUE, se = 2, se.smooth = c("smooth", 
        "mean", "none")[1], lty = rep(1, 4), lwd = 1, se.fill = TRUE, 
    ylim, xrange = 0.025, col = c("black", "grey80", "dodgerblue3", 
        "chartreuse4", "firebrick3", "purple3", "lightred", "grey40", 
        "grey60", "grey80"), pch = 21, pts = T, rug.add = TRUE, 
    leg.pos = "topright", leg.just = c(0, 1), ...) 
{
    require(lattice)
    if (one.abt <- any(class(x) == "single")) {
        n.trees <- x$n.trees
        se <- FALSE
    }
    else {
        tree.list <- x
        nlist <- length(x)
        x <- tree.list[[1]]
        if (missing(n.trees)) 
            n.trees <- sapply(tree.list, "[[", "n.trees")
        else if (length(n.trees) == 1) 
            n.trees <- rep(n.trees, nlist)
    }
    if (!is.numeric(i.var)) {
        i.var <- unique(match(i.var, x$var.names))
    }
    i.var <- i.var[order(x$var.type[i.var] > 0)]
    if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
        warning("i.var must be between 1 and ", length(x$var.names))
    }
    if (length(i.var) > 3) {
        warning("plot.abt creates up to 3-way interaction plots.\n        plot.abt will only return the plotting data structure.")
        return.grid = TRUE
    }
    se.sm <- c("smooth", "mean", "none")
    se.smooth <- se.sm[pmatch(se.smooth, c("smooth", "mean", 
        "none"))]
    if (all(is.na(se.smooth))) 
        warning("se.smooth not recognised!\n")
    pty <- c("contour", "persp", "lines", "trellis")
    ptype <- pty[pmatch(ptype, pty)]
    if (all(is.na(ptype))) 
        cat("ptype not recognised!\n")
    ty = c("link", "response")
    type <- ty[pmatch(type, ty)]
    if (all(is.na(type))) 
        warning("type not recognised!\n")
    if (length(lty) < 5) 
        lty <- rep(lty, 5)
    if (length(col) < 5) 
        col <- rep(col, 5)
    distrib <- x$distribution$name
    if (missing(npts)) {
        if ((ptype == "contour" | ptype == "persp") & (length(i.var) < 
            3)) 
            npts <- rep(50, length(i.var))
        else npts <- c(50, rep(5, length(i.var) - 1))
    }
    if (length(trim.grid) == 1) 
        trim.grid <- rep(trim.grid, length(i.var))
    grid.levels <- vector("list", length(i.var))
    for (i in 1:length(i.var)) {
        if (is.numeric(x$var.levels[[i.var[i]]])) {
            grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]), 
                max(x$var.levels[[i.var[i]]]), length = npts[i])
            if (trim.grid[i]) 
                grid.levels[[i]] <- grid.levels[[i]][-c(1, npts[i])]
        }
        else {
            grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]], 
                levels = x$var.levels[[i.var[i]]])) - 1
        }
    }
    X <- expand.grid(grid.levels)
    names(X) <- paste("x", 1:length(i.var), sep = "")
    x.num <- (x$var.type == 0)[i.var]
    if (clip && (length(i.var) > 1) && (sum(x.num) == 1)) {
        sub <- function(x, y) x[(x[, 1] > y[1] & x[, 1] < y[2]), 
            , drop = F]
        fv <- (x$var.type[i.var] > 0)
        nv <- (x$var.type[i.var] == 0)
        fact <- x$var.names[i.var][fv]
        numb <- x$var.names[i.var][nv]
        res <- lapply(split(x$mdata[, numb, drop = F], x$mdata[, 
            fact, drop = F]), range, na.rm=TRUE)
        Xs <- split(X, X[, !x.num, drop = F])
        for (i in 1:length(Xs)) Xs[[i]] <- sub(Xs[[i]], res[[i]])
        X <- do.call("rbind", Xs)
        row.names(X) <- 1:nrow(X)
    }
    else if (clip && (length(i.var) == 2) && all(x.num)) {
        xpdbox <- function(x, delta = 0.1) {
            x + delta * t(t(x) - apply(x, 2, mean))
        }
        poly.in <- function(xy, poly) {
            if (ncol(poly) == 2) 
                poly <- poly.tst(poly)
            n <- nrow(xy)
            np <- nrow(poly)
            nnp <- rep(n, np)
            check1 <- xor(xy[, 1] >= rep(poly[, 1], nnp), xy[, 
                1] > rep(poly[, 3], nnp))
            check2 <- rep(poly[, 2], nnp) + rep(poly[, 4], nnp) * 
                xy[, 1] > xy[, 2]
            as.vector(matrix(check1 & check2, n, np) %*% rep.int(1, 
                np)%%2 > 0)
        }
        poly.tst <- function(xy) {
            if (is.list(xy)) 
                xy <- as.data.frame(xy)
            if (!is.matrix(xy) || ncol(xy) != 2 || !is.numeric(xy[, 
                1]) || !is.numeric(xy[, 2])) 
                stop("xy must by nX2 numeric matrix or data.frame")
            xy <- as.matrix(xy)
            poly <- matrix(c(xy, xy[c(2:nrow(xy), 1), ]), ncol = 4)
            poly <- poly[poly[, 1] != poly[, 3], ]
            poly[, 4] <- (poly[, 4] - poly[, 2])/(poly[, 3] - 
                poly[, 1])
            poly[, 2] <- poly[, 2] - poly[, 1] * poly[, 4]
            poly
        }
        XX <- x$mdata[, x$var.names[i.var][1:2]]
        XX <- apply(XX, 2, function(x) {
            x[is.na(x)] <- mean(x, na.rm = T)
            x
        })
        hpts <- chull(XX)
        XXP <- xpdbox(XX[c(hpts, hpts[1]), ])
        ins <- poly.in(as.matrix(X), XXP)
    }
    else if (pts && (length(i.var) > 1) && (sum(x.num) == 2)) {
        XX <- x$mdata[, x$var.names[i.var][1:2]]
    }
    if (one.abt) {
        y <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
            cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
            n.class = as.integer(x$num.classes), i.var = as.integer(i.var - 
                1), n.trees = as.integer(n.trees), initF = as.double(x$initF), 
            trees = x$trees, c.splits = x$c.splits, var.type = as.integer(x$var.type), 
            PACKAGE = "abt")
    }
    else if (!se) {
        ytemp <- 0
        for (i in 1:nlist) {
            ytemp <- ytemp + .Call("gbm_plot", X = as.double(data.matrix(X)), 
                cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
                n.class = as.integer(x$num.classes), i.var = as.integer(i.var - 
                  1), n.trees = as.integer(n.trees[i]), initF = as.double(tree.list[[i]]$initF), 
                trees = tree.list[[i]]$trees, c.splits = tree.list[[i]]$c.splits, 
                var.type = as.integer(x$var.type), PACKAGE = "abt")
        }
        y <- ytemp/nlist
    }
    else {
        ytemp <- matrix(NA, nrow = nrow(X), ncol = nlist)
        for (i in 1:nlist) {
            ytemp[, i] <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
                cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
                n.class = as.integer(x$num.classes), i.var = as.integer(i.var - 
                  1), n.trees = as.integer(n.trees[i]), initF = as.double(tree.list[[i]]$initF), 
                trees = tree.list[[i]]$trees, c.splits = tree.list[[i]]$c.splits, 
                var.type = as.integer(x$var.type), PACKAGE = "abt")
        }
        y <- apply(ytemp, 1, mean, na.rm = T)
        y.se <- apply(ytemp, 1, var, na.rm = T)
        if (length(i.var) == 1 && is.numeric(x$var.levels[[i.var]])) {
            if (se.smooth == "smooth") 
                y.se <- runmed(y.se, k = 1 + 2 * min((npts - 
                  1)%/%2, ceiling(0.2 * npts)))
            else if (se.smooth == "mean") 
                y.se <- rep(mean(y.se), length(y.se))
            if (!one.abt) 
                y.se <- y.se * (nlist - 1)
        }
        y.se <- sqrt(y.se)
    }
    if (x$distribution$name == "multinomial") {
        X$y <- matrix(y, ncol = x$num.classes)
        colnames(X$y) <- x$classes
        if (se) {
            X$se <- matrix(y.se, ncol = x$num.classes)
            colnames(X$se) <- x$classes
            X$hi <- matrix(X$y + se * X$se, ncol = x$num.classes)
            colnames(X$hi) <- x$classes
            X$lo <- matrix(X$y - se * X$se, ncol = x$num.classes)
            colnames(X$lo) <- x$classes
        }
        if (type == "response") {
            X$y <- exp(X$y)
            X$y <- X$y/matrix(rowSums(X$y), ncol = ncol(X$y), 
                nrow = nrow(X$y))
        }
    }
    else {
        X$y <- y
        if (se) {
            X$se <- y.se
            X$hi <- X$y + se * X$se
            X$lo <- X$y - se * X$se
        }
        if (distrib == "bernoulli" && type == "response") {
            X$y <- 1/(1 + exp(-X$y))
            if (se) {
                X$hi <- 1/(1 + exp(-X$hi))
                X$lo <- 1/(1 + exp(-X$lo))
            }
        }
        else if (distrib == "poisson" && type == "response") {
            X$y <- exp(X$y)
            if (se) {
                X$hi <- exp(X$hi)
                X$lo <- exp(X$lo)
            }
        }
    }
    if (center) {
        if (missing(cvalue)) 
            cvalue <- cfun(X$y, na.rm = T)
        X$y <- X$y - cvalue
        if (se) {
            X$hi <- X$hi - cvalue
            X$lo <- X$lo - cvalue
        }
    }
    if (clip && (length(i.var) == 2) && all(x.num)) 
        X$y[!ins] <- NA
    f.factor <- rep(FALSE, length(i.var))
    for (i in 1:length(i.var)) {
        if (!is.numeric(x$var.levels[[i.var[i]]])) {
            X[, i] <- factor(x$var.levels[[i.var[i]]][X[, i] + 
                1], levels = x$var.levels[[i.var[i]]])
            f.factor[i] <- TRUE
        }
    }
    if (return.grid) {
        names(X)[1:length(i.var)] <- x$var.names[i.var]
        if (center) 
            attr(X, "cent") <- cvalue
        return(X)
    }
    nrows <- nrow(x$mdata)
    if (missing(ylim)) {
        if (se) 
            ylim <- c(min(X$lo, na.rm = T), max(X$hi, na.rm = T))
        else ylim <- range(X$y, na.rm = T)
        if (xrange) {
            rnge <- diff(ylim)
            ylim[1] <- ylim[1] - xrange * rnge
            ylim[2] <- ylim[2] + xrange * rnge
        }
    }
    if (length(i.var) == 1) {
        if (missing(xlab)) 
            xlab <- x$var.names[i.var]
        if (missing(ylab)) 
            ylab <- paste("f(", x$var.names[i.var], ")", sep = "")
        if (!f.factor) {
            if (x$distribution$name == "multinomial") {
                if (type == "response") {
                  ylabel <- "Predicted class probability"
                }
                else {
                  ylabel <- paste("f(", x$var.names[i.var], ")", 
                    sep = "")
                }
                plot(range(X$x1), range(X$y), type = "n", xlab = x$var.names[i.var], 
                  ylab = ylabel)
                for (ii in 1:x$num.classes) {
                  lines(X$x1, X$y[, ii], xlab = x$var.names[i.var], 
                    ylab = paste("f(", x$var.names[i.var], ")", 
                      sep = ""), col = ii, ...)
                }
            }
            else {
                plot(X$x1, X$y, type = "n", ylim = ylim, xlab = xlab, 
                  ylab = ylab, ...)
                if (se) {
                  if (!se.fill) {
                    lines(X$x1, X$hi, col = col[2], lty = lty[2], 
                      lwd = lwd)
                    lines(X$x1, X$lo, col = col[2], lty = lty[2], 
                      lwd = lwd)
                  }
                  else polygon(c(X$x1, rev(X$x1)), c(X$hi, rev(X$lo)), 
                    border = NA, col = col[2])
                }
                if (add.zero) 
                  abline(h = 0, lty = lty[2], col = col[2])
                lines(X$x1, X$y, col = col[1], ...)
                if (rug.add) 
                  rug(x$var.levels[[i.var]], ticksize = 0.05, 
                    side = 1, lwd = 1)
            }
        }
        else {
            if (x$distribution$name == "multinomial") {
                nX <- length(X$x1)
                dim.y <- dim(X$y)
                if (type == "response") {
                  ylabel <- "Predicted probability"
                }
                else {
                  ylabel <- paste("f(", x$var.names[i.var], ")", 
                    sep = "")
                }
                plot(c(0, nX), range(X$y), axes = FALSE, type = "n", 
                  xlab = x$var.names[i.var], ylab = ylabel)
                axis(side = 1, labels = FALSE, at = 0:nX)
                axis(side = 2)
                mtext(as.character(X$x1), side = 1, at = 1:nX - 
                  0.5, cex = par()$cex)
                segments(rep(1:nX - 0.75, each = dim.y[2]), as.vector(t(X$y)), 
                  rep(1:nX - 0.25, each = dim.y[2]), as.vector(t(X$y)), 
                  col = 1:dim.y[2])
                box()
            }
            else {
                plot(X$x1, X$y, ylim = ylim, xpd = FALSE, type = "n", 
                  axes = T, border = NA)
                box()
                if (se) {
                  len <- length(X$x1)
                  arrows(1:len, X$lo, 1:len, X$hi, length = par()$fin/len/4, 
                    angle = 90, code = 3, col = col[1])
                }
                points(X$x1, X$y, col = col[1], cex = 1.3 * par("cex"), 
                  bg = col[2], pch = pch, ...)
                mtext(x$var.names[i.var[1]], side = 1, line = par()$mgp[1], 
                  cex = par("cex") * par("cex.lab"))
                mtext(paste("f(", x$var.names[i.var[1]], ")", 
                  sep = ""), side = 2, las = 3, line = par()$mgp[1], 
                  cex = par("cex") * par("cex.lab"))
                box()
            }
        }
    }
    else if (length(i.var) == 2) {
        if (!f.factor[1] && !f.factor[2]) {
            cp <- (ptype[1] == "persp") | (ptype[1] == "contour")
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[1]]
            if (missing(ylab) & cp) 
                ylab <- x$var.names[i.var[2]]
            if (missing(ylab) & !cp) 
                ylab <- paste("f(", x$var.names[i.var[1]], ",", 
                  x$var.names[i.var[2]], ")", sep = "")
            if (ptype[1] == "contour") {
                z <- matrix(X$y, ncol = length(grid.levels[[2]]))
                image(grid.levels[[1]], grid.levels[[2]], z, 
                  xlab = xlab, ylab = ylab, zlim = range(z, na.rm = T), 
                  xaxs = "r", yaxs = "r", col = fpcol, ...)
                contour(grid.levels[[1]], grid.levels[[2]], z, 
                  zlim = range(z, na.rm = T), add = T, col = col[2], 
                  ...)
                if (pts) 
                  points(XX)
            }
            else if (ptype[1] == "persp") {
                z <- matrix(X$y, ncol = length(grid.levels[[2]]))
                if (pcol == "vary") 
                  pcol <- gcols(z, FUNCOL = fpcol)
                persp(grid.levels[[1]], grid.levels[[2]], z, 
                  zlab = zlab, xlab = xlab, ylab = ylab, zlim = range(z, 
                    na.rm = T), tick = "detailed", col = pcol, 
                  theta = theta, d = d, ...)
            }
            else if (ptype[1] == "lines") {
                X.lst <- split(X, X$x2)
                plot(X$x1, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x1, 
                  X.lst[[i]]$y, col = col[i], lty = lty[i], ...)
                if (add.legend) {
                  legend(leg.pos, as.character(signif(unique(sort(X$x2)), 
                    2)), lty = lty, xjust = leg.just[1], yjust = leg.just[2], 
                    col = col, cex = 0.8, bty = "n", title = x$var.names[i.var[2]])
                }
            }
            else if (ptype[1] == "trellis") 
                if (!se) {
                  print(xyplot(y ~ x1 | x2, data = X, xlab = xlab, 
                    ylab = ylab, type = "l", ylim = ylim, col = col[1], 
                    ...))
                }
                else {
                  X$yhi <- X$y + se * X$se
                  X$ylo <- X$y - se * X$se
                  print(xyplot(y + yhi + ylo ~ x1 | x2, data = X, 
                    type = "l", allow.mult = T, xlab = xlab, 
                    ylab = ylab, col = c(col[1], col[2], col[2]), 
                    lty = lty, ylim = ylim, ...))
                }
        }
        else if (f.factor[1] && !f.factor[2]) {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[2]]
            if (missing(ylab)) 
                ylab <- paste("f(", x$var.names[i.var[2]], "|", 
                  x$var.names[i.var[1]], ")", sep = "")
            if (ptype[1] == "lines") {
                X.lst <- split(X, X$x1)
                plot(X$x2, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x2, 
                  X.lst[[i]]$y, col = col[i], lty = lty[i], ...)
                if (add.legend) 
                  legend(leg.pos, levels(X$x1), lty = lty, col = col, 
                    cex = 0.8, xjust = leg.just[1], yjust = leg.just[2], 
                    bty = "n")
            }
            else if (!se) 
                print(xyplot(y ~ x2 | x1, data = X, xlab = xlab, 
                  ylab = ylab, type = "l", col = col[1], ylim = ylim, 
                  ...))
            else {
                X$yhi <- X$y + se * X$se
                X$ylo <- X$y - se * X$se
                print(xyplot(y + yhi + ylo ~ x2 | x1, data = X, 
                  type = "l", allow.mult = T, xlab = xlab, ylab = ylab, 
                  col = c(col[1], col[2], col[2]), lty = 1, ylim = ylim, 
                  ...))
            }
        }
        else if (!f.factor[1] && f.factor[2]) {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[1]]
            if (missing(ylab)) 
                ylab <- paste("f(", x$var.names[i.var[1]], "|", 
                  x$var.names[i.var[2]], ")", sep = "")
            if (ptype[1] == "lines") {
                X.lst <- split(X, X$x2)
                plot(X$x1, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x1, 
                  X.lst[[i]]$y, col = col[i], lty = lty[i])
                if (add.legend) 
                  legend(leg.pos, levels(X$x2), lty = lty, col = col, 
                    cex = 0.8, xjust = leg.just[1], yjust = leg.just[2], 
                    bty = "n")
            }
            else if (!se) 
                print(xyplot(y ~ x1 | x2, data = X, xlab = xlab, 
                  ylab = ylab, type = "l", col = col[1], ylim = ylim, 
                  ...))
            else {
                X$yhi <- X$y + se * X$se
                X$ylo <- X$y - se * X$se
                print(xyplot(y + yhi + ylo ~ x1 | x2, data = X, 
                  type = "l", allow.mult = T, xlab = xlab, ylab = ylab, 
                  col = c(col[1], col[2], col[2]), lty = lty, 
                  ylim = ylim, ...))
            }
        }
        else {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[i[1]]]
            if (missing(ylab)) 
                ylab <- paste("f(", paste(x$var.names[i.var[1:3]], 
                  collapse = ","), ")", sep = "")
            if (ptype[1] == "trellis") 
                barchart(x1 ~ y | x2, data = X, ylab = paste("f(~", 
                  x$var.names[i.var[1]], "*", x$var.names[i.var[2]], 
                  ")", sep = ""), orig = 0, ...)
            else {
                if (!se) {
                  barplot(X$y, names = outer(substring(x$var.levels[[i.var[1]]], 
                    1, 2), substring(x$var.levels[[i.var[2]]], 
                    1, 2), FUN = paste, sep = "\n"), ylim = ylim, 
                    xlab = xlab, ylab = ylab, col = col, xpd = FALSE, 
                    ...)
                  abline(h = 0)
                  box()
                }
                else {
                  yerr <- X$y + se * X$se * ifelse(X$y < 0, -1, 
                    1)
                  z <- barplot(X$y, names = outer(substring(x$var.levels[[i.var[1]]], 
                    1, 2), substring(x$var.levels[[i.var[2]]], 
                    1, 2), FUN = paste, sep = "\n"), ylim = ylim, 
                    xlab = xlab, ylab = ylab, col = rainbow(n = length(unlist(x$var.levels[i.var[1]]))), 
                    xpd = FALSE, ...)
                  abline(h = 0)
                  segments(z, X$y, z, yerr)
                  zerr <- 0.35 * (z[2] - z[1])
                  segments(z + zerr, yerr, z - zerr, yerr)
                  box()
                }
            }
        }
    }
    else if (length(i.var) == 3) {
        if (missing(xlab)) 
            xlab <- x$var.names[i.var[1]]
        if (missing(ylab)) 
            ylab <- paste("f(", paste(x$var.names[i.var[1:3]], 
                collapse = ","), ")", sep = "")
        if (sum(f.factor) < 3) {
            if (!se) {
                print(xyplot(y ~ x1 | x2 * x3, data = X, xlab = xlab, 
                  ylab = ylab, type = "l", ylim = ylim, ...))
            }
            else {
                X$yhi <- X$y + se * X$se
                X$ylo <- X$y - se * X$se
                print(xyplot(y + yhi + ylo ~ x1 | x2 * x3, data = X, 
                  type = "l", allow.mult = T, xlab = xlab, col = c(1, 
                    2, 2), lty = lty, ylab = ylab, ylim = ylim, 
                  ...))
            }
        }
        else if (sum(f.factor) == 3) {
            barchart(x1 ~ y | x2 * x3, data = X, xlab = xlab, 
                ylab = ylab, origin = 0, ...)
        }
    }
    if (center) 
        invisible(cvalue)
}

## ----end


## ---- findThresholds
findThreshold = function(x, y, tol=0.5, deriv=1, type=c('low','high','max','average')) {
  ffun=splinefun(x, y)
  fval=ffun(x, deriv=deriv)
  (rr = range(fval))
  wch=which(abs(rr)==max(abs(rr)))
  if (wch==1) {
    wch=which(fval<rr[1]*tol)
  } else {
    wch=which(fval>rr[2]*tol)
  }
  if (type=='high') {
    wch2 =which(y[wch]==max(y[wch]))
    return(x[wch][wch2])
  } else if (type=='low'){
    wch2 =which(y[wch]==min(y[wch]))
    return(x[wch][wch2])
  } else if (type=='max') {
    wch2 =which(fval[wch]==max(fval[wch]))
    return(x[wch][wch2])
  } else {
    return(mean(x[wch]))
  }
}
## ----end

## ---- commonLegend
common_legend = function(gg) {
    g_legend<-function(p1){
        tmp <- ggplot_gtable(ggplot_build(p1))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)}
    legend <- g_legend(gg[[2]] + theme(legend.position='right'))
    gg = lapply(gg, function(x) x + theme(legend.position='none'))
    gg[['Legend']] = legend
    gg
}
## ----end

## ---- plot.abts
plot.abts = function(mod, var.lookup,    center=FALSE, return.grid=TRUE,type='response', trans=NULL,groupby=NULL, ylab=NULL,manual.trans.y) {
    N=length(mod)
    mdata=mod[[1]]$mdata
    pdata=mod[[1]]$mdata[,-1]
    preds = colnames(pdata)
    p = list()
    thresholds=vector('list',length(preds))
    names(thresholds) <- preds
    Range = c(1e10,1e-05) # initialize range with too high min and too low max
    for (i in 1:length(preds)) {
        VAR = var.lookup[var.lookup$Abbreviation==preds[i],]
        pred=preds[i]
        xlab = as.character(VAR$pretty.name)
        
        grouped=FALSE
        if (is.null(groupby)) {
            newdata = plot.abt(mod, c(i), npts=100,center=center, type=type, return.grid=TRUE) %>%
                mutate(X = !!sym(preds[i]))
        } else if(preds[i]==groupby) {
            grouped=TRUE
            pretty.groupby=as.character(var.lookup[var.lookup$Abbreviation==groupby,'pretty.name'])
            newdata = plot.abt(mod, c(i), npts=100, center=center, type=type, return.grid=TRUE) %>%
                mutate_(.dots=setNames(preds[i], 'X')) %>%
                mutate_(.dots=setNames(groupby, 'group'))
        } else {
            grouped=TRUE
            pretty.groupby=as.character(var.lookup[var.lookup$Abbreviation==groupby,'pretty.name'])
            newdata = plot.abt(mod, c(which(preds==groupby), i), npts=100, center=center, type=type, return.grid=TRUE) %>%
                mutate_(.dots=setNames(preds[i], 'X')) %>%
                mutate_(.dots=setNames(groupby, 'group'))
        }

        resp = ifelse(length(mod[[1]]$response.name)==1, mod[[1]]$response.name, mod[[1]]$response.name[2])
        trans = as.character(var.lookup[var.lookup$Abbreviation==resp,'Transform'])
        trans = ifelse(trans=='log', 'exp', trans)
        newdata = newdata %>% mutate_at(vars(y,lo,hi), list(trans))
        
        if (is.numeric(newdata$X)) {
            if (grouped) {
                p[[i]] <- ggplot(newdata, aes(y=y, x=X, group=group, group=group, color=group, fill=group)) +
                    geom_ribbon(aes(ymin=lo, ymax=hi), color=NA,alpha=0.3) +
                    scale_fill_discrete(pretty.groupby) +
                    scale_color_discrete(pretty.groupby) + 
                    geom_line() +
                                        #scale_y_continuous(paste0('f(',parse(text=xlab),')')) +
                    scale_x_continuous(xlab) +
                    theme_classic()
                if (groupby=='NTR') {
                    p[[i]] = p[[i]] + scale_fill_manual(pretty.groupby, values=c('blue','lightgreen','darkgreen')) 
                    p[[i]] = p[[i]] + scale_color_manual(pretty.groupby, values=c('blue','lightgreen','darkgreen'))     
                }
                
                thresholds[[preds[i]]] = newdata %>% group_by_at(groupby) %>%
                    summarize(threshold.low=mean(findThreshold(x=X, y=y, deriv=1, type='low')),
                              threshold.high=mean(findThreshold(x=X, y=y, deriv=1, type='high')),
                              threshold.max=mean(findThreshold(x=X, y=y, deriv=1, type='max')),
                              threshold.av=mean(findThreshold(x=X, y=y, deriv=1, type='average'))) %>%
                    as.data.frame
                                        #print(paste('group=',unique(newdata$group)))
                                        #print(paste('groupby=',groupby))
                                        #print(paste('labels=',labels[[groupby]]))
            } else {
                thresholds[[preds[i]]] = c(threshold.low=mean(findThreshold(x=X, y=y, deriv=1, type='low')),
                                           threshold.high=mean(findThreshold(x=X, y=y, deriv=1, type='high')),
                                           threshold.max=mean(findThreshold(x=X, y=y, deriv=1, type='max')),
                                           threshold.av=mean(findThreshold(x=X, y=y, deriv=1, type='average')))
                p[[i]] <- ggplot(newdata, aes(y=y, x=X)) +
                    geom_ribbon(aes(ymin=lo, ymax=hi), fill='grey') +
                    geom_line() +
                                        #scale_y_continuous(paste0('f(',parse(text=xlab),')')) +
                    scale_x_continuous(xlab) +
                    theme_classic()
            }
        } else {
            thresholds[[preds[i]]]=NULL
            if (grouped) {
                p[[i]] <- ggplot(newdata, aes(y=y, x=X, group=group, color=group)) +
                    geom_blank() +
                    geom_linerange(aes(ymin=lo, ymax=hi), position=position_dodge(width=0.2)) +
                    geom_point(position=position_dodge(width=0.2)) +
                    scale_color_discrete(pretty.groupby) + 
                                        #scale_y_continuous(paste0('f(',parse(text=xlab),')')) +
                    scale_x_discrete(xlab) +
                    theme_classic()
                if (groupby=='NTR') {
                    p[[i]] = p[[i]] + scale_color_manual(pretty.groupby, values=c('blue','lightgreen', 'darkgreen')) 
                }
            } else  {
                p[[i]] <- ggplot(newdata, aes(y=y, x=X)) +
                    geom_blank() +
                    geom_linerange(aes(ymin=lo, ymax=hi)) +
                    geom_point() +
                                        #scale_y_continuous(paste0('f(',parse(text=xlab),')')) +
                    scale_x_discrete(xlab) +
                    theme_classic()
            }
        }
        Range[1] = ifelse(min(newdata$lo)<Range[1], min(newdata$lo), Range[1])
        Range[2] = ifelse(max(newdata$hi)>Range[2], max(newdata$hi), Range[2])
    }
    for (i in 1:length(preds)) {
        resp = ifelse(length(mod[[1]]$response.name)==1, mod[[1]]$response.name, mod[[1]]$response.name[2])
        pretty.resp = as.character(var.lookup$pretty.name[var.lookup$Abbreviation==resp])
        p[[i]] = p[[i]] + scale_y_continuous(stringr::str_wrap(pretty.resp,width=15), limits=Range) #+ scale_y_discrete(parse(text=ylab))
        #p[[i]] = p[[i]] + scale_y_continuous(ylab) #+ scale_y_discrete(parse(text=ylab)) 
    }
    r=NULL
    for (j in 1:length(mod)) {
        r = rbind(r, abt::relative.influence(mod[[j]]))
    }
    r=r[-which(rowSums(r)==0),]
    r=r %>% as.data.frame %>% mutate_all(function(x) 100*x/rowSums(.))
    rel.imps = apply(r, 2, function(x) c('Mean'=median(x),quantile(x, p=c(0.025,0.25,0.75,0.975)))) %>%
        as.data.frame %>%
        mutate(Stat=rownames(.)) %>%
        gather(key=Var, value=Value,-Stat) %>%
        left_join(var.lookup %>% dplyr::select(Var=Abbreviation, pretty.name)) %>%
        #left_join(xlabs %>% as.data.frame %>% gather(key=Var, value=var)) %>%
        spread(key=Stat, value=Value) %>% arrange(Mean) %>%
        mutate(Var=factor(Var, levels=unique(Var))) %>%
        dplyr::rename('Lower'=`2.5%`,'lower'=`25%`,'upper'=`75%`,'Upper'=`97.5%`) %>%
        mutate(Substantial=factor(ifelse(Mean>100/length(unique(Var)),1,0)), Substantial50=factor(ifelse(lower>100/length(unique(Var)),1,0)),Substantial95=factor(ifelse(Lower>100/length(unique(Var)),1,0))) 

    p2=ggplot(rel.imps, aes(y=Mean, x=Var)) +
            geom_linerange(aes(ymin=Lower,ymax=Upper, color=Substantial95), show.legend=FALSE) +
            geom_linerange(aes(ymin=lower,ymax=upper, color=Substantial50), show.legend=FALSE, size=1) +
            geom_point(aes(fill=Substantial, color=Substantial50),size=2, shape=21, show.legend=FALSE)+
            geom_hline(yintercept=100/nrow(rel.imps), linetype='dashed')+
            scale_y_continuous(expression(Relative~Importance)) +
            scale_x_discrete(labels=as.character(rel.imps$pretty.name)) + 
            scale_fill_manual(values=c('gray','black'))+
            scale_color_manual(values=c('gray','black'))+
            coord_flip() + theme_classic() +
            #ggtitle(paste0(letters[ii],')')) +
            theme(axis.title.y=element_blank(),plot.title=element_text(margin=ggplot2::margin(t=10,b=-10), hjust=0.01),, plot.margin=unit(c(0,1,0,0), 'lines'),
                  panel.spacing=unit(0,'lines'))
    p[['rel.imp']] = p2
    ri = rev(as.character((rel.imps %>% filter(Substantial==1))$Var))
    print(ri)
    #print(preds)
    wch=match(ri,preds)
    print(wch)
    #print(names(p))
    ps = p[c(length(p),wch)]

    
    list(p=p, ps=ps, thresholds=thresholds)
}

## ----end



## ---- stats.abt
## fitMethod=1: peffects.abt on newdata with NA values
## fitMethod=2; predict.abt on newdata
## fitMethod=3; predict.abt on data then averaged
stats.abt = function(mod, fitMethod=1, analysis) {
    N=length(mod)
    mdata=mod[[1]]$mdata
    pdata=mod[[1]]$mdata[,-1]
    preds = colnames(pdata)
    obs = mdata[[1]]
    fit = vector('list',N); names(fit) <- paste0('Boot',1:N);
    rel.imp = vector('list',N); names(rel.imp) <- paste0('Boot',1:N);
    R2.value = vector('list',N); names(R2.value) <- paste0('Boot',1:N);
    optim = vector('list',N); names(optim) <- paste0('Boot',1:N);

    m_gaussian <- function(y, f) sum((y - f)^2)
    m_binomial <- function(y, f) sum(2 * (y * log(ifelse(y == 0, 
        1, y/f)) + (1 - y) * log(ifelse(y == 0, 1, (1 - y)/(1 - 
        f)))))
    m_poisson <- function(y, f) sum(2 * (y * log(ifelse(y == 0, 
        1, y/f)) - (y - f)))
    m_laplace <- function(y, f) sum(abs(y - f))
    invlogit <- function(y) 1/(1 + exp(-y))
    dist <- mod[[1]]$distribution
    if (dist == "gaussian") {
        dev.fun <- m_gaussian
        fun <- I
        fun2=I
    } else if (dist == "bernoulli") {
        dev.fun <- m_binomial
        fun <- invlogit
        fun2<-binomial()$linkfun
    } else if (dist == "poisson") { # | names(mdata)[1]=='log(y)') {
        dev.fun <- m_poisson
        fun <- exp
        fun2<-log
    } else if (dist == "laplace") {
        dev.fun <- m_laplace
        fun <- I
        fun2<-I
    }
    
    for (i in 1:N) {
        #print(i)
        ri = relative.influence(mod[[i]])
        rel.imp[[i]] =  bind_rows(100*ri/sum(ri)) %>% gather(key=var, value=rel.inf) %>% as.data.frame
        fl=vector('list', length(preds)); names(fl) <- preds;
        pl=vector('list', length(preds)); names(pl) <- preds;
        a=cbind(pdata, y=predict.abt(mod[[i]], newdata=pdata)) #%>% mutate(y=binomial()$linkinv(y))
        ## Get the REGION averages on the response scale
        aa=a %>% mutate(y=fun(y))%>% lapply(FUN=function(x,y=.$y) data.frame(X1=x,y))
        for (j in preds) {
            #print(j)
            ## Make some prediction data.  Note, this will include means of non-focal predictors,
            ## Yet for some of the methods below, they will be replaced with NA values before predictions.
            newdata = prediction.data(Var=j, formula(mod[[i]]$Terms), data=mod[[1]]$mf, cyclic=NULL, groups=analysis[[j]], cat_levels=list(NTR.Pooled='NTR'))
            pred=j#preds[i]
            #xlab = xlabs[[pred]]
            ii = which(preds==pred)
            ## predicted values
            if (fitMethod==6) { #straight from plot.abt
                g=unique(c(j, ifelse(is.null(grps[[j]]), j, grps[[j]])))
                ad = newdata[[j]] %>% mutate_at(vars(-one_of(c(j,grps[[j]]))), funs(replace(.,TRUE,NA))) #%>%
                                        #dplyr::select_(.dots=g)
                ad = ad[names(ad) %in% g]
                fl[[j]] = peffects.abt(mod[[i]], data=ad, return.grid=TRUE, center=FALSE)
                }
            if (fitMethod==4) {
                g=ifelse(is.null(grps[[j]]), j, grps[[j]])
                group.means = a %>% mutate(y=binomial()$linkinv(y)) %>% group_by_(.dots=g) %>% summarize(ym=mean(y)) %>% ungroup %>% mutate(ym=binomial()$linkfun(ym))
                ad = newdata[[j]] %>% mutate_at(vars(-one_of(c(j,grps[[j]]))), funs(replace(.,TRUE,NA))) %>%
                    dplyr::select_(.dots=preds)
                fl[[j]] = peffects.abt(mod[[i]], data=ad, return.grid=TRUE, center=FALSE)  #%>% mutate(y=binomial()$linkinv(y))
                fl[[j]] = fl[[j]] %>% group_by_(.dots=g) %>% mutate(mu=mean(y)) %>% ungroup
                fl[[j]] = fl[[j]] %>% left_join(group.means) %>% mutate(y=y+(ym-mu))
            }
            if (fitMethod==2) fl[[j]] = cbind(newdata[[j]], y=predict.abt(mod[[i]], newdata=newdata[[j]])) #%>% mutate(y=binomial()$linkinv(y))
            if (fitMethod==1) {
                ad = newdata %>% mutate_at(vars(-one_of(c(j,analysis[[j]]))), funs(replace(.,TRUE,NA))) %>%
                    dplyr::select_(.dots=preds)
                fl[[j]] = peffects.abt(mod[[i]], data=ad, return.grid=TRUE, center=FALSE,n.trees=mod[[i]]$best) #%>% mutate(y=binomial()$linkinv(y))
            }
            if (fitMethod==5) {
                g=ifelse(is.null(grps[[j]]), j, grps[[j]])
                group.means = a %>% mutate(y=binomial()$linkinv(y)) %>% group_by_(.dots=g) %>% summarize(ym=mean(y)) %>% ungroup %>% mutate(ym=binomial()$linkfun(ym))
                ad = pdata %>% mutate_at(vars(-one_of(c(j,grps[[j]]))), funs(replace(.,TRUE,NA))) %>%
                    dplyr::select_(.dots=preds)
                fl[[j]] = peffects.abt(mod[[i]], data=ad, return.grid=TRUE, center=FALSE) #%>% mutate(y=binomial()$linkinv(y))
                fl[[j]] = fl[[j]] %>% group_by_(.dots=g) %>% mutate(mu=mean(y)) %>% ungroup
                fl[[j]] = fl[[j]] %>% left_join(group.means) %>% mutate(y=y+(ym-mu))
            }      
            if (fitMethod==3) fl[[j]] =a %>% group_by_(.dots=c(j,grps[[j]])) %>% summarize(y=mean(y))
            ## Generate data that is the obs for the focal predictor and means for all the rest
            dd = (mdata[,colnames(mdata)!=j] %>%
                mutate_all(.funs=function(x) ifelse(is.factor(x) | is.character(x), rev(sort(names(table(x))))[1] ,mean(x,na.rm=TRUE))) %>%
                #cbind(mdata[j]) %>% dplyr::select(-y)
                cbind(mdata[j]))[,-1]

            #dd = dd %>% mutate_at(vars(-one_of(c(j,analysis[[j]]))), funs(replace(.,TRUE,NA))) %>%
            #        dplyr::select_(.dots=preds)
       

            
            ## Create a data frame of predicted (y) and observed (y) on link scale for R2 calculations
            #pl[[j]]=data.frame(dd, y=predict.abt(object=mod[[i]], newdata=dd, n.trees=mod[[i]]$best), obs=fun2(mdata[[1]])) %>%
            #    dplyr::select(y,obs)
            pl[[j]]=data.frame(dd, y=fun(predict.abt(object=mod[[i]], newdata=dd, n.trees=mod[[i]]$best)), obs=(mdata[[1]])) %>%
                dplyr::select(y,obs)
            ## fl[[j]] = fl[[j]] %>% mutate_(.dots=setNames(pred,'X1'))
        }
        fit[[paste0('Boot',i)]] = fl
        R2.value[[paste0('Boot',i)]] = data.frame(lapply(pl,R2.gbm,method=1)) %>%
            mutate_all(.funs = function(x) ifelse(x<0,0,x))
       
        #optim[[paste0('Boot',i)]]=as.data.frame(lapply(aa ,getMax))
    }
    optim=abt.get.optim(mod) #as.data.frame(lapply(fl ,getMax))
    list(rel.imp=rel.imp, fit=fit,R2.value=R2.value, optim=optim)
}

## ----end


## ---- Prediction.data
## Generate a general prediction grid for a focal predictor (and optional group)
## This is a sequence for the focal predictor (optionally conditional on the levels of the grouping variable)
## and means for all other predictors
prediction.data = function(Var=NULL, formula, data, cyclic=NULL,groups=NULL, continuous.resolution=100, cat_levels=NULL) {
    terms=terms(formula)
    Var.names = attr(terms,'term.labels')
    num_vars = Var.names[unlist(lapply(data[,Var.names], 'is.numeric'))]
    cat_vars = Var.names[unlist(lapply(data[,Var.names], 'is.factor'))]
    dat = data
    if (!is.null(cyclic)) {
        for (i in length(cyclic)) {
            ii = which(Var.names %in% names(cyclic[[i]]))
            Var.names = Var.names[-1*ii]
            Var.names = c(Var.names, names(cyclic)[i])
            num_vars = num_vars[-1*ii]
            num_vars = c(num_vars, names(cyclic)[i])
        }
    }
    if (!is.null(groups)) if (is.na(groups)) groups=NULL; groupings=NULL
    if (!is.null(groups)) groupings=groups
    if (!is.null(groups)) {
        dat = data %>% dplyr::group_by_(.dots=groups)
    } else {
        dat= data
    }
    grid.levels = dat %>% do({
        z=. #%>% droplevels
        ll = list()
        focal = as.vector(as.data.frame(z[,Var])[,1])
        if (is.numeric(focal)) {
            if (all(is.na(focal))==TRUE) {
                ll[[Var]]=NA
            } else {
                ll[[Var]]=seq(min(focal,na.rm=TRUE), max(focal,na.rm=TRUE),len=continuous.resolution)
            }
        } else {
            ll[[Var]] = levels(as.data.frame(z)[,Var])    
        }
        k=0
        ## All other numeric predictors
        for (i in num_vars[num_vars!=Var]) {#x$var.names[x$var.names!=Var]) {
            focal = as.vector(as.data.frame(z[,i])[,1])
            if (all(is.na(focal))==TRUE) {
                ll[[i]]=NA
            } else {
                ll[[i]]=mean(focal,na.rm=TRUE)
            }
        }
        ## All other categorical predictors - just use a single level
        cat_var=cat_vars[cat_vars!=Var]
        if (!is.null(groupings)) cat_var=cat_var[cat_var!=groupings]
        for (i in cat_var) {#x$var.names[x$var.names!=Var]) {
            if (!is.null(cat_levels[[i]])) {
                ll[[i]]=factor(cat_levels[[i]], levels=levels(as.data.frame(z)[,i])) #levels(as.data.frame(z)[,i])
            } else {
                ll[[i]]=levels(as.data.frame(z)[,i])
            }
        }
        ll=expand.grid(ll)
        ## incase any of these are cyclic
        if (!is.null(cyclic)) {
            for (i in length(cyclic)) {
                ll = ll %>% mutate_(.dots=cyclic[[i]])
            }
        }
        ll 
    }) %>% na.omit %>% #droplevels %>%
    as.data.frame
    grid.levels
}

## ----end


## ---- Prediction.data.old
prediction_data.old = function(formula, data, group_by_cats=TRUE,groups=NULL,continuous.resolution=100, cyclic=NULL, cat_levels=NULL) {
    terms=terms(formula)
    Var.names = attr(terms,'term.labels')
    num_vars = Var.names[unlist(lapply(data[,Var.names], 'is.numeric'))]
    cat_vars = Var.names[unlist(lapply(data[,Var.names], 'is.factor'))]
    
    if (!is.null(cyclic)) {
        for (i in length(cyclic)) {
            ii = which(Var.names %in% names(cyclic[[i]]))
            Var.names = Var.names[-1*ii]
            Var.names = c(Var.names, names(cyclic)[i])
            num_vars = num_vars[-1*ii]
            num_vars = c(num_vars, names(cyclic)[i])
        }
    }
    grid.levels = list()
    for (ii in 1:length(Var.names)) {
        Var = Var.names[ii]
        if (group_by_cats) {
            groupings = cat_vars
            dat = data
            if (!is.null(groups)) {
                groupings = groups[[Var]]
                if (!is.null(groupings)) dat = data %>% dplyr::group_by_(.dots=groupings)
            }

        } else {
            dat= data
        }
        grid.levels[[Var]] = dat %>% do({
            z=.
            ll = list() #vector('list',length(num_vars))
            b = as.vector(as.data.frame(z[,Var])[,1])
            if (is.numeric(b)) {
                if (all(is.na(b))==TRUE) {
                    ll[[Var]]=NA
                } else {
                    ll[[Var]]=seq(min(b,na.rm=TRUE), max(b,na.rm=TRUE),len=continuous.resolution)
                }
            } else {
                ll[[Var]] = levels(as.data.frame(z)[,Var])    
            }
            k=0
            ## All other numeric predictors
            for (i in num_vars[num_vars!=Var]) {#x$var.names[x$var.names!=Var]) {
                b = as.vector(as.data.frame(z[,i])[,1])
                if (all(is.na(b))==TRUE) {
                    ll[[i]]=NA
                } else {
                    ll[[i]]=mean(b,na.rm=TRUE)
                }
            }
            ## All other categorical predictors - just use a single level
            cat_var=cat_vars[cat_vars!=Var]
            if (!is.null(groupings)) cat_var=cat_var[cat_var!=groupings]
            for (i in cat_var) {#x$var.names[x$var.names!=Var]) {
                if (!is.null(cat_levels[[i]])) {
                    ll[[i]]=factor(cat_levels[[i]], levels=levels(as.data.frame(z)[,i])) #levels(as.data.frame(z)[,i])
                } else {
                    ll[[i]]=levels(as.data.frame(z)[,i])
                }
            }
            ll=expand.grid(ll)
            ## incase any of these are cyclic
            if (!is.null(cyclic)) {
                for (i in length(cyclic)) {
                    ll = ll %>% mutate_(.dots=cyclic[[i]])
                }
            }
            ll 
        }) %>% na.omit %>% #droplevels %>%
        as.data.frame
    }
    grid.levels
    
}
## ----


## ---- R2.gbm
R2.gbm = function(x, method=1) {
    x=na.omit(x)
    y_i = x$obs
    u_i = x$y
    if (method==1) R2=cor(y_i,u_i)^2
    if (method==2) {
        r=y_i - u_i
        yadj = y_i - mean(y_i, na.rm=TRUE)
        R2=1 - sum(r^2)/sum(yadj^2)
    }
    if (method==3) {
        R2=(var(y_i, na.rm=TRUE) - mean((y_i - u_i)^2, na.rm=TRUE))/var(y_i, na.rm=TRUE)
    }
    if (method==4) {
        R2=var(u_i,na.rm=TRUE)/ (var(u_i,na.rm=TRUE) + var((y_i-u_i),na.rm=TRUE))
    }
    R2
}

## ----end


## ---- abt.get.optim
abt.get.optim = function(mod, center=FALSE, return.grid=TRUE,type='response', npts=100) {
    tree.list <- mod
    nlist=length(mod)
    x <- tree.list[[1]]
    n.trees <- sapply(tree.list, "[[", "n.trees")
    n.trees <- rep(n.trees, nlist)
    N=length(mod)
    mdata=mod[[1]]$mdata
    pdata=mod[[1]]$mdata[,-1]
    preds = colnames(pdata)
    optims = list()
    for (j in 1:length(preds)) {
        pred=preds[j]
        i.var <- j[order(x$var.type[j] > 0)]
        grid.levels <- vector("list", length(i.var))
        i = 1
        if (is.numeric(x$var.levels[[i.var[i]]])) {
            grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]), 
                                    max(x$var.levels[[i.var[i]]]), length = npts[i])
            XX = grid.levels[[i]]
        } else {
            grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]], 
                                                  levels = x$var.levels[[i.var[i]]])) - 1
            XX = factor(x$var.levels[[i.var[i]]], 
                                                  levels = x$var.levels[[i.var[i]]])
        }
        X <- expand.grid(grid.levels)
        names(X) <- paste("x", 1:length(i.var), sep = "")
        x.num <- (x$var.type == 0)[i.var]
        ytemp <- matrix(NA, nrow = nrow(X), ncol = nlist)
        for (i in 1:nlist) {
            ytemp[, i] <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
                cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
                n.class = as.integer(x$num.classes), i.var = as.integer(i.var - 
                  1), n.trees = as.integer(n.trees[i]), initF = as.double(tree.list[[i]]$initF), 
                trees = tree.list[[i]]$trees, c.splits = tree.list[[i]]$c.splits, 
                var.type = as.integer(x$var.type), PACKAGE = "abt")
        }
        wch=apply(ytemp, 2, which.max)
        optims[[pred]] = XX[wch]
    }
    optim = as.data.frame(optims) %>% mutate(Boot=paste0('Boot.',1:n()))
    optim = split(optim, optim$Boot)
    lapply(optim, function(x) x %>% dplyr::select(-Boot))
}

## ----end

## ---- summarize.values
summarize_values = function(val, type='optim',data.trans=data.trans, trans=trans) {
    optim = do.call('rbind',val) 
    optim1 = optim
    #if (type=='optim') {
    #    if (nrow(data.trans)==1) {
    #        optim1 = optim1 %>% cbind(data.trans) %>%
    #            mutate_(.dots=trans)
    #        optim1 = optim1 %>% dplyr:::select(-one_of(colnames(data.trans)))
    #    } else {
    #        optim1 = optim1 %>% left_join(data.trans) %>%
    #            mutate_(.dots=trans)
    #        optim1 = optim1 %>% dplyr:::select(-one_of(colnames(data.trans)))
    #    }
                                        #}
    if (type=='Rel.inf') {
        optim1 = optim1 %>% dplyr::rename('Var'='var', 'Value'='rel.inf') %>% group_by(Var)
    } else {
        optim1 = optim1 %>%
            dplyr::select_if(is.numeric) %>%
            gather(key=Var, value=Value) %>%
            group_by(Var)
    }
    
                                        #dplyr::summarize(Mean=sprintf('%2.2f (%2.2f,%2.2f)',mean(Value, na.rm=TRUE), lower=quantile(Value, p=0.025),upper=quantile(Value, p=0.975)))
    optim1 = optim1 %>%
        dplyr::summarize(Mean=sprintf('%2.2f',mean(Value, na.rm=TRUE)),
                         Median=sprintf('%2.2f',median(Value,na.rm=TRUE)),
                         lower=sprintf('%2.2f',quantile(Value, p=0.025, na.rm=TRUE)),
                         upper=sprintf('%2.2f',quantile(Value, p=0.975, na.rm=TRUE))) 
    optim2 = optim %>% dplyr::select_if(is.factor)
    if (ncol(optim2)>0) {
        optim2=optim2 %>% gather(key=Var, value=Value) %>%
        group_by(Var) %>%
                                        #dplyr::summarize(Mean=sprintf('%s',names(rev(sort(table(Value))))[1]))
        dplyr::summarize(Mean=names(rev(sort(table(Value))))[1], Median='', lower='', upper='')
        optim1=optim1 %>% full_join(optim2)
    }
    optim1
}

## ----end


## ---- get.response
get_response <- function(mod) {
    if (length(mod[[2]])==1) return(mod[[2]])
    if (length(mod[[2]])==2) return(mod[[2]][[2]])
}
## ----end


## ---- prettytables
pretty.stats = function(stats, var.lookup) {
    stats = stats %>% left_join(var.lookup %>% dplyr::select(pretty.name, Var=Abbreviation))
    stats %>%
        arrange(-as.numeric(as.character(Median.rel.inf))) %>%
        mutate(Optimum=sprintf('%s (%s-%s)', Mean.optim, lower.optim, upper.optim),
               `R-sq`=sprintf('%s (%s-%s)', Median.R2, lower.R2, upper.R2),
               `Rel inf`=sprintf('%s (%s-%s)', Median.rel.inf, lower.rel.inf, upper.rel.inf)) %>%
        select(Covariate=pretty.name, Optimum, `R-sq`, `Rel inf`)
}
pretty.thresholds = function(thresholds, stats, var.lookup) {
    thresholds =
        bind_rows(thresholds, .id='Var') %>%
        left_join(stats %>% dplyr::select(Var,Median.rel.inf)) %>%
        arrange(-as.numeric(as.character(Median.rel.inf))) %>%
        left_join(var.lookup %>% dplyr::select(pretty.name, Var=Abbreviation)) %>%
        select(Covariate=pretty.name, everything(), -Var,-Median.rel.inf)

    thresholds
}

## ----end

