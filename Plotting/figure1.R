

args = commandArgs(trailingOnly = TRUE)
args = as.vector(unlist(args))
library(dict)
library(gridDebug)
library(gplots) 
library(ape)

`%||%` <- function(a, b) {
  if (!is.null(a)) {
    a
  } else {
    b
  }
}

heatmapfix <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
    distfun = dist, hclustfun = hclust, reorderfun = function(d, 
        w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
        "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
    margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
    verbose = getOption("verbose"), ...) 
{
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv")) 
        doCdend <- FALSE
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc
    x <- x[rowInd, colInd]
    labRow <- labRow[rowInd] %||% rownames(x) %||% (1L:nr)[rowInd]
    labCol <- labCol[colInd] %||% colnames(x) %||% (1L:nc)[colInd]
    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    #lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 4)
    lhei <- c((0.00) + if (!is.null(main)) 0.2 else 0, 4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1L], 0.2, lhei[2L])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1L], 0.2, lwid[2L])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        print(lmat)
    }
    dev.hold()
    on.exit(dev.flush())
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1L], 0, 0, 0.5))
        image(rbind(if (revC) 
            nr:1L
        else 1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2L]))
        image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1L], 0, 0, margins[2L]))
    if (!symm || scale != "none") 
        x <- t(x)
    if (revC) {
        iy <- nr:1
        if (doRdend) 
            ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1L:nr
    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1L] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2L] - 1.25)
    if (!missing(add.expr)) 
        eval.parent(substitute(add.expr))
    par(mar = c(margins[1L], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
    if (doCdend) 
    	   wal=1
        #plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()
    if (!is.null(main)) {
        par(xpd = NA)
        title(main, cex.main = 1.5 * op[["cex.main"]])
    }
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}



normfile = "~/Documents/Marklab/SMN_group1_SMN1OOOSMN2OOORP11-974F13v3.fasta_annotate.fa_norm.txt"
treefile="~/Documents/Marklab/SMN_group1_SMN1OOOSMN2OOORP11-974F13v3.fasta_annotate.fa_tree.ph"
hclustaverage <-function(hclust)
{
	hclust(hclust,"average")
}


DistMat  = read.table(normfile, sep = ',', header = 0)
rownames(DistMat) = colnames(DistMat)
matrix = as.matrix(DistMat)
diags =  diag(matrix)
#DistMat= DistMat[which(diags>100),which(diags>100)]
#usenames = colnames(DistMat)[which (!(colnames(DistMat) %in% exclude))]
#matrix=DistMat
#matrix=DistMat[usenames,usenames]
#diags =  diag(as.matrix(matrix))
#matrix = matrix[which(diags>100),which(diags>100)]
#matrix=DistMat[colnames(DistMat)[which (DistMat[1,]>0)],colnames(DistMat)[which (DistMat[1,]>0)]]
#matrix=DistMat
#matrix = as.matrix(matrix)
#diags =  diag(matrix)

matrix_normed = matrix
for(row in 1:nrow(matrix)) {


    for(col in (row:ncol(matrix))) {
		association = (matrix_normed[row, col])/sqrt(diags[row]*diags[col])
        matrix[col, row] = association 
		matrix[row, col] = association 
    }
}

SMN_colors = c('lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightseagreen','lightsteelblue','lightsteelblue','lightsteelblue','lightsteelblue','lightsteelblue','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightcyan','lightpink','lightpink','lightpink','lightgreen','lightgrey','lightgrey','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightsalmon','lightsalmon','salmon','lightsalmon','lightsalmon','lightsalmon','lightsalmon','lightsalmon','lightslategray','lightslategray','lightslategray','lightslategray','lightslategray','lightslategray','slategray','lightslategray','lightyellow','lightyellow','lightyellow','lightyellow','lightyellow','lightyellow','lightyellow','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightskyblue','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightgray','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightcoral','lightseagreen')


#plot = heatmap(matrix,symm =T,margin = c(0,0),hclustfun = hclustaverage)
phylo_tree<-ape::read.tree(treefile)
dist_tree <- ape::cophenetic.phylo(phylo_tree)  
dist_tree <- as.dist(dist_tree)
hclust_tree <- hclust(dist_tree)
dendro_tree <- as.dendrogram(hclust_tree)
leaf_order <- rev(order.dendrogram(dendro_tree))
reversed_matrix <- matrix[order(rev(leaf_order)), ]

pdf("/Users/walfred/Desktop/SMN_heat.pdf", width = 10, height=10)

plot = heatmapfix(reversed_matrix,symm =T,margin = c(0,0),Colv = rev(dendro_tree), Rowv = rev(dendro_tree),RowSideColors =SMN_colors, labCol = FALSE)

dev.off()




plot = heatmap.2(reversed_matrix,symm =T,margin = c(0,0),Colv = rev(dendro_tree), Rowv = rev(dendro_tree),ColSideColors = rep("transparent", ncol(reversed_matrix)),RowSideColors =SMN_colors,density.info="none",trace="none",col = cols , scale = "none")

treenames = rownames(matrix)[plot$rowInd]

pheatmap(dm,symm =T,margin = c(0,0),clustering_method = "average",treeheight_row = 500,treeheight_col = 0,legend=FALSE)

regression = "~/Documents/Marklab/benchmark_out.txt"

reg_data = read.table(regression, sep = ',', header = 0)

reg_results = dict()

for (i in 1:length(rownames(reg_data)))
{
	reg_results[reg_data[i,1]] = reg_data[i,2]
}

reg_result = matrix()
for (i in 1:length(treenames) )
{
	name = treenames[i]
	if (reg_results[name] == 'NULL' )
	{
		next
	}
	else
	{
		matrix[i,which(treenames == reg_results[name])] = 
		
	}
}

par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
plot(reg_x,-reg_y,col='blue',lwd = 0.001, pch= 20,cex = 0.6,yaxt="n", ylab="",xlab="", axes=F,margin = c(0,0))

points(df2$x, df2$y, col='red')

points(0.03, 0.4, col='blue',lwd = 0.001, pch= 20,cex = 0.7)

X11()


par()$pin

par(new=T)


plot = heatmap.2(matrix,symm =T,density.info="none",trace="none",col=heat.colors(10) )


heatmap.2(matrix, Rowv=dend, Colv=dend)
phylo_tree<-ape::read.tree(treefile)
dist_tree <- ape::cophenetic.phylo(phylo_tree)  
dist_tree <- as.dist(dist_tree)
hclust_tree <- hclust(dist_tree)
dendro_tree <- as.dendrogram(hclust_tree)

library(RColorBrewer)
color_scheme <- colorRampPalette(c("blue", "white", "red"))

# Your color scale
colors = c("blue", "white", "red")

# The breakpoints you want. Adjust as needed.
breakpoints = c(min(matrix), 0, max(matrix))

# Create the color mapping function
color_mapping = colorRamp2(breakpoints, colors)

ht = Heatmap(
  matrix[, ncol(matrix):1],
  col = color_mapping,
  column_title = NULL,  # removes column title
  row_title = NULL,  # removes row title
  row_dend_reorder = FALSE,  # keep the order of dendrogram
  row_dend_width = unit(5, "cm"),  # adjust this to make tree larger or smaller
  show_row_dend = TRUE,  # shows the row dendrogram
  show_column_dend = FALSE,  # hides the column dendrogram
  show_row_names = FALSE,  # hides row labels
  show_column_names = FALSE,  # hides column labels
  cluster_columns = FALSE,  # don't cluster columns
  cluster_rows = hclust_tree,  # use your own dendrogram for row clustering
  heatmap_width = unit(15, "cm"),  # adjust the heatmap width
  heatmap_height = unit(10, "cm"),  # adjust the heatmap height
  show_heatmap_legend = FALSE
)

pdf("/Users/walfred/Desktop/SMN_heat.pdf", width = 10, height=10)

draw(heatmap)


dev.off()




pheatmap(matrix, cluster_rows = hclust_tree, cluster_cols = hclust_tree, treeheight_row = 200, treeheight_col = 200,legend=FALSE,margin = c(200,100), labels_row = NA, labels_col = NA,show_rownames = FALSE, show_colnames = FALSE)




plot = heatmap.2(matrix, Rowv = dendro_tree, Colv = dendro_tree, dendrogram = "none")



library(ComplexHeatmap)
library(circlize)










matrix = matrix/diag(matrix)

matrix = -log(1- matrix)
matrix[!is.finite(matrix) ] <- max(matrix[is.finite(matrix) ])
matrix[matrix > 3] <- 5 + 0.1*(matrix[matrix > 3]-3)








p$rowInd






