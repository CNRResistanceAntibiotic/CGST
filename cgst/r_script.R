#!/usr/bin/env Rscript



#Copyright 2020 Aurelien BIRER (abirer36@gmail.com)
#https://github.com/Nilad/CGST.git
#
#This script is the main step of CGST.
#
#It takes a single argument: a working directory containing the "similarity_matrix.csv" file produce by CGST.
#This R programm goal is to define few clusters of strains by the Ward Method with a similarty matrix.
#
#
#This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the Free Software Foundation,
#either version 3 of the License, or (at your option) any later version. CGST is distributed in
#the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#details. You should have received a copy of the GNU General Public License along with CGST. If
#not, see <http://www.gnu.org/licenses/>.


# Load Libraries
library(optparse)
library(questionr)
library(cluster)
library(fastcluster, warn.conflicts = FALSE)

################################ 

# Parse args
option_list = list(
make_option(c("-w", "--wd"), type = "character", default = "",
help = "Working Directory", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Read Data
data = read.csv(sprintf("%s/similarity_matrix.tsv", opt$wd), sep = "\t", row.names = 1)

# Change to matrix class
simi_matrix <- as.matrix(data)

# Get Disimilarity matrix with daisy
disi_matrix <- daisy(simi_matrix, metric = "euclidean", stand = TRUE)

# Get Hclust tree with fastcluster
arbre_ward <- hclust(disi_matrix, method = "ward.D2")

####
# get inertie results
inertie <- sort(arbre_ward$height, decreasing = TRUE)
# Get 3 Most Highest Inertie Break Point
l <- c()
n <- length(inertie) - 1
for (i in 1 : n) {
    l[i + 1] = inertie[i] - inertie[i + 1]
}
ndx <- sort(order(l, decreasing = T)[1 : 3])

#########
# Get Inertie Graph
pdf(sprintf("%s/inertial_graph.pdf", opt$wd))
plot(inertie[1 : n], type = "s", xlab = "Number of class", ylab = "Inertia", main = sprintf("Inertial Graph (3 best inertia)", ndx[1], ndx[2], ndx[3]),)
points(c(ndx[1], ndx[2], ndx[3]), inertie[c(ndx[1], ndx[2], ndx[3])], col = c("green3", "red3", "blue3"), cex = 2, lwd = 3)
legend("topright", bty = "n", legend = paste(ndx[1 : 3], " class"), text.col = c("green3", "red3", "blue3"), cex = .8)
invisible(dev.off())

#####
# Get Dendrogram Graph
pdf(sprintf("%s/dendrogram_with_class.pdf", opt$wd))
plot(arbre_ward, main = sprintf("Dendrogram partionned (3 best inertie)", ndx[1], ndx[2], ndx[3]), xlab = "", ylab = "", sub = "", axes = FALSE, hang = - 1)
legend("topright", bty = "n", legend = paste(ndx[1 : 3], " class"), text.col = c("green3", "red3", "blue3"), cex = .8)
rect.hclust(arbre_ward, ndx[1], border = "green3")
rect.hclust(arbre_ward, ndx[2], border = "red3")
rect.hclust(arbre_ward, ndx[3], border = "blue3")
invisible(dev.off())

##########
# Get Groups by Inertie Class
typo <- cutree(arbre_ward, k = c(ndx[1], ndx[2], ndx[3]))
write.table(typo, file = sprintf("%s/groups.tsv", opt$wd), quote = FALSE, sep = '\t', col.names = NA)
