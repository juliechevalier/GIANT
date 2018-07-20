# Copyright (c) 2011-2013 Trevor L. Davis <trevor.l.davis@stanford.edu>  
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  


##comment function to display message and optionnaly add it to log file

addComment <- function(text,addToFile=FALSE,fileName=NULL,append=TRUE,display=TRUE){
  if(display)cat(paste(c(text,"\n"),collapse = " ")) 
  if(addToFile)write(paste(text,collapse = " "),fileName,append=append)
}

##negative of a mathematical expression
negativeExpression <- function(expression){
  expression=gsub("\\+","_toMinus_",expression)
  expression=gsub("\\-","+",expression)
  expression=gsub("_toMinus_","-",expression)
  if(substr(expression,1,1)!="-" && substr(expression,1,1)!="+"){
    expression=paste(c("-",expression),collapse="")
  }

  return(expression)
}

#' Returns file name of calling Rscript
#'
#' \code{get_Rscript_filename} returns the file name of calling Rscript 
#' @return A string with the filename of the calling script.
#'      If not found (i.e. you are in a interactive session) returns NA.
#'
#' @export
get_Rscript_filename <- function() {
    prog <- sub("--file=", "", grep("--file=", commandArgs(), value=TRUE)[1])
    if( .Platform$OS.type == "windows") { 
        prog <- gsub("\\\\", "\\\\\\\\", prog)
    }
    prog
}

#' Recursively sorts a list
#'
#' \code{sort_list} returns a sorted list
#' @param unsorted_list A list.
#' @return A sorted list.
#' @export
sort_list <- function(unsorted_list) {
    for(ii in seq(along=unsorted_list)) {
        if(is.list(unsorted_list[[ii]])) {
            unsorted_list[[ii]] <- sort_list(unsorted_list[[ii]])
        }
    }
    unsorted_list[sort(names(unsorted_list))] 
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
