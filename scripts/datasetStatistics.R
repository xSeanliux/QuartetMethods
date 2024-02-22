

library(shiny)
library(plotly)
library(optparse)
library(dplyr)
library(stringr)
library(ape)
library(testit)
library(phangorn)
library(castor)

setwd("/Users/liusean/Desktop/Projects/Coding/Phylo/LingPhyloR")

source('inferenceUtils.R')

getCharCountsForMp1 <- function(numedge=0) {
  res = data.frame()
  for (c in list('low', 'mod', 'high', 'modhigh', 'veryhigh')) {
    print(c)
    vec = c()
    vec2 = c()
    for (t in 1:32) {
      for (r in 1:4) {
        if (numedge == 0) file = paste0('/Users/liusean/Desktop/Projects/Coding/Phylo/QuartetMethods/example/simulated_data_small-8.0/', c,'_noborrowing/', 'sim_tree',t,'_',r,'.csv')
        else {
          assert({FALSE})
        #   if (c %in% 1:5) ch = c+5
        #   else if (c==11) ch = 12
        #   else assert({FALSE})
        #   file = paste0("../SimulationPipeline/sim_outputs/config", ch, "/sim_net", numedge, '-', t, '_', r, '.csv')
        }
        print(c("file is ", file))
        df = read.csv(file)
        doubledf <- preprocessData(df, FALSE, 'Remove entire character', FALSE)$df
        uninf = findUninformativeCharacters(doubledf)$uninf
        vec <- c(vec, nrow(doubledf))
        vec2 <- c(vec2, nrow(doubledf)-uninf)
      }
    }
    print(c('The length of vec is ', length(vec)))
    # assert({length(vec) == 128})
    
    # if (c == 11) c=6
    res = rbind(res, data.frame(config=c, 
                           mp1=mean(vec), mp1std=sd(vec),
                           mp1inf=mean(vec2),mp1infstd=sd(vec2)/length(vec2)))
  }
  res
}


find_uninformative_characters <- function(df) {
  uninf <- c()
  numbigstates <- c()
  colstart <- 5
  df <- replace_qs_with_nums(df)
  for (char in seq_len(nrow(df))) {
    states <- as.character(df[char, colstart:ncol(df)])
    unique_states <- unique(states)
    states_with_atleast2 <- c()
    for (state in unique_states) {
      if (sum(states == state) >= 2) {
        states_with_atleast2 <- c(states_with_atleast2, state)
      }
    }
    if (length(states_with_atleast2) < 2) {
      uninf <- c(uninf, char)
    }
    numbigstates <- c(numbigstates, length(states_with_atleast2))
  }

  list(uninf = uninf, numbigstates = numbigstates)
}

