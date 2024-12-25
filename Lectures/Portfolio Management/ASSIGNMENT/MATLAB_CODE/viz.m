clear all
close all
clc

run("A_load_data.m")
load('portfolios.mat')

portfolios = {portfolioA, portfolioB, portfolioC, portfolioD, portfolioE, ...
    portfolioF, portfolioG, portfolioH, portfolioI, portfolioJ, portfolioK,...
    portfolioL, portfolioM, portfolioN, portfolioP, portfolioQ, portfolioCap};

PortfolioViz.plotEvolution(portfolios, dates_, Ret);