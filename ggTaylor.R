modObs2ggTaylor <- function(obs, mod) {
  ok  <- is.finite(obs) & is.finite(mod)
  obs <- obs[ok]
  mod <- mod[ok]
  sdObs <- sd(obs)
  sdModRel <- sd(mod)/sdObs
  correl <- cor(mod,obs)
  
  tibble(sdModRel=sdModRel, correl=correl, sdObs=sdObs, bias=mean(mod-obs))
  
}

##Taylor plot acc. to Taylor, K.E. (2001). "Summarizing multiple aspects of model performance in a single diagram". J. Geophys. Res. 106: 7183–7192. Bibcode:2001JGR...106.7183T. doi:10.1029/2000JD900719.
## cf. also https://en.wikipedia.org/wiki/Taylor_diagram
gg_taylor  <-
  function(sdModRel = NULL, ## The model SD relative to the  observed, i.e. sdMOD/sdObs
           correl = NULL,   ## Pearson correlation
           labels = NULL,   ##  Labels for the points, e.g. model names
           bias = NULL,     ## bias of the model mean(mod-obs)
           sdObs = 1.0,     ## the sd of the observations (everything is scaled with that)
           plotBias = TRUE, ## should the bias be indicated?
           plotNoVarErr = FALSE, ## should be model without a variance error be indicated?
           add = FALSE,     ## should only new points be plotted?
           sdMax = max(seq(0, max(dfTaylor$sdModRel, sdModRel, 1) +
                             0.099, 0.1)), ## What should be the maximum sd plotted (like xlim, but for the circle)
           dfTaylor = NULL, ## The statistics could also be in a tibble with respective column names
           ...) {
    
  if (is_tibble(dfTaylor)) {
    sdModRel <- dfTaylor$sdModRel
    correl <- dfTaylor$correl
    sdObs <- dfTaylor$sdObs
    labels <- dfTaylor$labels
    bias <- dfTaylor$bias
  }
  if (is.null(labels)) labels <- seq(length(correl)) %>% as.character()
  
  sdObs <- mean(sdObs)
  
  correl_sd2Taylor <- function(sd, correl) {
    tibble(x=correl*sd, y=sqrt(sd^2-x^2))
  }
  
  if (!add) { ## create a new Taylor coordinate system
    # Tibble with Grid of standard deviation (relative to obs) (=radius) and correlation (=sin(angle)
    taylorGrid <- cross_df(list(correls=c(-1, -0.99, -0.95, seq(-0.9,0.9,0.1), 0.95, 0.99, 1.0),
                         sds = seq(0,sdMax,0.1)))
    # Tibble with Isolines of equal RMSD (concentric around the perfect fit, correl=1, sd=1)
    isoRMSD <- cross_df(list( RMSD  = c( seq(0.2,2*sdMax+0.2, 0.2)),  
                              xVec = seq(-1.,1.,0.01))) %>%
      mutate( x2 = (xVec * RMSD + 1)*sdObs, y2 = sqrt(1-xVec^2)*RMSD*sdObs ) %>% 
      filter(sqrt(x2^2+y2^2)<=1.005*max(taylorGrid$sds)*sdObs)

    ## Plot the coordinate system first (project Standard devs and correls)
    p  <- taylorGrid %>% mutate(xv = correls*sds, yv = sqrt(sds^2-xv^2)*sdObs, xv=xv*sdObs) %>% 
      {
        ggplot(data=., mapping=aes(xv, yv)) + labs(x="Standard deviation", y="Standard deviation") + 
          # The grid
          geom_line(aes(group=sds), linetype=3) + geom_line(aes(group=correls), linetype=3) + 
          # The circle where sdMod = sdObs
          geom_line(data=filter(., sds==1), mapping=aes(xv,yv), linetype=2, size=1) +   
          # The RMSD isolines
          geom_line(data=isoRMSD, mapping = aes(x2, y2, group=RMSD), color="blue", size=0.5, alpha=0.5) +
          scale_x_continuous(expand = c(0,0) ) + scale_y_continuous(expand=c(0,0)) + theme_classic()  
          
      }
    ## Add labels for isolines and for circular ticks
    correlTicks <- c(seq(-0.9,0.9,0.1), 0.95, 0.99) %>% setdiff(0)
    p  <-  p +  geom_label(data=isoRMSD %>% filter(RMSD>0.05) %>% group_by(RMSD) %>% 
                             filter(y2>0.3^2*sdObs^2/abs(x2), x2>0) %>% 
                             summarise(across(.fns=first)), 
                           mapping = aes(x2,y2, label=paste0(RMSD,"\n", 1-RMSD^2), group=NULL), size=3) +
      annotate("text", x=correlTicks*sdMax*sdObs*1.04, 
               y=sqrt(sdMax^2-(correlTicks*sdMax)^2)*sdObs*1.04, label=correlTicks) +
      scale_x_continuous(expand=expansion(mult=c(0.03,0.05))) + scale_y_continuous(expand=expansion(mult=c(0,0.05)))
  } else { # Else just start from the last plot
    p <- last_plot()
  }
  
  ### Finally add the actual points
  
  ### First the classical Taylor plot, which "ignores" the bias
  modPoints <- correl_sd2Taylor(sdModRel,correl) * sdObs
   p  <- p + geom_text(data=modPoints %>% mutate(labels=labels), mapping=aes(x,y), fontface="bold", label=labels,...) +
     ## Add the reference observed point
    annotate("point", x=sdObs, y=0, color="blue", size=5, alpha=0.5)
  
   ## Where would the model be if there was no variance error?
  if (plotNoVarErr) {
    modPointsNoVarErr <- correl_sd2Taylor(1, correl) * sdObs
    p  <- p + geom_point(data=modPointsNoVarErr, mapping=aes(x,y), color="red", size=2, alpha=0.5)
  } 
  
   ## Add the bias into the plot: the distance to the blue dot is the RMSE incl. the bias
   ## MSE^2 = MSD^2 + bias^2 , so the bias is add perpendicular to the RMSD point
  xb <- NULL
  if (!is.null(bias) & plotBias) {
   
    RMSD <- sqrt(sdObs^2 + (sdModRel*sdObs)^2 - 2*sdObs^2*sdModRel*correl)
    RMSE <- sqrt(bias^2+RMSD^2)
    RMSDangle <- asin(modPoints$y/RMSD)
    RMSEangle  <- asin(bias/ RMSE)
    alpha <- ifelse(modPoints$x > sdObs, pi - RMSDangle - RMSEangle, RMSDangle - RMSEangle)
    xb <- sdObs - cos(alpha)* RMSE
    yb <- sin(alpha) * RMSE
    
    p <- p +     annotate("point", x=xb, y=yb, color="red", size=2, alpha=0.5) + 
      geom_segment(data=modPoints, aes(x=x, y=y, xend=xb, yend=yb), 
                   arrow=arrow(length=unit(0.1, "cm"), type="closed"))
  }

  p + coord_cartesian(xlim=c(min(0, modPoints$x*1.1, xb), NA))


}
