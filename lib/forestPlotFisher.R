forestPlotFisher <- function(fisherTable, SampleSheet, pVal = 0.05, fdr = NULL,
                             color=c('maroon','royalblue'),
                             geneFontSize = 0.8, titleSize = 1.2, lineWidth = 1,
                             m1Name = 'control', m2Name = 'case', m1.sampleSize = NULL, m2.sampleSize = NULL
){
  
  if(is.null(fdr)){
    m.sigs = fisherTable[fisherTable$pval < pVal,]
  } else{
    m.sigs = fisherTable[fisherTable$adjPval < fdr,]
  }
  
  m1.sampleSize <- length(unique(SampleSheet$patient_id[SampleSheet$case_control == m1Name]))
  m2.sampleSize <- length(unique(SampleSheet$patient_id[SampleSheet$case_control == m2Name]))
  
  # newly added, deal with parameter color -> vc_col for usage
  if (length(color)==1){
    vc_col = c(color,color)
  }else if(length(color==2)){
    vc_col = color
  }else{
    stop('colors length must be less equal than 2')
  }
  
  if (is.null(names(vc_col))){
    names(vc_col)=c(m2Name,m1Name)
  }else if (!(m1Name %in% names(vc_col)) && !(m2Name %in% names(vc_col))){
    stop(paste0('\ninput named vector [color] must contain both group name, \nwhich should be like c(',
                m1Name,', ', m2Name,') but now are ',list(names(vc_col)),'\n'))
  }else if ( !(m1Name %in% names(vc_col)) || !(m2Name %in% names(vc_col)) ){
    
    stop(paste0('\nif you pass [color] with named vector, then both of the names matching group names must be Explicitly declared,\n',
                'which should be like c(',
                m1Name,', ', m2Name,') but now are ',list(names(vc_col)),'\n'))
  }
  # end of newly added
  
  if(nrow(m.sigs) < 1){
    stop('No differetially mutated genes found !')
  }
  
  m.sigs$Aberration = factor(x = m.sigs$Aberration, levels = rev(m.sigs$Aberration))
  
  m.sigs$or_new = ifelse(test = m.sigs$or > 3, yes = 3, no = m.sigs$or)
  m.sigs$upper = ifelse(test = m.sigs$ci.up > 3, yes = 3, no = m.sigs$ci.up)
  m.sigs$lower = ifelse(test = m.sigs$ci.low > 3, yes = 3, no = m.sigs$ci.low)
  #Reverse by significance
  m.sigs$pos = rev(1:nrow(m.sigs))
  #m.sigs = m.sigs[order(m.sigs$pos)]
  m.sigs$pos = 1:nrow(m.sigs)
  xlims = c(0, 4)
  ylims = c(0.75, nrow(m.sigs))
  
  graphics::layout(mat = matrix(c(1, 2, 3, 4, 5, 6, 6, 6, 6, 6), byrow = TRUE, ncol = 5, nrow = 2), widths = c(4, 1, 1), heights = c(6, 1.2))
  par(mar = c(3, 1, 3, 5))
  plot(NA, xlim = xlims, ylim = ylims, axes= FALSE, xlab = NA, ylab = NA)
  
  apply(m.sigs[,c('or', 'ci.up', 'ci.low', 'ci.up', 'or_new', 'upper', 'lower', 'pos')], 1, function(x){
    p = x[5]; u = x[6]; u_orig = x[2]; l = x[7]; l_orig = x[3]; ypos = x[8]
    if (p<1){
      linecolor = vc_col[m2Name]
    }else if(p>1){
      linecolor = vc_col[m1Name]
    }else{
      linecolor = 'black'
    }
    points(x = p, y = ypos, pch = 16, cex = 1.1*(lineWidth))
    segments(x0 = l, y0 = ypos, x1 = u, y1 = ypos, lwd = lineWidth,col = linecolor)
    
    if(u_orig >3){
      segments(x0 = 3, y0 = ypos, x1 = 3.25, y1 = ypos, lwd = lineWidth,col=linecolor)
      points(x = 3.25, y = ypos, pch = ">", cex = 1.1*(lineWidth))
    }
    if(l_orig >3){
      segments(x0 = 3, y0 = ypos, x1 = 3.25, y1 = ypos, lwd = lineWidth,col=linecolor)
      points(x = 3.25, y = ypos, pch = ">", cex = 1.1*(lineWidth))
    }
    
  })
  
  geneFontSize = 0.8; titleSize = 1.2
  
  abline(v = 1, lty = 2, col = "gray", xpd = FALSE)
  axis(side = 1, at = 0:3, labels = c(0:3), font = 1, pos = 0.5, cex.axis = 1.3)
  
  mtext(text = m.sigs$Aberration, side = 4, line = 0.2, at = 1:nrow(m.sigs),
        font = 3, las= 2, cex = geneFontSize, adj = 0)
  mtitle = paste(m2Name, ' (n = ', m2.sampleSize, ')', ' v/s ' , m1Name, ' (n = ' ,m1.sampleSize, ')', sep='')
  title(main = mtitle, font = 1, adj = 0, cex.main = titleSize)
  #mtext(text = "Odds ratio", side = 1, line = 3, font = 1, cex = 0.7*(titleSize), adj = 0.25)
  
  # plot annotation columns of the graph c(group1,group2,OR,p-value) col (col 2 ~ col 5)
  # annotation columns group2
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = ylims)
  text(x = 0.5, y = 1:nrow(m.sigs), labels = as.numeric(unlist(m.sigs[,3])),
       adj = 0, font = 1, cex = 1.4*(geneFontSize))
  title(main = m2Name, cex.main = titleSize)
  # annotation columns group1
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = ylims)
  text(x = 0.5, y = 1:nrow(m.sigs), labels = as.numeric(unlist(m.sigs[,2])),
       adj = 0, font = 1, cex = 1.4*(geneFontSize))
  title(main = m1Name, cex.main = titleSize)
  # annotation columns OR
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = ylims)
  text(x = 0.5, y = 1:nrow(m.sigs), labels = round(m.sigs$or, digits = 3),
       adj = 0.5, font = 1, cex = 1.4*(geneFontSize))
  title(main = "OR", cex.main = titleSize)
  
  m.sigs$significance = ifelse(test =  as.numeric(m.sigs$adjPval) < 0.001, yes = "***", no =
                                 ifelse(test = as.numeric(m.sigs$adjPval) < 0.01, yes = "**", no =
                                          ifelse(test = as.numeric(m.sigs$adjPval) < 0.05, yes = "*", no = "NS")))
  # annotation columns P-value
  par(mar = c(3, 0, 3, 0))
  plot(rep(0, nrow(m.sigs)), 1:nrow(m.sigs), xlim = c(0, 1), axes = FALSE,
       pch = NA, xlab = "", ylab = "", ylim = ylims)
  text(x = 0.5, y = 1:nrow(m.sigs), labels = m.sigs$significance,
       adj = 0, font = 1, cex = 1.4*(geneFontSize))
  title(main = "AdjP-value", cex.main = titleSize)
  
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(0,1), ylim = c(0, 1), axes = FALSE, xlab = NA, ylab = NA)
  
  text(x = 0, labels = paste0(
    "Odds ratio with 95% CI\n(1 = no effect, < 1 ",
    m2Name,
    " has more mutants)"
  ), y = 0.6, adj = 0, xpd = TRUE, cex = 1.2)
  
}