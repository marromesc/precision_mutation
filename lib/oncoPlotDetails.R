col = c(missense = "violet", multi_hit = "black", nonsense = "red", UTR_variant = "yellow", splice_site = "green", inframe_indel = "orange",frameshift="blue",splicing="yellow")
alter_fun <- list(
  missense=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["missense"], col=NA))
  },
  UTR_variant=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["UTR_variant"], col=NA))
  },
  splice_site=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["splice_site"], col=NA))
  },
  inframe_indel=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["inframe_indel"], col=NA))
  },
  frameshift=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["frameshift"], col=NA))
  },
  nonsense=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["nonsense"], col=NA))
  },
  multi_hit=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["multi_hit"], col=NA))
  },
  splicing=function(x, y, w, h)
  {
    grid.rect(x, y, w*0.9, h*0.9, gp=gpar(fill=col["splicing"], col=NA))
  }
)

PastelOrange <- rgb(red=100/100, green=70/100, blue=28/100, alpha=1)
PastelRed <- rgb(red=100/100, green=41/100, blue=38/100, alpha=1)
PastelYellow <- rgb(red=99/100, green=99/100, blue=59/100, alpha=1)

PastelGreen1 <- rgb(red=218/255, green=241/255, blue=219/255, alpha=1)
PastelGreen2 <- rgb(red=143/255, green=214/255, blue=148/255, alpha=1)
PastelGreen3 <- rgb(red=60/255, green=165/255, blue=68/255, alpha=1)
PastelGreen4 <- rgb(red=14/255, green=37/255, blue=15/255, alpha=1)

PastelViolet1 <- rgb(red=229/255, green=204/255, blue=228/255, alpha=1)
PastelViolet2 <- rgb(red=177/255, green=102/255, blue=174/255, alpha=1)

PastelBlue1 <- rgb(red=22/255, green=232/255, blue=235/255, alpha=1)
PastelBlue2 <- rgb(red=126/255, green=164/255, blue=179/255, alpha=1)

DarkerTurquoise <- rgb(0/255,131/255,133/255)
Maroon <- rgb(128/255,0,0)
