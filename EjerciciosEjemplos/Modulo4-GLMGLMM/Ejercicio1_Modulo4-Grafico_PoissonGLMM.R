load('out.rda')

# PREDICHOS para graficar
vveg.pred <- seq(0,1.97, , 20000) # Covariate values (VOLVEG) for prediction ORGÁNICO (MAX = 1.97)

pred <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.pred 
                 
  ) 
}

cri<-function(x) quantile(x,prob=c(0.05,0.95))
cri1 <- apply(pred[,], 1, cri)

str(cri1[1,])
str(cri1[2,])
str(vveg.pred)


### GRÁFICO 
tiff(file = "M4-GLMMPoisson.tiff",                #Guardados como Tiff
     width = 120, height = 150,   #width=140 height=150
     units = "mm", res = 300)
par(mfrow = c(1,1),mar=c(4,4,1,1))

#
matplot(vveg.pred, cri1[1,], type = "l", lty = 0, lwd = 0, col = "blue", xlab = "",
        ylab = "Abundancia de artrópodos", frame.plot = F, ylim = c(0, 700)) 
polygon(x= c(vveg.pred, rev(vveg.pred)), y= c(cri1[1,], rev(cri1[2,])), 
        col =  adjustcolor("green", alpha.f = 0.15), border = NA)

lines(vveg.pred, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.pred
),
col = "darkgreen", lwd = 5)
mtext("Volumen vegetal",       font=1, side=1, cex=1, at=1, line = 1, padj= 1.3, )

#####
dev.off()
