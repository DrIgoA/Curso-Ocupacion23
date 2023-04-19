load('out.rda')

# PREDICHOS para graficar
vveg.predo <- seq(0,1.97, , 20000) # Covariate values (VOLVEG) for prediction ORGÁNICO (MAX = 1.97)
vveg.predc <- seq(0,1.55, , 20000) # Covariate values (VOLVEG) for prediction CONVENTIONAL (MAX = 1.55)

# ORGÁNICO - PRIMAVERA
pred1 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred1[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predo 
                  + out$mean$alpha[2] * 1
                  + out$mean$alpha[3] * 0
                  #+ out$mean$alpha[4] * vveg.predo
                  #+ out$mean$alpha[5] * vveg.predo * 0
  ) 
}
# CONVENCIONAL - PRIMAVERA
pred2 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred2[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predc 
                  + out$mean$alpha[2] * 0 
                  + out$mean$alpha[3] * 0
                  #+ out$mean$alpha[4] * vveg.predc * 0
                  #+ out$mean$alpha[5] * vveg.predc * 0
  ) 
}
# ORGÁNICO - VERANO
pred3 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred3[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predo 
                  + out$mean$alpha[2] * 1 
                  + out$mean$alpha[3] * 1
                  #+ out$mean$alpha[4] * vveg.predo * 1
                  #+ out$mean$alpha[5] * vveg.predo * 1
  ) 
}
# CONVENCIONAL - VERANO
pred4 <- array(NA, dim = c(20000, 32))
for(i in 1:32){
  pred4[,i]<- exp(out$mean$alpha0[i] 
                  + out$mean$alpha[1] * vveg.predc 
                  + out$mean$alpha[2] * 0 
                  + out$mean$alpha[3] * 1
                  #+ out$mean$alpha[4] * vveg.predc * 0
                  #+ out$mean$alpha[5] * vveg.predc * 1
  ) 
}

cri<-function(x) quantile(x,prob=c(0.05,0.95))
#cri1<-apply(out$q2.5$alpha,1, mean)
#cri2<-apply(out$q97.5$alpha,1, mean)
cri1 <- apply(pred1[,], 1, cri)
cri2 <- apply(pred2[,], 1, cri)
cri3 <- apply(pred3[,], 1, cri)
cri4 <- apply(pred4[,], 1, cri)
str(cri1[1,])
str(cri1[2,])
str(vveg.predo)


### GRÁFICO 
tiff(file = "M4-GLMMPoisson.tiff",                #Guardados como Tiff
     width = 190, height = 150,   #width=140 height=150
     units = "mm", res = 300)
par(mfrow = c(1,2),mar=c(4,4,1,1))

# ORGÁNICO - PRIMAVERA
matplot(vveg.predo, cri1[1,], type = "l", lty = 0, lwd = 0, col = "blue", xlab = "",
        ylab = "Abundancia de artrópodos", frame.plot = F, ylim = c(0, 700)) # Fig. 5.17 (b)
mtext("Primavera",          at=1,side=1, line = -12, cex = 1, padj= -20.5, font=1)
mtext("(a)",             at=0.07,       line=-2,    cex = 1)
polygon(x= c(vveg.predo, rev(vveg.predo)), y= c(cri1[1,], rev(cri1[2,])), 
        col =  adjustcolor("green", alpha.f = 0.15), border = NA)

lines(vveg.predo, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predo
                      + out$mean$alpha[2] * 1
                      + out$mean$alpha[3] * 0
                      #+ out$mean$alpha[4] * vveg.predo
                      #+ out$mean$alpha[5] * vveg.predo * 0
),
col = "darkgreen", lwd = 5)
mtext("Volumen vegetal",       font=1, side=1, cex=1, at=1, line = 1, padj= 1.3, )

# CONVENCIONAL - PRIMAVERA
polygon(x= c(vveg.predc, rev(vveg.predc)), y= c(cri2[1,], rev(cri2[2,])), 
        col =  adjustcolor("magenta", alpha.f = 0.15), border = NA)
lines(vveg.predc, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predc 
                      + out$mean$alpha[2] * 0 
                      + out$mean$alpha[3] * 0
                      #+ out$mean$alpha[4] * vveg.predc * 0
                      #+ out$mean$alpha[5] * vveg.predc * 0
),
col = "darkmagenta", lwd = 5)

# ORGÁNICO - VERANO
matplot(vveg.predo, cri3[1,], type = "l", lty = 0, lwd = 0, col = "darkgreen", xlab = "",
        ylab = "", frame.plot = F, ylim = c(0, 250)) 
#mtext("Organic farming", at=1  ,side=1, line = -4.5, cex = 1, padj= -20.5, font=2)
mtext("Verano", at=0.8,side=1, line = -12, cex = 1, padj= -20.5, font=1)
mtext("(b)",    at=0.07,       line=-2,    cex = 1)
polygon(x= c(vveg.predo, rev(vveg.predo)), y= c(cri3[1,], rev(cri3[2,])), 
        col =  adjustcolor("green", alpha.f = 0.15), border = NA)
#lines(vveg.predo, cri3[2,], lty = 2, lwd = 2, col = "darkgreen")
lines(vveg.predo, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predo 
                      + out$mean$alpha[2] * 1 
                      + out$mean$alpha[3] * 1
                      #+ out$mean$alpha[4] * vveg.predo * 1
                      #+ out$mean$alpha[5] * vveg.predo * 1
),
col = "darkgreen", lwd = 5)
mtext("Volumen vegetal",       font=1, side=1, cex=1, at=1, line = 1, padj= 1.3, )

# CONVENCIONAL - VERANO
polygon(x= c(vveg.predc, rev(vveg.predc)), y= c(cri4[1,], rev(cri4[2,])), 
        col =  adjustcolor("magenta", alpha.f = 0.15), border = NA)
#lines(vveg.predc, cri4[2,], lty = 2, lwd = 2, col = "red")
lines(vveg.predc, exp(out$mean$mu.alpha 
                      + out$mean$alpha[1] * vveg.predc 
                      + out$mean$alpha[2] * 0 
                      + out$mean$alpha[3] * 1
                      #+ out$mean$alpha[4] * vveg.predc * 0
                      #+ out$mean$alpha[5] * vveg.predc * 1
),
col = "darkmagenta", lwd = 5)
legend(0.8,250,legend=c( "Convencional","Orgánico" ), 
       lty=c(1,1),cex=1, col=c("darkmagenta","darkgreen"), bty="n", lwd=1.5,
       text.font=1)
#####
dev.off()
