
# a multiplier generating a decreasing trend after 10 years
trend<-c(rep(1,10),0.85,0.85^2,0.85^3,0.85^4,0.85^4)
# a random process, with the mean affected by the dreasing trend
TAC<-rnorm(15,trend*800,150)

banking     <-  - 0.1 * TAC # amont banked   each year if scenario banking
borrowing   <-    0.1 * TAC # amont borrowed each year if scenario borrowing


# scenario banking : TAC - amount banked the same year + amount banked the previous year
TACbanking    <-    TAC + banking   - c(0,banking[-length(TAC)])

# borrowing : TAC + amount borrow from next year - repay borrowed last year
TACborrowing  <-    TAC + borrowing - c(0,borrowing[-length(TAC)])


plot(TAC,type="l")
lines(TACbanking,col="red")
lines(TACborrowing,col="green")

