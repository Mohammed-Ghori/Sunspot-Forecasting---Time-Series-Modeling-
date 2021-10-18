#Part 1
install.packages("fpp")
library(fpp)
ssdata = sunspot.month # monthly sunspot data from 1749 to present
plot(ssdata)
head(ssdata)
summary(ssdata)
boxplot(ssdata)

# Decomposition
decomp = stl(ssdata, s.window='periodic')
plot(decomp)
plot(decomp$time.series[,'seasonal'])

# Initial Testing
ssts = ts(ssdata, frequency=12)
acf(ssts)
pacf(ssts)

# Differenced Testing
ssts_diff = diff(ssts, lag=1, differences=1)
plot(ssts_diff)
acf(ssts_diff)
pacf(ssts_diff)

# (3,1,0) Model
ssarima_310 = arima(ssts, order = c(0,1,0))
ssarima_310
acf(ssarima_310$residuals)
pacf(ssarima_310$residuals)
tsdiag(ssarima_310)
Arima(ssts, model=ssarima_310)
plot(ssts)
lines(fitted(Arima(ssts, model=ssarima_310)), col='red')

# (2,1,2) Model
ssarima = auto.arima(ssts,seasonal = TRUE, ic ='aic', trace =TRUE)
ssarima
tsdiag(ssarima)
Arima(sstest, model=sarima)
plot(sstest)
lines(fitted(Arima(sstest, model=sarima)), col='red')

# Seasonal ARIMA Model
sarima = arima(ssts, order = c(2,0,2), seasonal = c(0,1,2))
sarima
tsdiag(sarima)
acf(sarima$residuals)
pacf(sarima$residuals)
Arima(ssts, model=sarima)
plot(ssts)
lines(fitted(Arima(ssts, model=sarima)), col='red')

######### Rolling Forecast

sum1 = 0
for (i in 1589:3177){
  training2 = head(ssts,n=i)
  test2 = tail(ssts,n= (3177-i))
  m1 = sarima
  a = fitted(m1, h = 1) 
  b = test[i-1588]
  if ((a[i-1588]>0 & b>0) | (a[i-1588]<0 & b<0)) {
    sum1 = sum1 + 1
  }
}

print(sum1/3177)

####### Backtest

source('C:\\Users\\Mohammed\\Downloads\\backtest.R')
backtest(ssarima, ssts, 1589, 1, xre=NULL,fixed=NULL,inc.mean=TRUE) # 2,1,2 Model
backtest(ssarima_310, ssts, 1589, 1, xre=NULL,fixed=NULL,inc.mean=TRUE) # 3,1,0 Model
backtest(sarima, ssts, 1589, 1, xre=NULL,fixed=NULL,inc.mean=TRUE) # Seasonal Model

#Part 2
#viewing entire sunspot data
plot(sunspot.month)
z = sunspot.month
adj = log(sunspot.month) 
blah = ts(sunspot.month, start = 1940, end = 2019, frequency = 12) 
plot(stl(whole, s.window = "periodic", t.window = 5), main = "Decomposition of Sunspot Data")
plot(after)
#splitting into two time periods


whole = tail(sunspot.month, n = 885)
whole
plot(whole, main = "Time Plot of Sunspot Data")

before = head(whole, n = 492)
before
plot(before)
range(before)
plot(stl(before, s.window = "periodic", t.window = 301), main = "decomposition of before")

after = tail(whole, n = 393)
after
range(after)
plot(after)
plot(stl(after, s.window = "periodic", t.window = 301), main = "decomposition of after")


more_after = tail(after, n = 165)
more_after
range(more_after)
plot(more_after)



#looking ot ACF and PACF of before and after datasets
acf(before, main = "ACF of Before")
pacf(before, main = "PACF of Before")
acf(after, main = "ACF of After")
pacf(after, main = "PACF of After")


#modeling before dataset
c1 = auto.arima(before, seasonal = TRUE, ic = "aic", trace = TRUE)
c1
tsdiag(c1)
Box.test(r.c1, type = "Ljung-Box")
r.c1 = residuals(c1)
plot(residuals(c1))
acf(r.c1)


#testing other possible models *need to fix*

a1 = arima(before, order = c(1,0,2), seasonal = list(order = c(1,0,0), period=12)) 
a1
r.a1 = residuals(a1)
tsdiag(a1)
Box.test(r.a1, type = "Ljung-Box")


b1 = arima(before, order = c(2,0,2), seasonal = list(order = c(1,0,0), period=12)) 
b1
r.b1 = residuals(b1)
tsdiag(b1)
Box.test(r.b1, type = "Ljung-Box")

d1 = arima(before, order = c(2,0,3), seasonal = list(order = c(1,0,0), period=12)) 
d1
r.d1 = residuals(d1)
tsdiag(d1)
Box.test(r.d1, type = "Ljung-Box")





#modeling after dataset

g1 = auto.arima(after, seasonal = TRUE, ic = "aic", trace = TRUE)
g1
r.g1 = residuals(g1)
tsdiag(g1)
Box.test(r.g1, type = "Ljung-Box")


e1 = arima(after, order = c(1,1,2), seasonal = list(order = c(2,0,1), period=12)) 
e1
r.e1 = residuals(e1)
tsdiag(e1)
Box.test(r.e1, type = "Ljung-Box")

f1 = arima(after, order = c(2,1,2), seasonal = list(order = c(2,0,1), period=12)) 
f1
r.f1 = residuals(f1)
tsdiag(f1)
Box.test(r.f1, type = "Ljung-Box")

h1 = arima(after, order = c(2,1,3), seasonal = list(order = c(2,0,1), period=12)) 
h1
r.h1 = residuals(h1)
tsdiag(h1)
Box.test(r.h1, type = "Ljung-Box")


h2 = arima(after, order = c(2,0,3), seasonal = list(order = c(2,0,1), period=12)) 
h2
r.h2 = residuals(h2)
tsdiag(h2)
Box.test(r.h2, type = "Ljung-Box")


i1 = arima(after, order = c(1,0,3), seasonal = list(order = c(2,0,1), period=12)) 
i1
r.i1 = residuals(i1)
tsdiag(i1)
Box.test(r.i1, type = "Ljung-Box")

j1 = arima(after, order = c(1,1,3), seasonal = list(order = c(1,0,1), period=12)) 
j1
r.j1 = residuals(j1)
tsdiag(j1)
Box.test(r.j1, type = "Ljung-Box")

k1 = arima(after, order = c(2,0,3), seasonal = list(order = c(2,0,1), period=12)) 
k1
r.k1 = residuals(k1)
tsdiag(k1)
Box.test(r.k1, type = "Ljung-Box")

l1 = arima(after, order = c(2,0,3), seasonal = list(order = c(1,0,1), period=12)) 
l1
r.l1 = residuals(k1)
tsdiag(l1)
Box.test(r.l1, type = "Ljung-Box")


#Part 3 - States Data  
data.frame(YearlyData)
avg.spot = YearlyData$'Avg Sunspots'
before = ts(avg.spot, start=1940, end=1980, frequency=12)
plot(stl(before, s.window='periodic', t.window=301), main='Sunspot Yearly Data 1940-1980')
log.sun1 = log(before)
plot(stl(log.sun1, s.window='periodic', t.window=301), main='Log of Sunspot Yearly Data 1940-1980')

after = ts(avg.spot, start=1980, end=2019, frequency=12)
plot(stl(after, s.window='periodic', t.window=301), main='Sunspot Yearly Data 1980-2019')
log.sun2 = log(after)
plot(stl(log.sun2, s.window='periodic', t.window=301), main='Log of Sunspot Yearly Data 1980-2019')

#Minnesota Precipitation 
mn.prec = YearlyData$`MN Precipitation`
mn.prec.before = ts(mn.prec, start=1940, end=1980, frequency=12)
plot(stl(mn.prec.before, s.window='periodic', t.window=301), main='Precipitation in Minnesota 1940-1980')
log.mn.prec1 = log(mn.prec.before)
plot(stl(log.mn.prec1, s.window='periodic', t.window=301), main='Log of Precipitation in Minnesota 1940-1980') 

mn.prec.after = ts(mn.prec, start=1980, end=2019, frequency=12)
plot(stl(mn.prec.after, s.window='periodic', t.window=301), main='Precipitation in Minnesota 1980-2019')
log.mn.prec2 = log(mn.prec.after)
plot(stl(log.mn.prec2, s.window='periodic', t.window=301), main='Log of Precipitation in Minnesota 1980-2019')

#Minnesota Temperature 
mn.temp = YearlyData$`MN Avg Temperature`
mn.temp.before = ts(mn.temp, start=1940, end=1980, frequency=12)
plot(stl(mn.temp.before, s.window='periodic', t.window=301), main='Avg Temperature in Minnesota 1940-1980')
log.mn.temp1 = log(mn.temp.before)
plot(stl(log.mn.temp1, s.window='periodic', t.window=301), main='Log of Avg Temperature in Minnesota 1940-1980')

mn.temp.after = ts(mn.temp, start=1980, end=2019, frequency=12)
plot(stl(mn.temp.after, s.window='periodic', t.window=301), main='Avg Temperature in Minnesota 1980-2019')
log.mn.temp2 = log(mn.temp.after)
plot(stl(log.mn.temp2, s.window='periodic', t.window=301), main='Log of Avg Temperature in Minnesota 1980-2019')

#Alaska Precipitation 
ak.prec = YearlyData$`AK Precipitation`
ak.prec.before = ts(ak.prec, start=1940, end=1980, frequency=12)
plot(stl(ak.prec.before, s.window='periodic', t.window=301), main='Precipitation in Alaska 1940-1980')
log.ak.prec1 = log(ak.prec.before)
plot(stl(log.ak.prec1, s.window='periodic', t.window=301), main='Log of Precipitation in Alaska 1940-1980')

ak.prec.after = ts(ak.prec, start=1980, end=2019, frequency=12)
plot(stl(ak.prec.after, s.window='periodic', t.window=301), main='Precipitation in Alaska 1980-2019')
log.ak.prec2 = log(ak.prec.after)
plot(stl(log.ak.prec2, s.window='periodic', t.window=301), main='Log of Precipitation in Alaska 1980-2019')

#Alaska Temperature
ak.temp = YearlyData$`AK Avg Temperature`
ak.temp.before = ts(ak.temp, start=1940, end=1980, frequency=12)
plot(stl(ak.temp.before, s.window='periodic', t.window=301), main='Avg Temperature in Alaska 1940-1980')
log.ak.temp1 = log(ak.temp.before)
plot(stl(log.ak.temp1, s.window='periodic', t.window=301), main='Log of Avg Temperature in Alaska 1940-1980')

ak.temp.after = ts(ak.temp, start=1980, end=2019, frequency=12)
plot(stl(ak.temp.after, s.window='periodic', t.window=301), main='Avg Temperature in Alaska 1980-2019')
log.ak.temp2 = log(ak.temp.after)
plot(stl(log.ak.temp2, s.window='periodic', t.window=301), main='Log of Avg Temperature in Alaska 1980-2019')

#Louisiana Precipitation 
la.prec = YearlyData$`LA Precipitation`
la.prec.before = ts(la.prec, start=1940, end=1980, frequency=12)
plot(stl(la.prec.before, s.window='periodic', t.window=301), main='Precipitation in Louisiana 1940-1980')
log.la.prec1 = log(la.prec.before)
plot(stl(log.la.prec1, s.window='periodic', t.window=301), main='Log of Precipitation in Louisiana 1940-1980') 

la.prec.after = ts(la.prec, start=1980, end=2019, frequency=12)
plot(stl(la.prec.after, s.window='periodic', t.window=301), main='Precipitation in Louisiana 1980-2019')
log.la.prec2 = log(la.prec.after)
plot(stl(log.la.prec2, s.window='periodic', t.window=301), main='Log of Precipitation in Louisiana 1980-2019')

#Louisiana Temperature 
la.temp = YearlyData$`LA Avg Temperature`
la.temp.before = ts(la.temp, start=1940, end=1980, frequency=12)
plot(stl(la.temp.before, s.window='periodic', t.window=301), main='Avg Temperature in Louisiana 1940-1980')
log.la.temp1 = log(la.temp.before)
plot(stl(log.la.temp1, s.window='periodic', t.window=301), main='Log of Avg Temperature in Louisiana 1940-1980')

la.temp.after = ts(la.temp, start=1980, end=2019, frequency=12)
plot(stl(la.temp.after, s.window='periodic', t.window=301), main='Avg Temperature in Louisiana 1980-2019')
log.la.temp2 = log(la.temp.after)
plot(stl(log.la.temp2, s.window='periodic', t.window=301), main='Log of Avg Temperature in Lousiana 1980-2019')

#Linear Regression
sunprec.mn.1 =lm(before ~ mn.prec.before, data=YearlyData)
plot(sunprec.mn.1, main='Sunspots vs. Precipitation in MN 1940-1980')

#MN Precip. LM
log.sunprec.mn.1 = lm(log.sun1 ~ log.mn.prec1, data=YearlyData)
plot(log.sunprec.mn.1, main='Log of Sunspots vs. MN Precipitation 1940-1980')
log.sunprec.mn.2 = lm(log.sun2 ~ log.mn.prec2, data=YearlyData)
plot(log.sunprec.mn.2, main='Log of Sunspots vs. MN Precipitation 1980-2019')

summary(log.sunprec.mn.1)
summary(log.sunprec.mn.2)

#MN Temp. LM
log.suntemp.mn.1 = lm(log.sun1 ~ log.mn.temp1, data=YearlyData)
plot(log.suntemp.mn.1, main='Log of Sunspots vs. MN Avg Temperature 1940-1980')
log.suntemp.mn.2 = lm(log.sun2 ~ log.mn.temp2, data=YearlyData)
plot(log.suntemp.mn.1, main='Log of Sunspots vs. MN Avg Temperature 1980-2019')

summary(log.suntemp.mn.1)
summary(log.suntemp.mn.2)

#AK Precip. LM
log.sunprec.ak.1 = lm(log.sun1 ~ log.ak.prec1, data=YearlyData)
plot(log.sunprec.ak.1, main='Log of Sunspots vs. AK Precipitation 1940-1980')
log.sunprec.ak.2 = lm(log.sun2 ~ log.ak.prec2, data=YearlyData)
plot(log.sunprec.ak.2, main='Log of Sunspots vs. AK Precipitation 1980-2019')

summary(log.sunprec.ak.1)
summary(log.sunprec.ak.2)

#AK Temp. LM
log.suntemp.ak.1 = lm(log.sun1 ~ log.ak.temp1, data=YearlyData)
plot(log.suntemp.ak.1, main='Log of Sunspots vs. AK Avg Temperature 1940-1980')
log.suntemp.ak.2 = lm(log.sun2 ~ log.ak.temp2, data=YearlyData)
plot(log.suntemp.ak.1, main='Log of Sunspots vs. AK Avg Temperature 1980-2019')

summary(log.suntemp.ak.1)
summary(log.suntemp.ak.2)

#LA Precip. LM
log.sunprec.la.1 = lm(log.sun1 ~ log.la.prec1, data=YearlyData)
plot(log.sunprec.la.1, main='Log of Sunspots vs. LA Precipitation 1940-1980')
log.sunprec.la.2 = lm(log.sun2 ~ log.la.prec2, data=YearlyData)
plot(log.sunprec.la.2, main='Log of Sunspots vs. LA Precipitation 1980-2019')

summary(log.sunprec.la.1)
summary(log.sunprec.la.2)

#LA Temp. LM
log.suntemp.la.1 = lm(log.sun1 ~ log.la.temp1, data=YearlyData)
plot(log.suntemp.la.1, main='Log of Sunspots vs LA Avg Temperature 1940-1980')
log.suntemp.la.2 = lm(log.sun2 ~ log.la.temp2, data=YearlyData)
plot(log.suntemp.la.2, main='Log of Sunspots vs LA Avg Temperature 1980-2019')

summary(log.suntemp.la.1)
summary(log.suntemp.la.2)



