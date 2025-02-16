# Load required libraries
library(ggplot2)      # For data visualization
library(dplyr)        # For data manipulation
library(corrplot)     # For correlation matrix visualization
library(lmtest)       # For statistical tests
library(car)          # For multicollinearity check
library(tseries)      # For time series analysis
library(stargazer)    # For regression output formatting
library(forecast)     # For ARIMA modeling
library(nlme)         # For GLS
library(sandwich)     # For vcovHAC

# Load the dataset
data <- read.csv('C:/Users/ACER/Desktop/Econometrics/project/main.csv')

# Convert date column to Date type
data$date <- as.Date(data$date, format = "%Y-%m-%d")

# Inspect the dataset
summary(data)
str(data)

# Check for missing values in the entire dataset
print(paste("Number of missing values: ", sum(is.na(data))))
colSums(is.na(data))  # Check for missing values in each column

# Handle missing values if necessary
data$Interest_Rate[is.na(data$Interest_Rate)] <- mean(data$Interest_Rate, na.rm = TRUE)

# --------------------------
# View the Distribution of the data
# --------------------------


ggplot(data, aes_string(x = "GDP")) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.7) +
  ggtitle("Distribution of GDP") +
  theme_minimal()


ggplot(data, aes_string(x = "Inflation")) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.7) +
  ggtitle("Distribution of Inflation") +
  theme_minimal()


ggplot(data, aes_string(x = "sp500")) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.7) +
  ggtitle("Distribution of sp500") +
  theme_minimal()


ggplot(data, aes(x = GDP, y = sp500)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  ggtitle("Scatter Plot: GDP vs sp500") +
  xlab("GDP") +
  ylab("sp500") +
  theme_minimal()

ggplot(data, aes(x = GDP, y = Inflation)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  ggtitle("Scatter Plot: GDP vs Inflation") +
  xlab("GDP") +
  ylab("Inflation") +
  theme_minimal()

ggplot(data, aes(x = Inflation, y = sp500)) +
  geom_point(color = "steelblue", alpha = 0.7) +
  ggtitle("Scatter Plot: Inflation vs sp500") +
  xlab("Inflation") +
  ylab("sp500") +
  theme_minimal()

# --------------------------
# Check for Trends and Seasonality
# --------------------------
# Decompose time series data for S&P 500
sp500_ts <- ts(data$sp500, start = c(2000, 1), frequency = 12)
decomp <- decompose(sp500_ts)
plot(decomp)

# Check autocorrelation and partial autocorrelation
acf(sp500_ts, main = "ACF for S&P 500")
pacf(sp500_ts, main = "PACF for S&P 500")

# --------------------------
# Exploratory Data Analysis (EDA)
# --------------------------

# Correlation matrix
cor_matrix <- cor(data[, -1])  # Exclude the date column
print(cor_matrix)
corrplot(cor_matrix, method = "color", tl.cex = 0.6, number.cex = 0.9, addCoef.col = "black")

# Time series plot of S&P 500
ggplot(data, aes(x = date, y = sp500)) +
  geom_line(color = "blue") +
  labs(title = "S&P 500 Over Time", x = "Date", y = "S&P 500 Index")

# --------------------------
# Linear Regression Model
# --------------------------

# Fit the linear regression model
model <- lm(sp500 ~ Interest_Rate + Inflation + GDP + Unemployment + Ind_Prod, data = data)
summary(model)

# Multicollinearity check using Variance Inflation Factor (VIF)
vif_values <- vif(model)
print(vif_values)

# Remove variables with high VIF and refit the model
model_refit <- lm(sp500 ~ Interest_Rate + Inflation + Unemployment + Ind_Prod, data = data)  # Drop GDP if necessary
summary(model_refit)
anova(model, model_refit)

# Residual diagnostics
par(mfrow = c(2, 2))
plot(model)

# --------------------------
# Statistical Tests
# --------------------------

# Autocorrelation: Durbin-Watson test
dw_test <- dwtest(model)
print(dw_test)

gls_model <- gls(sp500 ~ Interest_Rate + Inflation + GDP + Unemployment + Ind_Prod, data = data, correlation = corAR1())
summary(gls_model)
#AIC,BIC->Used to compare model fit; lower values indicate a better model.
#phi->strong Autocorrelation
#log-likelihood->Higher values suggest a better fit

#Homoskedasticity: Breusch-Pagan test
bp_test <- bptest(model)
print(bp_test)

weights <- 1 / lm(abs(residuals(model)) ~ fitted(model))$fitted.values^2 # Calculate weights

wls_model <- lm(sp500 ~ Interest_Rate + Inflation + GDP + Unemployment + Ind_Prod, data = data, weights = weights) # Fit WLS model
summary(wls_model)
anova(model, wls_model)

# Shapiro-Wilk test for normality of residuals
shapiro.test(residuals(model))

# If non-normality exists, apply transformation
data$y_log <- log(data$sp500)  # Log transformation
model_log <- lm(y_log ~ Interest_Rate + log(Inflation) + log(GDP) + Unemployment + Ind_Prod, data = data)
summary(model_log)
anova(model, model_log)

#we can see that using log normal betters the model 

# --------------------------
# Lagged Variables
# --------------------------

# Create lagged variables
#Good idea is to automate the finding of the most optimal n
data$Interest_Rate_Lag <- dplyr::lag(data$Interest_Rate, n = 1)
data$Inflation_Lag <- dplyr::lag(data$Inflation, n = 1)

data$Interest_Rate_Lag[is.na(data$Interest_Rate_Lag)] <- mean(data$Interest_Rate_Lag, na.rm = TRUE)
data$Inflation_Lag[is.na(data$Inflation_Lag)] <- mean(data$Inflation_Lag, na.rm = TRUE)

# Fit the lagged model
model_lagged <- lm(sp500 ~ Interest_Rate_Lag + Inflation_Lag + GDP + Unemployment + Ind_Prod, data = data)
summary(model_lagged)

# Compare models
anova(model, model_lagged)

# --------------------------
# Variance-Covariance Robustness
# --------------------------

cov1 <- sqrt(diag(vcovHAC(model)))         # Robust standard errors for model
cov2 <- sqrt(diag(vcovHAC(model_lagged))) # Robust standard errors for lagged model

# Stargazer output
stargazer(model, model_lagged, align = TRUE, se = list(cov1, cov2),
          omit.stat = c("n"), report = "vc*", type = "text")

# --------------------------
# Time Series Modeling
# --------------------------

# Fit an ARIMA model
fitarima <- auto.arima(data$sp500)
summary(fitarima)

# Forecast
fsp500 <- forecast(fitarima, h = 5)
autoplot(fsp500)

#Оценени стойности на линейния модел
estim_y<-data.frame(estim_y=predict(model))
n<-nrow(estim_y)
ggplot(estim_y,aes(seq(1,n),estim_y))+geom_point()+
  geom_hline(yintercept=0,color="red")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

par(mfrow=c(1,1))
plot(predict(model))

# --------------------------
# Outlier Analysis
# --------------------------

# Cook's distance
cooksd <- data.frame(cook = cooks.distance(model))
cooksd$obs <- rownames(cooksd)

ggplot(cooksd) +
  geom_bar(aes(x = obs, y = cook), stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 4 * mean(cooksd$cook), color = "red") +
  geom_text(data = subset(cooksd, cook > 4 * mean(cooksd$cook)),
            aes(obs, cook, label = obs))


