
kalman_filter <- function(z, A, H, Q, R) {
  # Initialize state estimate and covariance matrix
  x <- matrix(0, ncol = 1)
  P <- matrix(0, ncol = 1)
  
  # Initialize result matrix to store estimates
  result <- matrix(0, ncol = length(z), nrow = length(x))
  
  for (i in 1:length(z)) {
    # Prediction step
    x <- A %*% x
    P <- A %*% P %*% t(A) + Q
    
    # Update step
    K <- P %*% t(H) %*% solve(H %*% P %*% t(H) + R)
    x <- x + K %*% (z[i] - H %*% x)
    I <- diag(ncol(P)) # Define the identity matrix explicitly
    P <- (I - K %*% H) %*% P
    
    # Store result
    result[, i] <- x
  }
  
  return(result)
}


# Load data
bitcoin_prices <- read.csv("CEX_BTCUSD_d.csv")

# Define parameters
A <- matrix(1, ncol = 1)
H <- matrix(1, ncol = 1)
Q <- matrix(0.01, ncol = 1)
R <- matrix(0.1, ncol = 1)

# Apply Kalman filter
filtered_prices <- kalman_filter(z = head(bitcoin_prices$open), A, H, Q, R)

# Plot filtered prices
plot(filtered_prices, xlab = "Time", ylab = "Price")
