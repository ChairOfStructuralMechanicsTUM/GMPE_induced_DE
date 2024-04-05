Function: ft_2024_DE
This MATLAB function provides ground-motion prediction equations for computing medians and standard deviations of different peak and spectral quantities for induced micro-earthquakes in Southern Germany.

Usage
matlab

[median, sigma, period1] = ft_2024_DE(response_choice, M_choice, M, T, rhyp, region)
Description
response_choice: A string with the requested quantity. Choose from: 'PGAH', 'PGAV', 'PGVH', 'PGVV', 'SPAH', 'SPAV', 'SPVH', 'SPVV', 'SVH', 'SVV'.

M_choice: A string with either 'MW' or 'ML' for the model with respect to Moment Magnitude or Local Magnitude, respectively.

M: Magnitude of the earthquake.

T: Period (in seconds) between 0.01 s and 1 s.

rhyp: Hypocentral distance (in kilometers).

region: A numerical value to specify the region:

0: General (fixed effect)
1: INS
2: G.MUC

Output
median: Median amplitude prediction in g or m/s.

sigma: Log of standard deviation (root mean square error, rmse).

period1: Natural periods in seconds.

Example
% Load coefficients

load R

% Calculate median and sigma for PGAH

[median, sigma, period1] = ft_2024_DE('PGAH', 'MW', 2, 0.1, 5, 1);


Note
Ensure that the necessary coefficients are loaded before calling this function.

Requirements
This function requires the following data:

Regression coefficients stored in the variable R.
Credits
This function was developed as part of the ground-motion prediction research for Southern Germany.

For detailed information on the function and its usage, refer to the source code.
