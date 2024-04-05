function [median, sigma, period1] = ft_2024_DE(response_choice,M_choice, M, T, rhyp, region)
% Provides ground-mtion prediction equations for computing medians and
% standard devisations of different peak and spectral quantities, 
% for periodos between 0.01 s and 1 s, for induced micro-earthquakes
% for different regions in Southern Germany
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
% response_choice = a string with the requested quantity, e.g. 'PGAH'
%                   select one of this:'PGAH','PGAV', 'PGVH','PGVV','SPAH',
%                   'SPAV','SPVH','SPVV', 'SVH','SVV';
%                   
% M_choice = a string with either 'MW' or 'ML' for model w.r.t to Moment
%           Magnitude or Local Magntude respectively
%
% M = Magnitude
% T = Period (sec); must be between 0.01 s and 1 s; can be a row vector
% rhyp = hypocentral distance (km)
% region        = 0 for general (fixed effect)
%               = 1 for INS
%               = 2 for G.MUC
% Output Variables
% median        = Median amplitude prediction in g or m/s
%
% sigma         = LOG of standard deviation (rmse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load coefficients
 load R
clear median period1 sigma

% Gravity constant
g=9.81; %m/s^2

% Lists of periods
period = logspace(-2,0,15); % nat. periods in s
N_period = length(period); % number of natural periods

% List of quantities
response_choices = {'PGAH','PGAV','PGVH','PGVV','SPAH','SPAV','SPVH','SPVV', 'SVH','SVV'};

% Original units of the regression coefficients
% stored in R
units = {'cm/s$^2$','cm/s$^2$','cm/s','cm/s','cm/s$^2$','cm/s$^2$','cm/s','cm/s','cm/s','cm/s' };

% Units of median output
units2 = {'g','g','m/s','m/s','g','g','m/s','m/s','m/s','m/s' };

% Check if the chosen quantity matches the possible choices
response_index = find(ismember(response_choices,response_choice));
if isempty(response_index)
    error('No matching quantity. Choices are: PGAH,PGAV,PGVH,PGVV,SPAH,SPAV,SVPH,SPVV,SVH,SVV')
end

% If chosen quantity matches --> store the name and units
response_choice = response_choices{response_index};
unit = units{response_index}; %cm/s^2, cm/s
unit_scales = [1/(100*g),1/(100*g),1/100,1/100,1/(100*g),1/(100*g),1/100,1/100,1/100,1/100]; %g, m/s
unit2 = units2{response_index}; %cm/s^2, cm/s
unit_scale = unit_scales(response_index);

% set region index
I = region;

% switch quantity case
    switch response_choice
       case {'PGAH','PGAV','PGVH','PGVV'}
           % Peak ground quantity
           % Ignore T and set the period index to 1
           period1 = 0;
           np = 1;
    
           % Print out information
           messg = ['Model coefficients for ',response_choice,' in ',unit,':'];
           disp(messg)
           disp(R.(M_choice).(response_choice){1,np}.PHI)
    
           % Evaluate median and standard deviation
           [median, sigma] =ft_2024_DE_sub(response_choice,M_choice, unit_scale,M, np, rhyp, I, R); %in g or m/s
    
           % Print out information
           messg2 = [response_choice,' in ',unit2,':'];
           disp(messg2)
           disp(median)
           
           % Printout to file
            %sprintf('%14.8f ,',real(SVV_sort))
    
       case {'SPAH','SPAV','SPVH','SPVV', 'SVH','SVV'}
           % Spectral quantity

            % Discharge T-values outside from the range of application
            % (0.01 s - 1 s)
            if isempty(find(T>0.01 & T<1))
                error(['All periods outside of the predefined range (0.01 s - 1s)']);
            else
             T_out = T(find(T<0.01 | T>1));   
             T = T(find(T>=0.01 & T<=1));
            
                % Print out discharged periods
                disp(['Discharged periods: ']);
                disp(T_out)
            end

            % Initialize vectors
            median=zeros(1, length(T));
            sigma=zeros(1, length(T));
            period1=T;
            
            for i=1:length(T)
            Ti = T(i);
               if (isempty(find(abs(period-Ti) < 0.0001))) % The user defined period requires interpolation
                    T_low = max(period(period < Ti));
                    T_high = min(period(period > Ti));
                    ip_low  = find(period==T_low);
                    ip_high = find(period==T_high);
                    
                    [Sa_low, sigma_low] = ft_2024_DE_sub(response_choice,M_choice, unit_scale, M, ip_low, rhyp, I, R);
                    [Sa_high, sigma_high] = ft_2024_DE_sub(response_choice,M_choice, unit_scale, M, ip_high, rhyp, I, R);
                    
                    x = [T_low T_high];
                    Y_sa = [Sa_low Sa_high];
                    Y_sigma = [sigma_low sigma_high];
                    median(i) = interp1(x, Y_sa, Ti);
                    sigma(i) = interp1(x, Y_sigma, Ti);
                else
                    ip_T = find(abs((period- Ti)) < 0.0001);
                    [median(i), sigma(i)] = ft_2024_DE_sub(response_choice,M_choice, unit_scale, M, ip_T, rhyp, I, R);
               end
            end

           % Print out information
           messg1 = 'Periods in s.:';
           messg2 = [response_choice,' in ',unit2,':'];
           disp(messg1)
           disp(T)
           disp(messg2)
           disp(median)
    
            % Plot spectrum against periods
            figure
            plot(period1, median,'.')
            xlabel('$T$ in s')
            ylabel([response_choice,' in ',unit2])
            % Plot spectrum against freq
            figure
            xlabel('$f$ in Hz')
            plot(1./period1, median,'.')
            ylabel([response_choice,' in ',unit2])
    end
end

%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%%%%%%%%%%

function [median, sigma] = ft_2024_DE_sub(response_choice,M_choice, unit_scale, M, np, rhyp, I, R)
       if I==0
           % General b_ij=0
           beta2 = R.(M_choice).(response_choice){1,np}.beta2;

           fittedmodel=@(pred)10.^(beta2(1)+beta2(2).*pred(:,1)+...
                    beta2(3).*log10(sqrt(pred(:,2).^2 + max(0.5,10.^(-0.28+0.19.*pred(:,1))).^2)));
       elseif (I == 1) || (I == 2)
           PHI = R.(M_choice).(response_choice){1,np}.PHI;

           fittedmodel=@(pred)10.^(PHI(1,I)+PHI(2,I).*pred(:,1)+...
                    PHI(3,I).*log10(sqrt(pred(:,2).^2 + max(0.5,10.^(-0.28+0.19.*pred(:,1))).^2)));
       else
           error('No matching region. Choices are: 0,1,2')
       end

       median = fittedmodel([M rhyp]).*unit_scale;
       sigma = R.(M_choice).(response_choice){1,np}.stats2.rmse; %The square root of the estimated error variance 
end
