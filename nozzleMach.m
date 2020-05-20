function Me = nozzleMach(eps,k,regime)
tol = 0.001;
unsolved = 1;
upFactor = 1.1;
downFactor = 1.2;
Me = zeros(1,length(eps));
for i = 1:length(eps)
    epsCurrent = eps(i);
    if epsCurrent == 1
        machE = 1;
    elseif epsCurrent < 1
        machE = 0;
    else
        if strcmp(regime,'sub')
            Mguess = .5;
            while unsolved
                epsRes = nozzleRatio(Mguess,k);
                if epsRes > (1+tol)*epsCurrent
                    Mguess = Mguess*upFactor;
                elseif epsRes < (1-tol)*epsCurrent
                    Mguess = Mguess/downFactor;
                else
                    machE = Mguess;
                    unsolved = 0;
                end 
            end
        elseif strcmp(regime,'sup')
            Mguess = 2;
            while unsolved
                epsRes = nozzleRatio(Mguess,k);
                if epsRes > (1+tol)*epsCurrent
                    Mguess = Mguess/downFactor;
                elseif epsRes < (1-tol)*epsCurrent
                    Mguess = Mguess*upFactor;
                else
                    machE = Mguess;
                    unsolved = 0;
                end
            end
        else
            error('The "regime" input must be either "sub" or "sup"!\n')
        end

    end
    Me(i) = machE;
    unsolved = 1;
end
end