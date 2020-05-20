function eps = nozzleRatio(Me,k)

eps = (1./Me).*sqrt((2/(k+1))*(1+0.5*(k-1)*Me.^2)).^((k+1)/(k-1));

end