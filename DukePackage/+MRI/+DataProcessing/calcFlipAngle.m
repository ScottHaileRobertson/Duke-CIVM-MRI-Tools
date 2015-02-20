function calcFlipAngle(pfile, dc_sample_idx)
fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1);   % cos theta decay
guess(1)=max(abs(pfile.data(dc_sample_idx,:)));
guess(2)=10*pi/180;       % just guess 10 degrees

xdata=1:size(pfile.data,2);
ydata = abs(pfile.data(dc_sample_idx,:));

[fitparams,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqcurvefit(fitfunct,guess,xdata,ydata);
ci = nlparci(fitparams,residual,jacobian);  % returns 95% conf intervals on fitparams by default
param_err=fitparams-ci(:,1)';

% focus on flip angle
flip_angle=fitparams(2)*180/pi;
flip_err=param_err(2)*180/pi;

disp(['Flip angle ~' num2str(flip_angle) ' (' num2str(flip_err) ' error)']);

figure();
plot(xdata,ydata,'-b');
hold on;
plot(xdata,fitfunct(fitparams,xdata),'-r');
legend('Acquired','Fit');
xlabel('Frame Number');
ylabel('Magnitude');
end