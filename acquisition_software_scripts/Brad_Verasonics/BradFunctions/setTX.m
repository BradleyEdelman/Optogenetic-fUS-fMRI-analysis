function TX=setTX(Resource,P,Trans)

TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', kaiser(Resource.Parameters.numTransmit,1)', ...
                   'focus', 0.0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.numAngle);
               
% - Set event specific TX attributes.
if fix(P.numAngle/2) == P.numAngle/2       % if na even
    P.startAngle = (-(fix(P.numAngle/2) - 1) - 0.5)*P.dtheta;
else
    P.startAngle = -fix(P.numAngle/2)*P.dtheta;
end
for n = 1:P.numAngle   % P.numAngle transmit events
    TX(n).Steer = [(P.startAngle+(n-1)*P.dtheta),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
end