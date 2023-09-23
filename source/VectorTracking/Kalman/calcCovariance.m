function updated_p_m = calcCovariance(state,p_m,Time_Step,Q,DCM,m,cg,imm)

PHI = calcJacobian(state,Time_Step,m,cg,imm);

% Rotating Jacobian to reflect ECEF
BLKDCM = blkdiag(DCM',DCM',DCM',DCM',[1 0;0 1]);
ROTPHI = BLKDCM'*PHI;
updated_p_m = ROTPHI*p_m*ROTPHI' + Q;

end

