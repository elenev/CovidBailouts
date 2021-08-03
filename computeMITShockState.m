function [res,simnext,transmat,qC]  = computeMITShockState( quantpoint , x, mobj, params, exogenv)

exst = quantpoint(1);
KB = quantpoint(2);
AB = quantpoint(3);
AI = quantpoint(4);
BI = quantpoint(5);
BG = quantpoint(6);

NSOL=mobj.NSOL;
exnpt=exogenv.exnpt;

mu_G = params.mu_G;
delta = params.delta;
deltaK = params.deltaK;
psi = params.psi;

solvec = x(1:NSOL);
KBtrans = x(NSOL+1);
LBtrans = x(NSOL+1+(1:exnpt));
WItrans = x(NSOL+1+exnpt+(1:exnpt));
BGtrans = x(NSOL+1+2*exnpt+(1:exnpt));
thisVI = x(NSOL+1+3*exnpt+1);
LB = x(NSOL+1+3*exnpt+2);
WI = x(NSOL+1+3*exnpt+3);

point = [exst,KB,LB,WI,BG];

%trans=[(1:exnpt)',KBtrans*ones(exnpt,1),LBtrans,WItrans,BGtrans*ones(exnpt,1)];
trans=[KBtrans*ones(exnpt,1);LBtrans;WItrans;BGtrans];
%trans=reshape(trans,exnpt,4);

[nextst,outstr]=mobj.calcStateTransition(point,solvec,0,trans,thisVI,params,exogenv);

qB = exp(solvec(1));
X = solvec(2);
p = 1 + psi*(X/KB-(mu_G-1+deltaK)); 
OmA = outstr.addvars.OmA;
MP = outstr.addvars.MP;

LB_check = qB * AB / ( p * KB );
WI_check = (MP + delta * OmA * qB) * AI + BI;

[fx,~,V]=mobj.calcEquations(point(1),nextst,solvec,outstr,1);
KBtrans_check = V{2}(1);
LBtrans_check = V{2}(exnpt+(1:exnpt))';
WItrans_check = V{2}(2*exnpt+(1:exnpt))';
BGtrans_check = V{2}(3*exnpt+(1:exnpt))';

res = zeros(size(x));
res(1:NSOL) = fx;
res(NSOL+1) = KBtrans - KBtrans_check;
res(NSOL+1+(1:exnpt)) = LBtrans - LBtrans_check;
res(NSOL+1+exnpt+(1:exnpt)) = WItrans - WItrans_check;
res(NSOL+1+2*exnpt+(1:exnpt)) = BGtrans - BGtrans_check;
res(NSOL+1+3*exnpt+1) = thisVI - V{1}(6);
res(NSOL+1+3*exnpt+2) = LB - LB_check;
res(NSOL+1+3*exnpt+3) = WI - WI_check;

if nargout>1
	outstr.addvars.sig2B_zeta_xi_next = outstr.addvars.sig2B_zeta_xi_next(5,2);
	addvec=model.DSGEModel.structToVec(outstr.addvars)';
	valvec=mobj.evaluateVal(point)'; 
	valvec( strcmp(mobj.V_names,'qB') ) = exp( solvec(strcmp(mobj.Sol_names,'qB')) );
	valvec( strcmp(mobj.V_names,'X') ) = solvec(strcmp(mobj.Sol_names,'X'));
	valvec( strcmp(mobj.V_names,'lamB') ) = solvec(strcmp(mobj.Sol_names,'lamB'));
	
	valvec( strcmp(mobj.V_names,'wbill') ) = outstr.addvars.Lscale * (...
		exp(solvec(strcmp(mobj.Sol_names,'wB'))) .* mobj.Params.Lbar(1) + ...
		exp(solvec(strcmp(mobj.Sol_names,'wS'))) .* mobj.Params.Lbar(2) );
	znsvec=[];
    if mobj.NZNS>0
	 	znsvec=mobj.Zfct.evaluateAt(point)';
    end
	% write different categories of variables in one row
	simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec,addvec,znsvec];
	
	transmat = zeros(exnpt,4);
	transmat(:,1) = KBtrans;
	transmat(:,2) = LBtrans;
	transmat(:,3) = WItrans;
	transmat(:,4) = BGtrans;
	
	% Compute consumption prices
	[~,~,V]=mobj.calcEquations(point(1),nextst,solvec,outstr,2);
	SDFB = V{3}.SDFB;
	SDFS = V{3}.SDFS;
	
	Qnext = mobj.Zfct.evaluateAt( nextst )';
	qCBnext = Qnext(:,2);
	qCSnext = Qnext(:,3);
	
	cB = exp( solvec(11) );
	cS = exp( solvec(6)  );
	
	prnext = exogenv.mtrans(exst,:);
	
	qC.qCB = cB + prnext * ( SDFB .* qCBnext );
	qC.qCS = cS + prnext * ( SDFS .* qCSnext );
end

end