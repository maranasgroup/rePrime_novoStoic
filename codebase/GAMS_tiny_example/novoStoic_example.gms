***************************************************************************************
* novoStoic example
* Author: Akhil Kumar, Lin Wang, Chiam Yu Ng, and Costas D. Maranas
* Pathway design using de novo steps through uncharted biochemical spaces
***************************************************************************************
OPTIONS
       decimals = 8
       solprint = on
       reslim = 20000000
       iterlim = 1000000000
       domlim = 20
       work = 10000000
       lp = cplex
       threads = 8
       sysout = off
       profile = 1
;
$offdigit
$inlinecom /* */
SETS

m moiety
/
$include "sets/moieties.txt"
/

i metabolite
/
$include "sets/metabolites.txt"
/

substrates(i)
/
'C00091'
'C00005'
'C00080'
/

products(i)
/
'ChEBI_41189'
'C00006'
'C00010'
'C00001'
/

EX_met(i)
/
'C00091'
'C00005'
'C00080'
'ChEBI_41189'
'C00006'
'C00010'
'C00001'
/
allow_imb(i)
/
'C00091'
'C00005'
'C00080'
'ChEBI_41189'
'C00006'
'C00010'
'C00001'
'C00232'
/
j rxn set
/
$include "sets/reactions.txt"
/

r rule set
/
$include "sets/rules.txt"
/
;

PARAMETERS
S(i,j) stoichiometry for each reaction j
/
$include "parameters/Sij.txt"
/

T(m,r) 'changes of moiety stoichiometry for each rule r'
/
$include "parameters/Tmr.txt"
/

C(m,i) 'moiety stoichiometry for each metabolite i'
/
$include "parameters/Cmi.txt"
/

LB(j)       lower bound of flux for reaction j

UB(j)       upper bound of flux for reaction j

LB_r(r)       lower bound of flux for rule r

UB_r(r)       upper bound of flux for rule r
;

SCALARS
dG_max        maximum Gibbs free energy /5/
bigM             big M value /1000/
vmax          maximum flux value /2/
;

VARIABLES
z
v(j)

vr(r)
vimb(i)
v_EX(i)
;

Binary VARIABLES
y_rxn(j)
y_rule(r)
;

Equations
obj             objective function for minimizing the number of rules
stoic           stoichiometric balance for metabolite i
moiety          moiety balance for moiety i
con1            lower bound constraints for reaction flux
con2            upper bound constraints for reaction flux
con3            lower bound constraints for rule flux
con4            upper bound constraints for rule flux
con5            exchange flux
con6		exchange flux
con7		exchange flux
con8		exchange flux
con9		exchange flux
con10		exchange flux
con11		exchange flux
con12		control total number of rxns
con13		control imbalance metabolites
;

* objective is to minimize the number of reaction rules
obj..           z =e= sum(r, y_rule(r));

* mass balance of known reactions
stoic(i)..      sum(j, S(i,j) * v(j) ) =e= vimb(i);

* moiety balance of reaction rules
moiety(m)..     sum(r, T(m,r) * vr(r)) + sum(i, C(m,i) * vimb(i)) =e= sum(i$EX_met(i), C(m,i) * v_EX(i));

* set LB and UB for reactiions and rules
con1(j)..       v(j) =g= LB(j)*y_rxn(j);
con2(j)..       v(j) =l= UB(j)*y_rxn(j);
con3(r)..       vr(r) =g= LB_r(r)*y_rule(r);
con4(r)..       vr(r) =l= UB_r(r)*y_rule(r);

*stoichiometric coefficients in the reaction equation
con5..		v_EX('C00091') =e= -1;
con6..		v_EX('ChEBI_41189') =e= 1;
con7..          v_EX('C00005') =e= -4;
con8..          v_EX('C00080') =e= -5;
con9..          v_EX('C00006') =e= 4;
con10..          v_EX('C00010') =e= 1;
con11..          v_EX('C00001') =e= 1;

*control the number of reactions
*con12..		sum(j, y_rxn(j)) =e= 0;

*set the imbalance metabolites
con13(i)..	 	   vimb(i)$(not allow_imb(i)) =e= 0;

*Setting and lower bound (LB) and upper bound (UB) of all reaction and rules
LB(j) = -vmax;
UB(j) = vmax;
LB_r(r) = -vmax;
UB_r(r) = vmax;


Model
novostoic
/
obj
Stoic
moiety
con1
con2
con3
con4
con5
con6
con7
con8
con9
con10
con11
*con12
con13
/
;

novostoic.optfile = 1;
novostoic.holdfixed = 1;

Solve novostoic Using mip Minimizing z;

File file1 /novoStoic_output.txt/;
Put file1;

Put "*** novoStoic output reactions and rules ***"//;
Put "known reactions and flux:"//;
loop(j$(v.l(j) > 0.01),
      put "'"j.tl"'", @50, v.l(j):0:8/;
);
Put //;

Put "reaction rule and flux:"//;
loop(r$(y_rule.l(r) = 1),
      put "'"r.tl"'", @50, vr.l(r):0:8/;
);