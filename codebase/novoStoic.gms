***************************************************************************************
* novoStoic 
* Author: Akhil Kumar
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

m atom set
/
$include "sets/atoms.index"
/

targ_m(m) only moieties related to target

i all metabolites
/
$include "sets/metabolite.index"
/

ec all ec numbers
/
$include "sets/ec.index"
/

ecFilter(ec) all ec numbers
/
*$include "sets/gutMeta.ec"
$include "sets/ec.index"
/

k(i) compound set with reactions
/
$include "sets/cpd.index"
/
k2(i) compound set with known price

k3(i) both k and k2
k3_cp(i) copy of k3

k4(i) in k2 but not in k
k5(i) in both k2 and k

substrates(i)
/
'br:mol:137282'
'br:mol:104413'
'br:mol:105482'
'br:mol:103930'
'br:mol:104555'
'br:mol:104080'
'br:mol:100484'
'ci:InChI:38649'
/

products(i)
/
'br:mol:32624'
'ci:InChI:79034'
'ci:InChI:79035'
'ci:InChI:79044'
'ci:InChI:50130'
'ci:InChI:79031'
'ci:InChI:79036'
/

sp(i) contains both substrate and product


j rxn set
/
$include "sets/rule.index"
/

r rule set
/
$include "sets/onlyRule.index"
/
***Linking sets*****
mi(m,i) linking moieties and sets
/
$include "parameters/atoms.formula"
/
ij(i,j)
/
$include "parameters/rule.scj"
/
mr(m,r)
/
$include "parameters/onlyRule.saj"
/

cut index for integer cuts
/1*10000/

l the number of rule combinations
/1*1/
trace(i,i) target start pairs
/**/

rxnEc(j,ec)
/
$include "sets/rule.ec"
/
ruleEc(r,ec)
/
$include "sets/onlyRule.ec"
/

rxnFiltered(j)
ruleFiltered(r)

;
Singleton set targ(i);

PARAMETERS

SIJ(i,j) contains cpd stoichiometry for each rule
/
$include "parameters/rule.scj"
/
SMR(m,r) contains atom stoichiometry for each rule
/
$include "parameters/onlyRule.saj"
/
DMI(m,i) contains atom cardinality for each metabolite
/
$include "parameters/atoms.formula"
/
PRICE(i)
/
$include "parameters/metabolite.price"
/
CARBON(i)
/
$include "parameters/metabolite.carbon"
/
FIX(r) contains atom cardinality for each metabolite
/
$include "parameters/fixMetabolitePerRule.fix"
/
TARGET(i)
******Total moieties*******
TOTALM_i(i)
TOTALM_targ
TOTALM_r(r)
*****Similarity/Dissimilarity*******
*SIM_ii(i,i)
SIM_ii(i)
SIM_ir(r)
SIM_ij(j)
************
ORD_i(i)

store_i(cut,i) store solutions
store_r(cut,r) store solutions
foundSoln(cut) solutionsFound
;

SCALAR
BOUNDS /100/
W_rxn /100/
W_rule /2/
W_met /4/
W_ex /4/
lambda /1/
price_threshold /5/
threshold /20/
ruleCombinations /0/
card_p /1/
card_s /0/
card_ecFilter /0/
;


card_s = card(substrates);
card_p = card(products);
card_ecFilter = card(ecFilter);



PRICE(k)$((Sum(m$mi(m,k),DMI(m,k)) <= 4) and (PRICE(k) = 0)) = eps;

if(card_p eq 0,
   k2(i)=yes$(PRICE(i) ne 0);
);

if(card_s >= 1,
   k2(i)=yes$substrates(i);
);


k3(i) = k(i) + k2(i) + products(i);
*k3(i) = k(i) + k2(i);

k3_cp(i) = k3(i);
k4(i) = k2(i)-k(i);
k5(i) = k2(i)*k(i);

if(card_ecFilter > 0,
   loop(rxnEc(j,ecFilter),
     rxnFiltered(j)=yes;
     ruleFiltered(r)=yes;
   );
);

sp(i) = substrates(i)+products(i);

/*two things left to try out, one is open up the exchanges to other molecules, similarity score to substrates, or use isoPath */

file file0 /k3.tsv/;
file0.pc=6;
put file0;

loop(k3,
     put k3.tl/;
);

$Offorder
ORD_i(i) = ord(i);
$Onorder

TOTALM_i(k3) = sum(m$mi(m,k3),DMI(m,k3));
*TOTALM_i(i) = sum(m$mi(m,i),DMI(m,i));
TOTALM_r(r) = sum(m$mr(m,r),abs(SMR(m,r)));

TARGET(i)=0;
SIM_ii(i)=0;
SIM_ir(r)=0;
SIM_ij(j)=0;

store_i(cut,i) = 0;
store_r(cut,r) = 0;
foundSoln(cut) = 0;

VARIABLES
z


POSITIVE VARIABLES
gm_p(m,l)
gm_n(m,l)

z1
z2
z3
z4
z5

BINARY VARIABLES
v_p(j)
v_n(j)

g_p(r,l)
g_n(r,l)

x_p(i,l)
x_n(i,l)

u_p(i)
u_n(i)

*y_rxn(j)
*y_rule(r,l)
*y_met(i,l)
*y_ex(i)
y_gm(m,l)
y_ctrl(i)
;

gm_p.up(m,l)=BOUNDS;
gm_n.up(m,l)=BOUNDS;

z1.up=BOUNDS*BOUNDS;
z2.up=BOUNDS*BOUNDS;
z3.up=BOUNDS*BOUNDS;
z4.up=BOUNDS*BOUNDS;

*******set all w_pi and w_n to zero where price is zero
*u_p.fx(i)$(not k2(i))=0;
*u_n.fx(i)$(not k2(i))=0;
*u_p.fx(i)$(not k3(i))=0;
*u_n.fx(i)$(not k3(i))=0;
****to force uptakes to be biological origin*****
*u_n.fx(i)$(not k(i))=0;
PRICE('br:mol:1')=eps;

*g_p.fx('KEGG_R00734#0',l)=0;
*g_n.fx('KEGG_R00734#0',l)=0;

y_ctrl.fx(i)$(not k3(i))=0;

if(card_ecFilter > 0,
   v_p.fx(j)$(not rxnFiltered(j))=0;
   v_n.fx(j)$(not rxnFiltered(j))=0;
   g_p.fx(r,l)$(not ruleFiltered(r))=0;
   g_n.fx(r,l)$(not ruleFiltered(r))=0;
);

EQUATIONS
***********************************step1************************************
eq1
eq2
eq2_1
eq3_b
eq3_b1
eq4_b
eq4_b1
eq4_b2
eq5_b
eq5_b1
eq6_b
eq6_b1
eq11
eq11_1
eq11_2
eq13
eq13_1
eq13_2
eq14_21
eq14_22
eq14_31
eq14_32
eq14_11
eq14_12
eq14_41
eq14_42
eq16_1
eq16_2
eq16_3
eq16_4
eq17_1
eq17_2
eq17_3
eq17_4
eq17_5
intcut
intcut_j
objeq_b1
objeq_b2
objeq_b3
objeq_b4
objeq_4
objeq_41
objective
objective2
;
***********************************step1************************************
eq1(k3)..                           Sum(j, SIJ(k3,j)*(v_p(j)-v_n(j))) + Sum(l,(x_p(k3,l)-x_n(k3,l))) + u_p(k3)-u_n(k3)=e=0;
eq2(m,l)..                          Sum(r, SMR(m,r)*(g_p(r,l)-g_n(r,l))) + Sum(k3, DMI(m,k3)*(x_p(k3,l)-x_n(k3,l)))=e=0;
*eq2_1(k3,r,l)$(combination(k3,r) = 1)..  g_p(r,l)+g_n(r,l)+x_p(k3,l)+x_n(k3,l)=l=1;
******to control the total number of rxns******
eq3_b(j)..                          v_p(j)+v_n(j)=l=1;
eq3_b1..                            Sum(j,v_p(j)+v_n(j))=l=W_rxn;
******atmost one rule will be active, for each set*****
eq4_b(r,l)..                        g_p(r,l)+g_n(r,l)=l=1;
eq4_b1(l)..                         Sum(r,g_p(r,l)+g_n(r,l))=l=W_rule;
eq4_b2..                            Sum((r,l),g_p(r,l)+g_n(r,l))=g=1;
******For each rule set, the number of metabolites that will be active in that set is controlled********
eq5_b(k3,l)..                       x_p(k3,l)+x_n(k3,l)=l=1;
eq5_b1(l)..                         Sum(k3,x_p(k3,l)+x_n(k3,l))=l=Sum(r,FIX(r)*g_p(r,l)+g_n(r,l));
******Control the total number of active exchanges******
eq6_b(k3)..                         u_p(k3)+u_n(k3)=l=1;
eq6_b1..                            Sum(k3,u_p(k3)+u_n(k3))=l=W_ex;
******Control based on carbon yield***********
eq11$(card_s eq 0)..               CARBON(targ)*u_p(targ)=g=lambda*Sum(k3$(not targ(k3)),CARBON(k3)*u_p(k3));
eq11_1$(card_s eq 0)..             PRICE(targ)*u_p(targ)=g=lambda*Sum(k3$(not targ(k3)),PRICE(k3)*u_p(k3));
eq11_2$(card_s eq 0)..             Sum(k3,CARBON(k3)*u_p(k3))=l=BOUNDS;
******Set the target > 1********
eq13_1$(card_s >= 1)..             Sum(substrates,u_n(substrates))=g=1;
eq13_2$(card_p >= 1)..             Sum(products,u_p(products))=g=1;
*******contraints mentioned in 14 will move down*****

eq14_21(k3)$(SIM_ii(k3) > threshold)..    u_p(k3)=e=0;
eq14_22(k3)$(SIM_ii(k3) > threshold)..    u_n(k3)=e=0;

eq14_31(j)$(SIM_ij(j) > threshold)..    v_p(j)=e=0;
eq14_32(j)$(SIM_ij(j) > threshold)..    v_n(j)=e=0;

eq14_11(k3,l)$(SIM_ii(k3) > threshold)..   x_p(k3,l)=e=0;
eq14_12(k3,l)$(SIM_ii(k3) > threshold)..   x_n(k3,l)=e=0;

eq14_41(r,l)$(SIM_ir(r) > threshold)..    g_p(r,l)=e=0;
eq14_42(r,l)$(SIM_ir(r) > threshold)..    g_n(r,l)=e=0;



*******Just for testing to see if these contraints speed up the execution*****
eq16_1(targ_m,l)..                       Sum(r, SMR(targ_m,r)*(g_p(r,l)-g_n(r,l)))=e=gm_p(targ_m,l)-gm_n(targ_m,l);
eq16_2(targ_m,l)..                       gm_p(targ_m,l)=l=BOUNDS*y_gm(targ_m,l);
eq16_3(targ_m,l)..                       -1*BOUNDS*(1-y_gm(targ_m,l))=l=-1*gm_n(targ_m,l);
eq16_4..                                  Sum((targ_m,l),gm_p(targ_m,l)+gm_n(targ_m,l))=g=1;
***********these contraints decide the way the network is shaped, this is experimental********
eq17_1(l)..                                Sum(k4,x_p(k4,l)+x_n(k4,l))=l=1;
*eq17_2(l)..                              x_p('ci:InChI:17572',l)+x_n('ci:InChI:17572',l)+g_p('KEGG_R06977#0',l)+g_n('KEGG_R06977#0',l)=l=1;
eq17_3(k3)..                              Sum(l,(x_p(k3,l)+x_n(k3,l)))+u_p(k3)+u_n(k3)=l=1 + y_ctrl(k3);
eq17_4..                                  Sum(k3,y_ctrl(k3))=l=W_met;
eq17_5(k3)$(SIM_ii(k3) > threshold)..     y_ctrl(k3)=e=0;
*************integer cuts, even this shapes the network, need to customize this more**********
intcut(cut)$(foundSoln(cut) = 1)..        Sum((r,l)$(store_r(cut,r) = 1),1-(g_p(r,l)+g_n(r,l))) + Sum((k3,l)$(store_i(cut,k3) = 1),1-(x_p(k3,l)+x_n(k3,l)))=g=1;
*************to prevent reproduction of existing reactions**********
intcut_j(j,l)$(Sum(k3,abs(SIJ(k3,j))) > 1).. Sum(k3$(ij(k3,j)),1-(x_p(k3,l)+x_n(k3,l)))=g=1;
*******objective*****

objeq_b1..                               z1=e=Sum((k3,l),SIM_ii(k3)*(x_p(k3,l)+x_n(k3,l)));
objeq_b2..                               z2=e=Sum(j,SIM_ij(j)*(v_p(j)+v_n(j)));
objeq_b3..                               z3=e=Sum((r,l),SIM_ir(r)*(g_p(r,l)+g_n(r,l)));
objeq_b4..                               z4=e=Sum(k3,SIM_ii(k3)*(u_p(k3)+u_n(k3)));

objeq_4$(card_s eq 0)..                   z5=e=Sum(k3,PRICE(k3)*(u_p(k3)-u_n(k3)));
objeq_41$(card_s eq 0)..                  z5=g=price_threshold;


objective..                            z=e= z1 + z2 + z3 + z4;
objective2..                             z=e= -z5 + eps*(z1 + z2 + z3 + z4);

Model quasiSynth
/
eq1
eq2
*eq2_1
eq3_b
eq3_b1
eq4_b
eq4_b1
eq4_b2
eq5_b
eq5_b1
eq6_b
eq6_b1
eq11
eq11_1
eq11_2
*eq13
eq13_1
eq13_2
eq14_21
eq14_22
eq14_31
eq14_32
eq14_11
eq14_12
eq14_41
eq14_42
*eq16_1
*eq16_2
*eq16_3
*eq16_4
*eq17_1
*eq17_2
eq17_3
eq17_4
eq17_5
intcut
intcut_j
objeq_b1
objeq_b2
objeq_b3
objeq_b4
objeq_4
objeq_41
objective
/;

scalar loopCount /0/;
scalar card_I /0/;
scalar proceedToSolve /0/;
scalar randomVal;
scalar flux /0/;
scalar metflux /0/;
scalar profit /0/;
scalar rules_rxns_total /0/;
scalar solution_found /1/;
scalar solution_count /0/;

*** if i know source and target, the card of card_I is changed here***
*set the scalar src_target = true;

card_I = card(i);



file file1 /optSynth_simvastatin_Ez.tsv/;
file1.pc=6;

put file1;
put "rid","flux","type"/;

while(loopCount < card_I,
      loopCount = loopCount + 1;
      proceedToSolve = 0;
      TARGET(i) = 0;
      targ(i) = no;
      targ_m(m)=no;

      if(card_p = 0,
          TARGET(i)$(k4(i) and (ORD_i(i) eq loopCount)) = 1;
      else
          TARGET(i)$(products(i) and (ORD_i(i) eq loopCount)) = 1;
      );

      proceedToSolve = sum(i,TARGET(i));

      price_threshold = 5;
      solution_found = 1;
      solution_count = 0;
      if(proceedToSolve = 1,

         SIM_ii(k3)=0;
         SIM_ir(r)=0;
         SIM_ij(j)=0;

         targ(i) = yes$(TARGET(i) eq 1);

         targ_m(m)=yes$(DMI(m,targ) ne 0);

         SIM_ii(k3)$(TOTALM_i(k3) <= 4) = 1;

         if(card_s eq 0,
             TOTALM_targ = TOTALM_i(targ);
             SIM_ii(k3) = sum(m$(mi(m,targ) and mi(m,k3)), 2*min(DMI(m,targ),DMI(m,k3)));
             SIM_ii(k3)$(SIM_ii(k3)=0)=10000;
             SIM_ii(k3)$((SIM_ii(k3) > 1) and (SIM_ii(k3) < 10000)) = TOTALM_i(k3) + TOTALM_targ - SIM_ii(k3) + 1;
         else
             SIM_ii(k3)$(SIM_ii(k3)=0) = smin(sp,TOTALM_i(k3) + TOTALM_i(sp) + sum(m$(mi(m,sp) and mi(m,k3)), 2*min(DMI(m,sp),DMI(m,k3)))) + 1;
         );



         if(card_s eq 0,
            TOTALM_targ = TOTALM_i(targ);
            SIM_ir(r) = sum(m$(mi(m,targ) and mr(m,r)), 2*min(DMI(m,targ),abs(SMR(m,r))));
            SIM_ir(r)$(SIM_ir(r)=0)=10000;
            SIM_ir(r)$((SIM_ir(r) > 1) and (SIM_ir(r) < 10000)) = TOTALM_targ + TOTALM_r(r) - SIM_ir(r) + 1;
         else
            SIM_ir(r) = smin(sp,TOTALM_r(r) + TOTALM_i(sp) + sum(m$(mi(m,sp) and mr(m,r)), 2*min(DMI(m,sp),abs(SMR(m,r))))) + 1;
         );



         /*a lot of them will be 1 since water and_ other small molecules will be there*/
         SIM_ij(j)=smax(ij(k3,j),SIM_ii(k3));
         SIM_ij(j)$(SIM_ij(j) > 1)=smin(ij(k3,j)$(SIM_ii(k3) > 1),SIM_ii(k3));

         SIM_ij(j)$(SIM_ij(j)=inf)=10000;
         SIM_ij(j)$(SIM_ij(j)=-inf)=10000;

         while(solution_found = 1,
               SOLVE quasiSynth USING MIP MINIMIZE z;
               if((quasiSynth.modelstat = 1 or quasiSynth.modelstat = 8),
                      if(abs(z.l) > eps,
                         solution_count = solution_count + 1;
                         rules_rxns_total = Sum(j,abs(v_p.l(j)-v_n.l(j)))+Sum((r,l),abs(g_p.l(r,l)-g_n.l(r,l)));
                         if(card_s = 0,
                            profit = Sum(k3,PRICE(k3)*(u_p.l(k3)-u_n.l(k3)));
                         else
                            profit = 0;
                         );


                           put file1;
                           file1.ap = 1;
                           put "solution for ",targ.tl,"below","****",solution_count,"","objective","profit","rxn+rules"/;
                           loop(j$(abs(v_p.l(j)-v_n.l(j)) > eps),
                                flux = v_p.l(j)-v_n.l(j);
                                loop(ij(i,j),
                                     metflux = SIJ(i,j)*flux;
                                     put "0",j.tl,flux,i.tl,metflux,"rxn",z.l,profit,rules_rxns_total/;
                                );
                           );
                           loop(l,
                              loop(r$(abs(g_p.l(r,l)-g_n.l(r,l)) > eps),
                                   flux = g_p.l(r,l)-g_n.l(r,l);
                                   ruleCombinations = ruleCombinations + 1;
                                   foundSoln(cut)$(ord(cut) = ruleCombinations) = 1;
                                   store_r(cut,r)$(ord(cut) = ruleCombinations) = g_p.l(r,l)+g_n.l(r,l);
                                   loop(k3$(abs(x_p.l(k3,l)-x_n.l(k3,l)) > eps),
                                        store_i(cut,k3)$(ord(cut) = ruleCombinations) = x_p.l(k3,l)+x_n.l(k3,l);
                                        metflux = x_p.l(k3,l)-x_n.l(k3,l);
                                        put l.tl,r.tl,flux,k3.tl,metflux,"rule",z.l,profit,rules_rxns_total/;
                                   );
                              );
                           );
                           loop(k3$(abs(u_p.l(k3)-u_n.l(k3)) > eps),
                                metflux = u_p.l(k3)-u_n.l(k3);
                                put "0","ex","0",k3.tl,metflux,"price = ",PRICE(k3),profit,rules_rxns_total/;
                           );
                           loop((k3,l)$(abs(x_p.l(k3,l)-x_n.l(k3,l)) > eps),
                                        store_i(cut,k3)$(ord(cut) = ruleCombinations) = x_p.l(k3,l)+x_n.l(k3,l);
                                        metflux = x_p.l(k3,l)-x_n.l(k3,l);
                                        put l.tl,k3.tl,metflux,"rule ex",z.l,profit,rules_rxns_total/;
                                   );
                           putclose file1;
*                           price_threshold = profit + 1;
                      );
             );
             if((quasiSynth.modelstat ne 1 and quasiSynth.modelstat ne 8),
                 threshold = threshold + 5;
                 if(threshold > 60,
                      solution_found = 0;
                    );
             );
         );

      );
);




