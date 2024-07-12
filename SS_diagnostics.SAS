title1 'Choosing Sample Size for Evaluating a Diagnostic Test';
data main;
label sn = 'sensitivity'
      sp = 'specificity'
      p  = 'prevalence of disease'
      w  = 'width of 95% CI'
      a_c = 'diseased subjects TP+FN'
      n1 = 'sample size for SN'
      n2 = 'sample size for SP'
      n  = 'total subjects N';

*step 1 specifications;
sn = 0.90; *substitute your value for sensitivity here;
sp = 0.85; *substitute your value for specificity here;
p  = 0.20; *substitute your value for prevalence here;
w  = 0.10; *substitute your value for width here;

*step 2 calculate TP+FN;
a_c = 1.96*1.96*sn*(1-sn)/(w*w);

*step 3 calculate N1;
n1 = a_c/p;
*round up to the next whole integer;
n1_int = int(n1);
if n1 ne n1_int then n1 = n1_int + 1;

*step 4 calculate FP+TN;
b_d = 1.96*1.96*sp*(1-sp)/(w*w);

*step 5 calculate N2;
n2 = b_d/(1-p);
*round up to the next whole integer;
n2_int = int(n2);
if n2 ne n2_int then n2 = n2_int + 1;

* step 6 get final sample size;
if n1 gt n2 then n=n1;
else if n2 gt n1 then n=n2;
else if n1=n2 then n=n1;

* print the sample size;
proc print label noobs;
var w sn sp p n1 n2 n;
run;
