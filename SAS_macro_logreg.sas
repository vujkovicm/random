* DEFINE new directory within SAS;
LIBNAME CCG 'F:\CCG1952\11_All_SNPs';

%macro survival; 
	%LET dsid=%SYSFUNC(OPEN(ccg.char_geno,i));
	%DO  i = 1 %to 144;
	%let var =%SYSFUNC(VARNAME(&dsid,&i*2));
	/*
	%let num = %input(&var, 3.0);
	*/

	/*
	PROC LIFETEST DATA = ccg.char_geno PLOTS=(survival) METHOD = KM cs=none NOTABLE NOPRINT;
	TIME year2event * all_replase(0) ;
	STRATA &var; 
	LABEL year2event = "Event-Free Survival (Years)";
	RUN;
	filename savep "D:\CCG1952\11_All_SNPs\KM curves\&var..gif"; 
	goptions reset=all gsfname= savep gsfmode=replace device=gif;
	*/ 
	
	PROC LIFETEST DATA = ccg.char_geno NOTABLE METHOD=KM;
	TIME year2event * all_replase(0);
	STRATA &var age_cat sex ethnic_origin immunoph /TEST = (logrank);
	LABEL year2event = "Event-Free Survival (Years)";
	RUN;
	
	/*
	proc logistic data=ccg.char_geno;
  	class &var (ref = '0')/ param=ref;
  	model all_replase = &var age_cat sex ethnic_origin immunoph;
	ODS SELECT ODDSRATIOS;
	ODS TRACE ON;
	ODS SHOW;
	run;
	*/

	/*
	proc freq data = ccg.char_geno;
	table &var * all_replase;
	run;
	*/

 	PROC PHREG DATA = ccg.char_geno NOSUMMARY;
	CLASS &var (ref = '0');
	MODEL year2event * all_replase(1) = &var age_cat sex ethnic_origin immunoph /ties = efron RL;
	ODS SELECT PARAMETERESTIMATES;
	ODS TRACE ON;
	ODS SHOW;
	RUN;
%END;
%LET close=%SYSFUNC(CLOSE(&dsid));
%mend;
%survival 


/* 4 SNPS */
/*
rs10925235
rs2842947
rs4712327
rs4819128
*/

/*
proc freq data = ccg.char_geno;
tables rs4819128_T * all_replase;
run;

proc logistic data=ccg.char_geno;
class rs4819128_T(ref = '0')/ param=ref ;
model all_replase = rs4819128_T;
ODS SELECT ODDSRATIOS;
ODS TRACE ON;
ODS SHOW;
run;

proc logistic data=ccg.char_geno;
class rs4819128_T(ref = '0')/ param=ref ;
model all_replase = rs4819128_T sex;
ODS SELECT ODDSRATIOS;
ODS TRACE ON;
ODS SHOW;
run;

proc logistic data=ccg.char_geno;
class rs4819128_T(ref = '0')/ param=ref ;
model all_replase = rs4819128_T age_cat sex ethnic_origin immunoph;
run;
*/
