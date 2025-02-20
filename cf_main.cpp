#include <iostream>
#include <stdio.h>

#include "TString.h"

#include "cf_fitter.h"

int main(int argc, char *argv[])
{
    TString data   = "240620";
    TString target = "test";

    VAR vars;
    set_vars(argc, argv, &vars);
    if (argc > 1) data	 = argv[1];
    if (argc > 2) target = argv[2];

    TString project_path = "/tank1/ge34zez/";
    TString ipath = project_path + "input/" + data + "/";
    TString opath = project_path + "output/" + data + "/" + target + "/";

    TString ifile = ipath + CHARGE_STR[vars.charge] + "/CATS_input_22all_" + CHARGE_STR[vars.charge] + Form("_mTBin_%i_Cent_%i.root", vars.mt, vars.mult);
    TString ifile_pp = ipath + "pp/CATS_input_22all_pp"     + Form("_mTBin_%i_Cent_%i.root", vars.mt, vars.mult);
    TString ifile_aa = ipath + "apap/CATS_input_22all_apap" + Form("_mTBin_%i_Cent_%i.root", vars.mt, vars.mult);

    if (vars.sample)
	printf("\n  \e[1;36m-->  \e[0;36mBootstrap: Sampling CF!  \e[1;36m<--\e[0m\n");
    else
	printf("\n  \e[1;36m-->  \e[0;36mCATS Fitting!  \e[1;36m<--\e[0m\n");

    printf("\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Input Path   │\e[0m   %s\e[0m\n", ipath.Data());
    printf("\e[1;34m  │  Output Path  │\e[0m   %s\e[0m\n", opath.Data());
    printf("\e[1;34m  └───────────────┘\e[0m\n");


    if (vars.system == PP)
    {
	if (vars.charge < BOTH)
	    cf_fitter(project_path, ifile, opath, &vars);
	else if (vars.charge == BOTH)
	    cf_combined_fitter(project_path, ipath, opath, ifile_pp, ifile_aa, &vars);
    }
    else if (vars.system = PL)
    {
	if (vars.charge < BOTH)
	{
	    ifile = "/tank1/ge34zez/input/pl_241129/analysis_pl.root";
	    cf_fitter_pl(project_path, ifile, opath, &vars);
	}
	else if (vars.charge == BOTH)
	{
	    ifile_pp = "/tank1/ge34zez/input/pl_241129/analysis_pl.root";
	    ifile_aa = "/tank1/ge34zez/input/pl_241129/analysis_pl.root";
	    cf_combined_fitter_pl(project_path, opath, ifile_pp, ifile_aa, &vars);
	}
    }

    printf("  \e[1;36m-->  \e[0;36mFitter finished!  \e[1;36m<--\e[0m\n\n");
    return 0;
}
