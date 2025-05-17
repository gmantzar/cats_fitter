#include <iostream>
#include <stdio.h>

#include "TString.h"

#include "cf_fitter.h"

int main(int argc, char *argv[])
{
    VAR vars;
    set_vars(argc, argv, &vars);

    if (vars.sample)
	printf("\n  \e[1;36m-->  \e[0;36mBootstrap: Sampling CF!  \e[1;36m<--\e[0m\n");
    else
	printf("\n  \e[1;36m-->  \e[0;36mCATS Fitting!  \e[1;36m<--\e[0m\n");

    string settings = (vars.system == PP)? "pp." : "pl.";
    TString ipath  = vars.get(settings + "input").data();
    TString opath  = vars.get(settings + "output").data();

    printf("\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Input Path   │\e[0m   %s\e[0m\n", ipath.Data());
    printf("\e[1;34m  │  Output Path  │\e[0m   %s\e[0m\n", opath.Data());
    printf("\e[1;34m  └───────────────┘\e[0m\n");


    if (vars.system == PP)
    {
	if (vars.charge < BOTH)
	    cf_fitter(&vars);
	else if (vars.charge == BOTH)
	    cf_combined_fitter(&vars);
    }
    else if (vars.system = PL)
    {
	if (vars.charge < BOTH)
	{
	    cf_fitter_pl(&vars);
	}
	else if (vars.charge == BOTH)
	{
	    cf_combined_fitter_pl(&vars);
	}
    }

    printf("  \e[1;36m-->  \e[0;36mFitter finished!  \e[1;36m<--\e[0m\n\n");
    return 0;
}
