/* includes */
#include <iostream>
#include <boost/property_tree/ptree.hpp>

#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Fitters.h"
#include "DLM_HistoAnalysis.h"
#include "DLM_RootWrapper.h"
#include "CommonAnaFunctions.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TString.h"
#include "TObjString.h"
#include "TGraph.h"

#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/Fitter.h"
#include "HFitInterface.h"

#include "cf_fitter.h"

/* defines  */
#define CUSP_WEIGHT 0.33
#define LAM_PL_FLAT 0.4969
#define CUTOFF_Q    200

#define POT_PP	    "AV18"
#define POT_PL	    "Chiral_Coupled_SPD"
#define POT_PS	    "DG_NLO19"
#define SOURCE	    "Gauss"
#define SOURCE_RSM  "McLevy_ResoTM"

#define out(X) cout << X << endl

using namespace std;

/* fit classes */
double bslfun(double *mom, double *par)
{
    double &k = *mom;

    double &N = par[0];
    double &a = par[1];
    double &b = par[2];
    double &c = par[3];
    double &d = par[4];

    return N*(1. + a*1. + b*k + c*pow(k, 2) + d*pow(k, 3.));
}

void baseline_fit(TH1F *histo, vector<double> &parameters, pair<double, double> xrange)
{
    // N = par[0]; a = par[1]; b = par[2]; c = par[3]; d = par[4];
    // pol = N * (1. + a*1. + b*x + c*x² + d*x³)
    auto npol = [](double *x, double par[], int degree)
    {
	double total = par[0];
	for (int order = 0; order <= degree; ++order) { total += par[0]*par[order + 1]*pow(*x, order); }
	return total;
    };

    // 3rd degree polynomial without 2nd degree term
    auto cpol3 = [&](double *x, double par[]) { par[2] = 0; return npol(x, par, 3); };

    auto fitter = make_unique<TF1>(TF1("prefitter", cpol3, xrange.first, xrange.second, 5));
    fitter->FixParameter(2, 0);
    histo->Fit(fitter.get(), "S, N, R, M");

    for (size_t npar = 0; npar < 5; ++npar) parameters.push_back(fitter->GetParameter(npar));
}

class plFitFunctionSimple
{
    public:
	DLM_CkDecomposition *mydecomp;
	TH1F *cf = nullptr;
	bool prefit = false;
	vector<double> prefit_pars;

	plFitFunctionSimple(DLM_CkDecomposition *comp) : mydecomp(comp) {}

	plFitFunctionSimple(DLM_CkDecomposition *comp, TH1F *hData, bool bsl_prefit)
	    : mydecomp(comp), prefit(bsl_prefit)
	{
	    cf = static_cast<TH1F*>(hData->Clone("correlation_function"));
	    if (prefit) baseline_fit(cf, prefit_pars, {240, 340});
	}

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    double &a = par[1];
	    double &b = par[2];
	    double &c = par[3];
	    double &d = par[4];
	    double &radius = par[5];

	    mydecomp->GetCk()->SetSourcePar(0, radius);
	    mydecomp->Update(false, false);

	    return N*(1. + a*1. + b*k + c*pow(k, 2) + d*pow(k, 3.)) * mydecomp->EvalCk(k);
	}
};

class FitFunctionSimplePLSB
{
    public:
	DLM_CkDecomposition *mydecomp;
	TF1 *sideband;
	double lam_flat;
	FitFunctionSimplePLSB(DLM_CkDecomposition *comp, double lam, TF1 *sb) :
	    mydecomp(comp), lam_flat(lam), sideband(sb) {}

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    double &a = par[1];
	    double &b = par[2];
	    double &c = par[3];
	    double &d = par[4];
	    double &radius = par[5];

	    mydecomp->GetCk()->SetSourcePar(0, radius);
	    mydecomp->Update(false, false);

	    return N*(1. + a*k + b*k*k + c*k*k*k + d*k*k*k*k) * (mydecomp->EvalCk(k) + lam_flat * (sideband->Eval(k) - 1));
	}
};

class FitFunSimpleAvg
{
    public:
	DLM_CkDecomposition *mydecomp = nullptr;
	TH1F *cf = nullptr;
	TH1F *me = nullptr;

	bool prefit = false;
	vector<double> prefit_pars;

	FitFunSimpleAvg(TH1F *hData, TH1F *hME, DLM_CkDecomposition *comp, bool bsl_prefit = false)
	    : mydecomp(comp), prefit(bsl_prefit)
	{
	    cf = static_cast<TH1F*>(hData->Clone("correlation_function"));
	    me = static_cast<TH1F*>(hME->Clone("mixed_event"));

	    if (prefit) baseline_fit(cf, prefit_pars, {300, 500});
	}

	FitFunSimpleAvg(TH1F **hData, TH1F **hME, DLM_CkDecomposition *comp, bool bsl_prefit = false)
	    : mydecomp(comp), prefit(bsl_prefit)
	{
	    cf = static_cast<TH1F*>((*hData)->Clone("correlation_function"));
	    me = static_cast<TH1F*>((*hME)->Clone("mixed_event"));

	    if (prefit) baseline_fit(cf, prefit_pars, {240, 400});
	}

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    double &a = par[1];
	    double &b = par[2];
	    double &c = par[3];
	    double &d = par[4];
	    double &radius = par[5];

	    mydecomp->GetCk()->SetSourcePar(0, radius);
	    mydecomp->Update(true, true);

	    double cf_avg = 0;
	    double me_sum = 0;
	    int bin = cf->GetXaxis()->FindBin(k);

	    double value_low = cf->GetXaxis()->GetBinLowEdge(bin);
	    double value_up  = cf->GetXaxis()->GetBinLowEdge(bin + 1);
	    int me_bin_low = me->GetXaxis()->FindBin(value_low);
	    int me_bin_up  = me->GetXaxis()->FindBin(value_up);

	    for (int nbin = me_bin_low; nbin < me_bin_up; ++nbin)
	    {
		cf_avg += me->GetBinContent(nbin) * mydecomp->EvalCk(me->GetBinCenter(nbin));
		me_sum += me->GetBinContent(nbin);
	    }
	    double ck_eval = (me_sum < 1e-6)? 0 : cf_avg / me_sum;

	    return N*(1. + a*1. + b*k + c*pow(k, 2) + d*pow(k, 3.)) * ck_eval;
	}
};

class GenuineFun
{
    public:
	DLM_CkDecomposition *mydecomp;

	GenuineFun(DLM_CkDecomposition *comp) { mydecomp = comp; }
	void set_decomp(DLM_CkDecomposition *comp) { mydecomp = comp; }

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &N = par[0];
	    return N*(mydecomp->GetCk()->Eval(k));
	}
};

class GenuineFunSmear
{
    public:
	DLM_CkDecomposition *mydecomp;

	GenuineFunSmear(DLM_CkDecomposition *comp) { mydecomp = comp; }
	void set_decomp(DLM_CkDecomposition *comp) { mydecomp = comp; }

	double operator() (double *mom, double *par)
	{
	    double &k = *mom;
	    double &lambda = par[0];

	    return lambda * (mydecomp->EvalSmearedMain(k)) + (1. - lambda);
	}
};

class GlobalChi2
{
    public:
	int *fit_pars_pp;
	int *fit_pars_aa;

	ROOT::Math::IMultiGenFunction *fChi2_1;
	ROOT::Math::IMultiGenFunction *fChi2_2;

	GlobalChi2(ROOT::Math::IMultiGenFunction &f1, ROOT::Math::IMultiGenFunction &f2, int *par_pp, int *par_aa)
	{
	    fit_pars_pp = par_pp;
	    fit_pars_aa = par_aa;

	    fChi2_1 = &f1;
	    fChi2_2 = &f2;
	}

	double operator() (const double *par) const
	{
	    double par1[6], par2[6];
	    for (int i = 0; i < 6; ++i)
	    {
		par1[i] = par[fit_pars_pp[i]];
		par2[i] = par[fit_pars_aa[i]];
	    }

	    return (*fChi2_1)(par1) + (*fChi2_2)(par2);
	}
};
/*-------------*/

void print_info_2(VAR *var, VAR_FMR *fmr, VAR_FR *fr)
{
    if (var->sample)
    {
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │ Random Sample │\e[0m   %i\n\e[0m", var->sample);
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");
    }

    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  System       │\e[0m   %i  (%s)\n\e[0m", var->system, SYSTEM_STR[var->system].Data());
    printf("\e[1;34m  │  Charge       │\e[0m   %i  (%s)\n\e[0m", var->charge, CHARGE_STR[var->charge].Data());
    printf("\e[1;34m  │  Bin mt       │\e[0m   %i\n\e[0m", var->mt);
    printf("\e[1;34m  │  Bin mult     │\e[0m   %i\n\e[0m", var->mult);
    printf("\e[1;34m  ├───────────────┤\e[0m\n");
    printf("\e[1;34m  │  LambdaPar    │\e[0m   %i  (%s)\n\e[0m", var->lam, (var->lam == 1)? "-10%" : (var->lam == 2)? "+10%" : "default");
    printf("\e[1;34m  │  Femtorange   │\e[0m   %i  (%.0f, %.0f)\n\e[0m", var->fmr, fmr->Min, fmr->Max);
    printf("\e[1;34m  │  Fitrange     │\e[0m   %i  (%.0f, %.0f)\n\e[0m", var->fr, fr->Min, fr->Max);
    printf("\e[1;34m  │  Baseline     │\e[0m   %i  (%s)\n\e[0m", var->bsl, (var->bsl)? "Pol2" : "Pol3");
    printf("\e[1;34m  │  Smearing     │\e[0m   %i  (%s)\n\e[0m", var->smear, (var->smear == 1)? "Run3" : (var->smear == 2)? "Run2" : "None");

    if (var->rsm)
    {
    printf("\e[1;34m  ├───────────────┤\e[0m\n");
    printf("\e[1;34m  │  RSM Frac     │\e[0m   %i  (%.4f)\n\e[0m", var->frac, RSM_FRAC[var->system][var->frac]);
    printf("\e[1;34m  │  RSM Mass     │\e[0m   %i  (%.1f)\n\e[0m", var->mass, RSM_MASS[var->system][var->mass]);
    }
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");
}

void gev_to_mev(unique_ptr<TH1F> &hist)
{
    unique_ptr<TH1F> hist_gev (static_cast<TH1F*>(hist->Clone("temporary")));

    int bin_ent = hist->GetNbinsX();
    int bin_width = 1; // MeV

    hist.reset(new TH1F(TString(hist->GetName()) + "_rescaled", hist->GetTitle(), bin_ent, 0, bin_ent*bin_width));

    for (size_t nbin = 0; nbin < bin_ent; ++nbin)
    {
	hist->SetBinContent(nbin, hist_gev->GetBinContent(nbin));
	hist->SetBinError(nbin, hist_gev->GetBinError(nbin));
    }
}

string settings_filename(VAR *var)
{
    string file = "";
    switch (var->system)
    {
	case PP: file = var->get("pp.input") + var->get("pp.file"); break;
	case PL: file = var->get("pl.input") + var->get("pl.file"); break;
	default:
	    cerr << "Error: Missing Input Histograms!" << endl;
    }
    return file;
}

string settings_potential(VAR *var)
{
    string potential = "";
    switch (var->system)
    {
	case PP: potential = var->get("pp.potential"); break;
	case PL: potential = var->get("pl.potential"); break;
    }
    return potential;
}

string transform_lower(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return tolower(c); });
    return s;
}

string transform_upper(string s)
{
    transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return toupper(c); });
    return s;
}

/* histogram getter functions */
void get_klambda_histos(TString filename, VAR *var, TH1F *lambda_pars[])
{
    // k* dependent λ parameters
    TString klambda_filename = filename;
    auto file = make_unique<TFile>(klambda_filename, "read");
    auto lambdas = make_unique<TDirectory*[]>(5);
    auto lambdas_th1 = make_unique<TH1F*[]>(5);
    lambdas[0] = static_cast<TDirectory*>(file->GetDirectory("genuine"));
    lambdas[1] = static_cast<TDirectory*>(file->GetDirectory("FD_lambda"));
    lambdas[2] = static_cast<TDirectory*>(file->GetDirectory("FD_sigma"));
    lambdas[3] = static_cast<TDirectory*>(file->GetDirectory("flat"));
    lambdas[4] = static_cast<TDirectory*>(file->GetDirectory("fake"));

    lambdas[1]->SetName("lambda");
    lambdas[2]->SetName("sigma");

    for (size_t nlam = 0; nlam < 5; ++nlam)
    {
	TString name = Form("Mult_%i/Mt_%i/Lambda", var->mult, var->mt + 1);
	lambdas_th1[nlam] = static_cast<TH1F*>(lambdas[nlam]->Get(name));
	lambdas_th1[nlam]->SetName(lambdas[nlam]->GetName());
	lambdas_th1[nlam]->SetDirectory(0);
    }

    //double scale = 0.1;
    double scale = 0.05;
    double total_gen_var = 0;
    double scale_non_gen = 0;
    double content = 0;

    for (int nbin = 1; nbin < lambdas_th1[0]->GetNbinsX() + 1; ++nbin)
    {
	total_gen_var = lambdas_th1[0]->GetBinContent(nbin) * scale;
	scale_non_gen = total_gen_var / (lambdas_th1[1]->GetBinContent(nbin) + lambdas_th1[2]->GetBinContent(nbin)
		    + lambdas_th1[3]->GetBinContent(nbin) + lambdas_th1[4]->GetBinContent(nbin));

	for (size_t nlam = 0; nlam < 5; ++nlam)
	{
	    if (nlam == 0)
	    {
		if  (var->lam == 1)
		{
		    content = lambdas_th1[nlam]->GetBinContent(nbin) * (1. - scale);
		    lambdas_th1[nlam]->SetBinContent(nbin, content);
		}
		else if (var->lam == 2)
		{
		    content = lambdas_th1[nlam]->GetBinContent(nbin) * (1. + scale);
		    lambdas_th1[nlam]->SetBinContent(nbin, content);
		}
	    }
	    else
	    {
		if  (var->lam == 1)
		{
		    content = lambdas_th1[nlam]->GetBinContent(nbin) * (1. + scale_non_gen);
		    lambdas_th1[nlam]->SetBinContent(nbin, content);
		}
		else if (var->lam == 2)
		{
		    content = lambdas_th1[nlam]->GetBinContent(nbin) * (1. - scale_non_gen);
		    lambdas_th1[nlam]->SetBinContent(nbin, content);
		}
	    }
	}
    }

    for (size_t nlam = 0; nlam < 5; ++nlam)
    {
	lambda_pars[nlam] = static_cast<TH1F*>(lambdas_th1[nlam]->Clone());
	lambda_pars[nlam]->SetDirectory(0);
    }

    file->Close();
}

void get_klambda_histos(TString filename, VAR *var, DLM_Histo<double> *lambda_pars[])
{
    auto lambdas_th1 = make_unique<TH1F*[]>(5);
    get_klambda_histos(filename, var, lambdas_th1.get());

    for (size_t nlam = 0; nlam < 5; ++nlam)
	lambda_pars[nlam] = Convert_TH1F_DoubleDlmHisto(static_cast<TH1F*>(lambdas_th1[nlam]->Clone(TString(lambdas_th1[nlam]->GetName()) + "_tmp")));
}

void get_input_histos(TString ifile, TH1F **input_cf, TH1F **input_me, TH1F **input_me_orig)
{
    auto input_file = make_unique<TFile>(ifile, "read");
    *input_cf = (TH1F*) input_file->Get("CF");
    *input_me = (TH1F*) input_file->Get("ME");
    *input_me_orig = (TH1F*) input_file->Get("ME_original");
    (*input_cf)->GetXaxis()->SetRangeUser(0., 700.);
    (*input_cf)->SetDirectory(0);
    (*input_me)->SetDirectory(0);
    (*input_me_orig)->SetDirectory(0);
    input_file->Close();
}

void get_input_histos(TString ifile, unique_ptr<TH1F> &input_cf,
	unique_ptr<TH1F> &input_me, unique_ptr<TH1F> &input_me_orig)
{
    auto input_file = make_unique<TFile>(ifile, "read");
    input_file->ls();
    input_cf.reset(static_cast<TH1F*>(input_file->Get("CF")));
    input_me.reset(static_cast<TH1F*>(input_file->Get("ME")));
    input_me_orig.reset(static_cast<TH1F*>(input_file->Get("ME_original")));
    input_cf->GetXaxis()->SetRangeUser(0., 700.);
    input_cf->SetDirectory(0);
    input_me->SetDirectory(0);
    input_me_orig->SetDirectory(0);
    input_file->Close();
}

void get_input_pp(TString ifile, unique_ptr<TH1F> &input_cf,
	unique_ptr<TH1F> &input_me, unique_ptr<TH1F> &input_me_orig, VAR *var, int charge)
{
    auto input_file = make_unique<TFile>(ifile, "read");

    auto mt = str_mt_bins_pp;
    auto mult = str_mult_bins;

    TString tdir_name = (charge == APAP) ? "apap/" : "pp/";
    tdir_name += "Rebin_8_Dim_1-" + mt[var->mt] + "_Dim_2-0-100_Dim_3-" + mult[var->mult];

    unique_ptr<TDirectory> tdir (static_cast<TDirectory*>(input_file->Get(tdir_name)));
    input_cf.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("CF/CF_Reweighted_rescaled"))->Clone("CF")));
    input_me.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("ME/ME_Reweighted_rescaled"))->Clone("ME")));
    input_me_orig.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("No_Rebin/ME"))->Clone("ME_original")));

    gev_to_mev(input_me_orig);

    input_cf->GetXaxis()->SetRangeUser(0., 700.);

    input_cf->SetDirectory(0);
    input_me->SetDirectory(0);
    input_me_orig->SetDirectory(0);
}

void get_input_pl(TString ifile, TH1F **input_cf, TH1F **input_me, TH1F **input_me_orig, VAR *var, int charge)
{
    auto input_file = make_unique<TFile>(ifile, "read");

    auto mt = str_mt_bins_pl;
    auto mult = str_mult_bins;

    TString tdir_name = (charge == APAP) ? "APAL/" : "PL/";
    tdir_name += "Rebin_8_Dim_1-" + mt[var->mt] + "_Dim_3-" + mult[var->mult];

    unique_ptr<TDirectory> tdir (static_cast<TDirectory*>(input_file->Get(tdir_name)));
    *input_cf = static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("CF/CF_Reweighted_rescaled"))->Clone("CF"));
    *input_me = static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("ME/ME_Reweighted_rescaled"))->Clone("ME"));
    *input_me_orig = static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("No_Rebin/ME/ME"))->Clone("ME_original"));

    (*input_cf)->GetXaxis()->SetRangeUser(0., 700.);

    (*input_cf)->SetDirectory(0);
    (*input_me)->SetDirectory(0);
    (*input_me_orig)->SetDirectory(0);
}

void get_input_pl(TString ifile, unique_ptr<TH1F> &input_cf,
	unique_ptr<TH1F> &input_me, unique_ptr<TH1F> &input_me_orig, VAR *var, int charge)
{
    auto input_file = make_unique<TFile>(ifile, "read");

    auto mt = str_mt_bins_pl;
    auto mult = str_mult_bins;

    //TString tdir_name = (charge == APAP) ? "apal/" : "pl/";
    TString tdir_name = (charge == APAP) ? "APAL/" : "PL/";
    //tdir_name += "Rebin_10_Dim_1-" + mt[var->mt] + "_Dim_2-0-100_Dim_3-" + mult[var->mult];
    tdir_name += "Rebin_10_Dim_1-" + mt[var->mt] + "_Dim_3-" + mult[var->mult];

    unique_ptr<TDirectory> tdir (static_cast<TDirectory*>(input_file->Get(tdir_name)));
    input_cf.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("CF/CF_Reweighted_rescaled"))->Clone("CF")));
    input_me.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("ME/ME_Reweighted_rescaled"))->Clone("ME")));
    input_me_orig.reset(static_cast<TH1F*>(static_cast<TH1F*>(tdir->Get("No_Rebin/ME"))->Clone("ME_original")));
    gev_to_mev(input_me_orig);

    input_cf->GetXaxis()->SetRangeUser(0., 700.);

    input_cf->SetDirectory(0);
    input_me->SetDirectory(0);
    input_me_orig->SetDirectory(0);
}

/*
void get_sideband_fit(string ifile,  unique_ptr<TF1> &sideband, string which_band, VAR_FR *fr, VAR *var, int charge)
{
    auto mt = str_mt_bins_pl;
    auto mult = str_mult_bins;

    auto input_file = make_unique<TFile>(ifile, "read");

    TString tdir_name = (charge == APAP) ? "APAL_sideBand" : "PL_sideBand";
    tdir_name += (!(transform_lower(which_band)).compare("high"))? "High" : "Low";
    tdir_name += "Rebin_10_Dim_1-" + mt[var->mt] + "_Dim_3-" + mult[var->mult];

    unique_ptr<TDirectory> tdir (static_cast<TDirectory*>(input_file->Get(tdir_name)));
    unique_ptr<TH1F> sideband_th1 (static_cast<TH1F*>(tdir->Get("CF/CF_Reweighted_rescaled")));
    sideband.reset(new TF1("sideband", "pol5", fr->Min, fr->Max));

    sideband_th1->Fit(sideband.get());
}
*/

void get_sample_histos(TString fname_cf, TString fname_syst, TH1F **input_cf, TH1F **input_me, TH1F **input_me_orig, VAR *var, int charge)
{
    get_input_histos(fname_cf, input_cf, input_me, input_me_orig);
    auto file_syst = make_unique<TFile>(fname_syst, "read");
    TString name_syst = Form("femto-dream-pair-task-track-track_std/mt_%i/mult_1/rebin_8/syst_th1", var->mt);

    unique_ptr<TH1F> input_syst (static_cast<TH1F*>(file_syst->Get(name_syst)));
    unique_ptr<TH1F> sample_cf (static_cast<TH1F*>((*input_cf)->Clone(Form("CF_sample_%s", CHARGE_STR[charge].Data()))));

    sample_cf->SetTitle(sample_cf->GetName());
    sample_cf->Reset();

    TRandom *random = new TRandom3(var->sample);
    double gaus, uniform, mean, sigma, syst;
    for (size_t n = 1; n < sample_cf->GetNbinsX(); ++n)
    {
	mean  = (*input_cf)->GetBinContent(n);
	sigma = (*input_cf)->GetBinError(n);
	syst  = input_syst->GetBinContent(n);

	gaus = random->Gaus(mean, sigma);
	uniform = random->Uniform(gaus - (syst / 2), gaus + (syst / 2));

	sample_cf->SetBinContent(n, uniform);
	if (var->stat)
	    sample_cf->SetBinContent(n, gaus);
	sample_cf->SetBinError(n, sigma);
    }
    *input_cf = static_cast<TH1F*>(sample_cf->Clone());
    (*input_cf)->SetDirectory(0);
}

void get_sample_histos(TString fname_cf, TString fname_syst, unique_ptr<TH1F> &input_cf,
	unique_ptr<TH1F> &input_me, unique_ptr<TH1F> &input_me_orig, VAR *var, int charge)
{
    get_input_pp(fname_cf, input_cf, input_me, input_me_orig, var, charge);
    auto file_syst = make_unique<TFile>(fname_syst, "read");
    // +2 for some reason the systematics file starts at 1 and the first bin is skipped
    TString name_syst = Form("femto-dream-pair-task-track-track_std/mt_%i/mult_1/rebin_8/syst_th1", var->mt + 2);

    unique_ptr<TH1F> input_syst (static_cast<TH1F*>(file_syst->Get(name_syst)));
    unique_ptr<TH1F> sample_cf (static_cast<TH1F*>(input_cf->Clone(Form("CF_sample_%s", CHARGE_STR[charge].Data()))));

    sample_cf->SetTitle(sample_cf->GetName());
    sample_cf->Reset();

    TRandom *random = new TRandom3(var->sample);
    double gaus, uniform, mean, sigma, syst;
    for (size_t n = 1; n < sample_cf->GetNbinsX(); ++n)
    {
	mean  = input_cf->GetBinContent(n);
	sigma = input_cf->GetBinError(n);
	syst  = input_syst->GetBinContent(n);

	gaus = random->Gaus(mean, sigma);
	uniform = random->Uniform(gaus - (syst / 2), gaus + (syst / 2));

	sample_cf->SetBinContent(n, uniform);
	if (var->stat)
	    sample_cf->SetBinContent(n, gaus);
	sample_cf->SetBinError(n, sigma);
    }
    input_cf.reset(static_cast<TH1F*>(sample_cf->Clone()));
    input_cf->SetDirectory(0);
}

template <typename T>
void get_input(TString fname_cf, T &input_cf, T &input_me, T &input_me_orig, VAR *var, int charge)
{
    if (var->sample)
    {
	string settings = (var->system == PP)? "pp.input" : "pl.input";
	TString name_syst = TString(var->get(settings)) + Form("UFFA_syst_%s_MultPercentile%i.root", CHARGE_STR[charge].Data(), var->mult);
	get_sample_histos(fname_cf, name_syst, input_cf, input_me, input_me_orig, var, charge);
    }
    else
    {
	get_input_pp(fname_cf, input_cf, input_me, input_me_orig, var, charge);
    }
}

void get_resolution_matrix(TString ipath, TH2F *&matrix, VAR *var, int charge)
{
    TString name_file, name_matrix;
    TFile *file;

    if (var->smear == 1)
    {
        name_file = "Resolution.root";
	name_matrix = "resolution_matrix_" + CHARGE_STR[charge];
    }
    else if (var->smear == 2)
    {
        name_file = "ALICE_pp_13TeV_MEpp.root";
        name_matrix = "h_RESO_pp_MeV";
    }

    if (var->smear)
    {
	auto file = make_unique<TFile>(ipath + "input/shared/" + name_file, "read");
	matrix = (TH2F*) file->Get(name_matrix);
        matrix->SetDirectory(0);

	printf("\e[1;34m  ┌───────────────┐\e[0m\n");
	printf("\e[1;34m  │  Res Matrix   │   \e[0m%s\n", name_file.Data());
	printf("\e[1;34m  │  Hist Name    │   \e[0m%s\n", name_matrix.Data());
	printf("\e[1;34m  └───────────────┘\e[0m\n\n");

	file->Close();
    }
    else
	matrix = NULL;
}

void get_corrected_cf(TH1F *cfs[], TH1F *lambdas[], double bsl_pars[], TH1F **corrected_cf)
{
    enum CF_ID { RAW = 0, LAM = 1, SIG = 2 };
    enum LAM_ID { GEN = 0, FLAT = 3, FAKE = 4 };

    double value_x, value_y, bsl, raw, lam, sig, flat;
    double lam_gen, lam_lam, lam_sig, lam_flat, lam_fake;
    for (size_t nbin = 1; nbin < (*corrected_cf)->GetNbinsX(); ++nbin)
    {
	value_x = cfs[RAW]->GetBinCenter(nbin);
	bsl = bsl_pars[0] * (1 + bsl_pars[2]*pow(value_x, 2) + bsl_pars[3]*pow(value_x, 3));

	raw = cfs[RAW]->GetBinContent(nbin);
	lam = cfs[LAM]->GetBinContent(nbin);
	sig = cfs[SIG]->GetBinContent(nbin);
	lam_gen = lambdas[GEN]->GetBinContent(nbin);
	lam_lam = lambdas[LAM]->GetBinContent(nbin);
	lam_sig = lambdas[SIG]->GetBinContent(nbin);
	lam_flat = lambdas[FLAT]->GetBinContent(nbin);
	lam_fake = lambdas[FAKE]->GetBinContent(nbin);

	value_y = (1/lam_gen)*((raw/bsl) - lam_lam*lam - lam_sig*sig - lam_flat - lam_fake);

	if (value_y > 10e-6 && value_y < 50)
	    (*corrected_cf)->SetBinContent(nbin, value_y);
    }
}
/*----------------------------*/

/* CATS getter functions for wave histos */
void get_cats_totalwave(CATS *cats, TH1F **container, double kstar, int channel)
{
    //unsigned bins = cats->GetNumMomBins();
    double kstar_bin = cats->GetMomBin(kstar);
    double xmin = 0, xmax = 30, bins = (xmax - xmin)*100;
    double radius;

    TString title = TString("hTotalWF_spin") + channel + "_" + static_cast<int>(kstar) + "mev";
    auto histo = make_unique<TH1F>(title, title, bins, xmin, xmax);
    for (size_t nbin = 0; nbin < bins; ++nbin)
    {
	radius = histo->GetBinCenter(nbin);
	histo->SetBinContent(nbin, cats->EvalWaveFun2(kstar_bin, radius, channel));
    }

    *container = static_cast<TH1F*>(histo->Clone());
    (*container)->SetDirectory(0);
}

void get_cats_totalwaves(CATS *cats, TH1F *container[], vector<double> kstars, int channel)
{
    for (size_t nkstar = 0; nkstar < kstars.size(); ++nkstar)
    {
	get_cats_totalwave(cats, &container[nkstar], kstars[nkstar], channel);
    }
}

void get_cats_radialwave(CATS *cats, TH1F **container, double kstar, TString target, TString dimension)
{
    int channel = 0, part_wave = 0;
    if	    (!strcmp("1S0", target)) { part_wave = 0; channel = 0; }
    else if (!strcmp("1P1", target)) { part_wave = 1; channel = 0; }
    else if (!strcmp("1D2", target)) { part_wave = 2; channel = 0; }
    else if (!strcmp("3P0", target)) { part_wave = 1; channel = 1; }
    else if (!strcmp("3P1", target)) { part_wave = 1; channel = 2; }
    else if (!strcmp("3P2", target)) { part_wave = 1; channel = 3; }
    else    out("Partial Wave \"" << target << "\" not supported!");

    int dim = 1;
    if (!strcmp("re",	dimension)) dim = 0;
    if (!strcmp("real",	dimension)) dim = 0;

    //unsigned bins = cats->GetNumMomBins();
    double kstar_bin = cats->GetMomBin(kstar);
    double bin_width = 0.1;
    double xmin = 0, xmax = 24;
    double bins = (xmax - xmin) / bin_width;
    double radius;
    bool divide_by_r = true;

    TString title = TString("hRadialWF_") + target + "_" + static_cast<int>(kstar) + "mev";
    auto histo = make_unique<TH1F>(title, title, bins, xmin, xmax);
    for (size_t nbin = 0; nbin < bins; ++nbin)
    {
	radius = histo->GetBinCenter(nbin);

	if (dim) histo->SetBinContent(nbin, cats->EvalRadialWaveFunction(kstar_bin, channel, part_wave, radius, divide_by_r).imag());
	else	 histo->SetBinContent(nbin, cats->EvalRadialWaveFunction(kstar_bin, channel, part_wave, radius, divide_by_r).real());
    }

    *container = static_cast<TH1F*>(histo->Clone());
    (*container)->SetDirectory(0);
}

void get_cats_radialwaves(CATS *cats, TH1F *container[], vector<double> kstars, TString target, TString dimension)
{
    for (size_t nkstar = 0; nkstar < kstars.size(); ++nkstar)
    {
	get_cats_radialwave(cats, &container[nkstar], kstars[nkstar], target, dimension);
    }
}

void get_cats_radialwave_epelbaum(CATS *cats, TH1F **container, double kstar, TString target, TString dimension)
{
    int channel = 0, part_wave = 0;
    if	    (!strcmp("1S0", target)) { part_wave = 0; channel = 0; }
    else if (!strcmp("1P1", target)) { part_wave = 1; channel = 0; }
    else if (!strcmp("1D2", target)) { part_wave = 2; channel = 0; }
    else if (!strcmp("3P0", target)) { part_wave = 1; channel = 1; }
    else if (!strcmp("3P1", target)) { part_wave = 1; channel = 4; }
    else if (!strcmp("3P2", target)) { part_wave = 1; channel = 7; }
    else    out("Partial Wave \"" << target << "\" not supported!");

    int dim = 1;
    if (!strcmp("re",	dimension)) dim = 0;
    if (!strcmp("real",	dimension)) dim = 0;

    //unsigned bins = cats->GetNumMomBins();
    double kstar_bin = cats->GetMomBin(kstar);
    double bin_width = 0.1;
    double xmin = 0, xmax = 200;
    double bins = (xmax - xmin) / bin_width;
    double radius;
    bool divide_by_r = true;

    TString title = TString("hRadialWF_") + target + "_" + static_cast<int>(kstar) + "mev";
    auto histo = make_unique<TH1F>(title, title, bins, xmin, xmax);
    for (size_t nbin = 0; nbin < bins; ++nbin)
    {
	radius = histo->GetBinCenter(nbin);

	//if (dim) histo->SetBinContent(nbin, cats->EvalRadialWaveFunction(kstar_bin, channel, part_wave, radius, divide_by_r).imag());
	//else	 histo->SetBinContent(nbin, cats->EvalRadialWaveFunction(kstar_bin, chanonel, part_wave, radius, divide_by_r).real());
	histo->SetBinContent(nbin, cats->EvalRadialWaveFunction(kstar_bin, channel, part_wave, radius, divide_by_r).real());
    }

    *container = static_cast<TH1F*>(histo->Clone());
    (*container)->SetDirectory(0);
}

void get_cats_radialwaves_epelbaum(CATS *cats, TH1F *container[], vector<double> kstars, TString target, TString dimension)
{
    for (size_t nkstar = 0; nkstar < kstars.size(); ++nkstar)
    {
	get_cats_radialwave_epelbaum(cats, &container[nkstar], kstars[nkstar], target, dimension);
    }
}

void get_cats_phaseshift_epelbaum(CATS *cats, TH1F **container, TString target, double xmin, double xmax)
{
    int channel = 0, part_wave = 0;
    if	    (!strcmp("1S0", target)) { part_wave = 0; channel = 0; }
    else if (!strcmp("1P1", target)) { part_wave = 1; channel = 0; }
    else if (!strcmp("1D2", target)) { part_wave = 2; channel = 0; }
    else if (!strcmp("3P0", target)) { part_wave = 1; channel = 1; }
    else if (!strcmp("3P1", target)) { part_wave = 1; channel = 4; }
    else if (!strcmp("3P2", target)) { part_wave = 1; channel = 7; }
    else    out("Partial Wave \"" << target << "\" not supported!");

    //unsigned bins = cats->GetNumMomBins();
    double bin_width = 0.1;
    double bins = (xmax - xmin) / bin_width;
    double momentum;

    TString title = TString("hPhase_") + target;
    auto histo = make_unique<TH1F>(title, title, bins, xmin, xmax);
    for (size_t nbin = 0; nbin < bins; ++nbin)
    {
	momentum = histo->GetBinCenter(nbin + 1);
	histo->SetBinContent(nbin + 1, 180. / TMath::Pi() * cats->EvalPhaseShift(momentum, channel, part_wave));
    }

    *container = static_cast<TH1F*>(histo->Clone());
    (*container)->SetDirectory(0);
}

template <typename T>
void get_cats_phaseshifts_epelbaum(CATS *cats, TH1F *container[], vector<T> targets, double xmin, double xmax)
{
    for (size_t ntarget = 0; ntarget < targets.size(); ++ntarget)
    {
	get_cats_phaseshift_epelbaum(cats, &container[ntarget], targets[ntarget], xmin, xmax);
    }
}

template <typename T>
void get_cats_radialwaves(CATS *cats, TH1F *container[], double kstar, vector<T> targets, TString dimension)
{
    for (size_t ntarget = 0; ntarget < targets.size(); ++ntarget)
    {
	get_cats_radialwave(cats, &container[ntarget], kstar, targets[ntarget], dimension);
    }
}

void get_cats_phaseshift(CATS *cats, TH1F **container, TString target, double xmin, double xmax)
{
    int channel = 0, part_wave = 0;
    if	    (!strcmp("1S0", target)) { part_wave = 0; channel = 0; }
    else if (!strcmp("1P1", target)) { part_wave = 1; channel = 0; }
    else if (!strcmp("1D2", target)) { part_wave = 2; channel = 0; }
    else if (!strcmp("3P0", target)) { part_wave = 1; channel = 1; }
    else if (!strcmp("3P1", target)) { part_wave = 1; channel = 2; }
    else if (!strcmp("3P2", target)) { part_wave = 1; channel = 3; }
    else    out("Partial Wave \"" << target << "\" not supported!");

    //unsigned bins = cats->GetNumMomBins();
    double bins = (xmax - xmin)*10;
    double momentum;

    TString title = TString("hPhase_") + target;
    auto histo = make_unique<TH1F>(title, title, bins, xmin, xmax);
    for (size_t nbin = 0; nbin < bins; ++nbin)
    {
	momentum = histo->GetBinCenter(nbin);
	histo->SetBinContent(nbin, 180. / TMath::Pi() * cats->EvalPhaseShift(momentum, channel, part_wave));
    }

    *container = static_cast<TH1F*>(histo->Clone());
    (*container)->SetDirectory(0);
}

template <typename T>
void get_cats_phaseshifts(CATS *cats, TH1F *container[], vector<T> targets, double xmin, double xmax)
{
    for (size_t ntarget = 0; ntarget < targets.size(); ++ntarget)
    {
	get_cats_phaseshift(cats, &container[ntarget], targets[ntarget], xmin, xmax);
    }
}

void get_cats_potential(CATS *cats, TH1F **pot_hist, double kstar, TString target)
{
    int channel = 0, part_wave = 0;
    if	    (!strcmp("1S0", target)) { part_wave = 0; channel = 0; }
    else if (!strcmp("1P1", target)) { part_wave = 1; channel = 0; }
    else if (!strcmp("1D2", target)) { part_wave = 2; channel = 0; }
    else if (!strcmp("3P0", target)) { part_wave = 1; channel = 1; }
    else if (!strcmp("3P1", target)) { part_wave = 1; channel = 2; }
    else if (!strcmp("3P2", target)) { part_wave = 1; channel = 3; }
    else    out("Partial Wave \"" << target << "\" not supported!");

    double xmin = 0, xmax = 8, bins = (xmax - xmin)*100 - 1;
    double radius;

    TString title = TString("hPot_") + target + "_" + static_cast<int>(kstar) + "mev";
    auto histo = make_unique<TH1F>(title, title, bins, xmin, xmax);
    for (size_t nbin = 1; nbin < bins; ++nbin)
    {
	radius = histo->GetBinCenter(nbin);
	histo->SetBinContent(nbin, cats->EvaluateThePotential(channel, part_wave, kstar, radius));
    }

    *pot_hist = static_cast<TH1F*>(histo->Clone());
    (*pot_hist)->SetDirectory(0);
}

void get_cats_potentials(CATS *cats, TH1F *pot_histos[], vector<double> kstars, TString target)
{
    for (size_t nkstar = 0; nkstar < kstars.size(); ++nkstar)
    {
	get_cats_potential(cats, &pot_histos[nkstar], kstars[nkstar], target);
    }
}

template <typename T>
void get_cats_potentials(CATS *cats, TH1F *pot_histos[], double kstar, vector<T> pot_names)
{
    for (size_t npot = 0; npot < pot_names.size(); ++npot)
    {
	get_cats_potential(cats, &pot_histos[npot], kstar, pot_names[npot]);
    }
}
/*---------------------------------------*/

/* CATS setup functions */
void setup_epos(DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm, int restype, const char *file)
{
    auto epos_file = make_unique<TFile>(file, "read");
    unique_ptr<TNtuple> tntuple (static_cast<TNtuple*>(epos_file->Get("InfoTuple_ClosePairs")));

    Float_t k_D, fP1, fP2, fM1, fM2, Tau1, Tau2;
    Float_t AngleRcP1, AngleRcP2, AngleP1P2;
    double RanVal1, RanVal2, RanVal3;
    DLM_Random RanGen(11);

    vector<pair<string, Float_t *>> values = {
	pair{"k_D", &k_D},
	pair{"P1", &fP1}, pair{"P2", &fP2},
	pair{"M1", &fM1}, pair{"M2", &fM2},
	pair{"Tau1", &Tau1}, pair{"Tau2", &Tau2},
	pair{"AngleRcP1", &AngleRcP1}, pair{"AngleRcP2", &AngleRcP2}, pair{"AngleP1P2", &AngleP1P2}
    };

    for (auto &branch : values)
    {
	tntuple->SetBranchAddress(branch.first.data(), branch.second);
    }

    for (size_t uEntry = 0; tntuple->GetEntry(uEntry); ++uEntry)
    {
	if (k_D > CUTOFF_Q) continue;

	RanVal1 = RanGen.Exponential(var_rsm->mass / (fP1 * var_rsm->tau));
	RanVal2 = RanGen.Exponential(var_rsm->mass2 / (fP2 * var_rsm->tau2));

	if	(restype == PR) MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2));
	else if (restype == RP) MagicSource->AddBGT_RP(RanVal1, -cos(AngleRcP1));
	else if (restype == RR) MagicSource->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
    }
}

void setup_rsm_pp(DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm)
{
    /* RSM Source Setup */
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, var_rsm->frac);
    MagicSource->SetUpReso(1, var_rsm->frac);
    MagicSource->InitNumMcIter(1000000);

    setup_epos(MagicSource, var_rsm, PR, var_rsm->epos + "EposDisto_p_pReso.root");
    setup_epos(MagicSource, var_rsm, RR, var_rsm->epos + "EposDisto_pReso_pReso.root");
}

void setup_rsm_pl(DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm)
{
    /* RSM Source Setup */
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, var_rsm->frac);
    MagicSource->SetUpReso(1, var_rsm->frac2);
    MagicSource->InitNumMcIter(1000000);

    setup_epos(MagicSource, var_rsm, PR, var_rsm->epos + "EposDisto_p_LamReso.root");
    setup_epos(MagicSource, var_rsm, RP, var_rsm->epos + "EposDisto_pReso_Lam.root");
    setup_epos(MagicSource, var_rsm, RR, var_rsm->epos + "EposDisto_pReso_LamReso.root");
}

void setup_rsm(DLM_CleverMcLevyResoTM *MagicSource, VAR *var, VAR_RSM *var_rsm)
{
    switch (var->system)
    {
	case PP: setup_rsm_pp(MagicSource, var_rsm); break;
	case PL: setup_rsm_pl(MagicSource, var_rsm); break;
    }
}

template <typename T>
void setup_cats(DLM_CommonAnaFunctions *setupper, CATS *cats, VAR_FMR *range, T t_target, VAR *var, DLM_CleverMcLevyResoTM *MagicSource, VAR_RSM *var_rsm)
{
    TString target = t_target;
    TString source = "Gauss";
    cout<<target<<endl;
    //target = "av18";
    cats->SetMomBins(range->Bins, range->Min, range->Max);
    if (var != NULL && var->rsm)
    {
	source = "";
	setup_rsm(MagicSource, var, var_rsm);
	cats->SetAnaSource(CatsSourceForwarder, MagicSource, 2);
	cats->SetAnaSource(0, 1.0);
	cats->SetAnaSource(1, 2.0);
	cats->SetUseAnalyticSource(true);
    }
    if	    (target == "pp")		{ setupper->SetUpCats_pp(  *cats, "AV18",   source, 0, 0); }
    else if (target == "reid93")	{ setupper->SetUpCats_pp(  *cats, "ReidV8", source, 0, 0); }
    else if (target == "reid68")	{ setupper->SetUpCats_pp(  *cats, "ReidSC", source, 0, 0); }
    else if (target == "epelbaum")	{ setupper->SetUpCats_pp(  *cats, "Applebaum", source, 0, 0); }
    else if (target == "psp")		{ setupper->SetUpCats_pSp( *cats, "DG_NLO19", source, 0, 0); }
    else if (target == "ps0")		{ setupper->SetUpCats_pS0( *cats, "Chiral", source); }
    else if (target == "pxm")		{ setupper->SetUpCats_pXim(*cats, "pXim_HALQCDPaper2020", source); }
    else if (target == "px0")		{ setupper->SetUpCats_pXi0(*cats, "pXim_HALQCDPaper2020", source); }
    else if (target == "px1530")	{ setupper->SetUpCats_pXim(*cats, "pXim1530", source); }
    else if (strstr(target.Data(), "av18"))
    {
	auto cPars = make_unique<CATSparameters>(CATSparameters::tSource, 1, true);
	cPars->SetParameter(0, 1.2);

	cats->SetMomentumDependentSource(false);
	cats->SetThetaDependentSource(false);
	cats->SetAnaSource(GaussSource, *cPars);
	cats->SetUseAnalyticSource(true);
	cats->SetExcludeFailedBins(false);

	// #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
	int POT_FLAG = v18_Coupled3P2;
	double PotPars1S0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 0, 0, 0};
	double PotPars3P0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 0};
	double PotPars3P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 1};
	double PotPars3P2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 2};
	double PotPars1D2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 0, 2, 2};

	auto cPotPars1S0 = make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars3P0 = make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars3P1 = make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars3P2 = make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);
	auto cPotPars1D2 = make_unique<CATSparameters>(CATSparameters::tPotential, 8, true);

	cPotPars1S0->SetParameters(PotPars1S0);
	cPotPars3P0->SetParameters(PotPars3P0);
	cPotPars3P1->SetParameters(PotPars3P1);
	cPotPars3P2->SetParameters(PotPars3P2);
	cPotPars1D2->SetParameters(PotPars1D2);

	cats->SetQ1Q2(1);
	cats->SetPdgId(2212, 2212);
	cats->SetRedMass(0.5 * Mass_p);

	if (!cats->GetNumChannels())
	{
	    cats->SetNumChannels(4);
	    cats->SetNumPW(0, 1);
	    cats->SetNumPW(1, 2);
	    cats->SetNumPW(2, 2);
	    cats->SetNumPW(3, 2);
	    cats->SetSpin(0, 0);
	    cats->SetSpin(1, 1);
	    cats->SetSpin(2, 1);
	    cats->SetSpin(3, 1);
	    cats->SetChannelWeight(0, 3. / 12.);
	    cats->SetChannelWeight(1, 1. / 12.);
	    cats->SetChannelWeight(2, 3. / 12.);
	    cats->SetChannelWeight(3, 5. / 12.);
	}

	if (strstr(target.Data(), "av18_s"))
	    cats->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
	if (strstr(target.Data(), "av18_sp"))
	{
	    cats->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
	    cats->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
	    cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
	}
    }
    else if (target == "reid68_b")
    {
	CATSparameters *cats_pars = NULL;

	CATSparameters *cPotPars1S0 = NULL;
	CATSparameters *cPotPars1P1 = NULL;

	CATSparameters *cPotPars3S1 = NULL;
	CATSparameters *cPotPars3P0 = NULL;
	CATSparameters *cPotPars3P1 = NULL;
	CATSparameters *cPotPars3P2 = NULL;
	CATSparameters *cPotPars1D2 = NULL;

	cats_pars = new CATSparameters(CATSparameters::tSource, 1, true);
	cats_pars->SetParameter(0, 1.2);
	cats->SetAnaSource(GaussSource, *cats_pars);
	cats->SetUseAnalyticSource(true);

	double PotPars1S0[8] = {pp_ReidOli, 0};
	double PotPars3P0[8] = {pp_ReidOli, 1};
	double PotPars3P1[8] = {pp_ReidOli, 1};
	double PotPars3P2[8] = {pp_ReidOli, 1};
	cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 2, true);
	cPotPars1S0->SetParameters(PotPars1S0);
	cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 2, true);
	cPotPars3P0->SetParameters(PotPars3P0);
	//cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 2, true);
	//cPotPars3P1->SetParameters(PotPars3P1);
	//cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 2, true);
	//cPotPars3P2->SetParameters(PotPars3P2);

	cats->SetMomentumDependentSource(false);
	cats->SetExcludeFailedBins(false);

	cats->SetQ1Q2(1);
	cats->SetPdgId(2212, 2212);
	cats->SetRedMass(0.5 * Mass_p);

	cats->SetNumChannels(4);
	cats->SetNumPW(0, 3);
	cats->SetNumPW(1, 2);
	cats->SetNumPW(2, 2);
	cats->SetNumPW(3, 2);
	cats->SetSpin(0, 0);
	cats->SetSpin(1, 1);
	cats->SetSpin(2, 1);
	cats->SetSpin(3, 1);
	cats->SetChannelWeight(0, 3. / 12.);
	cats->SetChannelWeight(1, 1. / 12.);
	cats->SetChannelWeight(2, 3. / 12.);
	cats->SetChannelWeight(3, 5. / 12.);

	if (cPotPars1S0) cats->SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
	if (cPotPars1P1) cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars1P1);
	if (cPotPars1D2) cats->SetShortRangePotential(0, 2, fDlmPot, *cPotPars1D2);
	if (cPotPars3S1) cats->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3S1);
	if (cPotPars3S1 && cPotPars3P1) cats->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3S1);
	if (cPotPars3S1 && cPotPars3P2) cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
	if (cPotPars3P0) cats->SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
	if (cPotPars3P1) cats->SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
	if (cPotPars3P2) cats->SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
    }
    else if (target == "pl")
    {
	setupper->SetUpCats_pL(*cats, "Chiral_Coupled_SPD", source, 11600, 0);
	cats->SetChannelWeight( 7,  1./4. * CUSP_WEIGHT);	//1S0 SN(s) -> LN(s)
	cats->SetChannelWeight( 8,  3./4. * CUSP_WEIGHT);	//3S1 SN(s) -> LN(s)
	cats->SetChannelWeight(10,  3./4. * CUSP_WEIGHT);	//3S1 SN(d) -> LN(s)
	cats->SetChannelWeight(13, 3./20. * CUSP_WEIGHT);	//3D1 SN(d) -> LN(d)
	cats->SetChannelWeight(15, 3./20. * CUSP_WEIGHT);	//3D1 SN(s) -> LN(d)
    }
    else if (target == "bonn")
    {
	auto cPars = make_unique<CATSparameters>(CATSparameters::tSource, 1, true);
	cPars->SetParameter(0, 1.2);

	cats->SetMomentumDependentSource(false);
	cats->SetThetaDependentSource(false);
	cats->SetAnaSource(GaussSource, *cPars);
	cats->SetUseAnalyticSource(true);
	cats->SetExcludeFailedBins(false);

	auto cPotPars1S0 = make_unique<CATSparameters>(CATSparameters::tPotential, 4, true);
	cPotPars1S0->SetParameter(0, 1304.76);
	cPotPars1S0->SetParameter(1, 2.607);
	cPotPars1S0->SetParameter(2, -113.946);
	cPotPars1S0->SetParameter(3, 0.9397);

	cats->SetQ1Q2(true);
	cats->SetPdgId(2212, 2212);
	cats->SetRedMass(0.5 * Mass_p);

	cats->SetNumChannels(4);
	cats->SetNumPW(0, 3);
	cats->SetNumPW(1, 2);
	cats->SetNumPW(2, 2);
	cats->SetNumPW(3, 2);
	cats->SetSpin(0, 0);
	cats->SetSpin(1, 1);
	cats->SetSpin(2, 1);
	cats->SetSpin(3, 1);
	cats->SetChannelWeight(0, 3. / 12.);
	cats->SetChannelWeight(1, 1. / 12.);
	cats->SetChannelWeight(2, 3. / 12.);
	cats->SetChannelWeight(3, 5. / 12.);

	auto double_gauss = [](double pars[]) { return pars[2]*exp(-pow(pars[0]*pars[3],2))+pars[4]*exp(-pow(pars[0]*pars[5],2)); };
	cats->SetShortRangePotential(0, 0, double_gauss, *cPotPars1S0);
    }

    cout<<"GG: Killing the Cat"<<endl;
    cats->KillTheCat();
}

template <typename T>
void setup_cats(DLM_CommonAnaFunctions *setupper, CATS *cats, VAR_FMR *range, T target, VAR *var)
{
    setup_cats(setupper, cats, range, target, var, nullptr, nullptr);
}

template <typename T>
void setup_cats(DLM_CommonAnaFunctions *setupper, CATS *cats, VAR_FMR *range, T target)
{
    setup_cats(setupper, cats, range, target, nullptr, nullptr, nullptr);
}

void setup_decomp(DLM_CkDecomposition *&decomp, DLM_Ck *&ck, CATS *cats, TH2F *res_matrix, VAR_FMR *range, const char *target)
{
    ck = new DLM_Ck(1, 0, *cats, range->Bins*2, range->Min, range->Max*2);
    ck->SetCutOff(range->Max, 700);
    ck->Update();

    const char *var1;
    int var2;
    if	    (!strcmp("pp",  target)) { var1 = "pp_Fitter";  var2 = 4; }
    else if (!strcmp("ppl", target)) { var1 = "pLambda";    var2 = 1; }
    else if (!strcmp("pl",  target)) { var1 = "pLambda";    var2 = 5; }
    else if (!strcmp("psp", target)) { var1 = "pSigmaPlus"; var2 = 1; }
    else if (!strcmp("ps0", target)) { var1 = "pSigma0";    var2 = 1; }
    else if (!strcmp("px0", target)) { var1 = "pXi0";	    var2 = 2; }
    else if (!strcmp("pxm", target)) { var1 = "pXim";	    var2 = 2; }
    else if (!strcmp("px1530", target)) { var1 = "pXim1530", var2 = 0; }

    decomp = new DLM_CkDecomposition(var1, var2, *ck, res_matrix);
}

void setup_decomp(DLM_CkDecomposition *&decomp, DLM_Ck *&ck, CATS *cats, unique_ptr<TH2F> &res_matrix, VAR_FMR *range, const char *target)
{
    setup_decomp(decomp, ck, cats, res_matrix.get(), range, target);
}
/*----------------------*/

void fix_fitter_pars(TF1 *fitter, vector<double> &parameters)
{
    for (size_t npar = 1; npar < parameters.size(); ++npar) fitter->FixParameter(npar, parameters[npar]);
}

void fix_fitter_radius(TF1 *fitter, double radius)
{
    fitter->FixParameter(6, radius);
}

void fix_global_fitter_pars(ROOT::Fit::Fitter *fitter, vector<double> &pars1, vector<double> &pars2)
{
    for (size_t npar = 1; npar < pars1.size(); ++npar)
    {
	fitter->Config().ParSettings(npar).SetValue(pars1[npar]);
	fitter->Config().ParSettings(npar).Fix();
    }

    for (size_t npar = 1; npar < pars2.size(); ++npar)
    {
	fitter->Config().ParSettings(pars1.size() + npar).SetValue(pars2[npar]);
	fitter->Config().ParSettings(pars1.size() + npar).Fix();
    }
}

void fix_global_fitter_radius(ROOT::Fit::Fitter *fitter, double radius)
{
    fitter->Config().ParSettings(10).SetValue(radius);
    fitter->Config().ParSettings(10).Fix();
}

void setup_fitter(TF1 *fitter, VAR *var, double radius)
{
    fitter->SetParameter(0, 1.);
    fitter->SetParameter(1, 1.);
    fitter->FixParameter(2, 0.);
    fitter->SetParameter(3, 1.);
    fitter->SetParameter(4, 1.);

    fitter->SetParameter(5, radius);
    fitter->SetParLimits(5, 0.4, 2);
}

void setup_global_fitter(ROOT::Fit::Fitter *fitter, VAR *var, double global_pars[])
{
    fitter->Config().SetParamsSettings(11, global_pars);
    fitter->Config().ParSettings(0).SetName("N pp");
    fitter->Config().ParSettings(1).SetName("a pp");
    fitter->Config().ParSettings(2).SetName("b pp");
    fitter->Config().ParSettings(3).SetName("c pp");
    fitter->Config().ParSettings(4).SetName("d pp");
    fitter->Config().ParSettings(5).SetName("N apap");
    fitter->Config().ParSettings(6).SetName("a apap");
    fitter->Config().ParSettings(7).SetName("b apap");
    fitter->Config().ParSettings(8).SetName("c apap");
    fitter->Config().ParSettings(9).SetName("d apap");
    fitter->Config().ParSettings(10).SetName("Radius");
    double radius = INIT_RADIUS[var->rsm][var->mult][var->mt];
    fitter->Config().ParSettings(10).SetValue(radius);
    //fitter->Config().ParSettings(10).Fix();

    fitter->Config().ParSettings(2).Fix();
    fitter->Config().ParSettings(7).Fix();

    fitter->Config().ParSettings(10).SetLimits(0.4, 2.0);

    for (size_t n = 1; n < 10; ++n)
    {
	if (n == 5) continue;
	fitter->Config().ParSettings(n).SetLimits(-1., 2.);
    }

    fitter->Config().MinimizerOptions().SetPrintLevel(0);
    fitter->Config().SetMinimizer("Minuit2", "Migrad");
}

void create_output_default(TFile *&file, TString iname, VAR *var)
{
    TString fname = iname + CHARGE_STR[var->charge];
    fname += (var->rsm)? "_rsm" : NULL;
    fname += Form("_mt%i_mult%i_fmr%i_fr%i_lam%i_bsl%i_smear%i",
	    var->mt, var->mult, var->fmr, var->fr, var->lam, var->bsl, var->smear);
    fname += (var->rsm)? Form("_frac%i_mass%i.root", var->frac, var->mass) : ".root";

    string settings = (var->system == PP)? "pp." : "pl.";
    file = new TFile(TString(var->get(settings + "output").data()) + fname, "recreate");
    printf("  \e[1;36m-->  \e[0;36mFile Created  \e[1;36m<--\e[0m\n\n");
    cout << "  " << var->get(settings + "output").data() << "\e[1;34m" << fname << "\e[0m\n\n";
}

void create_output_sample(TFile *&file, TString iname, VAR *var)
{
    TString fname = iname + CHARGE_STR[var->charge];
    fname += (var->rsm)? "_rsm" : NULL;
    if (var->sample && var->stat)
	fname += "_stat";
    fname += Form("_mt%i_mult%i_seed%i.root", var->mt, var->mult, var->sample);

    string settings = (var->system == PP)? "pp." : "pl.";
    TString opath = var->get(settings + "output").data();

    if (var->charge == 0) opath += "pp/";
    if (var->charge == 1) opath += "apap/";
    if (var->stat) opath += "stat/";
    if (var->rsm)  opath += "rsm/";

    file = new TFile(opath + fname, "recreate");
    printf("  \e[1;36m-->  \e[0;36mFile Created  \e[1;36m<--\e[0m\n\n");
    cout << "  " << opath << "\e[1;34m" << fname << "\e[0m\n\n";
}

void create_output(TFile *&file, TString iname, VAR *var)
{
    if (var->sample)
	create_output_sample(file, iname, var);
    else
	create_output_default(file, iname, var);
}

/* correlation fitter */
void cf_fitter(VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the fitter!  \e[1;36m<--\e[0m\n\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, PROTON, var->rsm, var->get("epos.path"));
    VAR_FMR *range_femto = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    if (var->mt > 4)
    {
	range_femto->Min = 8;
	range_fit->Min = 8;
    }
    if (var->mt == 2 && var->mult == 2)
    {
	range_femto->Min = 8;
	range_fit->Min = 8;
    }

    print_info_2(var, range_femto, range_fit);

    unique_ptr<TH1F> input_cf, input_me, input_me_orig;
    get_input(settings_filename(var), input_cf, input_me, input_me_orig, var, PP);

    unique_ptr<TH1F*> uptr_lambda_pars_th1 (new TH1F*[5]);
    auto lambda_pars_th1 = static_cast<TH1F**>(uptr_lambda_pars_th1.get());

    unique_ptr<DLM_Histo<double>*> uptr_lambda_pars (new DLM_Histo<double>*[5]);
    auto lambda_pars = static_cast<DLM_Histo<double>**>(uptr_lambda_pars.get());

    TString file_klambda = (var->charge)? var->get("pp.lampars_a") : var->get("pp.lampars_p");
    get_klambda_histos(file_klambda, var, lambda_pars_th1);
    get_klambda_histos(file_klambda, var, lambda_pars);

    TH2F *matrix_pp_res = nullptr;
    //get_resolution_matrix(project_path, matrix_pp_res, var, var->charge);

    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(var->get("cats.path").data());
    TH2F *matrix_pl_res = cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    TH2F *matrix_pl_dec = cats_setupper.GetResidualMatrix("pp", "pLambda");
    TH2F *matrix_ps_dec = cats_setupper.GetResidualMatrix("pp", "pSigmaPlus");

    DLM_CleverMcLevyResoTM MagicSource;
    DLM_CleverMcLevyResoTM *pMS = (var->rsm)? &MagicSource : NULL;
    VAR_RSM *pRSM = (var->rsm)? var_rsm : NULL;

    CATS cats_pp, cats_pl, cats_ps;
    setup_cats(&cats_setupper, &cats_pp, range_femto, settings_potential(var), var, pMS, pRSM);
    cout<<"Starting setup_cats"<<endl;
    setup_cats(&cats_setupper, &cats_pp, range_femto, target.Data(), pMS, pRSM);
    cout<<"GG: Done (1)"<<endl;
    setup_cats(&cats_setupper, &cats_pl, range_femto, "pl");
    cout<<"GG: Done (2)"<<endl;
    setup_cats(&cats_setupper, &cats_ps, range_femto, "psp");
    cout<<"GG: Done (3)"<<endl;

    DLM_Ck *ck_pp, *ck_pl, *ck_ps;
    cout<<"alive 0"<<endl;
    DLM_CkDecomposition *decomp_pp, *decomp_pl, *decomp_ps;
    setup_decomp(decomp_pp, ck_pp, &cats_pp, matrix_pp_res, range_femto, "pp");
    cout<<"GG: decomp Done (1)"<<endl;
    setup_decomp(decomp_pl, ck_pl, &cats_pl, matrix_pl_res, range_femto, "ppl");
    cout<<"GG: decomp Done (2)"<<endl;
    setup_decomp(decomp_ps, ck_ps, &cats_ps, matrix_pl_res, range_femto, "psp");
    cout<<"GG: decomp Done (3)"<<endl;

    // non-genuine contributions to pp
    decomp_pp->AddContribution(0, *lambda_pars[1], DLM_CkDecomp::cFeedDown, decomp_pl, matrix_pl_dec);
    decomp_pp->AddContribution(1, *lambda_pars[2], DLM_CkDecomp::cFeedDown, decomp_ps, matrix_ps_dec);
    decomp_pp->AddContribution(2, *lambda_pars[3], DLM_CkDecomp::cFeedDown);
    decomp_pp->AddContribution(3, *lambda_pars[4], DLM_CkDecomp::cFake);
    decomp_pp->AddPhaseSpace(input_me.get());
    decomp_pp->AddPhaseSpace(0, input_me.get());

    cout<<"alive 1"<<endl;
    // non-genuine contributions to pL
    decomp_pl->AddContribution(0, LAM_PL_FLAT, DLM_CkDecomp::cFeedDown);

    // update CkDec objects
    decomp_pp->Update(true, true);
    decomp_pl->Update(true, true);
    decomp_ps->Update(true, true);

    FitFunSimpleAvg fitfun(input_cf.get(), input_me_orig.get(), decomp_pp);
    TF1 *fitter = new TF1("Fit", fitfun, range_fit->Min, range_fit->Max, 6);
    TF1 *baseline = new TF1("Baseline", bslfun, range_fit->Min, range_fit->Max, 5);

    setup_fitter(fitter, var, 1.);
    input_cf->Fit(fitter, "S, N, R, M");

    double bsl_pars[5];
    for (int npar = 0; npar < 5; ++npar)
    {
	bsl_pars[npar] = fitter->GetParameter(npar);
	baseline->SetParameter(npar, fitter->GetParameter(npar));
    }

    GenuineFun genfun(decomp_pp);
    TF1* genuine = new TF1("Genuine", genfun, range_fit->Min, range_fit->Max, 1);
    genuine->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear(decomp_pp);
    TF1* genuine_smear = new TF1("Genuine Smeared", genfun_smear, range_fit->Min, range_fit->Max, 1);
    genuine_smear->SetParameter(0, 1.);

    decomp_pp->Update(true, true);
    TH1F *child_pl_smear = new TH1F("pl_smeared",   "pL smeared pp",    500, 0, 500);
    TH1F *child_ps_smear = new TH1F("ps_smeared",   "pS+ smeared pp",   500, 0, 500);
    for (short nbin = 0; nbin < 500; ++nbin)
    {
	child_pl_smear->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(0, child_pl_smear->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(1, child_ps_smear->GetBinCenter(nbin + 1) + 1) + 1);
    }

    TH1F **cfs_corr = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_corr[0] = static_cast<TH1F*>(input_cf->Clone());
    cfs_corr[1] = child_pl_smear;
    cfs_corr[2] = child_ps_smear;

    TH1F *corrected_cf = static_cast<TH1F*>(input_cf->Clone("CF_corrected"));
    corrected_cf->Reset();
    get_corrected_cf(cfs_corr, lambda_pars_th1, bsl_pars, &corrected_cf);

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", fitter->GetParameter(5), fitter->GetParError(5));
    printf("\e[1;34m  │  p-λ Size     │   \e[0m%.2f\n", decomp_pp->GetChild(0)->GetCk()->GetSourcePar(0));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", var);
    ofile->cd();

    TString str_charge = CHARGE_STR[var->charge];
    input_cf->Write("CF_" + str_charge);
    corrected_cf->Write("CF_corr_" + str_charge);
    input_me_orig->Write("ME_orig_" + str_charge);

    fitter->Write("fFitResult_" + str_charge);
    baseline->Write("fBl_" + str_charge);

    genuine->Write("fGenuine_" + str_charge);
    genuine_smear->Write("fGenuineSmeared_" + str_charge);

    child_pl_smear->Write("pl_smeared_" + str_charge);
    child_ps_smear->Write("ps_smeared_" + str_charge);

    TDirectory *tdir_lam = ofile->mkdir("lam_pars");
    tdir_lam->cd();

    lambda_pars_th1[0]->Write("genuine_" + str_charge);
    lambda_pars_th1[1]->Write("lambda_" + str_charge);
    lambda_pars_th1[2]->Write("sigma_" + str_charge);
    lambda_pars_th1[3]->Write("flat_" + str_charge);
    lambda_pars_th1[4]->Write("fake_" + str_charge);

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

void cf_fitter_pl(VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the fitter!  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Input File   │\e[0m   %s\e[0m\n", var->get("pl.file").data());
    printf("\e[1;34m  └───────────────┘\e[0m\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, LAMBDA, var->rsm, var->get("epos.path"));
    VAR_FMR *range_femto = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    range_femto->Max = 340;
    range_femto->Min = 0;

    range_fit->Max = 500;
    range_fit->Min = 0;

    print_info_2(var, range_femto, range_fit);

    TH1F *input_cf, *input_me, *input_me_orig;
    get_input_pl(string(var->get("pl.input") + var->get("pl.file")).data(), &input_cf, &input_me, &input_me_orig, var, var->charge);

    /* CATS Setup */
    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(var->get("cats.path").data());
    TH2F *matrix_pl_res = cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda");
    TH2F *matrix_pl_dec = cats_setupper.GetResidualMatrix("pp", "pLambda");
    TH2F *matrix_ls0_dec = cats_setupper.GetResidualMatrix("pLambda", "pSigma0");
    TH2F *matrix_lxm_dec = cats_setupper.GetResidualMatrix("pLambda", "pXim");
    TH2F *matrix_lx1530_dec = cats_setupper.GetResidualMatrix("pXim", "pXim1530");

    //vector<float> cusp_var {0.40, 0.27, 0.33}; // variations of the cusp

    /* CATS Objects */
    CATS cats_pl, cats_ps0, cats_px0, cats_pxm, cats_px1530;
    setup_cats(&cats_setupper, &cats_pl,  range_femto, "pl");
    setup_cats(&cats_setupper, &cats_ps0, range_femto, "ps0");
    setup_cats(&cats_setupper, &cats_px0, range_femto, "px0");
    setup_cats(&cats_setupper, &cats_pxm, range_femto, "pxm");
    setup_cats(&cats_setupper, &cats_px1530, range_femto, "px1530");

    /* Decomp Objects */
    DLM_Ck *ck_pl, *ck_ps0, *ck_px0, *ck_pxm, *ck_px1530;
    DLM_CkDecomposition *decomp_pl, *decomp_ps0, *decomp_px0, *decomp_pxm, *decomp_px1530;
    setup_decomp(decomp_pl,  ck_pl,  &cats_pl,  matrix_pl_res,	range_femto, "pl");
    setup_decomp(decomp_ps0, ck_ps0, &cats_ps0,	NULL,		range_femto, "ps0");
    setup_decomp(decomp_px0, ck_px0, &cats_px0, NULL,		range_femto, "px0");
    setup_decomp(decomp_pxm, ck_pxm, &cats_pxm, NULL,		range_femto, "pxm");
    setup_decomp(decomp_px1530, ck_px1530, &cats_px1530, NULL,	range_femto, "px1530");

    unique_ptr<double[]> lam_pars_pl  (new double[5]);
    unique_ptr<double[]> lam_pars_pxi (new double[5]);
    cats_setupper.SetUpLambdaPars_pL("pp13TeV_HM_Dec19", 0, 0, lam_pars_pl.get());
    cats_setupper.SetUpLambdaPars_pXim("pp13TeV_HM_Dec19", 0, 0, lam_pars_pxi.get());

    double lam_pXi_gen = lam_pars_pxi[0] / (1. - lam_pars_pxi[4]);
    double lam_pXi_pXi1530 = 0.5*(lam_pars_pxi[1]+lam_pars_pxi[2])/(1.-lam_pars_pxi[4]);
    double lam_pXi_flt = 1. - lam_pXi_gen - lam_pXi_pXi1530;

    double lam_gen = lam_pars_pl[0];
    double lam_pS0 = lam_pars_pl[1];
    double lam_pXim = lam_pars_pl[2];
    double lam_pXi0 = 0;
    double lam_flt = lam_pars_pl[3];
    double lam_mid = lam_pars_pl[4];

    double lam_pS0_flt = 0.18;

    printf(" lam_gen = %.2f%%\n",lam_gen*100.);
    printf(" lam_pS0 = %.2f%%\n",lam_pS0*100.);
    printf(" lam_pXim= %.2f%%\n",lam_pXim*100.);
    printf(" lam_pXi0= %.2f%%\n",lam_pXi0*100.);
    printf(" lam_flt = %.2f%%\n",lam_flt*100.);
    printf(" lam_mid = %.2f%%\n",lam_mid*100.);
    printf(" lam_pXi_gen = %.2f%%\n",lam_pXi_gen*100.);
    printf(" lam_pXi_pXi1530 = %.2f%%\n",lam_pXi_pXi1530*100.);
    printf(" lam_pXi_flt = %.2f%%\n",lam_pXi_flt*100.);
    printf(" lam_pS0_flt = %.2f%%\n",lam_pS0_flt*100.);

    decomp_pl->AddContribution(0, lam_pS0,  DLM_CkDecomp::cFeedDown, decomp_ps0, matrix_ls0_dec);
    decomp_pl->AddContribution(1, lam_pXim, DLM_CkDecomp::cFeedDown, decomp_pxm, matrix_lxm_dec);
    decomp_pl->AddContribution(2, lam_pXi0, DLM_CkDecomp::cFeedDown, decomp_px1530, matrix_lx1530_dec);
    decomp_pl->AddContribution(3, lam_flt,  DLM_CkDecomp::cFeedDown);
    decomp_pl->AddContribution(4, lam_mid,  DLM_CkDecomp::cFake);

    decomp_ps0->AddContribution(1, lam_pS0_flt, DLM_CkDecomposition::cFeedDown);

    decomp_pxm->AddContribution(0, lam_pXi_pXi1530, DLM_CkDecomposition::cFeedDown, decomp_px1530, matrix_lx1530_dec);
    decomp_pxm->AddContribution(1, lam_pXi_flt, DLM_CkDecomposition::cFeedDown);

    decomp_px0->AddContribution(0, lam_pXi_pXi1530, DLM_CkDecomposition::cFeedDown, decomp_px1530, matrix_lx1530_dec);
    decomp_px0->AddContribution(1, lam_pXi_flt, DLM_CkDecomposition::cFeedDown);

    ck_pl->Update();
    ck_ps0->Update();
    ck_px0->Update();
    ck_pxm->Update();
    ck_px1530->Update();

    //decomp_pl->AddContribution(0, LAM_PL["ps0"],  DLM_CkDecomp::cFeedDown, decomp_ps0, matrix_ls0_dec);
    //decomp_pl->AddContribution(1, LAM_PL["xim"],  DLM_CkDecomp::cFeedDown, decomp_pxm, matrix_lxm_dec);
    //decomp_pl->AddContribution(2, LAM_PL["flat"], DLM_CkDecomp::cFeedDown);
    //decomp_pl->AddContribution(3, LAM_PL["fake"], DLM_CkDecomp::cFake);
    //decomp_pl->AddPhaseSpace(input_me);
    //decomp_pl->Update(true, true);

    /* Fit Objects */
    plFitFunctionSimple fitfun(decomp_pl);
    TF1 *fitter = new TF1("Fit", fitfun, range_fit->Min, range_fit->Max, 6);
    TF1 *baseline = new TF1("Baseline", bslfun, range_fit->Min, range_fit->Max, 5);

    setup_fitter(fitter, var, 1.);
    input_cf->Fit(fitter, "S, N, R, M");

    double bsl_pars[5];
    for (int npar = 0; npar < 5; ++npar)
    {
	bsl_pars[npar] = fitter->GetParameter(npar);
	baseline->SetParameter(npar, fitter->GetParameter(npar));
    }

    GenuineFun genfun(decomp_pl);
    TF1* genuine = new TF1("Genuine", genfun, range_fit->Min, range_fit->Max, 1);
    genuine->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear(decomp_pl);
    TF1* genuine_smear = new TF1("Genuine Smeared", genfun_smear, range_fit->Min, range_fit->Max, 1);
    genuine_smear->SetParameter(0, 1.);

    decomp_pl->Update(true, true);

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", fitter->GetParameter(5), fitter->GetParError(5));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", var);
    ofile->cd();

    TString str_charge = CHARGE_STR[var->charge];
    input_cf->Write("CF_" + str_charge);
    input_me_orig->Write("ME_orig_" + str_charge);

    fitter->Write("fFitResult_" + str_charge);
    baseline->Write("fBl_" + str_charge);

    genuine->Write("fGenuine_" + str_charge);
    genuine_smear->Write("fGenuineSmeared_" + str_charge);

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

void cf_combined_fitter(VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the combined fitter!  \e[1;36m<--\e[0m\n\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, PROTON, var->rsm, var->get("epos.path"));
    VAR_FMR *range_femto_pp = new VAR_FMR(var->system, var->fmr);
    VAR_FMR *range_femto_aa = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit_pp = new VAR_FR(var->fr);
    VAR_FR  *range_fit_aa = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    // This setup is because we are missing the first cf bin for these mt and mult bins
    if (var->mt > 4)
    {
	range_femto_pp->Min = 8;
	range_femto_aa->Min = 8;
	range_fit_pp->Min = 8;
	range_fit_aa->Min = 8;
    }
    if (var->mt == 2 && var->mult == 2)
    {
	range_femto_aa->Min = 8;
	range_fit_aa->Min = 8;
    }

    print_info_2(var, range_femto_pp, range_fit_pp);

    // fitter parameter entry list, i.e. the last parameter is shared between both fitters
    int fit_pars_pp[6] = { 0, 1, 2, 3, 4, 10 };	    // N, a, b, c, d, r
    int fit_pars_aa[6] = { 5, 6, 7, 8, 9, 10 };
    double global_pars[11] = {0.98, 0, 0, 0, 0, 0.98, 0, 0, 0, 0, 1};	    // init parameters

    unique_ptr<TH1F> input_cf_pp, input_me_pp, input_me_orig_pp;
    unique_ptr<TH1F> input_cf_aa, input_me_aa, input_me_orig_aa;
    get_input(settings_filename(var), input_cf_pp, input_me_pp, input_me_orig_pp, var, PP);
    get_input(settings_filename(var), input_cf_aa, input_me_aa, input_me_orig_aa, var, APAP);

    unique_ptr<TH1F*[]> lambda_pars_pp_th1 (new TH1F*[5]);
    unique_ptr<TH1F*[]> lambda_pars_aa_th1 (new TH1F*[5]);
    get_klambda_histos(var->get("pp.lampars_p"), var, lambda_pars_pp_th1.get());
    get_klambda_histos(var->get("pp.lampars_a"), var, lambda_pars_aa_th1.get());

    unique_ptr<DLM_Histo<double>*[]> lambda_pars_pp (new DLM_Histo<double>*[5]);
    unique_ptr<DLM_Histo<double>*[]> lambda_pars_aa (new DLM_Histo<double>*[5]);
    get_klambda_histos(var->get("pp.lampars_p"), var, lambda_pars_pp.get());
    get_klambda_histos(var->get("pp.lampars_a"), var, lambda_pars_aa.get());

    TH2F *matrix_pp_res = nullptr;
    TH2F *matrix_aa_res = nullptr;
    //get_resolution_matrix(project_path, matrix_pp_res, var, PP);
    //get_resolution_matrix(project_path, matrix_aa_res, var, APAP);

    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(var->get("cats.path").data());

    unique_ptr<TH2F> matrix_pl_res (cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda"));
    unique_ptr<TH2F> matrix_pl_dec (cats_setupper.GetResidualMatrix("pp", "pLambda"));
    unique_ptr<TH2F> matrix_ps_dec (cats_setupper.GetResidualMatrix("pp", "pSigmaPlus"));

    DLM_CleverMcLevyResoTM MagicSource_pp, MagicSource_aa;
    DLM_CleverMcLevyResoTM *pMS_pp = (var->rsm)? &MagicSource_pp : nullptr;
    DLM_CleverMcLevyResoTM *pMS_aa = (var->rsm)? &MagicSource_aa : nullptr;
    VAR_RSM *pRSM = (var->rsm)? var_rsm : nullptr;
    TString target = settings_potential(var);
    target = "pp";
    //target = "av18_s";
    //target = "bonn";
    //target = "av18_s";
    //target = "reid93";
    //target = "reid68";

    cout<<"alive 0"<<endl;
    CATS cats_pp, cats_pl_pp, cats_ps_pp;
    setup_cats(&cats_setupper, &cats_pp,    range_femto_pp, target, var, pMS_pp, pRSM);
    setup_cats(&cats_setupper, &cats_pl_pp, range_femto_pp, "pl");
    setup_cats(&cats_setupper, &cats_ps_pp, range_femto_pp, "psp");

    cout<<"alive 1"<<endl;
    CATS cats_aa, cats_pl_aa, cats_ps_aa;
    setup_cats(&cats_setupper, &cats_aa,    range_femto_aa, target, var, pMS_aa, pRSM);
    setup_cats(&cats_setupper, &cats_pl_aa, range_femto_aa, "pl");
    setup_cats(&cats_setupper, &cats_ps_aa, range_femto_aa, "psp");

    cout<<"alive 2"<<endl;
    DLM_Ck *ck_pp, *ck_pl_pp, *ck_ps_pp;
    DLM_CkDecomposition *decomp_pp, *decomp_pl_pp, *decomp_ps_pp;
    setup_decomp(decomp_pp,    ck_pp,	 &cats_pp,    matrix_pp_res, range_femto_pp, "pp");
    setup_decomp(decomp_pl_pp, ck_pl_pp, &cats_pl_aa, matrix_pl_res, range_femto_pp, "ppl");
    setup_decomp(decomp_ps_pp, ck_ps_pp, &cats_ps_aa, matrix_pl_res, range_femto_pp, "psp");

    cout<<"alive 3"<<endl;
    DLM_Ck *ck_aa, *ck_pl_aa, *ck_ps_aa;
    DLM_CkDecomposition *decomp_aa, *decomp_pl_aa, *decomp_ps_aa;
    setup_decomp(decomp_aa,    ck_aa,	 &cats_aa,    matrix_aa_res, range_femto_aa, "pp");
    setup_decomp(decomp_pl_aa, ck_pl_aa, &cats_pl_aa, matrix_pl_res, range_femto_aa, "ppl");
    setup_decomp(decomp_ps_aa, ck_ps_aa, &cats_ps_aa, matrix_pl_res, range_femto_aa, "psp");

    // non-genuine contributions to pp
    decomp_pp->AddContribution(0, *lambda_pars_pp.get()[1], DLM_CkDecomp::cFeedDown, decomp_pl_pp, matrix_pl_dec.get());
    decomp_pp->AddContribution(1, *lambda_pars_pp.get()[2], DLM_CkDecomp::cFeedDown, decomp_ps_pp, matrix_ps_dec.get());
    decomp_pp->AddContribution(2, *lambda_pars_pp.get()[3], DLM_CkDecomp::cFeedDown);
    decomp_pp->AddContribution(3, *lambda_pars_pp.get()[4], DLM_CkDecomp::cFake);
    decomp_pp->AddPhaseSpace(input_me_pp.get());
    decomp_pp->AddPhaseSpace(0, input_me_pp.get());

    // non-genuine contributions to apap
    decomp_aa->AddContribution(0, *lambda_pars_aa.get()[1], DLM_CkDecomp::cFeedDown, decomp_pl_aa, matrix_pl_dec.get());
    decomp_aa->AddContribution(1, *lambda_pars_aa.get()[2], DLM_CkDecomp::cFeedDown, decomp_ps_aa, matrix_ps_dec.get());
    decomp_aa->AddContribution(2, *lambda_pars_aa.get()[3], DLM_CkDecomp::cFeedDown);
    decomp_aa->AddContribution(3, *lambda_pars_aa.get()[4], DLM_CkDecomp::cFake);
    decomp_aa->AddPhaseSpace(input_me_aa.get());
    decomp_aa->AddPhaseSpace(0, input_me_aa.get());

    // non-genuine contributions to pL
    decomp_pl_pp->AddContribution(0, LAM_PL_FLAT, DLM_CkDecomp::cFeedDown);
    decomp_pl_aa->AddContribution(0, LAM_PL_FLAT, DLM_CkDecomp::cFeedDown);

    // update CkDec objects
    decomp_pp->Update(true, true);
    decomp_pl_pp->Update(true, true);
    decomp_ps_pp->Update(true, true);

    decomp_aa->Update(true, true);
    decomp_pl_aa->Update(true, true);
    decomp_ps_aa->Update(true, true);

    // Fit function
    FitFunSimpleAvg fitfun_pp(input_cf_pp.get(), input_me_orig_pp.get(), decomp_pp, var->prefit);
    FitFunSimpleAvg fitfun_aa(input_cf_aa.get(), input_me_orig_aa.get(), decomp_aa, var->prefit);
    TF1 *fitter_pp = new TF1("Fit pp",	 fitfun_pp, range_fit_pp->Min, range_fit_pp->Max, 6);
    TF1 *fitter_aa = new TF1("Fit apap", fitfun_aa, range_fit_aa->Min, range_fit_aa->Max, 6);

    setup_fitter(fitter_pp, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);
    setup_fitter(fitter_aa, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);

    //vector<double> bsl_off {1, 0, 0, 0, 0};
    //fix_fitter_pars(fitter_pp, bsl_off);
    //fix_fitter_pars(fitter_aa, bsl_off);

    //fix_fitter_radius(fitter_pp, 1.10809);
    //fix_fitter_radius(fitter_aa, 1.10809);

    if (var->prefit)
    {
	fix_fitter_pars(fitter_pp, fitfun_pp.prefit_pars);
	fix_fitter_pars(fitter_aa, fitfun_aa.prefit_pars);
    }

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange opt_range_pp, opt_range_aa;
    opt_range_pp.SetRange(range_fit_pp->Min, range_fit_pp->Max);
    opt_range_aa.SetRange(range_fit_aa->Min, range_fit_aa->Max);

    ROOT::Fit::BinData wrap_data_pp(opt, opt_range_pp);
    ROOT::Fit::BinData wrap_data_aa(opt, opt_range_aa);

    ROOT::Math::WrappedMultiTF1 wrap_fit_pp(*fitter_pp, 1);
    ROOT::Math::WrappedMultiTF1 wrap_fit_aa(*fitter_aa, 1);

    ROOT::Fit::FillData(wrap_data_pp, input_cf_pp.get());
    ROOT::Fit::FillData(wrap_data_aa, input_cf_aa.get());

    ROOT::Fit::Chi2Function chi2_pp(wrap_data_pp, wrap_fit_pp);
    ROOT::Fit::Chi2Function chi2_aa(wrap_data_aa, wrap_fit_aa);

    GlobalChi2 global_chi2(chi2_pp, chi2_aa, fit_pars_pp, fit_pars_aa);

    ROOT::Fit::Fitter fitter;
    setup_global_fitter(&fitter, var, global_pars);
    //fix_global_fitter_pars(&fitter, bsl_off, bsl_off);
    //fix_global_fitter_radius(&fitter, 1.10809);

    if (var->prefit) fix_global_fitter_pars(&fitter, fitfun_pp.prefit_pars, fitfun_aa.prefit_pars);

    fitter.FitFCN(11, global_chi2, nullptr, wrap_data_pp.Size() + wrap_data_aa.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();

    fitter_pp->SetFitResult(result, fit_pars_pp);
    fitter_aa->SetFitResult(result, fit_pars_aa);

    double temp, total_chi2 = 0.;
    for (size_t nbin = input_cf_pp->FindBin(50); input_cf_pp->GetBinCenter(nbin) < 140; ++nbin)
    {
	temp = (pow(input_cf_pp->GetBinContent(nbin) - fitter_pp->Eval(input_cf_pp->GetBinCenter(nbin)), 2)) / pow(input_cf_pp->GetBinError(nbin), 2);
	total_chi2 += temp;
    }
    //cout << "pp fit Chi2 = " << total_chi2 << endl;
    out("pp fit Chi2 = " << total_chi2);

    fitter_pp->SetRange(opt_range_pp().first, opt_range_pp().second);
    fitter_aa->SetRange(opt_range_aa().first, opt_range_aa().second);

    double bsl_pars_pp[5], bsl_pars_aa[5];
    TF1 *baseline_pp = new TF1("Baseline pp",   bslfun, range_fit_pp->Min, range_fit_pp->Max, 5);
    TF1 *baseline_aa = new TF1("Baseline apap", bslfun, range_fit_aa->Min, range_fit_aa->Max, 5);
    for (size_t npar = 0; npar < 5; ++npar)
    {
	bsl_pars_pp[npar] = fitter_pp->GetParameter(npar);
	bsl_pars_aa[npar] = fitter_aa->GetParameter(npar);
	baseline_pp->SetParameter(npar, fitter_pp->GetParameter(npar));
	baseline_aa->SetParameter(npar, fitter_aa->GetParameter(npar));
    }

    GenuineFun genfun_pp(decomp_pp), genfun_aa(decomp_aa);
    TF1* genuine_pp = new TF1("Genuine pp", genfun_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_aa = new TF1("Genuine apap", genfun_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_aa->SetParameter(0, 1.);
    genuine_pp->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear_pp(decomp_pp), genfun_smear_aa(decomp_aa);
    TF1* genuine_smear_pp = new TF1("Genuine Smeared pp", genfun_smear_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_smear_aa = new TF1("Genuine Smeared apap", genfun_smear_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_smear_pp->SetParameter(0, 1.);
    genuine_smear_aa->SetParameter(0, 1.);

    decomp_pp->Update(true, true);
    decomp_aa->Update(true, true);
    TH1F *child_pl_smear_pp = new TH1F("pl_smeared_pp",   "pL smeared pp",    500, 0, 500);
    TH1F *child_pl_smear_aa = new TH1F("pl_smeared_apap", "pL smeared apap",  500, 0, 500);
    TH1F *child_ps_smear_pp = new TH1F("ps_smeared_pp",   "pS+ smeared pp",   500, 0, 500);
    TH1F *child_ps_smear_aa = new TH1F("ps_smeared_apap", "pS+ smeared apap", 500, 0, 500);
    for (short nbin = 0; nbin < 500; ++nbin)
    {
	child_pl_smear_pp->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(0, child_pl_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_pl_smear_aa->SetBinContent(nbin + 1, decomp_aa->EvalSignalSmearedChild(0, child_pl_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_pp->SetBinContent(nbin + 1, decomp_pp->EvalSignalSmearedChild(1, child_ps_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_aa->SetBinContent(nbin + 1, decomp_aa->EvalSignalSmearedChild(1, child_ps_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
    }

    /* create and fill histos of the total wave function for k* values */
    vector<double> wf_kstar {10, 20, 40, 80, 160};
    unique_ptr<TH1F*[]> twave_histos_s (new TH1F*[wf_kstar.size()]);
    unique_ptr<TH1F*[]> twave_histos_p (new TH1F*[wf_kstar.size()]);
    // channel 0 -> anti-symmetric wave function
    if (var->save_tw) get_cats_totalwaves(&cats_pp, twave_histos_s.get(), wf_kstar, 0);
    if (var->save_tw) get_cats_totalwaves(&cats_pp, twave_histos_p.get(), wf_kstar, 1);

    /* create and fill histos of partial waves for given k* values */
    vector<const char*> pwaves_names = {"1S0", "3P0", "3P1", "3P2", "1D2"};
    vector<double> pwaves_kstars {10, 20, 40, 80, 160};
    unique_ptr<unique_ptr<TH1F*[]>[]> pwaves_histos (new unique_ptr<TH1F*[]>[pwaves_names.size()]);
    for (size_t npwave = 0; npwave < pwaves_names.size(); ++npwave)
    {
	pwaves_histos[npwave] = make_unique<TH1F*[]>(pwaves_kstars.size());
	if (target == "epelbaum")
	    get_cats_radialwaves_epelbaum(&cats_pp, pwaves_histos[npwave].get(), pwaves_kstars, pwaves_names[npwave], "real");
	else
	    get_cats_radialwaves(&cats_pp, pwaves_histos[npwave].get(), pwaves_kstars, pwaves_names[npwave], "real");
    }

    /* create and fill histos of phaseshifts */
    vector<const char*> phase_names = {"1S0", "3P0", "3P1", "3P2", "1D2"};
    unique_ptr<TH1F*[]> phase_histos (new TH1F*[phase_names.size()]);
    if (var->save_ps)
    {
	if (target == "epelbaum")
	    get_cats_phaseshifts_epelbaum(&cats_pp, phase_histos.get(), phase_names, range_femto_pp->Min, range_femto_pp->Max);
	else
	    get_cats_phaseshifts(&cats_pp, phase_histos.get(), phase_names, range_femto_pp->Min, range_femto_pp->Max);
    }

    /* create and fill histos of potentials for partial waves */
    vector<const char*> pot_names = {"1S0", "3P0", "3P1", "3P2", "1D2"};
    unique_ptr<TH1F*[]> pot_histos (new TH1F*[pot_names.size()]);
    if (var->save_pot) get_cats_potentials(&cats_pp, pot_histos.get(), 40, pot_names); // k* = 40 MeV

    TH1F **cfs_pp = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_pp[0] = static_cast<TH1F*>(input_cf_pp->Clone());
    cfs_pp[1] = child_pl_smear_pp;
    cfs_pp[2] = child_ps_smear_pp;

    TH1F **cfs_aa = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_aa[0] = static_cast<TH1F*>(input_cf_aa->Clone());
    cfs_aa[1] = child_pl_smear_aa;
    cfs_aa[2] = child_ps_smear_aa;

    TH1F *corrected_cf_pp = (TH1F*) input_cf_pp->Clone("CF_corrected_pp");
    corrected_cf_pp->Reset();
    get_corrected_cf(cfs_pp, lambda_pars_pp_th1.get(), bsl_pars_pp, &corrected_cf_pp);

    TH1F *corrected_cf_aa = (TH1F*) input_cf_aa->Clone("CF_corrected_aa");
    corrected_cf_aa->Reset();
    get_corrected_cf(cfs_aa, lambda_pars_aa_th1.get(), bsl_pars_aa, &corrected_cf_aa);

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    result.Print(cout); printf("\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", result.Parameter(10), result.Error(10));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", var);
    ofile->cd();

    input_cf_pp->Write("CF_pp");
    input_cf_aa->Write("CF_apap");

    corrected_cf_pp->Write("CF_corr_pp");
    corrected_cf_aa->Write("CF_corr_apap");

    input_me_orig_pp->Write("ME_orig_pp");
    input_me_orig_aa->Write("ME_orig_apap");

    fitter_pp->Write("fFitResult_pp");
    fitter_aa->Write("fFitResult_apap");

    baseline_pp->Write("fBl_pp");
    baseline_aa->Write("fBl_apap");

    genuine_pp->Write("fGenuine_pp");
    genuine_aa->Write("fGenuine_apap");

    genuine_smear_pp->Write("fGenuineSmeared_pp");
    genuine_smear_aa->Write("fGenuineSmeared_apap");

    child_pl_smear_pp->Write("pl_smeared_pp");
    child_pl_smear_aa->Write("pl_smeared_apap");
    child_ps_smear_pp->Write("ps_smeared_pp");
    child_ps_smear_aa->Write("ps_smeared_apap");

    /* save channel wave */
    if (var->save_tw)
    {
	TDirectory *tdir_tw = ofile->mkdir("wave_function");
	tdir_tw->cd();

	for (size_t nwf = 0; nwf < wf_kstar.size(); ++nwf)
	{
	    twave_histos_s[nwf]->Write();
	    twave_histos_p[nwf]->Write();
	}
    }

    /* save radial waves */
    if (var->save_pw)
    {
	TDirectory *tdir_pw = ofile->mkdir("partial_waves");
	tdir_pw->cd();

	for (size_t nkstar = 0; nkstar < pwaves_kstars.size(); ++nkstar)
	{
	    for (size_t npwave = 0; npwave < pwaves_names.size(); ++npwave)
	    {
		pwaves_histos[npwave][nkstar]->Write();
	    }
	}
    }

    /* save phase shifts */
    if (var->save_ps)
    {
	TDirectory *tdir_ps = ofile->mkdir("phase_shifts");
	tdir_ps->cd();

	for (size_t nphase = 0; nphase < phase_names.size(); ++nphase)
	{
	    phase_histos[nphase]->Write();
	}
    }

    /* save potentials */
    if (var->save_pot)
    {
	TDirectory *tdir_pot = ofile->mkdir("potentials");
	tdir_pot->cd();

	for (size_t npot = 0; npot < pot_names.size(); ++npot)
	{
	    pot_histos[npot]->Write();
	}
    }

    TDirectory *tdir_lam = ofile->mkdir("lam_pars");
    tdir_lam->cd();

    lambda_pars_pp_th1.get()[0]->Write("genuine_pp");
    lambda_pars_pp_th1.get()[1]->Write("lambda_pp");
    lambda_pars_pp_th1.get()[2]->Write("sigma_pp");
    lambda_pars_pp_th1.get()[3]->Write("flat_pp");
    lambda_pars_pp_th1.get()[4]->Write("fake_pp");

    lambda_pars_aa_th1.get()[0]->Write("genuine_apap");
    lambda_pars_aa_th1.get()[1]->Write("lambda_apap");
    lambda_pars_aa_th1.get()[2]->Write("sigma_apap");
    lambda_pars_aa_th1.get()[3]->Write("flat_apap");
    lambda_pars_aa_th1.get()[4]->Write("fake_apap");

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

void cf_combined_fitter_pl(VAR *var)
{
    printf("\n  \e[1;36m-->  \e[0;36mEntering the combined fitter!  \e[1;36m<--\e[0m\n\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Input File   │\e[0m   %s\e[0m\n", var->get("pl.file").data());
    printf("\e[1;34m  └───────────────┘\e[0m\n");

    VAR_RSM *var_rsm = new VAR_RSM(PROTON, LAMBDA, var->rsm, var->get("epos.path"));
    VAR_FMR *range_femto_pp = new VAR_FMR(var->system, var->fmr);
    VAR_FMR *range_femto_aa = new VAR_FMR(var->system, var->fmr);
    VAR_FR  *range_fit_pp = new VAR_FR(var->fr);
    VAR_FR  *range_fit_aa = new VAR_FR(var->fr);
    VAR_QA  *qa_plots = new VAR_QA(var);

    //range_fit_pp->Max = 300;
    //range_fit_pp->Min = 0;
    //range_fit_aa->Max = 300;
    //range_fit_aa->Min = 0;

    //range_femto_pp->Max = 300;
    //range_femto_pp->Min = 0;
    //range_femto_aa->Max = 300;
    //range_femto_aa->Min = 0;

    range_femto_pp->Max = range_fit_pp->Max;
    range_femto_pp->Min = range_fit_pp->Min;
    range_femto_aa->Max = range_fit_aa->Max;
    range_femto_aa->Min = range_fit_aa->Min;

    print_info_2(var, range_femto_pp, range_fit_pp);

    // fitter parameter entry list, i.e. the last parameter is shared between both fitters
    int fit_pars_pp[6] = { 0, 1, 2, 3, 4, 10 };	    // N, a, b, c, d, r
    int fit_pars_aa[6] = { 5, 6, 7, 8, 9, 10 };
    double global_pars[11] = {0.98, 0, 0, 0, 0, 0.98, 0, 0, 0, 0, 1};	    // init parameters

    unique_ptr<TH1F> input_cf_pp, input_me_pp, input_me_orig_pp;
    unique_ptr<TH1F> input_cf_aa, input_me_aa, input_me_orig_aa;
    get_input_pl(string(var->get("pl.input") + var->get("pl.file")).data(), input_cf_pp, input_me_pp, input_me_orig_pp, var, PP);
    get_input_pl(string(var->get("pl.input") + var->get("pl.file")).data(), input_cf_aa, input_me_aa, input_me_orig_aa, var, APAP);

    DLM_CommonAnaFunctions cats_setupper;
    cats_setupper.SetCatsFilesFolder(var->get("cats.path").data());

    unique_ptr<TH2F> matrix_pl_res  (cats_setupper.GetResolutionMatrix("pp13TeV_HM_DimiJun20","pLambda"));
    unique_ptr<TH2F> matrix_pl_dec  (cats_setupper.GetResidualMatrix("pp", "pLambda"));
    unique_ptr<TH2F> matrix_ls0_dec (cats_setupper.GetResidualMatrix("pLambda", "pSigma0"));
    unique_ptr<TH2F> matrix_lxm_dec (cats_setupper.GetResidualMatrix("pLambda", "pXim"));

    DLM_CleverMcLevyResoTM MagicSource_pp, MagicSource_aa;
    DLM_CleverMcLevyResoTM *pMS_pp = (var->rsm)? &MagicSource_pp : NULL;
    DLM_CleverMcLevyResoTM *pMS_aa = (var->rsm)? &MagicSource_aa : NULL;
    VAR_RSM *pRSM = (var->rsm)? var_rsm : NULL;
    TString target = settings_potential(var);

    /* CATS Objects */
    CATS cats_pl_pp, cats_ps0_pp, cats_px0_pp, cats_pxm_pp;
    setup_cats(&cats_setupper, &cats_pl_pp,  range_femto_pp, target, var, pMS_pp, pRSM);
    setup_cats(&cats_setupper, &cats_ps0_pp, range_femto_pp, "ps0");
    setup_cats(&cats_setupper, &cats_px0_pp, range_femto_pp, "px0");
    setup_cats(&cats_setupper, &cats_pxm_pp, range_femto_pp, "pxm");

    CATS cats_pl_aa, cats_ps0_aa, cats_px0_aa, cats_pxm_aa;
    setup_cats(&cats_setupper, &cats_pl_aa,  range_femto_aa, target, var, pMS_pp, pRSM);
    setup_cats(&cats_setupper, &cats_ps0_aa, range_femto_aa, "ps0");
    setup_cats(&cats_setupper, &cats_px0_aa, range_femto_aa, "px0");
    setup_cats(&cats_setupper, &cats_pxm_aa, range_femto_aa, "pxm");

    /* Decomp Objects */
    DLM_Ck *ck_pl_pp, *ck_ps0_pp, *ck_px0_pp, *ck_pxm_pp;
    DLM_CkDecomposition *decomp_pl_pp, *decomp_ps0_pp, *decomp_px0_pp, *decomp_pxm_pp;
    setup_decomp(decomp_pl_pp,  ck_pl_pp,  &cats_pl_pp,  matrix_pl_res, range_femto_pp, "pl");
    setup_decomp(decomp_ps0_pp, ck_ps0_pp, &cats_ps0_pp, NULL,		range_femto_pp, "ps0");
    setup_decomp(decomp_px0_pp, ck_px0_pp, &cats_px0_pp, NULL,		range_femto_pp, "px0");
    setup_decomp(decomp_pxm_pp, ck_pxm_pp, &cats_pxm_pp, NULL,		range_femto_pp, "pxm");

    DLM_Ck *ck_pl_aa, *ck_ps0_aa, *ck_px0_aa, *ck_pxm_aa;
    DLM_CkDecomposition *decomp_pl_aa, *decomp_ps0_aa, *decomp_px0_aa, *decomp_pxm_aa;
    setup_decomp(decomp_pl_aa,  ck_pl_aa,  &cats_pl_aa,  matrix_pl_res, range_femto_aa, "pl");
    setup_decomp(decomp_ps0_aa, ck_ps0_aa, &cats_ps0_aa, NULL,		range_femto_aa, "ps0");
    setup_decomp(decomp_px0_aa, ck_px0_aa, &cats_px0_aa, NULL,		range_femto_aa, "px0");
    setup_decomp(decomp_pxm_aa, ck_pxm_aa, &cats_pxm_aa, NULL,		range_femto_aa, "pxm");

    decomp_pl_pp->AddContribution(0, LAM_PL["ps0"],  DLM_CkDecomp::cFeedDown, decomp_ps0_pp, matrix_ls0_dec.get());
    decomp_pl_pp->AddContribution(1, LAM_PL["xim"],  DLM_CkDecomp::cFeedDown, decomp_pxm_pp, matrix_lxm_dec.get());
    decomp_pl_pp->AddContribution(2, LAM_PL["flat"], DLM_CkDecomp::cFeedDown);
    decomp_pl_pp->AddContribution(3, LAM_PL["fake"], DLM_CkDecomp::cFake);

    decomp_pl_aa->AddContribution(0, LAM_PL["ps0"],  DLM_CkDecomp::cFeedDown, decomp_ps0_aa, matrix_ls0_dec.get());
    decomp_pl_aa->AddContribution(1, LAM_PL["xim"],  DLM_CkDecomp::cFeedDown, decomp_pxm_aa, matrix_lxm_dec.get());
    decomp_pl_aa->AddContribution(2, LAM_PL["flat"], DLM_CkDecomp::cFeedDown);
    decomp_pl_aa->AddContribution(3, LAM_PL["fake"], DLM_CkDecomp::cFake);

    // update CkDec objects
    decomp_pl_pp->Update(true, true);
    decomp_ps0_pp->Update(true, true);
    decomp_px0_pp->Update(true, true);
    decomp_pxm_pp->Update(true, true);

    decomp_pl_aa->Update(true, true);
    decomp_ps0_aa->Update(true, true);
    decomp_px0_aa->Update(true, true);
    decomp_pxm_aa->Update(true, true);

    // Fit function
    plFitFunctionSimple fitfun_pp(decomp_pl_pp);
    plFitFunctionSimple fitfun_aa(decomp_pl_aa);
    TF1 *fitter_pp = new TF1("Fit pp",	 fitfun_pp, range_fit_pp->Min, range_fit_pp->Max, 6);
    TF1 *fitter_aa = new TF1("Fit apap", fitfun_aa, range_fit_aa->Min, range_fit_aa->Max, 6);

    setup_fitter(fitter_pp, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);
    setup_fitter(fitter_aa, var, INIT_RADIUS[var->rsm][var->mult][var->mt]);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange opt_range_pp, opt_range_aa;
    opt_range_pp.SetRange(range_fit_pp->Min, range_fit_pp->Max);
    opt_range_aa.SetRange(range_fit_aa->Min, range_fit_aa->Max);

    ROOT::Fit::BinData wrap_data_pp(opt, opt_range_pp);
    ROOT::Fit::BinData wrap_data_aa(opt, opt_range_aa);

    ROOT::Math::WrappedMultiTF1 wrap_fit_pp(*fitter_pp, 1);
    ROOT::Math::WrappedMultiTF1 wrap_fit_aa(*fitter_aa, 1);

    ROOT::Fit::FillData(wrap_data_pp, input_cf_pp.get());
    ROOT::Fit::FillData(wrap_data_aa, input_cf_aa.get());

    ROOT::Fit::Chi2Function chi2_pp(wrap_data_pp, wrap_fit_pp);
    ROOT::Fit::Chi2Function chi2_aa(wrap_data_aa, wrap_fit_aa);

    GlobalChi2 global_chi2(chi2_pp, chi2_aa, fit_pars_pp, fit_pars_aa);

    ROOT::Fit::Fitter fitter;
    setup_global_fitter(&fitter, var, global_pars);

    fitter.FitFCN(11, global_chi2, nullptr, wrap_data_pp.Size() + wrap_data_aa.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();

    fitter_pp->SetFitResult(result, fit_pars_pp);
    fitter_aa->SetFitResult(result, fit_pars_aa);

    double temp, total_chi2 = 0.;
    for (size_t nbin = input_cf_pp->FindBin(50); input_cf_pp->GetBinCenter(nbin) < 140; ++nbin)
    {
	temp = (pow(input_cf_pp->GetBinContent(nbin) - fitter_pp->Eval(input_cf_pp->GetBinCenter(nbin)), 2)) / pow(input_cf_pp->GetBinError(nbin), 2);
	total_chi2 += temp;
    }
    //cout << "pp fit Chi2 = " << total_chi2 << endl;
    out("pp fit Chi2 = " << total_chi2);

    fitter_pp->SetRange(opt_range_pp().first, opt_range_pp().second);
    fitter_aa->SetRange(opt_range_aa().first, opt_range_aa().second);

    double bsl_pars_pp[5], bsl_pars_aa[5];
    TF1 *baseline_pp = new TF1("Baseline pp",   bslfun, range_fit_pp->Min, range_fit_pp->Max, 5);
    TF1 *baseline_aa = new TF1("Baseline apap", bslfun, range_fit_aa->Min, range_fit_aa->Max, 5);
    for (size_t npar = 0; npar < 5; ++npar)
    {
	bsl_pars_pp[npar] = fitter_pp->GetParameter(npar);
	bsl_pars_aa[npar] = fitter_aa->GetParameter(npar);
	baseline_pp->SetParameter(npar, fitter_pp->GetParameter(npar));
	baseline_aa->SetParameter(npar, fitter_aa->GetParameter(npar));
    }

    GenuineFun genfun_pp(decomp_pl_pp), genfun_aa(decomp_pl_aa);
    TF1* genuine_pp = new TF1("Genuine pp", genfun_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_aa = new TF1("Genuine apap", genfun_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_aa->SetParameter(0, 1.);
    genuine_pp->SetParameter(0, 1.);

    GenuineFunSmear genfun_smear_pp(decomp_pl_pp), genfun_smear_aa(decomp_pl_aa);
    TF1* genuine_smear_pp = new TF1("Genuine Smeared pp", genfun_smear_pp, range_fit_pp->Min, range_fit_pp->Max, 1);
    TF1* genuine_smear_aa = new TF1("Genuine Smeared apap", genfun_smear_aa, range_fit_aa->Min, range_fit_aa->Max, 1);
    genuine_smear_pp->SetParameter(0, 1.);
    genuine_smear_aa->SetParameter(0, 1.);

    decomp_pl_pp->Update(true, true);
    decomp_pl_aa->Update(true, true);
    TH1F *child_pl_smear_pp = new TH1F("ps0_smeared_pp",   "pS0 smeared pp",    500, 0, range_fit_pp->Max);
    TH1F *child_pl_smear_aa = new TH1F("ps0_smeared_apap", "pS0 smeared apap",  500, 0, range_fit_aa->Max);
    TH1F *child_ps_smear_pp = new TH1F("pxm_smeared_pp",   "pS- smeared pp",   500, 0, range_fit_pp->Max);
    TH1F *child_ps_smear_aa = new TH1F("pxm_smeared_apap", "pS- smeared apap", 500, 0, range_fit_aa->Max);
    for (short nbin = 0; nbin < 500; ++nbin)
    {
	child_pl_smear_pp->SetBinContent(nbin + 1, decomp_pl_pp->EvalSignalSmearedChild(0, child_pl_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_pl_smear_aa->SetBinContent(nbin + 1, decomp_pl_aa->EvalSignalSmearedChild(0, child_pl_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_pp->SetBinContent(nbin + 1, decomp_pl_pp->EvalSignalSmearedChild(1, child_ps_smear_pp->GetBinCenter(nbin + 1) + 1) + 1);
	child_ps_smear_aa->SetBinContent(nbin + 1, decomp_pl_aa->EvalSignalSmearedChild(1, child_ps_smear_aa->GetBinCenter(nbin + 1) + 1) + 1);
    }

    TH1F **cfs_pp = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_pp[0] = static_cast<TH1F*>(input_cf_pp->Clone());
    cfs_pp[1] = child_pl_smear_pp;
    cfs_pp[2] = child_ps_smear_pp;

    TH1F **cfs_aa = (TH1F**) malloc(sizeof(TH1F[3]));
    cfs_aa[0] = static_cast<TH1F*>(input_cf_aa->Clone());
    cfs_aa[1] = child_pl_smear_aa;
    cfs_aa[2] = child_ps_smear_aa;

    printf("\n\e[1;36m  -->  \e[0;36mData Fit  \e[1;36m<--\e[0m\n\n");
    result.Print(cout); printf("\n");
    printf("\e[1;34m  ┌───────────────┐\e[0m\n");
    printf("\e[1;34m  │  Source Size  │   \e[0m%.2f ± %.2f\n", result.Parameter(10), result.Error(10));
    printf("\e[1;34m  └───────────────┘\e[0m\n\n");

    TFile *ofile;
    create_output(ofile, "FitResults_", var);
    ofile->cd();

    input_cf_pp->Write("CF_pp");
    input_cf_aa->Write("CF_apap");

    input_me_orig_pp->Write("ME_orig_pp");
    input_me_orig_aa->Write("ME_orig_apap");

    fitter_pp->Write("fFitResult_pp");
    fitter_aa->Write("fFitResult_apap");

    baseline_pp->Write("fBl_pp");
    baseline_aa->Write("fBl_apap");

    genuine_pp->Write("fGenuine_pp");
    genuine_aa->Write("fGenuine_apap");

    genuine_smear_pp->Write("fGenuineSmeared_pp");
    genuine_smear_aa->Write("fGenuineSmeared_apap");

    child_pl_smear_pp->Write("ps0_smeared_pp");
    child_pl_smear_aa->Write("ps0_smeared_apap");
    child_ps_smear_pp->Write("pxm_smeared_pp");
    child_ps_smear_aa->Write("pxm_smeared_apap");

    TDirectory *tdir_qa = ofile->mkdir("qa_plots");
    tdir_qa->cd();

    qa_plots->lam->Write();
    qa_plots->fmr->Write();
    qa_plots->fr->Write();
    qa_plots->bsl->Write();
    qa_plots->smear->Write();

    if (var->rsm)
    {
	qa_plots->frac->Write();
	qa_plots->mass->Write();
    }

    ofile->Close();
}

